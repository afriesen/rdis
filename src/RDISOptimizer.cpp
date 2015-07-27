/*
 * OptimizerDG.cpp
 *
 *  Created on: Jul 25, 2013
 *      Author: afriesen
 */

#include "RDISOptimizer.h"
#include "util/numeric.h"
#include "util/utility.h"

#ifdef USE_HMETIS
//#include "external/include/hmetis.h"
#include "hmetis.h"
#endif // USE_HMETIS

#ifdef USE_PATOH
//#include "external/include/patoh/patoh.h"
#include "patoh.h"
#endif // USE_PATOH

#include BOOSTPATH/bind.hpp>

namespace rdis {

using namespace NumericVecOps;
namespace fs = BOOSTNS::filesystem;
namespace po = BOOSTNS::program_options;

RDISOptimizer::RDISOptimizer(
		OptimizableFunction & of,
		SubspaceOptimizer & ssopt_,
		po::variables_map options_,
		size_t cacheSize )
	: super( of )
	, ccpool( -1 )
	, rng( 1645333507 )
	, ssopt( ssopt_ )
	, options( options_ )
	, vardata( of.getNumVars() )
	, connGraph( of.getNumVars() )
	, epsilon( 0.5 )
	, cgNeedsInit( true )
	, nSimplifiedFactors( 0 )
	, bestFullFunctionEval( NumericMAX )
	, prevFullFunctionEval( NumericMAX )
	, approxFactors( true )
	, useCache( false )
	, useLinearCache( false )
	, useBounds( of.hasBounds() )
	, numassigned( 0 )
	, assignmentOrder( "" )
	, printedAO( false )
{
	m_timerCheckFrequency = 1;

	if ( m_func.semiring().Sum( 2, 6 ) != 2 ||
			m_func.semiring().Product( 2, 6 ) != 8 ) {
		throw "Unsupported semiring type in OptimizableFunction -- only the "
				"MinSum semiring is supported.";
	}

	processOptions();
	init( cacheSize );
}

RDISOptimizer::~RDISOptimizer() {}


void RDISOptimizer::addOptions( po::options_description & desc ) {
	desc.add_options()
			( "help,h", "Print help message" )
			( "epsilon,e", po::value< Numeric >()->default_value( 0.5 ),
					"error tolerance of optimizer" )
			( "timeout", po::value< Numeric >(),
					"timeout of entire optimization, in seconds" )
			( "doBnB", po::value< bool >()->default_value( true ),
					"true to perform branch and bound during the optimization" )
			( "approxfactors", po::value< bool >()->default_value( true ),
					"flag to enable/disable approximation of low-error factors" )
			( "timeAsSeed", po::value< bool >()->default_value( false ),
					"true to use current time as the random seed, false to use "
							"pre-defined constants" )
			( "staticDecomp", po::value< bool >()->default_value( false ),
					"true to decompose the function once and then re-use that "
							"decomposition each time" )
			( "steptol", po::value< Numeric >()->default_value( 1e-4 ),
					"tolerance to determine convergence between steps of our "
							"algorithm (at each level)" )
			( "varChoiceHeur", po::value< string >()->default_value( "patoh" ),
					"the heuristic to use for choosing variables (one of: "
							"patoh, hmetis, random, all)" )
			( "AVblksz", po::value< uint32 >(),
					"size of blocks of vars to assign at each level" )
			( "AVblkpct", po::value< Numeric >()->default_value( 0.2 ),
					"percentage of total #vars to assign at each level" )
			( "nRRperLvl", po::value< uint32 >()->default_value( 2 ),
					"number of random restarts to do per level of optimization" )
			( "nRRatTop", po::value< uint32 >(),
					"number of random restarts to do at the top level of "
							"optimization (overrides nRRperLvl if set)" )
			( "minRR", po::value< uint32 >()->default_value( 1 ),
					"min number of random restarts at each level of optimization" )
			( "maxNAtoRR", po::value< uint32 >()->default_value( 10 ),
					"max number of assignments (at each level) before forcing a "
							"random restart")
			( "SSmaxit", po::value< VariableCount >()->default_value( 25 ),
					"maximum number of iterations per subspace optimizer call" )
			( "SSftol", po::value< Numeric >(),
					"fractional tolerance to determine convergence of the "
							"subspace optimizer" )
			( "alwaysUseAllFactors", po::value< bool >()->default_value( false ),
					"always use all possible factors for g.d., even if it means "
							"setting some variables to arbitrary values" )
			( "noAssignLimitAtTop", po::value< bool >()->default_value( true ),
					"ignore the assignment limit for top level components (do "
							"each RR until convergence)" )
			( "printAllFullEvals", po::value< bool >()->default_value( false ),
					"true to enable printing of the current full function eval "
							"every time a component is optimized" )
			;
}


Numeric RDISOptimizer::optimize( const State * xinitial, Numeric finit,
		Clock::time_point timeout, bool printInfo,
		Clock::time_point actualStartTime ) {

	m_starttime = ( actualStartTime < Clock::now() ) ?
			actualStartTime : Clock::now();
	m_timeout = timeout;

	Numeric fevalSimpl = 0.0;
	Numeric ferrSimpl = 0.0;

	try {
		prepareToOptimize();

		if ( approxFactors ) fevalSimpl = simplifyFactors( ferrSimpl );

		// count the number of simplified factors
		nSimplifiedFactors = 0;

		// count empty factors (if there are any) towards the final output
		// (note that these won't get processed)
		for ( Factor * f : m_func.getFactors() ) {
			if ( f->areAllVarsAssigned() ) {
				assert( !f->isAssigned() );
				fevalSimpl += f->eval();
				++nSimplifiedFactors;
			} else if ( f->isAssigned() ) {
				++nSimplifiedFactors;
			}
		}

		if ( cgNeedsInit ) {
			connGraph.init( m_func, m_func.getVariables(), m_func.getFactors() );
		}

		if ( printInfo ) writeGraph();

		printedAO = !printInfo;
		if ( printInfo ) {
			std::cout <<
					"RDIS: bounds are " << ( useBounds ? "ON" : "OFF" ) <<
					", Cache is " << ( useCache ? "ON" : "OFF" ) << std::endl;
		}

		if ( xinitial != NULL ) setInitialState( *xinitial );
		topccomp->setFMin( finit - fevalSimpl );

		// decompose initial component if possible
		topccomp->decompose();
		SubdomainSP sd = sdpool.get();
		sd->set( fevalSimpl, ferrSimpl, false );
		topccomp->setCurSD( sd );

		if ( printInfo ) writeFactorGraph( *topccomp );

		std::cout << "initial feval from factor simplification (" <<
				nSimplifiedFactors  << " factors): " << fevalSimpl <<
				" (error: " << ferrSimpl << ")" << std::endl;

		for ( const ComponentIP & c : topccomp->getChildren() ) {
			c->computeFMin();

			// optimize this component
			doOptimization( *c );

			topccomp->onChildEvaluated( *c );
			productComponentVals( *topccomp, *c );

			assert( c->hasBeenEvald() || super::checkTimedOut() );
			if ( printInfo && c->getOptSD()->isBound ) {
				std::cout << "Component eval (" << c->getOptSD()->fx <<
						") is bound: " << *c << std::endl;
			}

			if ( c->haveVarsBeenAssigned() && c->getOptSD() ) {
				const SubdomainCSP & cmin( c->getOptSD() );
				copySubdomain( *cmin, m_optState, true, true );
				if ( m_func.hasBounds() ) {
					std::cout << *cmin << " has fmin: " << c->getFMin() <<
							", bounds: " << c->getUnassignedBounds() << std::endl;
				}
			}

			if ( !topccomp->getCurSD()->isBound ) {
				c->quickAssignOptVals( false );
			}
		}

		topccomp->quickUnassignChildren();

		m_wasBound = topccomp->getCurSD()->isBound;
		m_optVal = topccomp->getCurSD()->fx;
		m_optError = topccomp->getCurSD()->ferr;

		std::cout << "num assignments: "; printVector( numassign ); std::cout << std::endl;
//		std::cout << "num cache hits: "; printVector( numcachehits ); std::cout << std::endl;
//		std::cout << "num bound hits: "; printVector( numboundhits ); std::cout << std::endl;

	} catch ( std::exception & e ) {
		std::cerr << "Exception occurred during optimization: \n\t"
				<< e.what() << std::endl;
	}

	unsimplifyFactors();

	// clean up opt state -- set unknown var vals to 0
	for ( Numeric & x : m_optState ) {
		if ( x == NumericMAX ) x = 0;
	}

	m_optTime = ( Clock::now() - m_starttime );
	instrumentation::g_runTime = m_optTime;

	if ( !printedAO ) {
		std::cout << "partial assignment ordering: " << assignmentOrder << std::endl;
	}

	if ( printInfo ) {
		if ( m_stopped ) std::cerr << "*** Optimization was stopped ***" << std::endl;
		else if ( m_timedOut ) std::cerr << "*** Optimization timed out ***" << std::endl;

		if ( m_wasBound ) std::cerr << "*** Optimization result was a bound ***" << std::endl;
	}

	return m_optVal;
}


void RDISOptimizer::doOptimization( Component & ccomp ) {

	assert( !ccomp.haveVarsBeenAssigned() );

//	std::cout << "optimizing component " << ccomp.toString() <<
////			", fmin: " << ccomp.getFMin() <<
//		" at " << (Duration) ( Clock::now() - m_starttime ) << std::endl;

	// make sure there are factors in here to optimize
	if ( checkEmpty( ccomp ) ) return;

	// check if B&B or the cache allow us to bypass this computation
	if ( getFromCache( ccomp ) ) return;
	if ( checkUnassignedBound( ccomp ) ) return;

	if ( checkTimedOut( ccomp ) ) return;

	VarDataPtrVec & vds = ccomp.getVarsToAssign();
	vds.clear();

	// (heuristically) choose the variables to assign from this component
	getVarsToAssign( ccomp, vds );

	NumericVec & xvals = ccomp.getValToAssign();
	xvals.clear();

	while ( getValueFromDomain( ccomp, vds, xvals ) ) {

		// assign the variable, simplify factors, and update the conn. graph
		assign( ccomp, vds, xvals );

		// check "assigned" bound for this xval of this component and only
		// recurse if the bound wasn't hit
		if ( !checkAssignedBound( ccomp ) ) {

			// decompose into subcomponents based on the connectivity graph
			ccomp.decompose();

			size_t icc = 0, ncc = ccomp.getChildren().size();
			for ( const ComponentIP & c : ccomp.getChildren() ) {

				// compute and set fmin for this child
				c->computeFMin();

				if ( ccomp.getCurSD()->isBound ) {
					// if one child was already a bound, don't do wasted work
					c->setOptSD( sdpool.get(),
							c->getUnassignedBounds().upper(), true );
					c->setHasBeenEvald( true );
				} else {
					// otherwise, recurse
					doOptimization( *c );
				}

				// update the component based on the child's evaluation
				ccomp.onChildEvaluated( *c );

				// product the computed subcomponent values together
				productComponentVals( ccomp, *c );

				if ( ++icc < ncc && checkAssignedBound( ccomp ) ) break;
			}
		}

		updateDomain( ccomp );

		if ( checkTimedOut( ccomp ) ) {
			ccomp.setHasBeenEvald( false );
			break;
		}
		if ( ccomp.hasBeenEvald() ) break;
	}

	// restore simplified factors, unassign the var, and update the conn. graph
	unassign( ccomp, vds );

	// add this to the cache
	if ( useCache ) {
		addToCache( ccomp );
		if ( useLinearCache ) removeChildrenFromCache( ccomp );
	}
}


void RDISOptimizer::getVarsToAssign( Component & ccomp, VarDataPtrVec & vars ) {
	RDISOptionValues & rov = optionVals;

	if ( rov.useStaticDecomp && getVarsToAssignStatic( ccomp, vars ) ) return;

	if ( rov.heur == All || ( ccomp.getVars().size() <= rov.desiredNumVars ) ) {
		for ( VariableID vid : ccomp.getVars() ) {
			vars.push_back( &vardata.get( vid ) );
		}
		return;
	}

	std::vector< VariableCount > vidx;
	if ( rov.heur == Random ) {
		vidx.resize( ccomp.getVars().size() );
		for ( VariableCount i = 0; i < vidx.size(); ++i ) vidx[i] = i;
		permute( vidx, rng );
	}

	assert( vars.empty() );

	size_t nToAssign = 0;
	NumericVec val;

	for ( size_t counter = 0; counter <= rov.desiredNumVars; ++counter ) {
		// add the rest of the chosen variables' blocks to the list
		addBlockRemainders( ccomp, vars, nToAssign );

		// check if done
		if ( vars.size() >= rov.desiredNumVars ) break;

		// if not done, assign the chosen vars and get additional vars

		val.resize( vars.size(), 0 );
		for ( size_t i = nToAssign; i < vars.size(); ++i ) {
			val[i] = vars[i]->v->getDomain().min();
			assert( !vars[i]->v->isAssigned() );
		}
		quickAssignVals( vars, val );
		nToAssign = vars.size();

		size_t initsize = vars.size();
		if ( ( rov.heur == HMetis || rov.heur == PaToH ) &&
				( ccomp.getVars().size() - initsize ) > 4 ) {
			if ( rov.heur == HMetis ) getVarsToAssignHMetis( ccomp, vars );
			else if ( rov.heur == PaToH ) getVarsToAssignPaToH( ccomp, vars );

		} else if ( rov.heur == Random ) {
			uint32 numToGet = (uint32) ( vars.empty() ?
										 /*0.75 * */ rov.desiredNumVars :
										 ( rov.desiredNumVars - vars.size() ) );
			for ( uint32 i = 0; i < numToGet && !vidx.empty(); ++i ) {
				vars.push_back( &vardata.get( ccomp.getVars()[ vidx.back() ] ) );
				vidx.pop_back();
			}
		} else {
			getVarsToAssignMOMS( ccomp, vars, false );
		}

		// if made no progress, quit
		if ( vars.size() == initsize ) break;
	}

	// make sure to unassign any vars we assigned
	quickUnassignVals( vars );

	// make sure that at least one factor will be fully assigned
	ensureFactorWillBeAssigned( ccomp, vars );

	// if we still haven't found any vars, just assign them all (but this
	// shouldn't really happen)
	if ( vars.empty() ) {
		for ( VariableID vid : ccomp.getVars() ) {
			vars.push_back( &vardata.get( vid ) );
		}
	}

	// make sure the variable list is unique
	std::sort( vars.begin(), vars.end() );
	auto itunq = std::unique( vars.begin(), vars.end() );
	assert( itunq == vars.end() ); // should be unique already?
	vars.erase( itunq, vars.end() );
}


void RDISOptimizer::ensureFactorWillBeAssigned( Component & ccomp,
        VarDataPtrVec & vars ) {
	VariableIDVec & vids = vidstmp;
	vids.clear();
	vids.reserve( vars.size() );
	for ( const VarData * vd : vars ) vids.push_back( vd->v->getID() );

	// check if any factors will be assigned by this set of vars
	bool anyf = false;
	const Factor * fnlamin = NULL;
	VariableCount nlamin = this->vardata.numVars() + 1;
	for ( const Factor * f : ccomp.getFactors() ) {
		if ( f->isAssigned() && f->getAssignedKey() != vids[0] ) continue;
		if ( f->areAnyInFactor( vids ) ) {
			VariableCount nla = f->numLeftUnassigned( vids );
			if ( nla <= 0 ) {
				anyf = true;
//				break;
			} else if ( nla < nlamin ) {
				fnlamin = f;
				nlamin = nla;
			}
		}
	}

	// if none will be assigned, take the factor with the least number of vars
	// left to be assigned and assign those vars
	if ( fnlamin != NULL && ( !anyf || vars.size() <= 1 ) ) {
		size_t initpos = vars.size();
		for ( const Variable * v : fnlamin->getVariables() ) {
			if ( !this->vardata.isAssigned( v->getID() ) && std::find( vids.begin(),
					vids.end(), v->getID() ) == vids.end() &&
					ccomp.containsVar( v->getID() ) ) {
				vars.push_back( &this->vardata.get( v->getID() ) );
			}
		}

		addBlockRemainders( ccomp, vars, initpos );
	}
}


void RDISOptimizer::addBlockRemainders( Component & ccomp,
		VarDataPtrVec & vars, const size_t initpos ) {

	size_t initsize = vars.size();
	auto vitbeg = vars.begin() + initpos;
	auto vitend = vars.begin() + initsize;

	for ( VariableCount i = initpos; i < (VariableCount) initsize; ++i ) {
		VariableID vid = vars[i]->v->getID();
		VariableID vidlower = -1, vidupper = -1;
		m_func.getBlockRangeByVid( vid, vidlower, vidupper );

		for ( VariableID vidblk = vidlower; vidblk <= vidupper; ++vidblk ) {
			if ( vidblk == vid ) continue;
			VarData & vd = vardata.get( vidblk );

			// make sure this var isn't in the list already and is in the comp.
			if ( ccomp.containsVar( vidblk ) &&
				 std::find( vitbeg, vitend, &vd ) == vitend ) {
				vars.push_back( &vd );
			}
		}
	}
}


bool RDISOptimizer::getVarsToAssignStatic( Component & ccomp,
		VarDataPtrVec & vds ) {

	assert( !ccomp.haveVarsBeenAssigned() );
	bool gotVars = false;

	if ( !ccomp.getAssignedVids().empty() ) {
		// just reassign existing vars if they're already set
		for ( VariableID vid : ccomp.getAssignedVids() ) {
			vds.push_back( &vardata.get( vid ) );
		}

		gotVars = true;
	}

	return gotVars;
}


void RDISOptimizer::getVarsToAssignMOMS( Component & ccomp, VarDataPtrVec & vds,
		const bool weighted ) {

	// implementation of the MOMS heuristic

	VarData * bestvd( NULL );
	VariableCount minfsize = std::numeric_limits< VariableCount >::max();
	FactorPtrVec fsmall;
	fsmall.reserve( ccomp.getFactors().size() / 2 + 1 );

	NumericVec fweight;
	if ( weighted ) fweight.reserve( ccomp.getFactors().size() );

	for ( Factor * f : ccomp.getFactors() ) {
		const VariableCount fs = f->numUnassignedVars();

		if ( fs < minfsize ) {
			minfsize = fs;
			fsmall.clear();
			fsmall.push_back( f );
		} else if ( fs == minfsize ) {
			fsmall.push_back( f );
		}

		if ( weighted ) fweight.push_back( std::exp( (Numeric) -fs ) );
	}

	assert( !fsmall.empty() );

	VariableIDVec bestvids;
	bestvids.reserve( 5 );
	VariableCount bestval = 0;

	NumericVec vweights;
	vweights.reserve( ccomp.getVars().size() );

	for ( VariableID vid : ccomp.getVars() ) {
		if ( !weighted ) {
			// important to use v->isAssigned() instead of vardata.isAssigned() here
			if ( vardata.var( vid )->isAssigned() ) continue;

			VariableCount count = 0;
			for ( Factor * f : fsmall ) {
				if ( f->containsVar( vid ) ) ++count;
			}

			if ( count > bestval ) {
				bestvids.clear();
				bestvids.push_back( vid );
				bestval = count;
			} else if ( count == bestval ) {
				bestvids.push_back( vid );
			}
		} else {
			vweights.push_back( 0.0 );

			// important to use v->isAssigned() instead of vardata.isAssigned() here
			if ( vardata.var( vid )->isAssigned() ) continue;

			for ( size_t fidx = 0; fidx < ccomp.getFactors().size(); ++fidx ) {
				Factor * f = ccomp.getFactors()[fidx];
				if ( f->containsVar( vid ) ) vweights.back() += fweight[fidx];
			}
		}
	}

	if ( weighted ) {
		auto vit = std::max_element( vweights.begin(), vweights.end() );
		VariableCount vidx = ( vit - vweights.begin() );
		bestvd = &vardata.get( ccomp.getVars()[vidx] );
	} else if ( !bestvids.empty() ) {
		bestvd = &vardata.get( bestvids.front() );
	}

	if ( bestvd != NULL ) vds.push_back( bestvd );

	return;
}

void RDISOptimizer::getVarsToAssignHMetis( Component & ccomp,
		VarDataPtrVec & vds ) {
#ifdef USE_HMETIS
	const bool printdbg = true;
	const Clock::time_point start = Clock::now();

	typedef BOOSTNS::container::vector< int > IntVec;
	typedef BOOSTNS::container::flat_map< FactorID, int > FactorIndices;

	const FactorPtrVec & fctrs = ccomp.getFactors();
	const VariableIDVec & vars = ccomp.getVars();

	FactorIndices fctidx;
	fctidx.reserve( fctrs.size() );

	if ( printdbg ) std::cout << "calling HMetis on component " << ccomp << std::endl;

	const int nvtxs = fctrs.size();
	int nhedges = 0;

	// build the index for the factors in this component
	for ( int i = 0; i < (int) fctrs.size(); ++i ) {
		fctidx.emplace( fctrs[i]->getID(), i );
	}

	IntVec eptr;
	IntVec eind;
	eptr.reserve( vars.size() + 1 );
	eind.reserve( vars.size() * 2 );

	// create the hyperedges
	for ( VariableCount i = 0; i < (int) vars.size(); ++i ) {
		const VariableID vid = vars[i];
		const Variable * v = vardata.var( vid );

		if ( v->isAssigned() ) continue;

		int hesize = 0;
		eptr.push_back( eind.size() );

		for ( Factor * f : v->getFactors() ) {
			if ( f->containsVar( vid ) && !f->areAllVarsAssigned() ) {
				assert( fctidx.find( f->getID() ) != fctidx.end() );
				eind.push_back( fctidx[f->getID()] );
				++hesize;
			}
		}

		++nhedges;
	}

	if ( printdbg ) {
		std::cout << "HMetis with " << nhedges << " hyperedges and " << nvtxs <<
				" vertices from " << vars.size() << " vars and " <<
				fctrs.size() << " factors" << std::endl;
	}

	assert( nhedges == (int) eptr.size() );
	eptr.push_back( eind.size() );

	eptr.shrink_to_fit();
	eind.shrink_to_fit();

	const int nparts = 2;
	const int ubfactor = 40;

	IntVec options( 9, 0 );

	IntVec part( nvtxs, 0 );
	int edgecut = -1;

	HMETIS_PartRecursive( nvtxs, nhedges, NULL, eptr.data(), eind.data(), NULL,
			nparts, ubfactor, options.data(), part.data(), &edgecut );

	VariableIDVec sep;
	sep.reserve( edgecut );

	for ( VariableCount i = 0; i < (int) vars.size(); ++i ) {
		VariableID vid = vars[i];
		Variable * v = vardata.var( vid );
		if ( v->isAssigned() ) continue;
		int vpart = -1;
		for ( Factor * f : v->getFactors() ) {
			if ( !f->containsVar( vid ) || f->areAllVarsAssigned() ) continue;
			int fpart = part[fctidx[f->getID()]];

			if ( vpart < 0 ) {
				vpart = fpart;
			} else if ( vpart != fpart ) {
				sep.push_back( vid );
				break;
			}
		}

		// found all the separator vars
		if ( (int) sep.size() >= edgecut ) break;
	}

	if ( printdbg ) {
		std::cout << "vertex separator (hyperedges to cut) is: ";
		printVector( sep );
		std::cout << " in " << (Duration) (Clock::now() - start) << std::endl;
	}

	for ( VariableID vid : sep ) {
		vds.push_back( &vardata.get( vid ) );
	}
#endif // USE_HMETIS
}


void RDISOptimizer::getVarsToAssignPaToH( Component & ccomp, VarDataPtrVec & vds ) {

#ifdef USE_PATOH
	const bool printdbg = true;
	const Clock::time_point start = Clock::now();

	typedef BOOSTNS::container::vector<int> IntVec;
	typedef BOOSTNS::container::flat_map<FactorID, int> FactorIndices;

	const FactorPtrVec & fctrs = ccomp.getFactors();
	const VariableIDVec & vars = ccomp.getVars();

	FactorIndices fctidx;
	fctidx.reserve( fctrs.size());

	if ( printdbg )
		std::cout << "calling PaToH on component " << ccomp << std::endl;

	const int nvtxs = fctrs.size();
	int nhedges = 0;

	// build the index for the factors in this component
	for ( int i = 0; i < (int) fctrs.size(); ++i ) {
		fctidx.emplace( fctrs[i]->getID(), i );
	}

	// weights of vertices (factors) and hyperedges (variables)
	IntVec vwgts( nvtxs, 0 ), hewgts( nhedges, 1 );

	// eptr contains the index into eind that each hyperedge's set of factors
	// starts at, eind contains a list of factors for each hyperedge (in one
	// big list)
	IntVec eptr, eind;
	eptr.reserve( vars.size() + 1 );
	eind.reserve( vars.size() * 2 );

	// create the hyperedges
	for ( VariableCount i = 0; i < (int) vars.size(); ++i ) {
		const VariableID vid = vars[i];
		const Variable * v = vardata.var( vid );

		if ( v->isAssigned()) continue;

		eptr.push_back( eind.size());

		int hesize = 0;
		for ( Factor * f : v->getFactors()) {
			if ( f->containsVar( vid ) && !f->areAllVarsAssigned()) {
				assert( fctidx.find( f->getID()) != fctidx.end());
				eind.push_back( fctidx[f->getID()] );
				++hesize;
				++vwgts[fctidx[f->getID()]];
			}
		}

		++nhedges;
	}

	if ( printdbg )
		std::cout << "PaToH with " << nhedges << " hyperedges and " << nvtxs <<
				" vertices from " << vars.size() << " vars and " <<
				fctrs.size() <<
				" factors" << std::endl;

	assert( nhedges == (int) eptr.size());
	eptr.push_back( eind.size());    // make sure to put the last index in

	eptr.shrink_to_fit();
	eind.shrink_to_fit();

	// number of partitions
	const int nparts = 2;
	BOOSTNS::container::vector<float> targetWgts( nparts, 0.4 );
//	const int ubfactor = 30;

//	IntVec options( 9, 0 );

//	IntVec part( nvtxs, 0 );
//	int edgecut = -1;

	PaToH_Parameters args;
	PaToH_Initialize_Parameters( &args, PATOH_CONPART, PATOH_SUGPARAM_QUALITY );

	args._k = nparts;

	if ( optionVals.timeAsSeed ) {
		args.seed = std::time(NULL) + getpid();
	} else {
		args.seed = 823527793; // chosen randomly :)
	}

	int _c, _n, _nconst;
	int *cwghts, *nwghts, *xpins, *pins, cut;
	bool useFixCells = false;

	_c = nvtxs;				// number of cells (vertexes aka factors)
	_n = nhedges;			// number of nets (hyper edges aka variables)
	_nconst = 1; 			// number of constraints

	cwghts = vwgts.data();	// weight of each cell (aka vertexes == factors)
	nwghts = hewgts.data();	// weight of each net  (aka hyperedges == variables)
	xpins = eptr.data();	// index of first pins for each hyperedge
	pins = eind.data();		// IDs of pins in each hyperedge (factors

	PaToH_Alloc( &args, _c, _n, _nconst, cwghts, nwghts, xpins, pins );

	IntVec part( _c, 0 );
	IntVec partwgts( args._k * _nconst, 0 );

	PaToH_Part( &args, _c, _n, _nconst, useFixCells, cwghts, nwghts,
			xpins, pins, targetWgts.data(), part.data(), partwgts.data(), &cut );

	// VTX weights could be number of attached vars??

	if ( printdbg ) {
		std::cout << "called PaToH, resulting partition (" << cut << "): ";
//		for ( int i = 0; i < nvtxs; ++i ) {
//			std::cout << ( i == 0 ? "" : ", " ) << part[i];
//		}
		for ( int pw : partwgts ) {
			std::cout << pw << ", ";
		}
		std::cout << std::endl;

//		for ( int p = 0; p < nparts; ++p ) {
//			std::cout << "partition " << p << ": ";
//			for ( int i = 0; i < nvtxs; ++i ) {
//				if ( part[i] == p ) {
//					std::cout << fctrs[i]->getID() << ", ";
//				}
//			}
//			std::cout << std::endl;
//		}
	}

	VariableIDVec sep;
	sep.reserve( cut );

	for ( VariableCount i = 0; i < (int) vars.size(); ++i ) {
		VariableID vid = vars[i];
		Variable * v = vardata.var( vid );
		if ( v->isAssigned() ) continue;
		int vpart = -1;
		for ( Factor * f : v->getFactors() ) {
			if ( !f->containsVar( vid ) || f->areAllVarsAssigned() ) continue;
			int fpart = part[fctidx[f->getID()]];

			if ( vpart < 0 ) {
				vpart = fpart;
			} else if ( vpart != fpart ) {
				sep.push_back( vid );
				break;
			}
		}

		// found all the separator vars
		if ( (int) sep.size() >= cut ) break;
	}

	if ( printdbg ) {
		std::cout << "vertex separator (hyperedges to cut) is: ";
		printVector( sep );
		std::cout << " in " << (Duration) (Clock::now() - start) << std::endl;
	}

	// clean up
	PaToH_Free();

	for ( VariableID vid : sep ) {
		vds.push_back( &vardata.get( vid ) );
	}

#endif // USE_PATOH
}


bool RDISOptimizer::checkEmpty( Component & ccomp ) {
	if ( ccomp.getFactors().empty() ) {
		SubdomainSP newsd = sdpool.get();
		newsd->set( 0.0, 0.0, false );
		for ( VariableID vid : ccomp.getVars() ) {
			const VariableDomain & dom = vardata.var( vid )->getDomain();
			newsd->vids.push_back( vid );
			newsd->x.push_back( dom.closestVal( dom.median() ) );
		}

		ccomp.setCachedValue( newsd );

		return true;
	}

	return false;
}


bool RDISOptimizer::checkUnassignedBound( Component & ccomp ) {
	if ( !useBounds ) return false;
	bool result = false;

	// if this component's min value is larger than the current min, then don't
	// need to optimize this sub-function because it can't result in the optimum

	const Numeric lb( ccomp.getUnassignedBounds().lower() );
	const Numeric fmin( ccomp.getFMin() );

	if ( fmin < lb ) {
		std::cout << "UaB hit: fmin " << fmin << ", lb " << lb << std::endl;

		ccomp.setOptSD( sdpool.get(), lb, true );
		ccomp.setHasBeenEvald( true );
		result = true;
	}

	return result;
}


bool RDISOptimizer::checkAssignedBound( Component & ccomp ) {
	if ( !useBounds ) return false;
	bool result = false;

	const Numeric lb( ccomp.getAssignedBounds().lower() );
	const Numeric fmin( ccomp.getFMin() );

	if ( ccomp.getFMin() <= lb ) {
//		std::cout << "AB hit: fmin " << fmin << ", b " << lb << std::endl;

		ccomp.updateCurSD( lb, ccomp.getCurSD()->ferr, true );
		result = true;

		for ( VarData * vd : ccomp.getAssignedVDs() ) {
			++numboundhits[vd->v->getID()];
		}
	}

	return result;
}


bool RDISOptimizer::getFromCache( Component & ccomp ) {

	if ( !useCache || ccomp.getParent() == NULL ||
			ccomp.getCacheKey().factors.empty() ) {
		return false;
	}

	SubdomainCSP sdopt;
	bool result( cache.find( ccomp.getCacheKey(), sdopt ) );

	++instrumentation::g_cacheLookups;

	if ( result ) {
		// if the cached value is a bound then the fopt for which it was a bound
		// might not be the relevant fopt now (e.g., if the sibling components
		// of this ccomp have changed) so check if this bound is still valid
		if ( sdopt->isBound && ( sdopt->fx != ccomp.getFMin() ) &&
				( sdopt->fx <= ccomp.getFMin() ) ) {
			result = false;
		}

		if ( result ) result = ccomp.setCachedValue( sdopt );

		if ( result ) {
			for ( VariableID vid : sdopt->vids ) ++numcachehits[vid];
			++instrumentation::g_cacheHits;
		}
	}

	return result;
}


bool RDISOptimizer::getValueFromDomain( Component & ccomp,
		const VarDataPtrVec & vars, NumericVec & value ) {

	const bool printdbg = false;
	bool doAlternatingMin = true;
	RDISOptionValues & rov = optionVals;

	const bool isTopCComp = ( ccomp.getParent() == NULL ||
							  !ccomp.getParent()->haveVarsBeenAssigned() );
	const uint32 nRRthisLvl =
			( rov.nRRperLvl / ( ( (uint32) 1 ) << ( ccomp.getLevel() - 1 ) ) );
	const uint32 numRestarts = std::max( rov.minNumRestarts,
			isTopCComp ? rov.nRRatTop : nRRthisLvl );

	assert( !vars.empty() );
	if ( vars.empty() ) {
		std::cout << "vars vector empty in getValueFromDomain()" << std::endl;
		return false;
	}

	const VariableID repvid = vars.front()->v->getID();
	const size_t nrr = ccomp.getNumRandomRestarts();
	const bool forceRR = !( isTopCComp && rov.noAssignLimitAtTop ) &&
			( ccomp.getNumVAsinceLastRR() >= rov.maxNumAssignsUntilRR );

	// check if we're done with this component
	if ( ( ( forceRR || !ccomp.wasLastEvalOpt() ) && nrr > numRestarts ) ) {
		return false;
	}

	if ( printdbg ) {
		std::cout << "vars to assign: ";
		for ( const VarData * vd : vars ) std::cout << vd->v->getID() << ",";
		std::cout << std::endl;
		if ( m_func.hasBlockedVars() ) {
			std::cout << "-- blocks: ";
			for ( const VarData * vd : vars ) {
				std::cout << m_func.getBlockID( vd->v->getID() ) << ",";
			}
			std::cout << std::endl;
		}
	}

	NumericVec startvals;
	if ( vars.front()->v->isAssigned() ) {
		for ( VarData * vd : vars ) startvals.push_back( vd->v->eval() );
	}

	const bool wasAssigned = vardata.isAssigned( repvid );

	bool redoGD = false, success = false;
	FactorPtrVec & gdfs = factmp;
	VariablePtrVec & varptrs = varstmp;
	varptrs.reserve( vars.size() );

	SubdomainCSP sdprev;
	if ( ccomp.haveVarsBeenAssigned() ) sdprev = ccomp.getPrevSD();

	const bool useAllFactors = ccomp.hasNonmonotoneEvals() ||
			rov.alwaysUseAllFactors;

	do {
		redoGD = false;
		value.clear(); gdfs.clear(); varptrs.clear();

		// set a starting point, call the subspace opt. from there, and return
		// the point it gives

		ValSelType vst = getSSInitialVal( ccomp, vars, value, doAlternatingMin,
				forceRR, wasAssigned );

		quickAssignSSInitialVal( ccomp, vars, value, doAlternatingMin,
				sdprev, useAllFactors );

		if ( vst == Init || vst == Restart ) ccomp.incrementNumRandomRestarts();
		if ( vst == Restart && printdbg ) std::cout << "doing random restart" << std::endl;

		// create the list of factors for the subspace optimizer
		for ( Factor * f : ccomp.getFactors() ) {
			if ( f->areAllVarsAssigned() ||
				 ( f->isAssigned() && f->getAssignedKey() == repvid ) ) {
				for ( const VarData * vd : vars ) {
					if ( f->containsVar( vd->v->getID() ) ) {
						gdfs.push_back( f );
						break;
					}
				}
			}
		}

		// create the list of vars for the subspace optimizer
		for ( VarData * vd : vars ) varptrs.push_back( vd->v );

		Numeric deltafval = 0.0;

		// perform subspace optimization
		Numeric result = ssopt.optimize( varptrs, gdfs, value, deltafval,
										   printdbg || isTopCComp );

		if ( printdbg || isTopCComp ) {
//			std::cout << "SS opt. result (steps " << nsteps << "): " << result
//				<< ", diff: " << deltafval << " (" <<
			std::cout << "\t" << "step type = (" <<
					( vst == Init ? "initial values" :
					  ( vst == Restart ? "random restart" :
						"iterative improvement" )) << ")"
					<< ( isTopCComp ? " (is top)" : "" ) << ", nRR "
					<< ccomp.getNumRandomRestarts() << " / " << numRestarts;
			if ( printdbg ) std::cout << " -- vals: " << value;
			std::cout << std::endl;
		}

		quickUnassignSSInitialVal( ccomp, vars, value, doAlternatingMin,
				sdprev, useAllFactors );

		Numeric ftol = epsilon / 10.0;
		if ( options.count( "steptol" ) ) ftol = options[ "steptol" ].as< Numeric >();

		if ( sdprev && approxgeq( deltafval, 0.0, ftol ) ) {
			if ( ccomp.getNumRandomRestarts() < numRestarts ) {
				// try again from a random position
				doAlternatingMin = false;
				redoGD = true;
				if ( printdbg ) {
					std::cout << "subspace opt. made no progress -- doing RR ("
							<< deltafval << ", " << ftol << ")" << std::endl;
				}
			} else {
				// consider this a failure
				success = false;
				redoGD = false;
			}
		} else {
			success = true;
		}

	} while ( redoGD && !success );

	// return to initial assignment state before this func was called
	if ( !startvals.empty() ) quickAssignVals( vars, startvals );

	if ( success && doAlternatingMin && sdprev ) {
		setInitialValFromChildren( *sdprev );
	}

	return success;
}


RDISOptimizer::ValSelType RDISOptimizer::getSSInitialVal( Component & ccomp,
        const VarDataPtrVec & vars, NumericVec & value,
		const bool doAlternatingMin, const bool forceRR,
		const bool wasAssigned ) {

	RDISOptimizer::ValSelType vst = Init;

	// randomly choose values for each var from their domains
	for ( VarData * vd : vars ) {
		if ( doAlternatingMin && !wasAssigned && vd->xvinitIsSet ) {
			value.push_back( vd->xvalinit );
			vst = Init;

			// set to false so initial vals don't propagate across restarts
			vd->xvinitIsSet = false;
		} else if ( doAlternatingMin && wasAssigned && !forceRR ) {
			value.push_back( vd->v->eval() );
			vst = AlternatingMin;
		} else {
			value.push_back( sampleRandomState( *vd ) );
			vst = Restart;
		}
	}

	return vst;
}


void RDISOptimizer::quickAssignSSInitialVal( Component & ccomp,
		const VarDataPtrVec & vars, const NumericVec & value,
		const bool doAlternatingMin, const SubdomainCSP & sdprev,
		const bool useAllFactors ) {

	// assign vals to all the other vars so that their factors have evals (even
	// if those evals are not very good)
	if ( doAlternatingMin && sdprev ) {
		for ( const auto & c : sdprev->children ) {
			ccomp.quickAssignFromSD( *c );
		}
	} else if ( useAllFactors ) {
		// have no info for these vars, so set them to their initial val
		// or a random val
		for ( VariableID vid : ccomp.getVars() ) {
			VarData & vd = vardata.get( vid );
			Numeric xval = 0;
			if ( vd.xvinitIsSet ) {
				xval = vd.xvalinit;
			} else {
				xval = sampleRandomState( vd );
				vd.xvalinit = xval;
				vd.xvinitIsSet = true;
			}

			vd.v->assign( xval );
			m_func.onVarAssigned( vid, xval );
		}
	}

	// (quick) assignment of vars to their values (don't go through all the
	// other assignment machinery)
	quickAssignVals( vars, value );
}


void RDISOptimizer::quickUnassignSSInitialVal( Component & ccomp,
		const VarDataPtrVec & vars, const NumericVec & value,
		const bool doAlternatingMin, const SubdomainCSP & sdprev,
		const bool useAllFactors ) {
	if ( doAlternatingMin && sdprev ) {
		//		std::cout << "unassigning from SD" << std::endl;
		for ( const SubdomainCSP & c : sdprev->children ) {
			ccomp.quickUnassignFromSD( *c );
		}
	} else if ( useAllFactors ) {
		for ( VariableID vid : ccomp.getVars() ) {
			vardata.var( vid )->unassign();
			m_func.onVarUnassigned( vid );
		}
		quickAssignVals( vars, value );
	}

	// (quick) unassignment of vars
	if ( !ccomp.haveVarsBeenAssigned() ) quickUnassignVals( vars );
}


Numeric RDISOptimizer::sampleRandomState( const VarData & vd ) {

	const VariableDomain & dom = vd.v->getDomain();

	NumericInterval sampleItvl = dom.getSamplingInterval();
	NumericInterval fullItvl = dom.interval();

	// if the sampling interval has been specified, then use it directly
	if ( sampleItvl.lower() != fullItvl.lower() ||
		 sampleItvl.upper() != fullItvl.upper() ) {
		const Numeric w = BOOSTNS::numeric::width( sampleItvl );
		if ( approxeq( w, 0, 1e-12 ) ) {
			std::cout << "sampling interval of variable " << vd.v->getID() <<
			" has width 0 " << sampleItvl << " (domain: " << dom << ")"
			<< std::endl;
			return BOOSTNS::numeric::median( sampleItvl );
		} else {
			BOOSTNS::random::uniform_real_distribution<> unif(
					sampleItvl.lower(), sampleItvl.upper() );
			return dom.closestVal( unif( rng ) );
		}
	}

	const Numeric w = dom.realWidth();
	if ( approxeq( w, 0, 1e-12 ) ) return dom.median();

	BOOSTNS::random::uniform_real_distribution<> unif( 0, w );
	Numeric wval = unif( rng );

	const auto subitvls = dom.subintervals();
	Numeric realval = dom.subintervals()[0].lower();
	long long int sii = 0;

//	std::cout << "var " << vd.v->getID() << ", width: " << w <<
//			", wval: " << wval << " -- ";
//	printContainer( subitvls ); std::cout << std::endl;

	while ( wval > 0.0 && sii < (long long int) subitvls.size() ) {
		const NumericInterval & si = subitvls[sii];

		if ( realval + wval > si.upper() ) {
			wval -= ( si.upper() - realval );
			realval = ( sii + 1 < (int) subitvls.size() ?
						subitvls[sii+1].lower() : si.upper() );
		} else {
			realval += wval;
			wval = 0.0;
		}

		++sii;
	}

	assert( approxeq( wval, 0.0, 1e-8 ) );

	return realval;
}


void RDISOptimizer::quickAssignVals( const VarDataPtrVec & vars,
									 const NumericVec & xval ) {
	assert( vars.size() == xval.size() );
	for ( size_t i = 0; i < vars.size(); ++i ) {
		vars[i]->v->assign( xval[i] );
		m_func.onVarAssigned( vars[i]->v->getID(), xval[i] );
	}
}


void RDISOptimizer::quickUnassignVals( const VarDataPtrVec & vars ) {
	for ( VarData * vd : vars ) {
		if ( vd->v->isAssigned()) {
			vd->v->unassign();
			m_func.onVarUnassigned( vd->v->getID());
		}
	}
}


void RDISOptimizer::assign( Component & ccomp, const VarDataPtrVec & vds,
		const NumericVec & xval ) {

	const bool wasUnassigned( !vardata.isAssigned( vds.front()->v->getID() ) );

	if ( wasUnassigned && approxFactors ) {
		std::cout << ((Duration) ( Clock::now() - m_starttime )).count() <<
				": assigning vars (";
		printVars( vds );
		std::cout << ") blocks (";
		for ( const VarData * vd : vds ) {
			std::cout << ( vd == vds.front() ? "" : "," );
			std::cout << m_func.getBlockID( vd->v->getID() );
		}
		std::cout << ") val (" << xval <<
				( wasUnassigned ? ") (was unassigned)" : "" ) <<
				" for ccomp " << ccomp << std::endl;
	}

	if ( numassign[vds.front()->v->getID()] == 0 ) {
		assignmentOrder += "(";
		BOOSTNS::format fmt1( "%1%, " ), fmt2( "%1%" );
		for ( VarData * vd : vds ) {
			assignmentOrder += BOOSTNS::str( ( vd == vds.back() ? fmt2 : fmt1 )
					% vd->v->getID() );
			++numassigned;
		}
		assignmentOrder += "), ";
	}

	for ( VarData * vd : vds ) ++numassign[vd->v->getID()];

	if ( !printedAO && numassigned >= (VariableCount) numassign.size() ) {
		std::cout << "assignment ordering: " << assignmentOrder << std::endl;
		printedAO = true;
	}

	// perform the assignment and update the connectivity graph
	for ( size_t i = 0; i < vds.size(); ++i ) {
		vardata.assign( vds[i]->v->getID(), xval[i] );
		m_func.onVarAssigned( vds[i]->v->getID(), xval[i] );

		if ( wasUnassigned ) connGraph.onVariableAssigned( vds[i]->v );
	}

	ccomp.onVarsAssigned( vds, xval );

	// simplify all the factors that these variables are in
	if ( approxFactors && simplifyFactors( ccomp, vds ) ) {
		std::cout << "num simplified factors: " << nSimplifiedFactors <<
				" / " << m_func.getFactors().size() << ", num assigned vars: "
				<< vardata.numAssigned() << " after " <<
				(Duration) ( Clock::now() - m_starttime ) << std::endl;
	}

    SubdomainSP sd = sdpool.get();
    for ( VarData * vd : vds ) {
        sd->vids.push_back( vd->v->getID() );
        sd->x.push_back( vd->v->eval() );
    }
    ccomp.setCurSD( sd );

	// compute the partial evaluation of f(x) caused by assigning these vars
	Numeric ferr( 0.0 );
	Numeric pfeval = m_func.evalFactors( ccomp.getFactors(), ferr );

	ccomp.setPartialEval( pfeval );
	ccomp.updateCurSD( pfeval, ferr, false );
}


void RDISOptimizer::unassign( Component & ccomp, const VarDataPtrVec & vds ) {

	if ( vardata.isAssigned( vds.front()->v->getID() ) ) {
//		std::cout << "unassigning vars "; printVars( vds ); std::cout << std::endl;

		// unassign all of the factors that were assigned
		if ( approxFactors && unsimplifyFactors( ccomp ) ) {
			std::cout << "num simplified factors: " << nSimplifiedFactors <<
					" / " << m_func.getFactors().size() << " after " <<
					(Duration) ( Clock::now() - m_starttime ) << std::endl;
		}

		// unassign the variables -- in reverse assigned order -- and reconnect
		// them to their neighbours
		for ( auto rit = vds.rbegin(), rend = vds.rend(); rit != rend; ++rit ) {
			vardata.unassign( (*rit)->v->getID() );
			m_func.onVarUnassigned( (*rit)->v->getID() );
			connGraph.onVariableUnassigned( (*rit)->v );
		}

		// finish up processing of this component
		ccomp.onVarUnassigned();
	}
}


void RDISOptimizer::unassign( VariableID vid ) {
	if ( vardata.isAssigned( vid ) ) vardata.unassign( vid );
}


bool RDISOptimizer::simplifyFactors( Component & ccomp,
		const VarDataPtrVec & vds ) {
	assert( ccomp.haveVarsBeenAssigned() );

	bool changed = false;

	SortedFactorIDs & smplfact = ccomp.getSimplifiedFactors();
	const FactorPtrVec & fctrs( ccomp.getFactors() );

	for ( Factor * f : fctrs ) {
		bool skip = true;
		for ( VarData * vd : vds ) {
			if ( f->containsVar( vd->v->getID() ) ) {
				skip = false;
				break;
			}
		}
		if ( skip ) continue;

		const bool wasSimplified = f->isSimplified();

		bool needsUnsimp = f->updateSimplification(
//				10.0 * epsilon / (Numeric) factors.size(),
				epsilon,
				ccomp.getAssignedVids(), &connGraph );

		if ( f->isSimplified() ) {
			if ( !wasSimplified || needsUnsimp ) {
				++nSimplifiedFactors;
				smplfact.insert( f->getID() );
				changed = true;
			}
		} else {
			if ( wasSimplified ) {
				smplfact.erase( f->getID() );
				--nSimplifiedFactors;
				assert( nSimplifiedFactors >= 0 );
				changed = true;
			}
		}
	}

	return changed;
}


bool RDISOptimizer::unsimplifyFactors( Component & ccomp ) {
	SortedFactorIDs & simplfactors = ccomp.getSimplifiedFactors();
	bool changed = false;

	if ( !simplfactors.empty() ) {
		const FactorPtrVec & factors = m_func.getFactors();
		for ( auto it = simplfactors.rbegin(); it != simplfactors.rend(); ++it ) {
			assert( factors[*it]->isSimplified() );
			factors[*it]->unsimplify( ccomp.getAssignedVids(), &connGraph );
			--nSimplifiedFactors;
		}
		changed = true;
	}

	simplfactors.clear();
	return changed;
}


Numeric RDISOptimizer::simplifyFactors( Numeric & ferr ) {
	Numeric feval = 0.0;
	for ( Factor * f : m_func.getFactors() ) {
		f->updateSimplification(
//				10.0 * epsilon / (Numeric) factors.size(),
				epsilon,
				VariableIDVec( 1, -1 ), NULL );

		if ( f->isAssigned() ) {
			feval += f->eval();
			++nSimplifiedFactors;
		}

		ferr += f->getSimplificationError();
	}
	return feval;
}


void RDISOptimizer::unsimplifyFactors() {
	for ( Factor * f : m_func.getFactors() ) {
		if ( f->isSimplified() && f->getAssignedKey() == -1 ) {
			f->unsimplify( VariableIDVec( 1, -1 ), NULL );
			--nSimplifiedFactors;
		}
	}
}


void RDISOptimizer::productComponentVals( Component & parent,
		Component & child ) {

	assert( parent.getEpsilon() == epsilon && child.getEpsilon() == epsilon );
	assert( child.hasBeenEvald() );

	const SubdomainCSP pcursd = parent.getCurSD();
	const SubdomainCSP & coptsd = child.getOptSD();

	// perform the summation and keep track of whether or not this is a bound or
	// is (partially) from the cache
	parent.updateCurSD( ( pcursd->fx + coptsd->fx ),
						( pcursd->ferr + coptsd->ferr ),
			            ( pcursd->isBound || coptsd->isBound ),
			            coptsd );

	assert( coptsd->isBound || child.haveVarsBeenAssigned() ||
			child.getFactors().empty() );

//	std::cout << "component " << child << " eval is " << coptsd->fx
//			<< " (" << ( coptsd->isBound ? "bound)" :
//					"not bound) " ) << " at " << *coptsd << std::endl;

	++instrumentation::g_numProducts;
}


void RDISOptimizer::updateDomain( Component & ccomp ) {

	assert( ccomp.haveVarsBeenAssigned() );

	VarData & vd( *ccomp.getRepresentativeVD() );
	assert( vd.v->isAssigned() );

	SubdomainCSP sdopt = ccomp.getOptSD();
	SubdomainCSP newsd = ccomp.getCurSD();

	Numeric ftol = epsilon / 10.0;
	if ( options.count( "steptol" ) ) ftol = options[ "steptol" ].as< Numeric >();

	bool isNewMin = false;
	if ( !sdopt || ( newsd->fx < sdopt->fx &&
			( !newsd->isBound || sdopt->isBound ) ) ) {
		isNewMin = !sdopt || !approxeq( newsd->fx, sdopt->fx, ftol );
		ccomp.setOptSD( newsd );
	}

	// print this for top ccomps
	if ( ccomp.getParent() && !ccomp.getParent()->haveVarsBeenAssigned() ) {
//		std::cout << "NEW EVAL: " << newsd->fx << " (" <<
//			( sdopt ? sdopt->fx : nan( "" ) ) << ", progress: " <<
//			( newsd->fx - ( ccomp.getPrevSD() ? ccomp.getPrevSD()->fx : 0 ) )
//			<< ") after " <<
//			(Duration) ( Clock::now() - m_starttime ) <<
////				" for ccomp " << ccomp /*<< "\n\t" << *newsd*/ <<
//			std::endl;

		if ( ccomp.getNumVAsinceLastRR() > 2 && sdopt && !newsd->isBound &&
			 !approxleq( newsd->fx, ccomp.getPrevSD()->fx, epsilon )
//				newsd->fx >= ccomp.getPrevSD()->fx
				) {
			std::cout << "NEW EVAL NOT CONVERGING PROPERLY -- " <<
				ccomp.getNumVAsinceLastRR() << ", nRRs " <<
				ccomp.getNumRandomRestarts() << " -- " << newsd->fx << " > "
				<< ccomp.getPrevSD()->fx << ", " <<
				newsd->fx - ccomp.getPrevSD()->fx << ", is bound? "
				<< newsd->isBound << std::endl;
			ccomp.setNonmonotoneEvals( true );
		}
	}

	// update fmin/fmax with a better value if one was found
	if ( newsd->fx < ccomp.getFMin() ) ccomp.setFMin( newsd->fx );

	ccomp.setLastEvalWasOpt( isNewMin );

	++instrumentation::g_numSums;
	const bool alwaysPrint = options[ "printAllFullEvals" ].as< bool >();

	if ( alwaysPrint || !ccomp.getChildren().empty() ||
		 ( ccomp.getParent() && !ccomp.getParent()->haveVarsBeenAssigned() ) ) {
		Numeric cureval = 0;
		bool res = computeCurrentFunctionEval( ccomp, newsd, cureval );
		if ( res ) {
			std::cout << "new full eval: " << cureval << ", best " <<
				bestFullFunctionEval << " (diff: " <<
				( cureval - bestFullFunctionEval ) <<
				"), prev " << prevFullFunctionEval << " (diff: " <<
				( cureval - prevFullFunctionEval ) <<
				") after " << (Duration) ( Clock::now() - m_starttime ) <<
				std::endl;
			bestFullFunctionEval = std::min( cureval, bestFullFunctionEval );
			prevFullFunctionEval = cureval;
		}
	}

	ccomp.setPrevSD( newsd );
}


bool RDISOptimizer::computeCurrentFunctionEval( Component & ccompIn,
		const SubdomainCSP & newsd, Numeric & curFuncEval ) {

	State & x = ccompIn.getCurFuncEvalState();

	if ( !ccompIn.getCurFEStateIsSet() ) {
		x.assign( vardata.numVars(), NumericMAX );

		// get the top component
		const Component * topcc = &ccompIn;
		while ( topcc->getParent() != NULL ) topcc = topcc->getParent();

		assert( topcc->getParent() == NULL && !topcc->haveVarsBeenAssigned() );

		const bool ccinIsTop = ( ccompIn.getParent() == topcc );
		const Subcomponents & topccomps = topcc->getChildren();

		std::queue< const Component * > ccs;

		for ( const ComponentIP & cc : topccomps ) {
//			if ( !cc->haveVarsBeenAssigned() ) {
//				if ( ccinIsTop ) {
//					std::cout << "can't compute current function eval because uncomputed component ("
//							<< *cc << std::endl;
//				}
//				return false;
//			} else if ( !cc->getOptSD() ) {
//				if ( ccinIsTop ) {
//					std::cout << "can't compute current function eval because component has no eval ("
//							<< *cc << std::endl;
//				}
//				return false;
//			}

			ccs.push( cc.get() );
		}

		for ( VariableID vid = 0; vid < vardata.numVars(); ++vid ) {
			VarData & vd = vardata.get( vid );
			if ( vd.xvinitIsSet ) x[vid] = vd.xvalinit;
		}

		// make sure we copy from the top down, because the bottom ccomps have
		// more up-to-date information so they'll overwrite the out-of-date info
		while ( !ccs.empty() ) {
			const Component * cc = ccs.front();
			ccs.pop();

			if ( cc->haveVarsBeenAssigned() ) {
				if ( !cc->hasBeenEvald() ) {
					for ( const ComponentIP & c : cc->getChildren() ) {
						ccs.push( c.get() );
					}

					if ( cc->getPrevSD() ) {
						copySubdomain( *cc->getPrevSD(), x, true, false );
					}

				} else if ( cc->getOptSD() ) {
					copySubdomain( *cc->getOptSD(), x,
								   true, false );
				}
			}
		}

		for ( VariableID vid = 0; vid < vardata.numVars(); ++vid ) {
			// even more up-to-date are those variables that are currently assigned
			if ( vardata.isAssigned( vid ) ) {
				x[vid] = vardata.var( vid )->eval();
			}

			if ( x[vid] == std::numeric_limits< Numeric >::max() ) {
//				std::cout << "can't compute current function eval because "
//						"unknown variable state " << vid << std::endl;
				return false;
			}
		}

	}

	// do this last so it overwrites the static values from before
	if ( newsd ) {
		copySubdomain( *newsd, x, true, false );
	} else if ( ccompIn.getPrevSD() ) {
		copySubdomain( *ccompIn.getPrevSD(), x, true, false );
	}

	curFuncEval = m_func.eval( x );

	ccompIn.setCurFEStateIsSet( true );

//	std::cout << "at state: " << x << std::endl;

	return true;
}


void RDISOptimizer::addToCache( Component & ccomp ) {

	assert( useCache );
	if ( ccomp.dontCache() || ccomp.getParent() == NULL
			|| !ccomp.getParent()->haveVarsBeenAssigned()
			|| ccomp.getCacheKey().factors.empty() ) {
		return;
	}

	assert( ccomp.getRepresentativeVD() != NULL );
	const SubdomainCSP & sdopt( ccomp.getOptSD() );
	if ( !sdopt ) return;

	assert( !ccomp.getFullEvalWasCached() );

	cache.insert( ccomp.getCacheKey(), sdopt );
	++instrumentation::g_cacheInserts;
}


void RDISOptimizer::removeChildrenFromCache( Component & /*ccomp*/ ) {
	assert( false ); // not supported yet
}


void RDISOptimizer::setInitialState( const State & xinit ) {
	if ( xinit.empty() ) return;
	assert( xinit.size() == vardata.numVars() );

	for ( VariableID vid = 0; vid < vardata.numVars(); ++vid ) {
		vardata.get( vid ).xvalinit = xinit[vid];
		vardata.get( vid ).xvinitIsSet = true;
	}
}


void RDISOptimizer::setInitialValFromChildren( const Subdomain & sd ) {
	for ( const SubdomainCSP & c : sd.children ) {
		for ( size_t i = 0; i < c->vids.size(); ++i ) {
			if ( vardata.isUnassigned( c->vids[i] ) ) {
				VarData & vd = vardata.get( c->vids[i] );
				vd.xvalinit = c->x[i];
				vd.xvinitIsSet = true;
			}
		}
		setInitialValFromChildren( *c );
	}
}


void RDISOptimizer::init( size_t cachesize ) {
	cache.init( cachesize );
	vardata.init( m_func.getVariables() );
}


void RDISOptimizer::processOptions() {

	RDISOptionValues & rov = optionVals;

	if ( options.count( "epsilon" ) ) {
		epsilon = options[ "epsilon" ].as< Numeric >();
	}

	if ( options.count( "doBnB" ) && !options[ "doBnB" ].as< bool >() ) {
		useBounds = false;
	}

	if ( options.count( "timeAsSeed" ) &&
			options[ "timeAsSeed" ].as< bool >() ) {
		rov.timeAsSeed = true;
		rng.seed( std::time(NULL) + getpid() );
	}

	rov.useStaticDecomp = options.count( "staticDecomp" )
			&& options[ "staticDecomp" ].as< bool >();


	uint32 & dnv = optionVals.desiredNumVars;
	dnv = (uint32) ( (Numeric) vardata.numVars() * 0.2 );

	if ( options.count( "AVblksz" ) ) {
		dnv = options[ "AVblksz" ].as< uint32 >();
	} else if ( options.count( "AVblkpct" ) ) {
		dnv = (uint32) round( options[ "AVblkpct" ].as< Numeric >() *
				(Numeric) vardata.numVars() );
	}
	dnv = std::max( std::min( dnv, (uint32) vardata.numVars() ), 1u );

	optionVals.heur = PaToH;
	if ( options.count( "varChoiceHeur" ) ) {
		string heur = options[ "varChoiceHeur" ].as< string >();
		if ( heur == "hmetis" ) {
			optionVals.heur = HMetis;
		} else if ( heur == "random" ) {
			optionVals.heur = Random;
		} else if ( heur == "all" ) {
			optionVals.heur = All;
		}
	}

	rov.minNumRestarts = options.count( "minRR" ) ?
						 options[ "minRR" ].as< uint32 >() : 1;

	rov.nRRperLvl = options.count( "nRRperLvl" ) ?
					options[ "nRRperLvl" ].as< uint32 >() : 2;

	rov.nRRatTop = options.count( "nRRatTop" ) ?
				   options[ "nRRatTop" ].as< uint32 >() : rov.nRRperLvl;

	rov.maxNumAssignsUntilRR = options.count( "maxNAtoRR" ) ?
							   options[ "maxNAtoRR" ].as< uint32 >() : 10;

	rov.alwaysUseAllFactors = ( options.count( "alwaysUseAllFactors" ) &&
			options[ "alwaysUseAllFactors" ].as< bool >() );

	rov.noAssignLimitAtTop = options.count( "noAssignLimitAtTop" ) ?
						 options[ "noAssignLimitAtTop" ].as< bool >() : true;
}


void RDISOptimizer::prepareToOptimize() {
	for ( VariableID vid = 0; vid < vardata.numVars(); ++vid ) {
		if ( vardata.isAssigned( vid ) ) unassign( vid );
	}
	assert( vardata.numAssigned() == 0 );
	assert( vardata.numUnassigned() == vardata.numVars() );

	numassigned = 0;
	numassign.assign( vardata.numVars(), 0 );
	numcachehits.assign( vardata.numVars(), 0 );
	numboundhits.assign( vardata.numVars(), 0 );

	cache.clear();

	// reset the connectivity graph if the factors got reset
	vardata.reset();
	if ( m_func.resetFactorsOnChange() || cgNeedsInit ) {
		connGraph.clear();
		cgNeedsInit = true;
	}

	m_optState.clear();
	m_optState.resize( vardata.numVars(), 0 );
	m_optTime = BOOSTNS::chrono::seconds( 0 );

	m_optVal = std::numeric_limits< Numeric >::max();
	m_optError = std::numeric_limits< Numeric >::max();

	instrumentation::clearInstrumentation();

	m_timedOut = false;
	m_timerCheckCount = 0;
	m_wasBound = false;

	nSimplifiedFactors = 0;
	assignmentOrder = "";
	printedAO = false;

	bestFullFunctionEval = NumericMAX;
	prevFullFunctionEval = NumericMAX;

	m_optState.assign( vardata.numVars(), NumericMAX );

	for ( const Variable * v : m_func.getVariables() ) {
		if ( v->isAssigned() ) m_optState[v->getID()] = v->eval();
	}

	bool isNewTopCC = false;
	if ( !topccomp || !optionVals.useStaticDecomp ) {
		topccomp = ccpool.get( m_func, m_func.getFactors(), &vardata,
				&connGraph, optionVals.useStaticDecomp );
		isNewTopCC = true;
	}

	if ( isNewTopCC ) {
		topccomp->init( epsilon, 0, m_func.getVariables() );
	} else {
		topccomp->reinit( epsilon, NULL );
	}
}


bool RDISOptimizer::checkTimedOut( Component & ccomp ) {
	if ( super::checkTimedOut() ) {
		if ( ! ( ccomp.haveVarsBeenAssigned() && ccomp.getOptSD() ) ) {
			NumericInterval b = ( ccomp.haveVarsBeenAssigned() ?
					  ccomp.getAssignedBounds() : ccomp.getUnassignedBounds() );
			ccomp.setOptSD( sdpool.get(), b.upper(), true );
		}
		ccomp.setHasBeenEvald( true );
		return true;
	}
	return false;
}


void RDISOptimizer::printStats( std::ostream & os, bool printInfo ) const {
	using namespace instrumentation;
	const bool pi( printInfo );
	os << ( pi ? "instrumentation: " : "" ) <<
		g_runTime.count() << ( pi ? " secs, " : ", " ) <<
		g_numVarOps << ( pi ? " var ops, " : ", " ) <<
		g_numFactEvals << ( pi ? " fevals, " : ", " ) <<
		g_numProducts << ( pi ? " prods, " : ", " ) <<
		g_numSums << ( pi ? " sums, " : ", " ) <<
		g_cacheInserts << ( pi ? " (" : ", " ) <<
		g_cacheMidInserts << ( pi ? ")" : "" ) << ( pi ? " inserts, " : ", " ) <<
		g_cacheHits << ( pi ? " / " : ", " ) << g_cacheLookups <<
		( pi ? " hits, " : ", " ) <<
		g_numDownhillSteps << ( pi ? " / " : ", " ) << g_numSteps <<
		( pi ? " downhill, " : ", " ) <<
		"(" << g_counter1 << ", " << g_counter2 << ")";
	if ( pi ) os << std::endl;
}


void RDISOptimizer::writeGraph() const {

	fs::fstream ofs( "variable_block_graph.gv", std::ios_base::out );
	if ( ofs.is_open() ) {
		ofs << "strict graph G {" << std::endl;
		ofs << "layout=sfdp ;" << std::endl;
		ofs << "overlap = prism ;" << std::endl;

		const VariableCount nblocks = m_func.getNumBlocks();
		for ( VariableCount blkid = 0; blkid < nblocks; ++blkid ) {
			ofs << "b" << blkid << " [color=green, style=filled] ; " << std::endl;
		}

		BOOSTNS::container::vector< BOOSTNS::dynamic_bitset<> > edges(
				nblocks, BOOSTNS::dynamic_bitset<>( nblocks, false ) );

		for ( VariableID vid = 0; vid < vardata.numVars(); ++vid ) {
			if ( vardata.isAssigned( vid ) ) continue;
			VariableCount blkid = m_func.getBlockID( vid );

			for ( const Factor * f : vardata.var( vid )->getFactors() ) {
				if ( f->isAssigned() || !f->containsVar( vid ) ) continue;
				for ( const Variable * v2 : f->getVariables() ) {
					VariableID vid2 = v2->getID();
					VariableCount blkid2 = m_func.getBlockID( vid2 );
					if ( edges[blkid][blkid2] || vardata.isAssigned( vid2 ) ) continue;
					if ( blkid != blkid2 && f->containsVar( vid2 ) ) {
						ofs << "b" << blkid << " -- " << "b" << blkid2 << " ; "
								<< std::endl;
						edges[blkid][blkid2] = true;
						edges[blkid2][blkid] = true;
					}
				}
			}
		}


		ofs << " } ";
		ofs.close();
	}
}


void RDISOptimizer::writeFactorGraph( const Component & ccomp ) const {

	if ( ccomp.getVars().empty() ) return;

	string sgid = BOOSTNS::str( BOOSTNS::format( "%1%_%2%.%3%" )
			% ccomp.getVars().front()
			% ccomp.getVars().size()
			% ccomp.getFactors().size() );

	BOOSTNS::filesystem::fstream ofs( "logs/factor_graph" + sgid + ".gv",
			std::ios_base::out );
	if ( ofs.is_open() ) {
		ofs << "strict graph G {" << std::endl;
		ofs << "layout=sfdp ;" << std::endl;
		ofs << "overlap=prism ;" << std::endl;

//		for ( VariableID vid = 0; vid < vardata.numVars(); ++vid ) {
		for ( VariableID vid : ccomp.getVars() ) {
			if ( vardata.isAssigned( vid ) ) continue;
			ofs << "x" << vid << " [color=blue, style=filled] ; " << std::endl;
		}

//		for ( FactorID fid = 0; fid < factors.size(); ++fid ) {
////		for ( const Factor * f : ccomp.getFactors() ) {
//			if ( factors[fid]->isAssigned() ) continue;
////			if ( factors[fid]->getVariables().size() == 1 ) continue;
//			ofs << "f" << fid << " [shape=box, color=gray, style=filled] ; " << std::endl;
//		}

//		for ( FactorID fid = 0; fid < (FactorID) m_func.getFactors().size(); ++fid ) {
		for ( const Factor * f : ccomp.getFactors() ) {
//			Factor * f = m_func.getFactors()[fid];
			if ( f->isAssigned() ) continue;
			bool addedFactor = false;
//			if ( f->getVariables().size() == 1 ) continue;
			for ( const Variable * v : f->getVariables() ) {
				if ( v->isAssigned() || !f->containsVar( v->getID() ) ) continue;
				ofs << "x" << v->getID() << " -- " << "f" << f->getID() << " ; "
						<< std::endl;
				if ( !addedFactor ) {
					ofs << "f" << f->getID()
							<< " [shape=box, color=gray, style=filled] ; "
							<< std::endl;
					addedFactor = true;
				}
			}
		}

		ofs << " } ";
		ofs.close();
	}
}


void RDISOptimizer::printVars( const VarDataPtrVec & vds ) const {
	for ( const VarData * vd : vds ) {
		std::cout << ( vd == vds.front() ? "" : ", " ) << vd->v->getID();
	}
}


std::ostream & operator<<( std::ostream & out, RDISOptimizer const & o ) {
	o.printStats( out, true );
	return out;
}

} // namespace rdis
