/*
 * optBA.cpp
 *
 *  Created on: Sep 10, 2014
 *      Author: afriesen
 */

#include "common.h"

#include "OptimizableFunction.h"
#include "optimizers/CGDSubspaceOptimizer.h"
#include "optimizers/LMSubspaceOptimizer.h"

#include "bundleadjust/BundleAdjustmentFactor.h"
#include "bundleadjust/BundleAdjustmentFunction.h"

#include "Optimizer.h"
#include "RDISOptimizer.h"

#include "optimizers/BCDOptimizer.h"

#include "util/utility.h"
#include "util/signal_handler.h"

#include BOOSTPATH/foreach.hpp>
#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/filesystem.hpp>

#include BOOSTPATH/assign/list_of.hpp> // for 'list_of()'
#include BOOSTPATH/assign/std/vector.hpp> // for 'operator+=()'

#include BOOSTPATH/math/distributions/normal.hpp>


namespace fs = BOOSTNS::filesystem;
namespace bc = BOOSTNS::chrono;
namespace psi = rdis::instrumentation;
namespace po = BOOSTNS::program_options;

using namespace rdis;



void optimize_BA( BOOSTNS::random::mt19937 & rng,
		const po::variables_map & options, const string & pdbfile );

bool doBCDOpt( OptimizableFunction & func, SubspaceOptimizer & ssopt,
		const po::variables_map & options, const State & xinit,
		Clock::time_point timeout, Clock::time_point starttime,
		State & xopt, Numeric & fopt );

void addOptions( po::options_description & desc );


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main( int argc, const char * const argv[] )
{
	BOOSTNS::random::mt19937 rng( 834725927 );
	const Clock::time_point start = Clock::now();

	std::cout << "called with: " << std::endl << "\t";
	for ( int i = 0; i < argc; ++i ) {
		std::cout << argv[i] << " ";
	}
	std::cout << std::endl;

	po::options_description desc( "Bundle Adjustment options" );

	try {
		addOptions( desc );

		po::variables_map options;
		po::store( po::parse_command_line( argc, argv, desc ), options );

//		po::parsed_options parsed =
//				po::command_line_parser( argc, argv ).options( desc ).allow_unregistered().run();
//		po::store( parsed, vm );
//		std::vector< string > to_pass_further =
//				po::collect_unrecognized( parsed.options, po::exclude_positional );

		po::notify( options );

		if ( options.count( "help" ) ) {
			std::cout << desc << std::endl;

		} else if ( !options.count( "file" ) ) {
			std::cout << "\nat least one input BA file must be specified\n" << std::endl;
			std::cout << desc << std::endl;

		} else {
			// get the set of files to optimize
			const std::vector< string > files =
					options[ "file" ].as< std::vector< string > >();

			// optimize each one
			for ( string file : files ) {
				optimize_BA( rng, options, file );
			}
		}

	} catch ( po::error & e) {
		  std::cerr << "\nCommand line error: " << e.what() << std::endl << std::endl;
		  std::cerr << desc << std::endl;
	} catch ( std::exception & e ) {
		std::cerr << "\nException occurred in main (B.A. opt): \n\t" << e.what()
				<< ", application will now exit. Goodbye!" << std::endl;
	}

	const Duration elapsed = Clock::now() - start;
	std::cout << "total elapsed time: " << elapsed.count() << " seconds" << std::endl;
}



void optimize_BA( BOOSTNS::random::mt19937 & rng,
		const po::variables_map & options, const string & bafile ) {

	// get optimization type from command line
	const string opttype = options[ "opttype" ].as< string >();

	const bool isRDIS = ( opttype == "rdis" );
	const bool isBCD = ( opttype == "bcd" );
//	const bool isDisc = ( opttype == "discrete" );
//	const bool isGrid = ( opttype == "grid" );

	Numeric regionwidth = 0.0;
	Numeric eps = 1.0, lipsch = 20;
	RDISOptimizer * popt = NULL;

//	regionwidth = ( isRDIS || isBCD ? 2500 : ( isDisc ? 0.0 : 5.0 ) );
	regionwidth = 2500;

	assert( isRDIS || isBCD /*|| isDisc || isGrid*/ );

	eps = options.count( "epsilon" ) ? options[ "epsilon" ].as< Numeric >() : eps;
	lipsch = options.count( "lipschitz" ) ? options[ "lipschitz" ].as< Numeric >() : lipsch;

	if ( options.count( "timeAsSeed" ) && options[ "timeAsSeed" ].as< bool >() ) {
		rng.seed( std::time(NULL) + getpid() );
		if ( options[ "sampleOffset"].as< size_t >() > 0 ) {
			std::cout << "WARNING: sampleOffset and timeAsSeed should not both be set" << std::endl;
		}
	}

	VariableCount ncams = options.count( "ncams" ) ?
						  options["ncams"].as< size_t >() : 0;
	VariableCount npts = options.count( "npts" ) ?
						 options["npts"].as< size_t >() : 0;

	// load and initialize the BA problem we'll optimize
	BundleAdjustmentFunction bafunc;
	if ( !bafunc.load( bafile, ncams, npts ) ) {
		std::cout << "loading file " << bafile << " failed -- exiting" << std::endl;
		return;
	}

	SubspaceOptimizer * pssopt;
	if ( options["useCGD"].as< bool >() ) {
		pssopt = new CGDSubspaceOptimizer( bafunc );
	} else {
		pssopt = new LMSubspaceOptimizer( bafunc );
	}
	SubspaceOptimizer & ssopt = *pssopt;
////	CGDSubspaceOptimizer ssopt( bafunc );
//	LMSubspaceOptimizer ssopt( bafunc );
	ssopt.setParameters( options );

	// create the optimizer and initialize it properly
//	popt = ( isRDIS ?
//			new OptimizerISGD( bafunc, ssopt, options ) :
//			new OptimizerDG( bafunc, options ) );
	popt = new RDISOptimizer( bafunc, ssopt, options );

	RDISOptimizer & opt = *popt;
	setSignalHandler( &opt );

//	opt.setOptType( 3 );
//	opt.setUseFiniteDiffs( true );
//	opt.setUseLinearLB( false );

	opt.setApproxFactors( options[ "approxfactors" ].as< bool >() );
	opt.setUseCache( !isRDIS );
//	opt.setUseExactCaching( isDisc );

	const Numeric MaxVal = std::numeric_limits< Numeric >::max();

	size_t sampleOffset = options[ "sampleOffset" ].as< size_t >();
	size_t nsamples = options[ "nsamples" ].as< size_t >() + sampleOffset;

	std::vector< State > initialStates( nsamples );

	if ( options.count( "randinit" ) && !options["randinit"].as< bool >() ) {
		sampleOffset = 0;
		nsamples = 1;
		initialStates.clear();
		initialStates.push_back( bafunc.getInitialState() );

	} else {
		// create all the random initial states up front (so they're consistent
		// across runs if the seed doesn't change)
		for ( size_t i = 0; i < nsamples; ++i ) {
			for ( VariableID vid = 0; vid < bafunc.getNumVars(); ++vid ) {
				NumericInterval si =
						bafunc.getVariables()[vid]->getDomain().getSamplingInterval();
				BOOSTNS::random::uniform_real_distribution<> unif(
						si.lower(), si.upper() );
				initialStates[i].push_back( unif( rng ) );
			}

//			if ( i >= sampleOffset ) {
//				std::cout << "state " << i << ": " << initialStates[i] << std::endl;
//			}
		}
	}

	const Clock::time_point start = Clock::now();
	Clock::time_point timeout = Clock::time_point::max();

	if ( options.count( "timeout" ) ) {
		size_t ms = 1000.0 * options[ "timeout" ].as< Numeric >();
		timeout = Clock::now() + bc::milliseconds( ms );
	}

	Numeric bestfopt = NumericMAX;
	State bestx;

	size_t i = sampleOffset;
	for ( ; i < nsamples; ) {

		Numeric finit = MaxVal;
		const State & xinit = initialStates[i];

		std::cout << "[" << i << "] optimization starting -- bounds: " <<
				bafunc.computeBounds() << std::endl;
		std::cout << "eps: " << eps << ", lipschitz: " << lipsch <<
				", step: " << eps / lipsch <<
				", finit: " << finit << std::endl << "xinit: " << xinit << std::endl;

		std::cout << "initial eval: " << bafunc.eval( xinit ) << std::endl;

		if ( options.count( "useDA" ) ) {
			throw "deterministic annealing not yet supported";
		} else if ( options.count( "useCP" ) ) {
			throw "cutting plane not yet supported";
		}

		++i;

		State xopt;
		Numeric result = MaxVal;
		bool timedOut = false;

		if ( !isBCD ) {
			if ( options["ignoreCamVars"].as< bool >() ) {
				for ( Variable * v : bafunc.getVariables() ) {
					if ( bafunc.getVarType( v->getID() ) == BundleAdjust::CAM_FOCAL ||
							bafunc.getVarType( v->getID() ) == BundleAdjust::CAM_RDL_K1 ||
							bafunc.getVarType( v->getID() ) == BundleAdjust::CAM_RDL_K2 ) {
						bafunc.getVariables()[i]->assign( bafunc.getInitialState()[i] );
					}
				}
			}
			result = opt.optimize( &xinit, finit, timeout, true, start );
			xopt = opt.getOptState();
			opt.printStats();
			timedOut = ( opt.wasStopped() || opt.timedOut() );
		} else {
			if ( options["ignoreCamVars"].as< bool >() ) {
				std::cout << "IGNORECAMVARS OPTION NOT SUPPORTED FOR BCD" << std::endl;
			}
			timedOut = doBCDOpt( bafunc, ssopt, options, xinit, timeout, start,
					xopt, result );
		}

		assert( !xopt.empty() );

		const Duration elapsed = Clock::now() - start;

		std::cout << "bundle adjustment opt. result: " << result << std::endl
				<< "\tat " << xopt << std::endl;
		Numeric actual = bafunc.eval( xopt );
		std::cout << "actual eval: " << actual << " after "
				<< elapsed.count() << " seconds" << std::endl;

		if ( actual < bestfopt ) {
			bestfopt = actual;
			bestx = xopt;
		}

		if ( timedOut ) {
			std::cout << "timeout detected..." << std::endl;
			break;
		}
	}

	const Duration elapsed = Clock::now() - start;
	std::cout << "----------------" << std::endl;
	std::cout << "best: " << bestfopt << " in " << ( i - sampleOffset ) <<
			" samples after " << elapsed.count() << " seconds at state " <<
			std::endl << bestx << std::endl;
}


bool doBCDOpt( OptimizableFunction & func, SubspaceOptimizer & ssopt,
		const po::variables_map & options, const State & xinit,
		Clock::time_point timeout, Clock::time_point starttime,
		State & xopt, Numeric & fopt ) {

	size_t blksize = std::ceil( (Numeric) func.getNumVars() / 5.0 );
	if ( options.count( "AVblksz" ) ) {
		blksize = options["AVblksz"].as< size_t >();
	} else if ( options.count( "AVblkpct" ) ) {
		blksize = std::ceil( options["AVblkpct"].as< Numeric >() *
				( (Numeric) func.getNumVars() ) );
	}

	Numeric tol = 1e-6;
	if ( options.count( "steptol" ) ) tol = options["steptol"].as< Numeric >();

	size_t maxit = 100;
	if ( options.count( "SSmaxit" ) ) maxit = options["SSmaxit"].as< VariableCount >();

	std::cout << "BCD OPTIMIZATION (block size: " << blksize << ", numvars: " <<
			func.getNumVars() << ")" << std::endl;

	fopt = NumericMAX;
	BCDOptimizer bcdopt( func, ssopt, blksize, tol );

	fopt = bcdopt.optimize( maxit, 100000, &xinit, true, timeout, 0, starttime );

	xopt = bcdopt.getOptState();
	std::cout << "BCD opt eval " << fopt << " in " << bcdopt.getRuntime() << std::endl;

	return bcdopt.getTimedOut();
}


// add all options to boost program_options
void addOptions( po::options_description & desc ) {

	desc.add_options()
			( "file,f", po::value< std::vector< string > >(),
					"input file(s) containing bundle adjustment problems to optimize (can be multiple)" )
			( "ncams", po::value< size_t >(),
				  "number of cameras to load from the file" )
			( "npts", po::value< size_t >(),
				  "number of points to load from the file" )
			( "opttype", po::value< string >()->default_value( "rdis" ),
					"type of optimization to perform (can be: rdis, bcd, discrete, grid" )
			( "nsamples", po::value< size_t >()->default_value( 1 ),
					"number of times to optimize each specified PDB file" )
			( "sampleOffset", po::value< size_t >()->default_value( 0 ),
					"sample number to start at (0-indexed)" )
			( "randinit", po::value< bool >()->default_value( true ),
					"BA with random initialization (don't use existing values)" )
			( "useCGD", po::value< bool >()->default_value( true ),
					"optimize BA with CGD instead of LM" )
			( "ignoreCamVars", po::value< bool >()->default_value( false ),
					"true to ignore the focal length and radial distortion params during the optmization" )
//			( "perruntimeout", po::value< Numeric >(),
//					"per-run optimization timeout, in seconds (applies to det. annealing and cutting plane)" )
//			( "useDA", po::value< bool>()->default_value( false ),
//					"flag to enable/disable use of deterministic annealing" )
//			( "useCP", po::value< bool >()->default_value( false ),
//					"flag to enable/disable use of cutting plane optimization" )
//			( "CPftol", po::value< Numeric >(),
//					"the initial tolerance for unsetting factors when using cutting plane" )
					;

	RDISOptimizer::addOptions( desc );
}
