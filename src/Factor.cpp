/*
 * Factor.cpp
 *
 *  Created on: Jan 12, 2013
 *      Author: afriesen
 */

#include "Factor.h"
#include "ConnectivityGraph.h"

#include BOOSTPATH/assign/ptr_map_inserter.hpp>

namespace rdis {

Factor::Factor( FactorID factorID )
	: id( factorID )
	, isPartlySimplified( false )
	, isAssignedConstant( false )
	, vidAssigned( -1 )
	, assignmentError( 0.0 )
	, simplificationError( 0.0 )
	, recomputeFactorEval( true )
	, numVarsAssigned( 0 )
	, factorEval( 0.0 )
	, recomputeBounds( true )
	, cachedBounds( 0.0 )
{
}

Factor::~Factor() {
	for ( VariableData * vd : variableData ) {
		if ( vd != NULL ) delete vd;
	}
	variableData.clear();
}


void Factor::addVariable( Variable * v ) {

	const VariableID vid( v->getID() );

	if ( vid < (VariableCount) variableData.size() && variableData[vid] != NULL ) {
		std::cout << "trying to re-add variable " << v->getID() << " to factor "
				<< id << ": " << *this << std::endl;
		return;
	}

	assert( vid >= (VariableCount) variableData.size() || variableData[vid] == NULL );

	if ( vid >= (VariableCount) variableData.size() ) {
		variableData.resize( vid+1, NULL );
	}

	if ( variableData[vid] == NULL ) {
		variableData[vid] = new VariableData( v, variables.size() );
		variables.push_back( v );
		v->addFactor( this );
	}

	if ( (VariableCount) isVarAssigned.size() <= vid ) {
		isVarAssigned.resize( vid + 1 );
		m_isVarInFactor.resize( vid + 1 );
	}

	m_isVarInFactor.set( vid );

	recomputeFactorEval = true;
	recomputeBounds = true;
}


void Factor::init() {
	isPartlySimplified = false;
	isAssignedConstant = false;
	vidAssigned = -1;
	assignmentError = 0.0;
	simplificationError = 0.0;
	recomputeFactorEval = true;
	numVarsAssigned = 0;
	factorEval = 0.0;
	recomputeBounds = true;
	cachedBounds = 0.0;

	isVarAssigned.reset();
	m_isVarInFactor.reset();

	for ( Variable * v : variables ) {
		assert( variableData[v->getID()] != NULL );
		variableData[v->getID()]->wasAssignedBeforeFactorWasSimplified = true;
		m_isVarInFactor.set( v->getID() );
	}
}


void Factor::clear() {

	variables.clear();
	isVarAssigned.clear();
	m_isVarInFactor.clear();

	for ( VariableData * vd : variableData ) {
		if ( vd != NULL ) delete vd;
	}
	variableData.clear();

	init();
}


Numeric Factor::eval() const {
	assert( isAssignedConstant || areAllVarsAssigned() );
	if ( !isAssignedConstant ) {
		factorEval = evalFactorCached();
		++instrumentation::g_numFactEvals;
	}

	assert( !std::isnan( factorEval ) );
	return factorEval;
}


NumericInterval Factor::computeBounds( const VariableIDVec & vidsToIgnore ) const {

	if ( areAnyInFactor( vidsToIgnore ) && areAnyAssigned( vidsToIgnore ) ) {
		return computeFactorBounds( vidsToIgnore );
	}

	if ( isAssignedConstant || areAllVarsAssigned() ) return eval();

	if ( recomputeBounds ) {
		cachedBounds = computeFactorBounds( vidsToIgnore );
		recomputeBounds = false;
	}

	assert( !std::isnan( cachedBounds.lower() ) &&
			!std::isnan( cachedBounds.upper() ) );

	return cachedBounds;
}


void Factor::computeGradient( PartialGradient & g, bool doGradCheck ) const {
	if ( isAssignedConstant ) return;

	g.reserve( variables.size() );
	for ( const auto & v : variables ) {
		g[v->getID()] = getDerivative( v->getID() );
	}

	if ( doGradCheck ) checkGradient( g );
}


void Factor::onVarAssigned( VariableID varID, Numeric /*newVal*/ ) {
	assert( !isVarAssigned.test( varID ));
	if ( !isVarAssigned.test( varID ) ) {
		++numVarsAssigned;
		isVarAssigned.set( varID );
		recomputeFactorEval = true;
		recomputeBounds = true;
	}
}


void Factor::onVarChanged( VariableID varID, Numeric /*oldVal*/,
		Numeric /*newVal*/ ) {
	assert( isVarAssigned.test( varID ) );
	recomputeFactorEval = true;
	recomputeBounds = true;
}


void Factor::onVarUnassigned( VariableID varID, Numeric /*oldVal*/ ) {
	assert( isVarAssigned.test( varID ) );
	if ( isVarAssigned.test( varID ) ) {
		isVarAssigned.reset( varID );
		recomputeFactorEval = true;
		recomputeBounds = true;
		--numVarsAssigned;
	}
}


void Factor::onVarDomainChanged( VariableID vid ) {
	assert( !isVarAssigned.test( vid ) );
	recomputeFactorEval = true;
	recomputeBounds = true;
}


void Factor::checkGradient( const PartialGradient & g, const Numeric h ) const {

//	const Numeric h = 1e-8;
	Numeric origeval = 0.0;
	PartialGradient fd;
	fd.reserve( g.size() );

	for ( Variable * v : variables ) {
		origeval = v->eval();

		v->assign( v->eval() - h/2.0 );
		Numeric fl = evalFactor();

		// TODO: GET RID OF THIS BECAUSE IT DOESN'T CALL OPTFUNC.ONVARASSIGNED()
		// 		 SO IS WRONG FOR CERTAIN FUNCTIONS THAT RELY ON THAT

		v->assign( v->eval() + h );
		Numeric fh = evalFactor();

		// restore v
		v->assign( origeval );

		fd[v->getID()] = ( fh - fl ) / h;
	}

	assert( fd.size() == g.size() );

	bool diff = false;
	for ( Variable * v : variables ) {
		assert( g.find( v->getID() ) != g.end() );
		if ( !approxeq( fd.at( v->getID() ), g.at( v->getID() ), 1e-4 ) ) {
			diff = true;

			std::cout << "fac " << id << " var " << v->getID() << ": df/dv = " <<
					g.at( v->getID() ) << ", fd = " << fd.at( v->getID() ) <<
					" -- diff = " << fd.at( v->getID() ) - g.at( v->getID() ) <<
					std::endl;
		}
	}

//	assert( !diff );
}


VariableCount Factor::numLeftUnassigned( const VariableIDVec & toassign ) const {
	if ( areAllVarsAssigned() ) return 0;

	VariableCount nva = numVarsAssigned;
	for ( VariableID vid : toassign ) {
		if ( isVarInFactor( vid ) && isVarAssigned.test( vid ) ) ++nva;
	}

	assert( nva <= (VariableCount) variables.size() );
	return ( variables.size() - nva );
}


bool Factor::updateSimplification( Numeric epsilon,
		const VariableIDVec & simplificationKeys, ConnectivityGraph * cg ) {

	assert( !simplificationKeys.empty() );
	assert( simplificationKeys.front() >= 0 || simplificationKeys.size() == 1 );

	// need this to know which var(s) to ignore
	assert( simplificationKeys.front() < 0 || areAnyInFactor( simplificationKeys ) );

	if ( isAssignedConstant && vidAssigned != simplificationKeys.front() ) {
		return false;
	} else if ( areAllVarsAssigned() && !areAnyInFactor( simplificationKeys) ) {
		return false;
	}

	using namespace BOOSTNS::numeric;
	bool result = false;
	const bool wasAssigned = isAssignedConstant;

	NumericInterval bounds = computeBounds( simplificationKeys );

//	std::cout << "bounds of " << id << ": " << bounds << " (simpl keys: ";
//	printContainer( simplificationKeys ); std::cout << std::endl;

	// if the error is less than epsilon in each direction, assign this
	// factor -- updating the connectivity graph only if it changed

	if ( width( bounds ) <= 2*epsilon ) {
		std::cout << "assigning factor " << id << " (width " << width(bounds) <<
				", eps " << epsilon << ") to " << median(bounds) << " (" <<
				bounds << ") at vid ";
		printContainer( simplificationKeys );
		std::cout << std::endl;

		assign( median( bounds ), simplificationKeys.front(), cg );
		assignmentError = width( bounds ) / 2.0;
		result = true;

	} else {
		if ( wasAssigned ) {
			std::cout << "unassigning factor " << id << " at vid ";
			printContainer( simplificationKeys ); std::cout << std::endl;

			unassign( simplificationKeys.front(), cg );
			assignmentError = 0.0;
		}

		result = updatePartialSimplification( epsilon, simplificationKeys,
				bounds, cg );
	}

	return result;
}


void Factor::unsimplify( const VariableIDVec & simplKeys,
		ConnectivityGraph * cg ) {
	assert( isAssignedConstant || isPartlySimplified );
	assert( !isAssignedConstant || vidAssigned == simplKeys.front() );

	if ( isAssignedConstant && vidAssigned == simplKeys.front() ) {
		std::cout << "unassigning factor " << id << " at vid ";
		printContainer( simplKeys ); std::cout << std::endl;

		unassign( simplKeys.front(), cg );
		assignmentError = 0.0;
	}

	if ( isPartlySimplified ) undoPartialSimplification( simplKeys, cg );
}


void Factor::assign( Numeric fval, VariableID assignmentKey,
		ConnectivityGraph * cg ) {

	assert( !isAssignedConstant || assignmentKey == vidAssigned );

	factorEval = fval;
	vidAssigned = assignmentKey;

	assert( !std::isnan( fval ) );

	if ( !isAssignedConstant && cg != NULL ) cg->onFactorAssigned( this );

	isAssignedConstant = true;
}


void Factor::unassign( VariableID assignmentKey, ConnectivityGraph * cg ) {
	assert( isAssigned() );
	assert( vidAssigned == assignmentKey );

	isAssignedConstant = false;
	vidAssigned = -1;
	factorEval = 0;

	if ( cg != NULL ) cg->onFactorUnassigned( this );
}


void Factor::print( std::ostream & os ) const {
	os << "[ ";
	for ( VariablePtrVec::const_iterator it( variables.begin() );
			it != variables.end(); ++it ) {
		os << ( it != variables.begin() ? " + " : "" ) << (*it)->getName();
	}
	os << " ]";
}


std::ostream & operator<<( std::ostream & out, Factor const & f ) {
	f.print( out );
	return out;
}

} // namespace rdis
