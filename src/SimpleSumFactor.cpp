/*
 * PolynomialFactor.cpp
 *
 *  Created on: Apr 6, 2013
 *      Author: afriesen
 */

#include "SimpleSumFactor.h"
#include "util/numeric.h"

namespace rdis {

SimpleSumFactor::SimpleSumFactor( FactorID factorID, Numeric constant,
		Numeric exponent, Numeric coefficient )
	: super( factorID )
	, constant( constant )
	, exponent( exponent )
	, coefficient( coefficient )
	, flags( NUMFLAGS ) {

	flags[HAS_CONST] = ( constant != 0 );
	flags[HAS_EXP] = ( exponent != 1 );
	flags[HAS_COEFF] = ( coefficient != 1 );
}

SimpleSumFactor::~SimpleSumFactor() {}


void SimpleSumFactor::addVariable( Variable * vp, Numeric exponent,
		Numeric coefficient ) {

	VariableID vid( vp->getID() );

	assert( vid >= (VariableCount) variableData.size() || variableData[vid] != NULL );
	if ( vid >= (VariableCount) variableData.size() ) variableData.resize( vid+1, NULL );

	if ( variableData[vid] == NULL ) {
		variableData[vid] = new PolynomialVariableData( vp, variables.size(),
				exponent, coefficient );
		variables.push_back( vp );
		vp->addFactor( this );
	}

	if ( (VariableCount) isVarAssigned.size() <= vid ) {
		isVarAssigned.resize( vid + 1 );
		m_isVarInFactor.resize( vid + 1 );
	}

	m_isVarInFactor.set( vid );

	recomputeFactorEval = true;
	recomputeBounds = true;
}


void SimpleSumFactor::init() {

	super::init();

	assert( flags[HAS_EXP] || exponent == 1 );
	assert( flags[HAS_COEFF] || coefficient == 1 );
	assert( flags[HAS_CONST] || constant == 0 );
}


NumericInterval SimpleSumFactor::computeFactorBounds(
		const VariableIDVec & vidsToIgnore ) const {
	NumericInterval veval( partialVarEval );

	for ( VariableID vid : vidsToIgnore ) {
		if ( isVarInFactor( vid ) && isVarAssigned[vid] ) {
			const auto & vd = getVarData( vid );
			veval -= evalVariable( vid, vd.vp->eval() );
			veval += evalVariable( vid, vd.vp->getDomain().interval() );
		}
	}

	if ( !areAllVarsAssigned() ) {
		for ( VariablePtrVec::const_iterator it( variables.begin() );
				it != variables.end(); ++it ) {

			VariableID vid( (*it)->getID() );

			if ( !isVarAssigned.test( vid ) ) {
				veval += evalVariable( vid, ( *it )->getDomain().interval() );
			}
		}
	}

	veval = evalFactor( veval );
	return veval;
}


Numeric SimpleSumFactor::getDerivative( const VariableID & vid ) const {
	assert( areAllVarsAssigned() );

	Numeric dfdv = 1.0;

	if ( flags[HAS_COEFF] ) dfdv *= coefficient;
	if ( !approxeq( exponent, 0.0 ) ) {
		Numeric feval = partialVarEval + constant;
		feval = power( feval, exponent - 1.0 );
		dfdv *= exponent * feval;
	}

	const PolynomialVariableData & vd( getVarData( vid ) );

	if ( vd.hasCoeff() ) dfdv *= vd.coeff;
	if ( !approxeq( vd.exponent, 0.0 ) ) {
		Numeric veval = vd.vp->eval();
		veval = power( veval, vd.exponent - 1.0 );
		dfdv *= vd.exponent * veval;
	}

	return dfdv;
}


Numeric SimpleSumFactor::evalVariable( VariableID varID, Numeric val ) const {
	const PolynomialVariableData & vfd( getVarData( varID ) );
	if ( vfd.hasExp() ) val = power( val, vfd.exponent );
	if ( vfd.hasCoeff() ) val *= vfd.coeff;

	return val;
}


NumericInterval SimpleSumFactor::evalVariable( VariableID varID,
		NumericInterval val ) const {

	const PolynomialVariableData & vfd( getVarData( varID ) );
	if ( vfd.hasExp() ) val = power( val, vfd.exponent );
	if ( vfd.hasCoeff() ) val *= vfd.coeff;

	return val;
}


Numeric SimpleSumFactor::evalFactor() const {
	Numeric varEval = partialVarEval;
	if ( flags[HAS_CONST] ) varEval += constant;
	if ( flags[HAS_EXP] ) varEval = power( varEval, exponent );
	if ( flags[HAS_COEFF] ) varEval *= coefficient;
	return varEval;
}


NumericInterval SimpleSumFactor::evalFactor( NumericInterval varEval ) const {
	if ( flags[HAS_CONST] ) varEval += constant;
	if ( flags[HAS_EXP] ) varEval = power( varEval, exponent );
	if ( flags[HAS_COEFF] ) varEval *= coefficient;
	return varEval;
}


void SimpleSumFactor::onVarAssigned( VariableID vid, Numeric newval ) {

	super::onVarAssigned( vid, newval );

//	if ( !isVarAssigned.test( vid ) ) {
		partialVarEval += evalVariable( vid, newval );
		++instrumentation::g_numVarOps;
//	}
}


void SimpleSumFactor::onVarChanged( VariableID varID, Numeric oldVal,
		Numeric newVal ) {

	super::onVarChanged( varID, oldVal, newVal );

	partialVarEval -= evalVariable( varID, oldVal );
	partialVarEval += evalVariable( varID, newVal );
	instrumentation::g_numVarOps += 2;
}


void SimpleSumFactor::onVarUnassigned( VariableID varID, Numeric oldVal ) {
	super::onVarUnassigned( varID, oldVal );

//	if ( isVarAssigned.test( varID ) ) {
		partialVarEval -= evalVariable( varID, oldVal );
		++instrumentation::g_numVarOps;
//	}
}


void SimpleSumFactor::print( std::ostream & os ) const {
	if ( flags[HAS_COEFF] ) os << "(" << coefficient << ")";
	os << "[ ";
	if ( flags[HAS_CONST] ) os << constant << " + ";

	for ( VariablePtrVec::const_iterator it( variables.begin() );
			it != variables.end(); ++ it ) {
		const PolynomialVariableData & vfd( getVarData( (*it)->getID() ) );

		os << ( it != variables.begin() ? " + " : "" );
		if ( vfd.hasCoeff() ) os << "(" << vfd.coeff << ")";
		os << (*it)->getName();
		if ( vfd.hasExp() ) os << "^(" << vfd.exponent << ")";
	}

	os << " ]";
	if ( flags[HAS_EXP] ) os << "^(" << exponent << ")";
}

} // namespace rdis
