/*
 * SimpleProductFactor.cpp
 *
 *  Created on: May 8, 2013
 *      Author: afriesen
 */

#include "NonlinearProductFactor.h"
#include "util/numeric.h"

namespace rdis {

NonlinearProductFactor::NonlinearProductFactor( FactorID factorID,
		Numeric coefficient, /*bool useSigmoid_,*/ bool useExponential_,
		VariableCount totalNVars )
	: super( factorID )
	, coefficient( coefficient )
//	, useSigmoid( useSigmoid_ )
	, useExponential( useExponential_ )
{
	if ( totalNVars > 0 ) variableData.reserve( totalNVars );
}

NonlinearProductFactor::~NonlinearProductFactor() {}


void NonlinearProductFactor::addVariable( Variable * vp, Numeric exponent,
		Numeric constant, bool useSine ) {

	// don't bother adding a var with exponent zero
	if ( exponent == 0.0 ) return;

	VariableID vid( vp->getID() );

	assert( vid >= (VariableCount) variableData.size() || variableData[vid] != NULL );

	if ( vid >= (VariableCount) variableData.size() ) {
		variableData.resize( vid+1, NULL );
	}

	if ( variableData[vid] == NULL ) {
		variableData[vid] = new SimpleProductVariableData( vp, variables.size(),
				exponent, constant, useSine );

		variables.push_back( vp );

		vp->addFactor( this );
	}

	if ( (VariableCount) isVarAssigned.size() <= vid ) {
		isVarAssigned.resize( vid + 1 );
		m_isVarInFactor.resize( vid + 1 );
	}
}


void NonlinearProductFactor::computeGradient( PartialGradient & g,
		bool doGradCheck ) const {

	g.clear();
	g.reserve( variables.size() );

//	Numeric feval( evalFactor( 0.0 ) );
//
//	Numeric k = 1.0;
//	if ( useSigmoid ) {
////		Numeric fe = 1.0 / ( 1.0 + std::exp( -feval ) );
////		k = ( 1 - fe ) * fe;
//		// the old switcheroo
//		feval /= coefficient;
//		k = ( 1.0 - feval ) * feval * coefficient;
//		feval = evalFactor( false ) / coefficient;
//	}
//
//	std::cout << "init feval " << feval << " at " << std::endl;
//	for ( const Variable * v : variables ) {
//		std::cout << v->getID() << ", " << v->eval() << "; ";
//	}
//	std::cout << std::endl;
//
//	for ( const Variable * v : variables ) {
//		const SimpleProductVariableData & vfd = getVarData( v->getID() );
//		Numeric inner = ( vfd.vp->eval() - vfd.constant );
//		std::cout << "var " << v->getID() << " eval " << v->eval() <<std::endl;
//		if ( !vfd.takeSine() ) {
//			if ( !approxeq( inner, 0.0, 1e-20 ) ) {
//				g[v->getID()] = k * feval * vfd.exponent / inner;
//			} else {
//				g[v->getID()] = 0.0;
//			}
//		} else {
//			Numeric innerexp = vfd.hasExp() ? pow( inner, vfd.exponent ) : inner;
//			Numeric siniex = sin( innerexp );
//			Numeric cosiex = cos( innerexp );
//			Numeric k2 = vfd.hasExp() ? pow( inner, vfd.exponent - 1.0 ) : 1.0;
//
//			std::cout << "var " << v->getID() << ": " << inner << ", " <<
//					innerexp << ", " << siniex << ", " << cosiex << std::endl;
//
//			if ( !approxeq( siniex, 0.0, 1e-20 ) ) {
//				g[v->getID()] = k * feval * vfd.exponent * cosiex / siniex *
//						vfd.exponent * k2;
//			} else {
//				g[v->getID()] = 0.0;
//			}
//		}
//	}

//	assert( !useSigmoid );
	assert( !useExponential );

	for ( const Variable * v : variables ) {
		g[v->getID()] = getDerivative( v->getID() );
	}

	if ( doGradCheck ) checkGradient( g );
}


NumericInterval NonlinearProductFactor::computeFactorBounds(
		const VariableIDVec & vidsToIgnore ) const {

	NumericInterval feval( 1 );
	for ( const Variable * v : variables ) {
		bool toignore = ( std::find( vidsToIgnore.begin(),
				vidsToIgnore.end(), v->getID() ) != vidsToIgnore.end() );
		NumericInterval val( v->isAssigned() && !toignore ?
				v->eval() : v->getDomain().interval() );

		const SimpleProductVariableData & vfd = getVarData( v->getID() );

		if ( vfd.hasConstant() ) val -= vfd.constant;
		if ( vfd.hasExp() ) val = power( val, vfd.exponent );
		if ( vfd.takeSine() ) val = BOOSTNS::numeric::sin( val );

		feval *= val;
	}

//	if ( useSigmoid ) feval = 1.0 / ( 1.0 + BOOSTNS::numeric::exp( -feval ) );
	if ( useExponential ) feval = BOOSTNS::numeric::exp( -feval );

	feval *= coefficient;

	return feval;
}



Numeric NonlinearProductFactor::getDerivative( const VariableID & vid ) const {
	assert( areAllVarsAssigned() );
	Numeric varEval( 1 );

//	assert( !useSigmoid );

	for ( Variable * v : variables ) {
		const SimpleProductVariableData & vfd = getVarData( v->getID() );
		Numeric val( v->eval() );

		if ( v->getID() == vid ) {
			if ( !vfd.hasExp() && !vfd.takeSine() ) continue;
			val -= vfd.constant;
			Numeric inner = val;
			Numeric innerexp = power( inner, vfd.exponent );
			val = power( val, vfd.exponent - 1.0 );
			val *= vfd.exponent;
			if ( vfd.takeSine() ) val *= cos( innerexp );
			varEval *= val;

		} else {
			if ( vfd.hasConstant() ) val -= vfd.constant;
			if ( vfd.hasExp() ) val = power( val, vfd.exponent );
			if ( vfd.takeSine() ) val = sin( val );
			varEval *= val;
		}
	}

	return varEval * coefficient;
}


Numeric NonlinearProductFactor::evalFactor() const {
	return evalFactor( /*useSigmoid,*/ useExponential );
}


Numeric NonlinearProductFactor::evalFactor( /*bool usesigm,*/ bool useexp ) const {
	assert( areAllVarsAssigned() );
//	if ( usesigm && useexp ) {
//		throw "NonlinearProductFactor doesn't support both a sigmoid and exponential"
//				" at the same time";
//	}

	Numeric varEval( 1 );
	for ( const Variable * v : variables ) {
		Numeric val = v->eval();
		const SimpleProductVariableData & vfd = getVarData( v->getID() );
		if ( vfd.hasConstant() ) val -= vfd.constant;
		if ( vfd.hasExp() ) val = power( val, vfd.exponent );
		if ( vfd.takeSine() ) val = sin( val );
		varEval *= val;
	}

	Numeric feval = varEval;
//	if ( usesigm ) feval = 1.0 / ( 1.0 + std::exp( -feval ) );
	if ( useexp ) feval = std::exp( -feval );
	feval *= coefficient;

	return feval;
}


void NonlinearProductFactor::print( std::ostream & os ) const {
	const bool printMatlab( false );

	if ( coefficient != 1 ) {
		os << "(" << coefficient << ")";
		if ( !variables.empty() ) os << " * ";
	}
//	if ( useSigmoid ) os << "sigm";
	if ( useExponential ) os << "exp";
	if ( !variables.empty() ) os << "( ";
	if ( useExponential ) os << "-";

	for ( VariablePtrVec::const_iterator it( variables.begin() );
			it != variables.end(); ++ it ) {

		os << ( it != variables.begin() ? " * " : "" );
		const SimpleProductVariableData & vfd( getVarData( (*it)->getID() ) );

		if ( vfd.takeSine() ) os << "sin";
		if ( vfd.hasConstant() || vfd.takeSine() ) os << "(";

		if ( printMatlab ) os << "x(" << (*it)->getID()+1 << ")";
		else os << (*it)->getName();

		if ( vfd.hasExp() ) os << "^" << vfd.exponent << "";
		if ( vfd.hasConstant() ) os << " - " << vfd.constant;
		if ( vfd.hasConstant() || vfd.takeSine() ) os << ")";
	}

	if ( !variables.empty() ) os << " )";
	if ( printMatlab ) os << " + ...";
}

} // namespace rdis
