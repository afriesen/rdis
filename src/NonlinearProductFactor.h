/*
 * SimpleProductFactor.h
 *
 *  Created on: May 8, 2013
 *      Author: afriesen
 */

#ifndef NONLINEARPRODUCTFACTOR_H_
#define NONLINEARPRODUCTFACTOR_H_

#include "Factor.h"

namespace rdis {

// this factor is of the form:
//		fev = c * [ \prod_i( sin( c_i * X_i^e_i ) ) ] if useExponential = false
// 	OR	fev = c * exp{ -1 * [ \prod_i( sin( c_i * X_i^e_i ) ) ] } otherwise
// where c and c_i are (constant) coefficients, X_i are the variables,
// e_i are the exponents, and the sin(...) is optional (could
// just be \prod_i( c_i * X_i^e_i )) -- on a per variable basis
class NonlinearProductFactor
		: public Factor {

	typedef rdis::Factor super;

public:
	struct SimpleProductVariableData
			: public Factor::VariableData {

		typedef Factor::VariableData super;

		SimpleProductVariableData()
			: super( NULL, -1 )
			, exponent( -1 )
			, constant( -1 )
			, useSine( false )
		{}

		SimpleProductVariableData( Variable * vp_, size_t index_,
				Numeric exponent_, Numeric constant_, bool useSine_ )
			: super( vp_, index_ )
			, exponent( exponent_ )
			, constant( constant_ )
			, useSine( useSine_ )
		{}

		virtual ~SimpleProductVariableData() {}

		bool hasExp() const { return exponent != 1; }
		bool hasConstant() const { return constant != 0; }
		bool takeSine() const { return useSine; }
		Numeric exponent;
		Numeric constant;
		bool useSine;
	};


public:
	// construct this NLP factor. totalNVars is an optional param to pre-allocate
	// memory for the specified number of variables
	NonlinearProductFactor( FactorID factorID, Numeric coefficient = 1,
			/*bool useSigmoid = false,*/ bool useExponential = false,
			VariableCount totalNVars = 0 );
	virtual ~NonlinearProductFactor();


	// don't support this override
	void addVariable( Variable * /*vp*/ ) { assert( false ); }

	// add the specified variable (x) to this factor with the specified exponent (e)
	// and constant (k) and whether or not we should take sin(k*x^e) when computing
	// the value of this variable
	void addVariable( Variable * vp, Numeric exponent = 1,
			Numeric constant = 0, bool useSine = false );

public:
	void setCoeff( Numeric coeff ) {
		this->coefficient = coeff;
	}

public:
	void computeGradient( PartialGradient & g, bool doGradCheck ) const;
	Numeric getDerivative( const VariableID & vid ) const;

private:
	NumericInterval computeFactorBounds( const VariableIDVec & vidsToIgnore ) const;

	// evaluate the factor
	Numeric evalFactor() const;
	Numeric evalFactor( /*bool useSigmoid,*/ bool useExponential ) const;

	// evaluate this factor without altering anything inside this function
	// (especially w.r.t. cached values)
	virtual Numeric evalNoCache() const { return evalFactor(); }

public:
	void print( std::ostream & os ) const;

private:
	const SimpleProductVariableData & getVarData( VariableID vid ) const {
		assert( variableData[vid] != NULL );
		return (const SimpleProductVariableData &) *variableData[vid];
	}

private:
	// the coefficient of this factor
	Numeric coefficient;

//	bool useSigmoid;

	// true to take exp( -feval ) after multiplying all of the variable values but
	// before multiplying with the factor coefficient
	bool useExponential;
};

} // namespace rdis

#endif // NONLINEARPRODUCTFACTOR_H_
