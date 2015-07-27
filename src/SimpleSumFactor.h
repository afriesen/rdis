/*
 * PolynomialFactor.h
 *
 *  Created on: Apr 6, 2013
 *      Author: afriesen
 */

#ifndef RDIS_POLYNOMIALFACTOR_H_
#define RDIS_POLYNOMIALFACTOR_H_

#include "Factor.h"

namespace rdis {

class SimpleSumFactor
		: public Factor {

protected:
	typedef rdis::Factor super;

	enum PolynomialFactorFlags {
		HAS_CONST = 0,
		HAS_EXP = 1,
		HAS_COEFF = 2,
		NUMFLAGS = 3 // keep at end of list
	};

	struct PolynomialVariableData
			: public SimpleSumFactor::VariableData {

		typedef SimpleSumFactor::VariableData super;

		PolynomialVariableData()
			: super( NULL, -1 )
			, exponent( -999999 )
			, coeff( -999999 )
		{}

		PolynomialVariableData( Variable * vp_, size_t index_,
				Numeric exponent_, Numeric coeff_ )
			: super( vp_, index_ )
			, exponent( exponent_ )
			, coeff( coeff_ )
		{}

		virtual ~PolynomialVariableData() {}

		bool hasExp() const { return exponent != 1; }
		bool hasCoeff() const { return coeff != 1; }
		Numeric exponent;
		Numeric coeff;
	};

public:
	SimpleSumFactor( FactorID factorID, Numeric constant = 0,
			Numeric exponent = 1, Numeric coefficient = 1 );
	virtual ~SimpleSumFactor();


	// don't support this override
	void addVariable( Variable * /*vp*/ ) { assert( false ); }

	void addVariable( Variable * vp, Numeric exponent = 1,
			Numeric coefficient = 1 );

	void init();

public:
	void setCoeff( Numeric coeff ) {
		this->coefficient = coeff;
		flags[HAS_COEFF] = ( coeff != 1 );
	}

	void setConstant( Numeric constant ) {
		this->constant = constant;
		flags[HAS_CONST] = ( constant != 0 );
	}

	void setExponent( Numeric exp ) {
		this->exponent = exp;
		flags[HAS_EXP] = ( exp != 1 );
	}

protected:
	virtual NumericInterval computeFactorBounds(
			const VariableIDVec & vidsToIgnore ) const;

	virtual Numeric getDerivative( const VariableID & vid ) const;

protected:
	// combine each variable with its coefficient and exponent
	Numeric evalVariable( VariableID varID, Numeric val ) const;
	NumericInterval evalVariable( VariableID varID,
			NumericInterval interval ) const;

	// evaluate the factor given the evaluation of the internal variables
	Numeric evalFactor() const;
	NumericInterval evalFactor( NumericInterval interval ) const;

	// evaluate this factor without altering anything inside this function
	// (especially w.r.t. cached values)
	virtual Numeric evalNoCache() const { return evalFactor(); }

public:
	// called when a variable is assigned/changed/unassigned a new value
	virtual void onVarAssigned( VariableID vid, Numeric newVal);
	virtual void onVarChanged( VariableID vid, Numeric oldVal, Numeric newVal );
	virtual void onVarUnassigned( VariableID vid, Numeric oldVal );

	void print( std::ostream & os ) const;

protected:
	const PolynomialVariableData & getVarData( VariableID vid ) const {
		assert( isVarInFactor( vid ) );
		assert( variableData[vid] != NULL );
		return (const PolynomialVariableData &) *variableData[vid];
	}

protected:
	// the constant, exponent, and coefficient of this factor
	Numeric constant;
	Numeric exponent;
	Numeric coefficient;

	// current evaluation of those vars that have been assigned
	Numeric partialVarEval;

	// flags indicating the components of this factor (has a coeff, has a
	// constant, has an exponent)
	BOOSTNS::dynamic_bitset<> flags;

};

typedef BOOSTNS::shared_ptr<SimpleSumFactor> PolyFactorSP;
typedef std::vector< BOOSTNS::shared_ptr<SimpleSumFactor> > PolyFactorSPVec;
typedef std::vector< SimpleSumFactor * > PolyFactorPtrVec;

} // namespace rdis

#endif // RDIS_POLYNOMIALFACTOR_H_
