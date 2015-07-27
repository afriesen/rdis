/*
 * OptimizableFunction.h
 *
 *  Created on: Jan 11, 2013
 *      Author: afriesen
 */

#ifndef RDIS_OPTIMIZABLEFUNCTION_H_
#define RDIS_OPTIMIZABLEFUNCTION_H_

#include "common.h"
#include "State.h"
#include "Variable.h"
#include "Factor.h"
#include "VariableDomain.h"

#include BOOSTPATH/filesystem/fstream.hpp>

#include <vector>

namespace rdis {

class OptimizableFunctionGenerator;

// class representing an optimizable function
class OptimizableFunction {
protected:
	typedef std::vector< string > SplitVector;

public:
	// The semirings for function and bounds evaluation, respectively.
	// (note: these must be defined by all subclasses (see Semiring.h))
	virtual const Semiring< Numeric > & semiring() const = 0;
	virtual const Semiring< NumericInterval > & semiringB() const = 0;

public:
	OptimizableFunction(
			bool hasDerivatives = false,
			bool hasBounds = false,
			bool resetFactorsOnChange = false,
			bool hasBlockedVars = false );

	explicit OptimizableFunction(
			const VariableDomain & defaultDomain,
			bool hasDerivatives = false,
			bool hasBounds = false,
			bool resetFactorsOnChange = false,
			bool hasBlockedVars = false );
//	explicit OptimizableFunction( const string & filename );

	virtual ~OptimizableFunction();

public:
	// initialize all of the factors and variables of this function (call only
	// once, after all factors and variables have been added)
	virtual void init();

	// reset the memberships and connections of all the different factors --
	// call when hyperrects change -- NEVER CALL DURING OPTIMIZATION
	virtual void resetFactors() {}

public:
	// called after a variable has been assigned, in case this function needs to
	// update something internally
	inline virtual void onVarAssigned( const VariableID /*vid*/,
			const Numeric /*val*/ ) const {}

	// called after a variable has been unassigned
	inline virtual void onVarUnassigned( const VariableID /*vid*/ ) const {}

protected:
	// add a variable to this function
	virtual Variable * addVariable( string const & name, string const & domain,
			VariableID & id );

public:
	virtual VariableCount getNumVars() const { return numVars; }
	virtual const VariableDomain & getDefaultDomain() const { return defaultDomain; }
	virtual VariablePtrVec & getVariables()  { return variables; }
	virtual FactorPtrVec const & getFactors() const { return factors; }

public:
	// evaluate the function (or its bounds) at the current values of the
	// variables -- note that all variables must be assigned for eval(), but not
	// for computing bounds
	virtual Numeric eval() const;
	virtual Numeric evalFactors( const FactorPtrVec & fctrs,
			Numeric & ferr, bool useCached = true ) const;

	// evaluate the function at the specified values of the variables (leaves
	// this function unaltered)
	virtual Numeric eval( const State & x );

public:
	// compute the bounds of this function (or the subfunction defined by the
	// specified set of factors) based on the values currently assigned to the
	// variables
	virtual NumericInterval computeBounds() const;
	virtual NumericInterval computeBounds( const FactorPtrVec & fctrs,
			VariableID assignedVID ) const;

	// compute the gradient of this function at the current assigned point (all
	// vars must be assigned)
	virtual void computeGradient( NumericVec & gradient,
			bool checkGrad = false ) const;

	// compute the gradient of this function over the specified set of factors
	// at the currently assigned point
	virtual void computeGradient( const FactorPtrVec & facs,
			PartialGradient & gradient, bool checkGrad = false ) const;

protected:
	// compute the gradient of this function at the current point assuming that
	// it is a sum of terms
	virtual void computeGradientOfSum( const FactorPtrVec & facs,
			PartialGradient & gradient ) const;

	// compute the gradient of this function at the current point assuming that
	// it is a product of factors
	virtual void computeGradientOfProduct( const FactorPtrVec & facs,
			PartialGradient & gradient ) const;

public:
	// perform one iteration of expectation-maximization (if this function
	// supports it) on the subset of variables and factors specified. If
	// sanitizeValues is true, ensure that vals of vars never exceed the var's
	// hypperrect. Returns false on failure.
	// NOTE: this changes the values currently assigned to the variables
	virtual bool EMStep( const VariablePtrVec & /*vars*/,
			const FactorPtrVec & /*facs*/, bool /*sanitizeValues*/,
			bool /*useFactorCaching*/ = true ) {
		std::cout << "EM not supported for target function" << std::endl;
		assert( false );
		return false;
	}

public:
	// true if this function's factors support derivative computation
	inline virtual bool hasDerivatives() const { return m_hasDerivatives; }

	// true if this function's factors support estimation of their upper and
	// lower bounds
	inline virtual bool hasBounds() const { return m_hasBounds; }

	// true if this function needs to have its factors reset whenever the
	// hyperrects / variables / etc change
	inline virtual bool resetFactorsOnChange() const {
		return m_resetFactorsOnChange;
	}

	// true if the variables in this function are grouped into tightly connected
	// blocks that should be optimized together
	inline virtual bool hasBlockedVars() const { return m_hasBlockedVars; }

public:
	// return the number of blocks
	inline virtual VariableCount getNumBlocks() const { return numVars; }

	// returns the range of vars that comprises the block specified by the
	// <blockid> parameter (inclusive) -- if vidlower / vidupper are negative,
	// this means the block ID is invalid
	inline virtual void getBlockRangeByBlkId( VariableCount blockid,
			VariableID & vidlower, VariableID & vidupper ) const {
		vidlower = blockid;
		vidupper = blockid;
	}

	// returns the range of vars that comprises the block containing the
	// specified variable (note that the default is that each var is a block)
	inline virtual void getBlockRangeByVid( VariableID vid,
			VariableID & vidlower, VariableID & vidupper ) const {
		vidlower = vid;
		vidupper = vid;
	}

	// return the ID of the block that contains this variable
	inline virtual VariableCount getBlockID( VariableID vid ) const {
		return vid;
	}

public:
	// load/save this function from/to a file
	virtual bool load( const string &file ) = 0;
	virtual bool save( const string &file, State * state = NULL ) = 0;

	// print a representation of this function to an output stream
	virtual void print( std::ostream & os ) const;

	// write a graphviz representation of this function -- does nothing here
	virtual void writeGraph() const {}

	friend std::ostream & operator<<( std::ostream & out,
			OptimizableFunction const & p );

	friend bool operator==( const OptimizableFunction & lhs,
			const OptimizableFunction & rhs );
	friend bool operator!=( const OptimizableFunction & lhs,
			const OptimizableFunction & rhs );

protected:
	bool readAndSplit( BOOSTNS::filesystem::fstream & ifs, SplitVector & sv );

	virtual void checkGradient( const FactorPtrVec & facs,
			const PartialGradient & g, const Numeric h = 1e-8 ) const;

protected:
	friend class OptimizableFunctionGenerator;

	// the number of variables in this function
	VariableCount numVars;

	// the default domain for each variable
	VariableDomain defaultDomain;

	// flags indicating whether this function (and its factors) support
	// derivative and bound evaluations, or whether the variables are blocked or
	// factors need a reset when the connectivity changes
	const bool m_hasDerivatives;
	const bool m_hasBounds;
	const bool m_resetFactorsOnChange;
	const bool m_hasBlockedVars;

	// vector containing the variables in this polynomial (indexed by their IDs)
	VariablePtrVec variables;

	// vector containing the factors in this polynomial (indexed by their IDs)
	FactorPtrVec factors;

protected: // stored for memory efficiency
	PartialGradient pgtemp1;
};

typedef BOOSTNS::shared_ptr< OptimizableFunction > OptimizableFunctionSP;

typedef std::vector< FactorID > FactorIDVec;
typedef std::vector< FactorIDVec > FactorIDVecVec;

} // namespace rdis

#endif // RDIS_OPTIMIZABLEFUNCTION_H_
