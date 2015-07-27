/*
 * Factor.h
 *
 *  Created on: Jan 12, 2013
 *      Author: afriesen
 */

#ifndef RDIS_FACTOR_H_
#define RDIS_FACTOR_H_

#include "common.h"
#include "State.h"
#include "Variable.h"

#include BOOSTPATH/dynamic_bitset.hpp>
#include BOOSTPATH/container/vector.hpp>

namespace rdis
{

// forward declaration to deal with circular dependency
class Variable;
class ConnectivityGraph;

class Factor {

protected:
	struct VariableData {
		VariableData()
			: vp( NULL )
			, index( -1 )
			, wasAssignedBeforeFactorWasSimplified( true )
		{}
		VariableData( Variable * vp_, size_t index_ )
			: vp( vp_ )
			, index( index_ )
			, wasAssignedBeforeFactorWasSimplified( true )
		{}
		virtual ~VariableData() {}

		// pointer to the variable itself
		Variable * vp;

		// index into the variables vector
		size_t index;

		// true if this variable was assigned before the factor was simplified
		bool wasAssignedBeforeFactorWasSimplified;
	};

	typedef BOOSTNS::container::vector< VariableData * > VariableDataVec;

public:
	// constructor that takes the unique ID of this factor
	Factor( FactorID factorID );
	Factor( const Factor & f ); // make sure no copy constructor is generated
	virtual ~Factor();

	// add the specified variable to the list of variables in this factor
	virtual void addVariable( Variable * vp );

	// initialize this factor (call this before using the factor)
	virtual void init();

	// clear everything in this factor (all var info, etc)
	virtual void clear();

public:

	// evaluate this factor -- all variables must currently be assigned values
	// or this factor must have been assigned a constant value
	virtual Numeric eval() const;

	// compute the bounds for this factor -- calls computeFactorBounds() after
	// checking cached value (derived classes should override
	// computeFactorBounds()). Can ignore the fact that some variables in this
	// factor are currently assigned (useful for determining if this factor
	// should remain simplified)
	virtual NumericInterval computeBounds(
			const VariableIDVec & vidsToIgnore = VariableIDVec() ) const;


	// evaluate this factor without altering anything inside this function
	// (especially w.r.t. cached values)
	virtual Numeric evalNoCache() const {
		assert( false );
		return 0;
	}

public:
	// puts the current partial derivatives into PartialGradient
	virtual void computeGradient( PartialGradient & g,
			bool checkGradient = false ) const;


protected:
	// compute dF/dv (partial derivative of this factor w.r.t. this var)
	virtual Numeric getDerivative( const VariableID & /*vid*/ ) const {
		assert( false );
		return 0;
	}

public:
	// called when a variable is assigned/changed/unassigned a new value
	virtual void onVarAssigned( VariableID vid, Numeric newVal);
	virtual void onVarChanged( VariableID vid, Numeric oldVal, Numeric newVal );
	virtual void onVarUnassigned( VariableID vid, Numeric oldVal );

	// called when a variable's domain is changed
	virtual void onVarDomainChanged( VariableID vid );

protected:
	virtual void checkGradient( const PartialGradient & g,
			const Numeric h = 1e-8 ) const;

public:
	// return the unique ID of this factor
	FactorID getID() const { return id; }

	// return a const reference to the list of variables
	const std::vector< Variable * > & getVariables() const {
		return variables;
	}

	bool isSimplified() const {
		return ( isPartlySimplified || isAssignedConstant );
	}

	// return true if this factor has been partially / fully simplified
	bool isPartiallySimplified() const { return isPartlySimplified; }

	// check if this factor has been assigned a constant value
	bool isAssigned() const { return isAssignedConstant; }

	// return the key (variable ID) for when the full assignment happened
	VariableID getAssignedKey() const { return vidAssigned; }


	// return the current error caused by simplification of this factor
	Numeric getSimplificationError() const {
		assert( isAssignedConstant || ( assignmentError == 0.0 ) );
		assert( isPartlySimplified || approxeq( simplificationError, 0.0, 1e-8 ) );
		Numeric error( 0.0 );
		if ( isAssignedConstant ) error = assignmentError;
		else if ( isPartlySimplified ) error = simplificationError;
		return error;
	}


	 // check if this factor's variables are all assigned values
	bool areAllVarsAssigned() const {
		assert( numVarsAssigned <= (VariableCount) variables.size() );
		return ( numVarsAssigned == (VariableCount) variables.size() );
	}

	// number of vars left unassigned after assigning those in <toassign>
	VariableCount numLeftUnassigned( const VariableIDVec & toassign ) const;


	// true if this variable is currently considered part of this factor (based
	// on the existing assignments and simplifications)
	virtual bool containsVar( VariableID vid ) const {
		assert( !isPartlySimplified ); // not supported here
		return ( !isAssignedConstant && isVarInFactor( vid ) );
	}

	bool isVarInFactor( VariableID vid ) const {
		return ( vid >= 0 && vid < (VariableID) m_isVarInFactor.size() &&
				m_isVarInFactor[vid] );
	}

	bool areAnyInFactor( const VariableIDVec & vids ) const {
		for ( VariableID vid : vids ) {
			if ( isVarInFactor( vid ) ) return true;
		}
		return false;
	}

	bool wasVarAssignedBeforeSimpl( VariableID vid ) const {
		const VariableData & vd = getVarData( vid );
		return vd.wasAssignedBeforeFactorWasSimplified;
	}

	// return the number of variables in this factor
	size_t numVars() const { return variables.size(); }

	// return the number of (un)assigned variables in this factor
	size_t numAssignedVars() const { return numVarsAssigned; }
	size_t numUnassignedVars() const { return ( numVars() - numVarsAssigned ); }

	// overload the implicit conversion operator
	operator FactorID() const { return id; }

	// print this factor to an output stream
	virtual void print( std::ostream & os ) const;
	friend std::ostream & operator<<( std::ostream & out, Factor const & f );

public:
	// update the simplification of this factor -- return true if
	// unsimplification will be needed later with this simplificationKey
	virtual bool updateSimplification( Numeric epsilon,
			const VariableIDVec & simplificationKeys, ConnectivityGraph * cg );

	virtual void unsimplify( const VariableIDVec & simplKeys,
			ConnectivityGraph * cg );

//protected:
	// assign this factor to a constant value, despite its variable's values
	virtual void assign( Numeric fval, VariableID assignmentKey,
			ConnectivityGraph * cg = NULL );
	virtual void unassign( VariableID assignmentKey,
			ConnectivityGraph * cg = NULL );
	// NOTE: ASSIGN/UNASSIGN ARE DEPRECATED -- DON'T CALL THEM DIRECTLY -- USE
	// SIMPLIFICATION FUNCTIONS ABOVE

protected:
	virtual bool updatePartialSimplification( Numeric /*epsilon*/,
			const VariableIDVec & /*simplificationKey*/,
			NumericInterval /*factorBounds*/, ConnectivityGraph * /*cg*/ ) {
		return false;
	}
	virtual void undoPartialSimplification(
			const VariableIDVec & /*simplificationKeys*/,
			ConnectivityGraph * /*cg*/ ) {}


protected:
	virtual Numeric evalFactorCached() const {
		if ( recomputeFactorEval ) {
			factorEval = evalFactor();
			recomputeFactorEval = false;
		}
		return factorEval;
	}

	// evaluate the factor (assumes all vars are assigned)
	virtual Numeric evalFactor() const = 0;

	// returns the current bounds of this factor (based on the current assigned
	// values of the variables) -- if vidsToIgnore is non-empty, ignore any
	// assignments of the specified variables
	virtual NumericInterval computeFactorBounds(
			const VariableIDVec & /*vidsToIgnore*/ = VariableIDVec() ) const {
		assert( false );
		Numeric maxv( std::numeric_limits< Numeric >::max() );
		return NumericInterval( -maxv, maxv );
	}


protected:
	virtual VariableData & getVarData( VariableID vid ) {
		assert( variableData[vid] != NULL );
		return *variableData[vid];
	}

	virtual const VariableData & getVarData( VariableID vid ) const {
		assert( variableData[vid] != NULL );
		return *variableData[vid];
	}

	bool areAnyAssigned( const VariableIDVec & vids ) const {
		for ( VariableID vid : vids ) {
			if ( vid < (VariableID) isVarAssigned.size() &&
				 isVarAssigned.test( vid )) {
				return true;
			}
		}
		return false;
	}


	inline friend bool operator==( const Factor & f1, const Factor & f2 ) {
		return ( f1.id == f2.id );
	}

	inline friend bool operator<( const Factor & f1, const Factor & f2 ) {
		return ( f1.id < f2.id );
	}


protected:
	// the unique ID of this factor
	FactorID id;

	// a list of the variables in this factor (not indexed by vid)
	std::vector< Variable * > variables;

	// a vector with the data for the variables in this factor (indexed by vid)
	VariableDataVec variableData;

	// true if this factor is currently partially simplified in some way
	bool isPartlySimplified;

	// true when this factor has been set to a constant value
	bool isAssignedConstant;

	// the variable that was being processed when this factor was assigned
	VariableID vidAssigned;

	// the possible error created by assigning / simplifying this factor
	Numeric assignmentError;
	Numeric simplificationError;

	// true to cause a recomputation of the cached factor eval
	mutable bool recomputeFactorEval;

	// the number of variables
	VariableCount numVarsAssigned;

	// computed full evaluation of this factor (only valid when isAssigned() is
	// true) -- mutable because it's a cached value
	mutable Numeric factorEval;

	// cache the current bounds
	mutable bool recomputeBounds;
	mutable NumericInterval cachedBounds;

	// flags indicating whether the corresponding variable is in this factor
	// (indexed by vid)
	BOOSTNS::dynamic_bitset<> isVarAssigned;

	// flag indicating whether the corresponding variable is in this factor
	// (indexed by vid)
	BOOSTNS::dynamic_bitset<> m_isVarInFactor;

}; // class Factor


typedef std::vector< Factor * > FactorPtrVec;
typedef BOOSTNS::container::flat_set< Factor * > SortedFactorPtrs;


struct FactorPtrComparator {
	bool operator()( const Factor * lhs, const Factor * rhs ) {
		return ( lhs->getID() < rhs->getID() );
	}
};


} // namespace rdis

#endif // RDIS_FACTOR_H_
