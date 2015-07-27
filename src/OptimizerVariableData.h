/*
 * OptimizerVariableData.h
 *
 *  Created on: May 3, 2013
 *      Author: afriesen
 */

#ifndef RDIS_OPTIMIZERVARIABLEDATA_H_
#define RDIS_OPTIMIZERVARIABLEDATA_H_

#include "common.h"
#include "Variable.h"
#include "ExclusionRegion.h"
#include "ComponentCache.h"

#include BOOSTPATH/dynamic_bitset.hpp>
#include BOOSTPATH/container/vector.hpp>

namespace rdis {

class VarData;
class OptimizerVariableData;

typedef std::vector< Numeric > NumericVec;


class VarData {

public:
	// constructor that takes a pointer to the variable, the lispchitz constant
	// of this variable, and the global epsilon to use
	VarData( Variable * v_ );
	~VarData() {}

	// clear the data in here, and prepare for it to be used again
	void reset();

public:
	// pointer to the variable itself
	Variable * v;

	// the index in the assignment stack that this variable is at (tells us the
	// order of assignment of variables)
	VariableCount assignmentIndex;

	// initial value to set this variable to when first assigning it
	Numeric xvalinit;
	bool xvinitIsSet;
};


std::ostream & operator<<( std::ostream & os, const VarData & vd );

typedef BOOSTNS::container::vector< VarData > VarDataVec;
typedef BOOSTNS::container::vector< VarData * > VarDataPtrVec;


class OptimizerVariableData {
public:
	OptimizerVariableData( const VariableCount & numVars );
	~OptimizerVariableData();

	// initialize memory, etc
	void init( const VariablePtrVec & vars );

	// reset (and prepare to use again)
	void reset();

public:
	// return the number of variables
	inline VariableCount numVars() const { return m_numVars; }

	// return the corresponding variable
	inline Variable * var( const VariableID & vid ) {
		return m_vardata[vid].v;
	}

	inline const Variable * var( const VariableID & vid ) const {
		return m_vardata[vid].v;
	}


	// return the corresponding variable data struct
	inline VarData & get( const VariableID & vid ) {
		assert( vid >= 0 && vid < (int) m_vardata.size() );
		return m_vardata[vid];
	}

	inline const VarData & get( const VariableID & vid ) const {
		assert( vid >= 0 && vid < (int) m_vardata.size() );
		return m_vardata[vid];
	}


	// return the count of assigned (unassigned) variables
	inline VariableCount numAssigned() const { return m_numAssigned; }

	inline VariableCount numUnassigned() const {
		return ( m_numVars - m_numAssigned );
	}


	// return true if a var is assigned / unassigned
	inline bool isAssigned( const VariableID & vid ) const {
		return m_isAssigned[vid];
	}

	inline bool isUnassigned( const VariableID & vid ) const {
		return !m_isAssigned[vid];
	}

public:
	void assign( VariableID vid, const Numeric & val );
	void unassign( VariableID vid );

private:
	inline void setAssigned( const VariableID & vid ) {
		m_isAssigned.set( vid );
		++m_numAssigned;
	}

	inline void setUnassigned( const VariableID & vid ) {
		m_isAssigned.reset( vid );
		--m_numAssigned;
	}

private:
	friend class VarData;

private:
	const VariableCount m_numVars;
	VarDataVec m_vardata;

	// number of assigned variables and flags indicating whether each var is
	// currently assigned
	VariableCount m_numAssigned;
	BOOSTNS::dynamic_bitset<> m_isAssigned;
};

} // namespace rdis

#endif // RDIS_OPTIMIZERVARIABLEDATA_H_
