/*
 * ComponentCache.h
 *
 *  Created on: May 19, 2013
 *      Author: afriesen
 */

#ifndef RDIS_COMPONENTCACHE_H_
#define RDIS_COMPONENTCACHE_H_

#include "common.h"
#include "State.h"
#include "Variable.h"

#include BOOSTPATH/tuple/tuple.hpp>
#include BOOSTPATH/unordered_map.hpp>
#include BOOSTPATH/container/deque.hpp>
#include BOOSTPATH/container/flat_set.hpp>

#include BOOSTPATH/geometry.hpp>

#include <functional>
#include <map>
#include <deque>


namespace rdis {

// -------- definitions --------

struct AssignedVarLookup;

struct VariableAssignment {

	VariableAssignment()
		: priority( -1 )
		, v( NULL )
	{}

	VariableAssignment( const VariableCount priority_, const VariableID & vid_,
			const Variable * v_ )
		: priority( priority_ )
		, v( v_ )
	{}

	// allow cast to a variable ID
	inline operator VariableID() { return v->getID(); }

	bool operator==( const VariableAssignment & rhs ) const {
		return ( v->getID() == rhs.v->getID() );
	}

	VariableCount priority;

	// store the variable so we don't have to always update the key each time
	// it gets assigned, it just knows the assigned value from this
	const Variable * v;
};


typedef BOOSTNS::container::flat_set< FactorID > SortedFactorIDs;
typedef BOOSTNS::container::flat_set< VariableID > SortedVariableIDs;
typedef BOOSTNS::container::flat_set< VariableAssignment > SortedVarAssignments;



// the set of factors and set of unassigned vars determines the key of the hash
// table we're using as the first part of the cache
struct CacheKey {

	CacheKey()
		: storedHashValid( false )
		, storedHash( 0 )
	{}

	inline void clear() {
		factors.clear();
		unassignedVars.clear();
		storedHashValid = false;
		storedHash = 0;
		assignedVars.clear();
	}

	SortedFactorIDs factors;
	SortedVariableIDs unassignedVars;

	SortedVarAssignments assignedVars;

	mutable bool storedHashValid;
	mutable size_t storedHash;
};


// useful comparison operations
inline bool operator<( const VariableAssignment & va1,
		const VariableAssignment & va2 ) {
	assert( va1.priority >= 0 && va2.priority >= 0 );
	assert( va1.v->getID() != va2.v->getID() || va1.priority == va2.priority );
	return ( va1.priority < va2.priority );
}


// hash function and equality operator for the cache key
std::size_t hash_value( const CacheKey & k );
bool operator==( const CacheKey & c1, const CacheKey & c2 );



// -------- first part of cache --------
// for a specific component/subfunction (set of factor IDs and set of unassigned
// variables), look up the second part of the cache

// data structure to cache components and their values
class ComponentCache {

public:
	ComponentCache();
	virtual ~ComponentCache();

	// initialize the cache with the specified size
	void init( size_t size );

	void clear();
	void clear( size_t newsize );

public:
	void insert( const CacheKey & ck, const SubdomainCSP & sd );
	bool find( const CacheKey & c, SubdomainCSP & sd ) const;
	void remove( const CacheKey & ck );

private:
	bool printDebug;
};


std::ostream & operator<<( std::ostream & os, const CacheKey & ck );


} // namespace rdis

#endif // RDIS_COMPONENTCACHE_H_
