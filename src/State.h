/*
 * State.h
 *
 *  Created on: May 1, 2013
 *      Author: afriesen
 */

#ifndef RDIS_STATE_H_
#define RDIS_STATE_H_

#include "common.h"

#include "IntrusivePtrPool.h"

#include BOOSTPATH/shared_ptr.hpp>
#include BOOSTPATH/intrusive_ptr.hpp>
#include BOOSTPATH/container/flat_map.hpp>
#include BOOSTPATH/container/flat_set.hpp>

#include <map>
#include <vector>
#include <iostream>

namespace rdis {

typedef std::vector< Numeric > State;
//typedef std::map< VariableID, Numeric > PartialState;
//typedef BOOSTNS::container::flat_map< VariableID, Numeric > PartialState;

typedef BOOSTNS::container::flat_map< VariableID, Numeric > PartialGradient;
//typedef std::vector< PartialGradient > PartialGradientVec;

class Subdomain;
typedef IntrusivePtrPool< Subdomain > SubdomainPool;

typedef BOOSTNS::intrusive_ptr< Subdomain > SubdomainSP;
typedef BOOSTNS::intrusive_ptr< const Subdomain > SubdomainCSP;


// tree structure that (recursively) represents a point in a subdomain of the
// entire state space
class Subdomain
	: public IntrusivePtrPoolObj< Subdomain, SubdomainPool > {
public:
	typedef BOOSTNS::container::vector< SubdomainCSP > Children;

	Subdomain();
	Subdomain( const VariableIDVec & vids_ );
	Subdomain( const VariableID & vid_ );
	~Subdomain();

public:
	// set this subdomain's <x, fx> values with no sub-subdomains
	void set( Numeric x_, Numeric fx_, Numeric ferr_, bool isBound_ );
	void set( const NumericVec & x_, Numeric fx_, Numeric ferr_, bool isBound_ );

	// assume x is set elsewhere
	void set( Numeric fx_, Numeric ferr_, bool isBound_ );

	void setAsCopy( const Subdomain & rhs );

	void addChild( const SubdomainCSP & csd );

	void clear();

public:
	inline Subdomain & operator=( const Subdomain & rhs ) {
		if ( this != &rhs ) {
			vids = rhs.vids;
			x = rhs.x;
			fx = rhs.fx;
			ferr = rhs.ferr;
			children = rhs.children;
		}
		return *this;
	}

	friend inline bool operator<( const Subdomain & lhs, const Subdomain & rhs ) {
		return lhs.vids < rhs.vids;
	}

public:
	// the "top" variable(s) of this subdomain
	VariableIDVec vids;

	// the value(s) of the variable(s)
	NumericVec x;

	// the evaluation of this subdomain: f(x, {y_i \in descendants(x)})
	Numeric fx;

	// the possible error in the evaluation of this subdomain
	Numeric ferr;

	// true if the fx in here is the result of a bounding operation, instead of
	// a full eval of the subdomain (must be true if any children have this set)
	bool isBound;

	// the sub-subdomains comprising this subdomain
	Children children;
};


// copy the values in this subdomain into the specified state -- return the
// number of values copied
VariableCount copySubdomain( const Subdomain & sd, State & x, bool copyChildren,
		bool checkUnique = true );

template< class T >
void productGradient( PartialGradient & g1, const PartialGradient & g2,
		const Semiring< T > & sr );


inline State operator-( const State & x1, const State & x2 ) {
	assert( x1.size() == x2.size() );
	State xout( x1.size(), 0 );
	for ( size_t i = 0; i < x1.size(); ++i ) {
		xout[i] = std::abs( x1[i] - x2[i] );
	}
	return xout;
}

inline Numeric distL2( const State & x1, const State & x2 ) {
	assert( x1.size() == x2.size() );
	Numeric dist( 0.0 );
	for ( size_t i = 0; i < x1.size(); ++i ) {
		dist += ( x1[i] - x2[i] ) * ( x1[i] - x2[i] );
	}
	return std::sqrt( dist );
}

inline Numeric distL1( const State & x1, const State & x2 ) {
	assert( x1.size() == x2.size() );
	Numeric dist( 0.0 );
	for ( size_t i = 0; i < x1.size(); ++i ) {
		dist += std::abs( x1[i] - x2[i] );
	}
	return dist;
}

inline bool operator==( const State & x1, const State & x2 ) {
	if ( x1.size() != x2.size() ) return false;
	for ( size_t i = 0; i < x1.size(); ++i ) {
		if ( x1[i] != x2[i] ) return false;
	}
	return true;
}

void printState( const State & x, std::ostream & os = std::cout );
void printSubdomain( const Subdomain & sd, std::ostream & os = std::cout );
std::ostream & operator<<( std::ostream & os, const Subdomain & sd );
std::ostream & operator<<( std::ostream & os, const State & x );



template< class T >
void productGradient( PartialGradient & g1, const PartialGradient & g2,
		const Semiring< T > & sr ) {
	if ( g1.empty() ) {
		g1 = g2;
		return;
	}

	// necessary to keep iterators valid as we insert (and it's more efficient)
	g1.reserve( g1.size() + g2.size() - 1 );

//#ifdef DEBUG
//	size_t g1s( g1.size() ), g2s( g2.size() );
//	PartialGradient u1, u2;
//	std::set_intersection( g1.begin(), g1.end(), g2.begin(), g2.end(), std::inserter( u1, u1.begin() ) );
//	std::set_intersection( g2.begin(), g2.end(), g1.begin(), g1.end(), std::inserter( u2, u2.begin() ) );
//#endif // DEBUG

	auto it1( g1.begin() );
	auto it2( g2.begin() );

	// merge these two PartialGradients, taking the product if there is overlap
	while ( true ) {
		if ( it2 == g2.end() ) break;

		if ( it1 == g1.end() || it1->first > it2->first ) {
			// insert if g2.vid wasn't found in g1
			g1.insert( it1, *it2 );
			++it1; ++it2;
		} else if ( it2->first == it1->first ) {
			// product if the same variable is in both
			it1->second = sr.Product( it1->second, it2->second );
			++it1; ++it2;
		} else {
			// move g1 iterator until it reaches
			assert( it1->first < it2->first );
			++it1;
		}
	}

//#ifdef DEBUG
//	assert( g1.size() < g1s + g2s );
//
//	for ( auto it( u1.begin() ); it != u1.end(); ++it ) {
//		assert( g1[it->first] != it->second );
//	}
//	for ( auto it( u2.begin() ); it != u2.end(); ++it ) {
//		assert( g1[it->first] != it->second );
//	}
//
//	if ( g1.size() <= g1s+g2s-2 || g1.size() > 10 ) {
//		g1s += 0;
//	}
//#endif // DEBUG
}

} // namespace rdis


#endif // RDIS_STATE_H_
