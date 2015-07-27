/*
 * ComponentCache.cpp
 *
 *  Created on: May 19, 2013
 *      Author: afriesen
 */

#include "ComponentCache.h"

#include BOOSTPATH/functional/hash.hpp>


namespace rdis {

std::size_t hash_value( const CacheKey & ck ) {
	const uint32_t separator( 0xffffffff );
	if ( !ck.storedHashValid ) {
		std::size_t seed( 17 );
		BOOSTNS::hash_range( seed, ck.factors.begin(), ck.factors.end() );
		BOOSTNS::hash_combine( seed, separator );
		BOOSTNS::hash_range( seed, ck.unassignedVars.begin(),
				ck.unassignedVars.end() );

		BOOSTNS::hash_combine( seed, separator );
		for ( const auto & va : ck.assignedVars ) {
			BOOSTNS::hash_combine( seed, va.v->getID() );
		}

		ck.storedHash = seed;
		ck.storedHashValid = true;
	}

	return ck.storedHash;
}


bool operator==( const CacheKey & c1, const CacheKey & c2 ) {
	if ( c1.factors.size() != c2.factors.size() ) return false;
	if ( c1.unassignedVars.size() != c2.unassignedVars.size() ) return false;
	if ( c1.assignedVars.size() != c2.assignedVars.size() ) {
		return false;
	}

	for ( auto it1( c1.factors.begin() ), end1( c1.factors.end() ),
			it2( c2.factors.begin() ); it1 != end1; ++it1, ++it2 ) {
		if ( *it1 != *it2 ) return false;
	}

	for ( auto it1( c1.unassignedVars.begin() ),
			end1( c1.unassignedVars.end() ), it2( c2.unassignedVars.begin() );
			it1 != end1; ++it1, ++it2 ) {
		if ( *it1 != *it2 ) return false;
	}

	for ( auto it1( c1.assignedVars.begin() ), end1( c1.assignedVars.end() ),
			it2( c2.assignedVars.begin() ); it1 != end1; ++it1, ++it2 ) {
		if ( it1->v->getID() != it2->v->getID() ) return false;
	}

	return true;
}




// -------- main access point of cache --------

ComponentCache::ComponentCache()
		: printDebug( false ) {
}

ComponentCache::~ComponentCache() {}



void ComponentCache::init( const size_t /* size */ ) {}


void ComponentCache::clear() {}


void ComponentCache::clear( const size_t /* size */ ) {
	clear();
}


void ComponentCache::insert( const CacheKey & ck, const SubdomainCSP & sd ) {

	assert( !ck.factors.empty() );
	assert( !ck.unassignedVars.empty() );

	if ( printDebug ) {
		std::cout << "inserting " << *sd << " into cache at " << ck << std::endl;
	}

	// perform insertion into cache here
	assert( false ); // not currently supported
}


bool ComponentCache::find( const CacheKey & ck, SubdomainCSP & sd ) const {

	assert( !ck.factors.empty() );
	assert( !ck.unassignedVars.empty() );

	if ( printDebug ) std::cout << "looking in cache for key: " << ck << " -- (";

	bool result = false;

	// perform cache lookup here
	assert( false ); // not currently supported

	if ( printDebug ) {
		if ( result ) std::cout << "found: " << *sd << ")" << std::endl;
		else std::cout << "not found)" << std::endl;
	}

	return result;
}


void ComponentCache::remove( const CacheKey & ck ) {

	assert( false ); // not currently supported

	assert( !ck.factors.empty() );
	assert( !ck.unassignedVars.empty() );

	if ( printDebug ) std::cout << "removing from cache at " << ck << std::endl;
}


std::ostream & operator<<( std::ostream & os, const CacheKey & ck ) {
	assert( !ck.factors.empty() );
	assert( !ck.unassignedVars.empty() );

	os << "f: {";

	for ( auto it( ck.factors.begin() ); it != ck.factors.end(); ++it ) {
		os << (it == ck.factors.begin() ? "" : ", " ) << *it;
	}

	os << "}, uv: {";

	for ( auto it( ck.unassignedVars.begin() ); it != ck.unassignedVars.end(); ++it ) {
		os << ( it == ck.unassignedVars.begin() ? "" : ", " ) << *it;
	}

	os << "}, av: {";

	for ( auto it( ck.assignedVars.begin() ); it != ck.assignedVars.end(); ++it ) {
		os << (it == ck.assignedVars.begin() ? "" : ", " ) << "("
				<< it->v->getID() << "=" << it->v->eval() << ")";
	}

	os << "}";
	return os;
}

} // namespace rdis
