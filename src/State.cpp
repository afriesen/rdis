/*
 * State.cpp
 *
 *  Created on: May 5, 2013
 *      Author: afriesen
 */


#include "State.h"

#include <iterator>


namespace rdis {

Subdomain::Subdomain()
	: x( 0 )
	, fx( 0 )
	, ferr( 0 )
	, isBound( false )
{}

Subdomain::Subdomain( const VariableIDVec & vids_ )
	: vids( vids_ )
	, x( 0 )
	, fx( 0.0 )
	, ferr( 0.0 )
	, isBound( false )
{}

Subdomain::Subdomain( const VariableID & vid_ )
	: x( 0 )
	, fx( 0 )
	, ferr( 0 )
	, isBound( false )
{
	assert( vid_ >= 0 );
	vids.push_back( vid_ );
}


Subdomain::~Subdomain() {}


void Subdomain::set( Numeric x_, Numeric fx_, Numeric ferr_, bool isBound_ ) {
	assert( vids.size() == 1 );
	x.clear(); x.push_back( x_ );
	fx = fx_;
	ferr = ferr_;
	isBound = isBound_;
	children.clear();
}


void Subdomain::set( const NumericVec & x_, Numeric fx_, Numeric ferr_,
		bool isBound_ ) {
	assert( vids.size() == x_.size() );
	x = x_;
	fx = fx_;
	ferr = ferr_;
	isBound = isBound_;
	children.clear();
}


void Subdomain::set( Numeric fx_, Numeric ferr_, bool isBound_ ) {
	fx = fx_;
	ferr = ferr_;
	isBound = isBound_;
	children.clear();
}


void Subdomain::setAsCopy( const Subdomain & rhs ) {
	clear();

	assert( rhs.vids.size() == rhs.x.size() );

	vids.insert( vids.end(), rhs.vids.begin(), rhs.vids.end() );
	x.insert( x.end(), rhs.x.begin(), rhs.x.end() );

	set( rhs.fx, rhs.ferr, rhs.isBound );

	children.insert( children.end(), rhs.children.begin(), rhs.children.end() );
}


void Subdomain::clear() {
	vids.clear();
	x.clear();
	fx = std::numeric_limits< Numeric >::max();
	children.clear();
}


void Subdomain::addChild( const SubdomainCSP & csd ) {
	if ( csd->isBound ) isBound = true;
//	children.push_back( csd );
	children.emplace_back( csd );
}


VariableCount copySubdomain( const Subdomain & sd, State & x, bool copyChildren,
		bool checkUnique ) {
	if ( checkUnique ) {
#ifdef DEBUG
		const Numeric maxv( std::numeric_limits< Numeric >::max() );
		for ( VariableID vid : sd.vids ) {
			assert( x[vid] == maxv );
		}
#endif // DEBUG
//		if ( sd.isBound ) {
//			std::cout << "copySubdomain: var " << sd.vid << " is a bound" << std::endl;
//		}
	}

	VariableCount ncopied = sd.vids.size();
	for ( size_t i = 0; i < sd.vids.size(); ++i ) {
		x[sd.vids[i]] = sd.x[i];
	}

	if ( copyChildren ) {
		for ( Subdomain::Children::const_iterator it( sd.children.begin() ),
				end( sd.children.end() ); it != end; ++it ) {
			ncopied += copySubdomain( **it, x, copyChildren, checkUnique );
		}
	}

	return ncopied;
}


void printState( const State & x, std::ostream & os ) {
	os << "(";
	for( size_t i = 0; i < x.size(); ++i ) {
		os << ( i == 0 ? "" : ", " ) << x[i];
	}
	os << ")";
}

void printSubdomain( const Subdomain & sd, std::ostream & os ) {

	std::cout << "(";
	printContainer( sd.vids );
	std::cout << ") : {" << sd.x << ", " << sd.fx << ( sd.isBound ? ", B" : "" ) << "}";

	if ( !sd.children.empty() ) {
		std::cout << " [";
		for ( Subdomain::Children::const_iterator it( sd.children.begin() ),
				end( sd.children.end() ); it != end; ++it ) {
			if ( it != sd.children.begin() ) std::cout << ", ";
			printSubdomain( **it, os );
		}
		std::cout << "] ";
	}
}


std::ostream & operator<<( std::ostream & os, const Subdomain & sd ) {
	printSubdomain( sd );
	return os;
}

std::ostream & operator<<( std::ostream & os, const State & x ) {
	printState( x, os );
	return os;
}

} // namespace rdis

