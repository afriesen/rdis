/*
 * VariableDomain.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: afriesen
 */

#include "VariableDomain.h"

#include BOOSTPATH/format.hpp>
#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/lexical_cast.hpp>
#include BOOSTPATH/algorithm/string.hpp>

namespace rdis {


bool VariableDomain::Comp::operator()( const NumericInterval & ni,
		const Numeric & val ) const {
	return BOOSTNS::numeric::interval_lib::cerlt( ni, val );
}

bool compareIntervalsLower( const NumericInterval & lhs,
		const NumericInterval & rhs ) {
	return ( lhs.lower() < rhs.lower() );
}


VariableDomain::VariableDomain()
		: realwidthCached( -1.0 )
{}

VariableDomain::VariableDomain( Numeric value )
		: m_interval( value, value )
		, m_samplingInterval( value, value )
		, realwidthCached( -1.0 ) {
	m_subintervals.push_back( m_interval );
}

VariableDomain::VariableDomain( Numeric lower, Numeric upper )
		: m_interval( lower, upper )
		, m_samplingInterval( lower, upper )
		, realwidthCached( -1.0 ) {
	m_subintervals.push_back( m_interval );
}

VariableDomain::VariableDomain( const string & domain )
		: realwidthCached( -1.0 ) {
	parse( domain );
	m_samplingInterval = m_interval;
}

VariableDomain::VariableDomain( const IntervalVec & subitvls )
	: m_subintervals( subitvls.begin(), subitvls.end() ) {
	assert( !subitvls.empty() );

	cleanSubintervals();

	m_samplingInterval = m_interval;
}


void VariableDomain::parse( string const & domain ) {

	using namespace BOOSTNS;
	SplitVector sv;
	split( sv, domain, is_any_of( " ~:," ), token_compress_on );

	assert( sv.size() % 2 == 0 || sv.size() == 3 );
	if ( sv.size() == 2 || sv.size() == 3 ) {
		m_interval.assign( lexical_cast< Numeric >( sv.front() ),
				lexical_cast< Numeric >( sv.back() ) );
		m_subintervals.push_back( m_interval );

	} else {
		m_subintervals.reserve( sv.size() / 2 );

		for ( unsigned int i = 0; i < sv.size()-1; i += 2 ) {
			NumericInterval intvl( lexical_cast< Numeric >( sv[i] ),
					lexical_cast< Numeric >( sv[i+1] ) );

			m_subintervals.push_back( intvl );
		}

		cleanSubintervals();
	}

	assert( !m_subintervals.empty() );
}


void VariableDomain::cleanSubintervals() {
	using namespace BOOSTNS::numeric;

	// make sure the subintervals are ordered
	std::sort( m_subintervals.begin(), m_subintervals.end(),
			compareIntervalsLower );

	for ( size_t i = 1; i < m_subintervals.size(); ++i ) {
		NumericInterval & psi = m_subintervals[i - 1];
		NumericInterval & si = m_subintervals[i];

		if ( overlap( psi, si ) ) {
//			std::cout << i << " replacing " << psi << " and " << si << " with ";
			psi.assign( psi.lower(), std::max( psi.upper(), si.upper() ) );
//			std::cout << psi << " (" << m_subintervals.size() << ")" << std::endl;
			m_subintervals.erase( m_subintervals.begin() + i );
			--i;
		}
	}

	m_interval.assign( m_subintervals.front().lower(),
			m_subintervals.back().upper() );
}


Numeric VariableDomain::realWidth() const {
	if ( realwidthCached < 0.0 ) {
		Numeric w = 0.0;
		for ( const auto & si : m_subintervals ) {
			w += BOOSTNS::numeric::width( si );
		}
		realwidthCached = w;
	}

	return realwidthCached;
}


size_t VariableDomain::closestSubinterval( const Numeric val ) const {
		// TODO: make this more efficient if there are ever lots of subintervals

	assert( !m_subintervals.empty() );
	if ( m_subintervals.size() == 1 ) return 0;

	IntervalVec::const_iterator lb( std::lower_bound( m_subintervals.begin(),
			m_subintervals.end(), val, VariableDomain::Comp() ) );

	if ( lb == m_subintervals.end() ) {
		assert( ( lb - 1 )->upper() < val );
		return ( (lb - 1) - m_subintervals.begin() );
	}

	assert( val <= lb->upper() );

	if ( BOOSTNS::numeric::in( val, *lb ) || lb == m_subintervals.begin() ) {
		return ( lb - m_subintervals.begin() );
	}

	if ( ( val - (lb-1)->upper() ) < ( lb->lower() - val ) ) {
		return ( (lb-1) - m_subintervals.begin() );
	}

	return ( lb - m_subintervals.begin() );
}


Numeric VariableDomain::closestVal( const Numeric val ) const {
	const NumericInterval & si = m_subintervals[closestSubinterval( val )];
	if ( BOOSTNS::numeric::in( val, si ) ) return val;
	else if ( val < si.lower() ) return si.lower();
	else return si.upper();
}


std::ostream & operator<<( std::ostream & out, const VariableDomain & vd ) {

	if ( vd.subintervals().size() == 1 ) {
		out << "[" << vd.min() << "," << vd.max() << "]";
	} else {
		out << "{ ";
		for ( unsigned long int i = 0; i < vd.subintervals().size(); ++i ) {
			const NumericInterval & si( vd.subintervals()[i] );
			out << (i == 0 ? "" : ", ") << "[" << si.lower() << "," <<
					si.upper() << "]";
		}
		out << " }";
	}

	return out;
}

} // namespace rdis
