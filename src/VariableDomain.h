/*
 * VariableDomain.h
 *
 *  Created on: Jan 14, 2013
 *      Author: afriesen
 */

#ifndef RDIS_VARIABLEDOMAIN_H_
#define RDIS_VARIABLEDOMAIN_H_

#include "common.h"

namespace rdis {

typedef std::vector< NumericInterval > IntervalVec;

// class representing a custom domain for a variable in a polynomial
class VariableDomain {

	struct Comp {
		bool operator()( const NumericInterval & ni, const Numeric & val ) const;
	};

	typedef std::vector< string > SplitVector;

public:
	VariableDomain();
	VariableDomain( Numeric value );
	VariableDomain( Numeric lower, Numeric upper );
	VariableDomain( const string & domain );
	VariableDomain( const IntervalVec & subintvls );
	~VariableDomain() {};

	void setSamplingInterval( const NumericInterval & samplingIntvl ) {
		assert( BOOSTNS::numeric::subset( samplingIntvl, m_interval ) );
		m_samplingInterval = samplingIntvl;
	}
	NumericInterval getSamplingInterval() const { return m_samplingInterval; }

	void operator=( const NumericInterval & iv ) {
		m_interval = iv;
		m_samplingInterval = iv;
		m_subintervals.clear();
		m_subintervals.push_back( m_interval );
		realwidthCached = -1.0;
	}

	void operator=( const string & domain ) {
		realwidthCached = -1.0;
		parse( domain );
	}

	const NumericInterval & interval() const { return m_interval; }
	const IntervalVec & subintervals() const { return m_subintervals; }

	Numeric max() const { return m_interval.upper(); }
	Numeric min() const { return m_interval.lower(); }

	// note: this could be outside one of the subintervals...
	Numeric median() const { return BOOSTNS::numeric::median( m_interval ); }

	// return the true width of this domain (i.e., the sum of the widths of the
	// subintervals)
	Numeric realWidth() const;

	// return the index of the subinterval that is closest to the specified val
	size_t closestSubinterval( const Numeric val ) const;

	// return the value (in the domain) that is closest to the specified val
	Numeric closestVal( const Numeric val ) const;

protected:
	void parse( string const & domain );

	void cleanSubintervals();

protected:
	// the full, top-level interval representing this domain (the convex hull of
	// all of the subintervals)
	NumericInterval m_interval;

	// if drawing samples from this var's domain, draw them from this interval
	// instead of the full one (useful for rotations, etc.)
	NumericInterval m_samplingInterval;

	// the set of sub-intervals within this domain (if empty, only
	IntervalVec m_subintervals;

	mutable Numeric realwidthCached;
};


std::ostream & operator<<( std::ostream & out, const VariableDomain & vd );

typedef std::vector< VariableDomain > VarDomainVec;

} // namespace rdis

#endif // RDIS_VARIABLEDOMAIN_H_
