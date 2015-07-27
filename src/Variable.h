/*
 * Variable.h
 *
 *  Created on: Jan 12, 2013
 *      Author: afriesen
 */

#ifndef RDIS_VARIABLE_H_
#define RDIS_VARIABLE_H_

#include "common.h"
#include "VariableDomain.h"
#include "Factor.h"

namespace rdis
{

// forward declaration to deal with circular dependency
class Factor;


// a variable contained in a polynomial
class Variable {
protected:
public:
	// constructor for unspecified domain
	Variable( VariableID id, const string & varname );

	// constructor for custom domain
	Variable( VariableID id, const string & varname, const string & vdomain );

	// constructor for default or custom domain
	Variable( VariableID id, const string & varname,
			const VariableDomain & vdomain );

	// destructor
	~Variable();

public:
	// add the specified factor to the list of factors that this variable is in
	void addFactor( Factor * f );

	// clear the factor info stored in this variable so factors can be reset
	void clearFactors();

public:
	void assign( Numeric newval, bool notifyFactors = true );
	void unassign();

	inline Numeric eval() const {
		assert( m_isAssigned );
		return m_value;
	}

	inline NumericInterval curBounds() const {
		return ( m_isAssigned ? NumericInterval( m_value ) : m_domain.interval() );
	}

public:
	inline const VariableID & getID() const { return m_id; }

	inline string const & getName() const { return m_name; }
	inline void setName( const string & name ) { m_name = name; }

	// get / set the full domain of this variable
	inline const VariableDomain & getDomain() const { return m_domain; }
	inline void setDomain( string const & dom ) { setDomain( VariableDomain( dom )); }
	void setDomain( const VariableDomain & dom );

	// get / set the max domain of this variable (the full domain must be a subdomain
	// of the max domain -- defaults to [Numeric::minval, Numeric::maxval] )
	inline const VariableDomain & getMaxDomain() const { return m_maxDomain; }
	void setMaxDomain( string const & dom ) { m_maxDomain = dom; }
	void setMaxDomain( const VariableDomain & dom ) { m_maxDomain = dom; }

	inline void setSamplingInterval( const NumericInterval & intvl ) {
		m_domain.setSamplingInterval( intvl );
	}

	std::vector< Factor * > & getFactors() { return m_factors; }
	const std::vector< Factor * > & getFactors() const { return m_factors; }

	inline bool isAssigned() const { return m_isAssigned; }

	friend bool operator==( const Variable & v1, const Variable & v2 ) {
		return( v1.m_id == v2.m_id );
	}

private:
	// the unique ID of this variable (these should be 0-indexed)
	VariableID m_id;

	// the name of this variable
	string m_name;

	// the maximum domain bounds of this variable (defaults to [-Inf, Inf])
	// -- the domain must be a subrange of this
	VariableDomain m_maxDomain;

	// the current (full) domain of this variable
	VariableDomain m_domain;

	// the factors that this variable is in
	std::vector< Factor * > m_factors;

	// true if this variable is currently assigned a value
	bool m_isAssigned;

	// the current domain interval assigned to this variable
	Numeric m_value;
};


typedef std::vector< Variable * > VariablePtrVec;

} // namespace rdis

#endif // RDIS_VARIABLE_H_
