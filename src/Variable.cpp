/*
 * Variable.cpp
 *
 *  Created on: Jan 12, 2013
 *      Author: afriesen
 */

#include "Variable.h"

namespace rdis
{

Variable::Variable( VariableID id, string const & varname )
	: m_id( id )
	, m_name( varname )
	, m_maxDomain( NumericMIN, NumericMAX )
	, m_domain( 0 )
	, m_isAssigned( false )
	, m_value( 0 )
{
	m_factors.reserve( 10 );
}

Variable::Variable( VariableID id, string const & varname,
		string const & vdomain )
	: m_id( id )
	, m_name( varname )
	, m_maxDomain( NumericMIN, NumericMAX )
	, m_domain( vdomain )
	, m_isAssigned( false )
	, m_value( 0 )
{
	m_factors.reserve( 10 );
}

Variable::Variable( VariableID id, string const & varname,
		VariableDomain const & vdomain )
	: m_id( id )
	, m_name( varname )
	, m_maxDomain( NumericMIN, NumericMAX )
	, m_domain( vdomain )
	, m_isAssigned( false )
	, m_value( 0 )
{
	m_factors.reserve( 10 );
}

Variable::~Variable() {}


void Variable::addFactor( Factor * f ) {

	// require factors to be sorted -- useful for merging factors later
	assert( m_factors.empty() || f->getID() > m_factors.back()->getID() );

	m_factors.push_back( f );
}


void Variable::clearFactors() {
	assert( !m_isAssigned );
	m_factors.clear();
}


void Variable::assign( Numeric newval, bool notifyFactors ) {

	if ( m_isAssigned ) {
		const Numeric tol = 1e-12;

		// don't do notifications if nothing changed
		if ( notifyFactors && !approxeq( newval, m_value, tol ) ) {
			for ( FactorID fid = 0; fid < (FactorID) m_factors.size(); ++fid ) {
				m_factors[fid]->onVarChanged( m_id, m_value, newval );
			}
		}
	} else {
		m_isAssigned = true;

		if ( notifyFactors ) {
			for ( FactorID fid = 0; fid < (FactorID) m_factors.size(); ++fid ) {
				m_factors[fid]->onVarAssigned( m_id, newval );
			}
		}
	}

	m_value = newval;
}

void Variable::unassign() {
	assert( m_isAssigned );

	m_isAssigned = false;

	// tell associated factors that this var has been unassigned
	for ( FactorID fid = 0; fid < (FactorID) m_factors.size(); ++fid ) {
		m_factors[fid]->onVarUnassigned( m_id, m_value );
	}

	m_value = 0;
}


void Variable::setDomain( const VariableDomain & dom ) {
	if ( !BOOSTNS::numeric::overlap( m_maxDomain.interval(), dom.interval() ) ) {
		std::cout << "failed to set domain for var " << m_id << " to " <<
				dom.interval() << " (max: " << m_maxDomain.interval() << ")" << std::endl;
	}
	assert( BOOSTNS::numeric::overlap( m_maxDomain.interval(), dom.interval() ) );
	m_domain = dom;

	for ( const auto & f : m_factors ) {
		f->onVarDomainChanged( m_id );
	}
}

} // namespace rdis
