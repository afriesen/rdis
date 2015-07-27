/*
 * OptimizerVariableData.cpp
 *
 *  Created on: May 3, 2013
 *      Author: afriesen
 */

#include "OptimizerVariableData.h"

namespace rdis {


VarData::VarData( Variable * v_ )
	: v( v_ )
	, assignmentIndex( -1 )
	, xvalinit( 0 )
	, xvinitIsSet( false )
{
}


void VarData::reset() {
	assignmentIndex = -1;

	xvalinit = 0;
	xvinitIsSet = false;
}


std::ostream & operator<<( std::ostream & os, const VarData & vd ) {
	os << vd.v->getID();
	return os;
}


OptimizerVariableData::OptimizerVariableData( const VariableCount & numVars )
	: m_numVars( numVars )
	, m_numAssigned( 0 )
{}

OptimizerVariableData::~OptimizerVariableData() {}


void OptimizerVariableData::init( const VariablePtrVec & vars ) {

	m_vardata.reserve( m_numVars );
	m_isAssigned.resize( m_numVars, false );

	for ( VariableCount i = 0; i < m_numVars; ++i ) {
		m_vardata.emplace_back( vars[i] );
	}
}


void OptimizerVariableData::reset() {

	for ( VariableID vid = 0; vid < m_numVars; ++vid ) {
		get( vid ).reset();
	}

	assert( m_numAssigned == 0 );
}


void OptimizerVariableData::assign( VariableID vid, const Numeric & val ) {

	VarData & vd( get( vid ) );
	bool wasUnassigned( isUnassigned( vid ) );

	if ( wasUnassigned ) {
		vd.assignmentIndex = m_numAssigned;
		setAssigned( vid );
	}

	// assign the variable
	vd.v->assign( val );
}


void OptimizerVariableData::unassign( VariableID vid ) {
	assert( isAssigned( vid ) );

	VarData & vd( get( vid ) );
	vd.v->unassign();

	setUnassigned( vid );
}

} // namespace rdis
