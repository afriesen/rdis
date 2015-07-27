/*
 * SubspaceOptimizer.cpp
 *
 *  Created on: Mar 25, 2014
 *      Author: afriesen
 */

#include "SubspaceOptimizer.h"

namespace rdis {

SubspaceOptimizer::SubspaceOptimizer( OptimizableFunction & f_ )
	: f( f_ )
	, doAscent( f.semiring().Sum( 2, 6 ) == 6 )
	, maxiters( 50 )
	, ftol( 3.0e-8 )
{
}

SubspaceOptimizer::~SubspaceOptimizer() {}


void SubspaceOptimizer::setParameters(
		const BOOSTNS::program_options::variables_map & options )
{
	if ( options.count( "SSmaxit" ) ) {
		maxiters = options[ "SSmaxit" ].as< VariableCount >();
	}

	if ( options.count( "SSftol" ) ) {
		ftol = options[ "SSftol" ].as< Numeric >();
	}

	assert( maxiters > 0 );
}


void SubspaceOptimizer::quickAssignVals( const VariablePtrVec & vars,
		const NumericVec & xval, bool sanitizeVals ) {
	assert( vars.size() == xval.size() );
	for ( size_t i = 0; i < vars.size(); ++i ) {
		Numeric val = sanitizeVals ?
					  vars[i]->getDomain().closestVal( xval[i] ) : xval[i];

		if ( val != xval[i] ) {
			std::cout << "var " << vars[i]->getID() << " xval " << xval[i] <<
					" changed to " << val << std::endl;
		}

		vars[i]->assign( val );
		f.onVarAssigned( vars[i]->getID(), val );
	}
}

} // namespace rdis
