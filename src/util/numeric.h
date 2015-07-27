/*
 * numeric.h
 *
 *  Created on: Oct 2, 2013
 *      Author: afriesen
 */

#ifndef RDIS_NUMERIC_H_
#define RDIS_NUMERIC_H_

#include "common.h"

namespace rdis {

// somewhat efficiently compute the power of a value
Numeric power( Numeric val, Numeric exp );

// somewhat efficiently compute the power of an interval
NumericInterval power( NumericInterval ival, Numeric exp );


namespace NumericVecOps {

NumericVec & operator+=( NumericVec & lhs, const NumericVec & rhs );
NumericVec & operator-=( NumericVec & lhs, const NumericVec & rhs );

NumericVec & operator*=( NumericVec & lhs, Numeric rhs );

NumericVec operator*( const NumericVec & lhs, Numeric rhs );
NumericVec operator/( const NumericVec & lhs, Numeric rhs );

NumericVec & abs( NumericVec & vec );
NumericVec & log( NumericVec & vec );

Numeric dot( const NumericVec & lhs, const NumericVec & rhs );

bool normalize( NumericVec & v );

} // namespace NumericVecOps

} // namespace rdis

#endif // RDIS_NUMERIC_H_
