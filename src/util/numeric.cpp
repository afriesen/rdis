/*
 * numeric.cpp
 *
 *  Created on: Oct 6, 2013
 *      Author: afriesen
 */

#include "util/numeric.h"

namespace rdis {

Numeric power( Numeric val, Numeric exp ) {
	if ( exp == 0. ) {
		val = 1.;
	} else if ( exp == 1. ) {
		// do nothing
	} else if ( exp == 2. ) {
		val *= val;
	} else {
		val = std::pow( val, exp );
	}
	return val;
}


NumericInterval power( NumericInterval ival, Numeric exp ) {
	if ( exp == 0. ) {
		ival.assign( 1., 1. );
	} else if ( exp == 1. ) {
		// do nothing
	} else if ( exp == 2. ) {
		ival = BOOSTNS::numeric::square( ival );
	} else if ( exp < 0 ) {
		power( ival, -exp );
		ival = BOOSTNS::numeric::interval_lib::multiplicative_inverse(
				ival );
	} else {
		ival = BOOSTNS::numeric::pow( ival, exp );
	}

	return ival;
}


namespace NumericVecOps {

NumericVec & operator+=( NumericVec & lhs, const NumericVec & rhs ) {
	assert( lhs.size() == rhs.size() );
	for ( size_t i = 0; i < lhs.size(); ++i ) {
		lhs[i] += rhs[i];
	}
	return lhs;
}

NumericVec & operator-=( NumericVec & lhs, const NumericVec & rhs ) {
	assert( lhs.size() == rhs.size() );
	for ( size_t i = 0; i < lhs.size(); ++i ) {
		lhs[i] -= rhs[i];
	}
	return lhs;
}


NumericVec & operator*=( NumericVec & lhs, Numeric rhs ) {
	for ( size_t i = 0; i < lhs.size(); ++i ) {
		lhs[i] *= rhs;
	}
	return lhs;
}


NumericVec operator*( const NumericVec & lhs, Numeric rhs ) {
	NumericVec newlhs = lhs;
	for ( size_t i = 0; i < lhs.size(); ++i ) {
		newlhs[i] = lhs[i] * rhs;
	}
	return newlhs;
}


NumericVec operator/( const NumericVec & lhs, Numeric rhs ) {
	NumericVec newlhs = lhs;
	for ( size_t i = 0; i < lhs.size(); ++i ) {
		newlhs[i] = lhs[i] / rhs;
	}
	return newlhs;
}


NumericVec & abs( NumericVec & vec ) {
	for ( Numeric & v : vec ) {
		v = std::abs( v );
	}
	return vec;
}


NumericVec & log( NumericVec & vec ) {
	for ( Numeric & v : vec ) {
		v = std::log( v );
	}
	return vec;
}


Numeric dot( const NumericVec & lhs, const NumericVec & rhs ) {
	assert( lhs.size() == rhs.size() );
	Numeric eval = 0;
	for ( size_t i = 0; i < lhs.size(); ++i ) {
		eval += lhs[i] * rhs[i];
	}
	return eval;
}


bool normalize( NumericVec & v ) {
	Numeric sum = 0;
	for ( Numeric & val : v ) {
		val = std::abs( val );
		sum += val;
	}

	if ( sum > 0.0 ) {
		for ( Numeric & val : v ) {
			val /= sum;
		}
	} else {
		return false;
	}

	return true;
}

} // namespace NumericVecOps

} // namespace rdis
