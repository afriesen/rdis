/*
 * BundleAdjustmentCommon.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: afriesen
 */

#include "BundleAdjustmentCommon.h"

namespace rdis {

NumericInterval normalize( const NumericInterval vun[3], NumericInterval v[3] ) {
	using namespace BOOSTNS::numeric;
	const int idx[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
	const NumericInterval one( 1 );

	// for each dimension,
	// 	compute x / (x^2 + y^2 + z^2) = sign(x) * 1 / sqrt( 1 + (y^2 + z^2) / x^2 )

	for ( int i = 0; i < 3; ++i ) {
		bool doSecond( false );
		NumericInterval first( vun[i] ), second( 0 );
		if ( zero_in( first ) ) {
			doSecond = true;
			second.assign( 0, first.upper() );
			first.assign( first.lower(), 0 );
		}

		NumericInterval val( square( vun[idx[i][0]] ) + square( vun[idx[i][1]] ) );
		v[i] = one / sqrt( one + val / square( first ) );
		if ( std::signbit( first.lower() ) ) v[i] *= -1;

		if ( doSecond ) {
			v[i] = hull( v[i], one / sqrt( one + val / square( second ) ) );
		}
	}

	// return the normalizer
	return sqrt( square( vun[0] ) + square( vun[1] ) + square( vun[2] ) );
}


NumericInterval simpleNormalize( const NumericInterval vun[3],
		NumericInterval v[3] ) {
	namespace bn = BOOSTNS::numeric;
	NumericInterval norm( bn::sqrt( bn::square( vun[0] ) + bn::square( vun[1] ) +
			bn::square( vun[2] ) ) );
	for ( int i = 0; i < 3; ++i ) {
		v[i] = vun[i] / norm;
	}

	return norm;
}


//Numeric normalize( const Numeric vun[3], Numeric v[3] ) {
//	Numeric norm( std::sqrt( vun[0]*vun[0] + vun[1]*vun[1] + vun[2]*vun[2] ) );
//	v[0] = vun[0] / norm;
//	v[1] = vun[1] / norm;
//	v[2] = vun[2] / norm;
//	return norm;
//}

} // namespace rdis

