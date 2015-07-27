/*
 * BundleAdjustmentCommon.h
 *
 *  Created on: Apr 6, 2013
 *      Author: afriesen
 */

#ifndef BUNDLEADJUSTMENTCOMMON_H_
#define BUNDLEADJUSTMENTCOMMON_H_

#include "common.h"

#include BOOSTPATH/numeric/interval.hpp>

//#include BOOSTPATH/numeric/ublas/vector.hpp>

namespace rdis {


//#define SEPARATE_THETA_VAR


typedef long long int CameraID;
typedef long long int PointID;

//namespace ObservationDimension {
//enum ObservationDim {
//	x = 0,
//	y = 1,
//	ndims
//};
//} // namespace ObservationDimension

namespace BundleAdjust {

enum VariableType {
	CAM_ROT_X = 0,	// camera rotation params for rotating into camera
	CAM_ROT_Y,		// coord. frame (can be angle+vector, or quaternion)
	CAM_ROT_Z,

#ifdef SEPARATE_THETA_VAR
	CAM_THETA,
#endif

	CAM_TRANS_X,	// translation to camera coord. frame
	CAM_TRANS_Y,
	CAM_TRANS_Z,

	CAM_FOCAL,		// focal length

	CAM_RDL_K1, 	// radial distortion params
	CAM_RDL_K2,

	POINT_X,		// the point that will be projected
	POINT_Y,
	POINT_Z,

	NUMPARAMS 		// keep at end of list
};

} // namespace BundleAdjust


template< class T >
void cross( const T a[3], const T b[3], T a_cross_b[3] ) {
	a_cross_b[0] = a[1]*b[2] - a[2]*b[1];
	a_cross_b[1] = a[2]*b[0] - a[0]*b[2];
	a_cross_b[2] = a[0]*b[1] - a[1]*b[0];
}

template< class T >
T dot( const T a[3], const T b[3] ) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// normalize the vector vun and put it into v, return the normalizing constant
NumericInterval normalize( const NumericInterval vun[3], NumericInterval v[3] );
NumericInterval simpleNormalize( const NumericInterval vun[3],
		NumericInterval v[3] );

inline Numeric normalize( const Numeric vun[3], Numeric v[3] ) {
	Numeric norm( std::sqrt( vun[0]*vun[0] + vun[1]*vun[1] + vun[2]*vun[2] ) );
	if ( norm != 0.0 ) {
		v[0] = vun[0] / norm;
		v[1] = vun[1] / norm;
		v[2] = vun[2] / norm;
	} else {
		v[0] = vun[0];
		v[1] = vun[1];
		v[2] = vun[2];
	}
	return norm;
}

} // namespace rdis

#endif // BUNDLEADJUSTMENTCOMMON_H_
