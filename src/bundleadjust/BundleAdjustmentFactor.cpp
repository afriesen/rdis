/*
 * BundleAdjustmentFactor.cpp
 *
 *  Created on: Apr 6, 2013
 *      Author: afriesen
 */

#include "bundleadjust/BundleAdjustmentFactor.h"

namespace rdis {

using namespace BundleAdjust;

BundleAdjustmentFactor::BundleAdjustmentFactor( FactorID factorID,
			CameraID cameraID, PointID pointID, Numeric obsX, Numeric obsY )
		: super( factorID )
		, m_cameraID( cameraID )
		, m_pointID( pointID )
		, m_obsX( obsX )
		, m_obsY( obsY ) {
	bundleAdjVars.resize( NUMPARAMS );
	m_tempVarIntervals.resize( NUMPARAMS, 0 );
	m_tempVarVals.resize( NUMPARAMS, 0 );
	m_tempGradient.resize( NUMPARAMS, 0 );
}

BundleAdjustmentFactor::~BundleAdjustmentFactor() {}


void BundleAdjustmentFactor::addVariable( Variable * vp,
		VariableType index ) {

//	std::cout << "adding var " << vp->getID() << " (" << vp->getName() <<
//			") to factor " << id << ", cam: " << m_cameraID << ", point: " <<
//			m_pointID << " at index " << index << std::endl;

	super::addVariable( vp );
	bundleAdjVars[index] = vp;
}


Numeric BundleAdjustmentFactor::evalNoCache() const {
	return evalFactor();
}


NumericInterval BundleAdjustmentFactor::computeFactorBounds(
		const VariableIDVec & vidsToIgnore ) const {
	getVarVals( m_tempVarIntervals, false, !vidsToIgnore.empty(),
			vidsToIgnore, true );
	return evalFactor( m_tempVarIntervals );
}


Numeric BundleAdjustmentFactor::evalFactor() const {

	assert( areAllVarsAssigned() );
//	checkDerivatives();
	getVarVals( m_tempVarVals );
	const Numeric fev = evalFactor( m_tempVarVals );
//	std::cout << "f " << id << " vars: " << m_tempVarVals << std::endl;
//	std::cout << "f " << id << " eval: " << fev << std::endl;
	return fev;
}


void BundleAdjustmentFactor::getVarVals( IntervalVec & vals,
		bool useFullDomains, bool ignoreSpecificVarAssignments,
		const VariableIDVec & vidsToIgnore, bool useAssigned ) const {

	assert( vals.size() == NUMPARAMS );

	for ( int i = 0; i < NUMPARAMS; ++i ) {
		Variable & v( *bundleAdjVars[i] );

		if ( useFullDomains ) {
			vals[i] = v.getDomain().interval();

		} else if ( ignoreSpecificVarAssignments &&
				std::find( vidsToIgnore.begin(), vidsToIgnore.end(),
						v.getID() ) != vidsToIgnore.end() ) {
//				v.getID() == vidToIgnore ) {
			vals[i] = v.getDomain().interval();

		} else if ( useAssigned && v.isAssigned() ) {
			vals[i] = v.eval();

		} else {
			vals[i] = v.getDomain().interval();
		}
	}
}


void BundleAdjustmentFactor::getVarVals( NumericVec & vals ) const {

	assert( vals.size() == NUMPARAMS );
	for ( int i = 0; i < NUMPARAMS; ++i ) {
		vals[i] = bundleAdjVars[i]->eval();
	}
}


NumericInterval BundleAdjustmentFactor::evalFactor(
		const IntervalVec & vals ) const {
	namespace bn = BOOSTNS::numeric;

	NumericInterval p[3];
	angleAxisRotatePoint( vals, p );

////    std::cout << "\tv: (" << camera[0] << ", " << camera[1] << ", " << camera[2]
////        << ")" << std::endl;
//std::cout << "\trotated: (" << p[0].lower() << ", " << p[1].lower() << ", "
//		<< p[2].lower() << ")" << std::endl;

	// translate the points to the camera's coordinate frame
	p[0] += vals[CAM_TRANS_X];
	p[1] += vals[CAM_TRANS_Y];
	p[2] += vals[CAM_TRANS_Z];

	// perspective divison: -[px, py] /= pz
	p[0] = -p[0] / p[2];
	p[1] = -p[1] / p[2];

	// r^2 = px^2 + py^2
	NumericInterval r2( bn::square( p[0] ) + bn::square( p[1] ) );

	const NumericInterval & f( vals[CAM_FOCAL] );
	const NumericInterval & k1( vals[CAM_RDL_K1] );
	const NumericInterval & k2( vals[CAM_RDL_K2] );

	// distortion = ( 1 + k1*r^2 + k2*r^4 )
	NumericInterval dstn( NumericInterval( 1 ) + k1 * r2 + k2 * bn::square( r2 ) );

	// final projection: p[d] = f * distortion * p[d]
	p[0] = f * dstn * p[0];
	p[1] = f * dstn * p[1];

//std::cout << "\tdistortion: " << dstn.lower() << std::endl;
//std::cout << "\tfocal: " << f.lower() << std::endl;
//std::cout << "\tfinal: (" << p[0] << ", " << p[1] << ")"
//		<< std::endl;


//std::cout << "\tobs: (" << m_obsX << ", " << m_obsY << ")" << std::endl;
//std::cout << "\tresiduals: (" << p[0] - m_obsX << ", " <<
//		p[1] - m_obsY << ")" << std::endl;

	// error = (px - obsx)^2 + (py - obsy)^2 -- have to do computation
	// separately because was getting weird results otherwise
	NumericInterval errorX( bn::square( p[0] - m_obsX ) );
	NumericInterval errorY( bn::square( p[1] - m_obsY ) );

	NumericInterval error( ( errorX + errorY ) / 2.0 );
//	std::cout << "factor " << id << ": " << error << std::endl;
	return error;
}


Numeric BundleAdjustmentFactor::evalFactor( const NumericVec & vals ) const {

	Numeric v[3], theta( 0 ), P[3], pp[2], r2, dstn, px[2], res[2];
	evalPixelVals( vals, v, theta, P, pp, r2, dstn, px );
	return error( px, res );
}


void BundleAdjustmentFactor::evalPixelVals( const NumericVec & vals,
		Numeric v[3], Numeric & theta, Numeric P[3], Numeric pp[2],
		Numeric & r2, Numeric & dstn, Numeric pxl[2] ) const {

	Numeric vxp[3], vdp;

	// compute the rotated points p, then translate them into the camera's
	// reference frame and do perspective division
	angleAxisRotateTranslatePoint( vals, v, theta, vxp, vdp, P );
	perspectiveDivide( P, pp );

	// final projection: p[d] = f * distortion * p[d]
	distort( vals, pp, r2, dstn, pxl );

	//std::cout << "\tobs: (" << m_obsX << ", " << m_obsY << ")" << std::endl;
	//std::cout << "\tresiduals: (" << pOut[0] - m_obsX << ", " <<
	//		pOut[1] - m_obsY << ")" << std::endl;
}


void BundleAdjustmentFactor::angleAxisRotatePoint( const IntervalVec & vals,
		NumericInterval p[3] ) const {

	namespace bn = BOOSTNS::numeric;
	NumericInterval v[3];

	// unnormalized rotation vector (vun)
	const NumericInterval vun[3] = { vals[CAM_ROT_X], vals[CAM_ROT_Y],
			vals[CAM_ROT_Z] };

	// get the point and theta
	p[0] = vals[POINT_X]; p[1] = vals[POINT_Y]; p[2] = vals[POINT_Z];


#ifndef SEPARATE_THETA_VAR
	NumericInterval theta( bn::sqrt( bn::square( vun[0] ) + bn::square( vun[1] ) +
			bn::square( vun[2] ) ) );
#else
	NumericInterval theta( vals[CAM_THETA] );
#endif // SEPARATE_THETA_VAR

	// normalized rotation vector (v)
//	normalize( vun, v );
	simpleNormalize( vun, v );


//std::cout << "\ttheta: " << theta << std::endl;
//std::cout << "\tv: (" << vun[0] << ", " << vun[1] << ", " << vun[2] << ")" << std::endl;
//std::cout << "\tvnorm: (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;

//	if ( theta > 0.0 ) {
		NumericInterval costheta, sintheta;
		if ( bn::width( theta ) < 1e-6 ) {
			costheta = std::cos( bn::median( theta ) );
			sintheta = std::sin( bn::median( theta ) );
		} else {
			costheta = bn::cos( theta );
			sintheta = bn::sin( theta );
		}
		NumericInterval oneminusctheta( NumericInterval( 1 ) - costheta );

		// compute cross(v, p)
		NumericInterval v_cross_p[3];
		cross( v, p, v_cross_p );

		// compute the dot(v, p)
		NumericInterval v_dot_p( dot( v, p ) );

		// rotate the points into the camera's coordinate frame
		p[0] = p[0]*costheta + v_cross_p[0]*sintheta + v[0]*oneminusctheta*v_dot_p;
		p[1] = p[1]*costheta + v_cross_p[1]*sintheta + v[1]*oneminusctheta*v_dot_p;
		p[2] = p[2]*costheta + v_cross_p[2]*sintheta + v[2]*oneminusctheta*v_dot_p;
//	} else {
//	    // ---- the following is taken from Sameer Agarwal's Ceres code ----
//
//		// Near zero, the first order Taylor approximation of the rotation
//	    // matrix R corresponding to a vector w and angle theta is
//	    //
//	    //   R = I + hat(w) * sin(theta)
//	    //
//	    // But sintheta ~ theta and theta * w = angle_axis, which gives us
//	    //
//	    //  R = I + hat(w)
//	    //
//	    // and actually performing multiplication with the point pt, gives us
//	    // R * pt = pt + w x pt.
//	    //
//	    // Switching to the Taylor expansion at zero helps avoid all sorts
//	    // of numerical nastiness.
//	    NumericInterval v_cross_p[3];
//	    cross( v, p, v_cross_p );
//	    for (int i = 0; i < 3; ++i) {
//	    	p[i] = p[i] + v_cross_p[i];
//	    }
//	}
}


void BundleAdjustmentFactor::angleAxisRotateTranslatePoint(
		const NumericVec & vals, Numeric v[3], Numeric & theta,
		Numeric v_cross_p[3], Numeric & v_dot_p,
		Numeric P[3] ) const {

	// unnormalized rotation vector (vun)
	const Numeric vun[3] = { vals[CAM_ROT_X], vals[CAM_ROT_Y], vals[CAM_ROT_Z] };

	// get the point and theta
	P[0] = vals[POINT_X]; P[1] = vals[POINT_Y]; P[2] = vals[POINT_Z];


	// get theta and normalized v
#ifndef SEPARATE_THETA_VAR
	theta = normalize( vun, v );
#else
	theta = vals[CAM_THETA];
	normalize( vun, v );
#endif // SEPARATE_THETA_VAR


//std::cout << "\ttheta: " << theta << std::endl;
//std::cout << "\tv: (" << vun[0] << ", " << vun[1] << ", " << vun[2] << ")" << std::endl;
//std::cout << "\tvnorm: (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;

	if ( theta > 0.0 ) {
		Numeric costheta( std::cos( theta ) );
		Numeric sintheta( std::sin( theta ) );
		Numeric oneminusctheta( 1 - costheta );

		// compute the cross product, v x p
		cross( v, P, v_cross_p );

		// compute the dot product, v.p
		v_dot_p = dot( v, P );

		// rotate the points into the camera's coordinate frame
		P[0] = P[0]*costheta + v_cross_p[0]*sintheta + v[0]*oneminusctheta*v_dot_p;
		P[1] = P[1]*costheta + v_cross_p[1]*sintheta + v[1]*oneminusctheta*v_dot_p;
		P[2] = P[2]*costheta + v_cross_p[2]*sintheta + v[2]*oneminusctheta*v_dot_p;
	} else {
	    // ---- the following is taken from Sameer Agarwal's Ceres code ----

		// Near zero, the first order Taylor approximation of the rotation
	    // matrix R corresponding to a vector w and angle theta is
	    //
	    //   R = I + hat(w) * sin(theta)
	    //
	    // But sintheta ~ theta and theta * w = angle_axis, which gives us
	    //
	    //  R = I + hat(w)
	    //
	    // and actually performing multiplication with the point pt, gives us
	    // R * pt = pt + w x pt.
	    //
	    // Switching to the Taylor expansion at zero helps avoid all sorts
	    // of numerical nastiness.
//	    Numeric v_cross_p[3];
	    cross( v, P, v_cross_p );
	    v_dot_p = 0;
	    for (int i = 0; i < 3; ++i) {
	    	P[i] = P[i] + v_cross_p[i];
	    }
	}

	// translate the points to the camera's coordinate frame
	P[0] += vals[CAM_TRANS_X];
	P[1] += vals[CAM_TRANS_Y];
	P[2] += vals[CAM_TRANS_Z];
}


void BundleAdjustmentFactor::computeGradient( PartialGradient & g,
		bool doGradCheck ) const {
	getVarVals( m_tempVarVals );
	/*Numeric feval =*/ computeGradient( m_tempVarVals, m_tempGradient );

	for ( VariableCount i = 0; i < NUMPARAMS; ++i ) {
		g[bundleAdjVars[i]->getID()] = m_tempGradient[i];
	}

	if ( doGradCheck ) checkGradient( g, 1e-8 );
}


Numeric BundleAdjustmentFactor::computeGradient( const NumericVec & vals,
		NumericVec & grad ) const {

	// this is the full factor eval, giving us all the intermediates we'll need
	// for gradient computations
	Numeric v[3], theta( 0 ), vxp[3], vdp, P[3], pp[2], r2, dstn, pix[2], res[2];
	angleAxisRotateTranslatePoint( vals, v, theta, vxp, vdp, P );
	perspectiveDivide( P, pp );
	distort( vals, pp, r2, dstn, pix );
	const Numeric feval( error( pix, res ) );

	const Numeric & f( vals[CAM_FOCAL] );
	const Numeric & k1( vals[CAM_RDL_K1] ), & k2( vals[CAM_RDL_K2] );

	const Numeric tmp1( 2.0*(k1 + 2.0*k2*r2) );

	const Numeric q[3] = { vals[POINT_X], vals[POINT_Y], vals[POINT_Z] };
	const Numeric sinth( std::sin( theta ) ), costh( std::cos( theta ) );

#ifdef SEPARATE_THETA_VAR
	const Numeric vun[3] = { vals[CAM_ROT_X], vals[CAM_ROT_Y], vals[CAM_ROT_Z] };
	const Numeric dthdvx( 0 ), dthdvy( 0 ), dthdvz( 0 );
	const Numeric vnorm( std::sqrt( vun[0]*vun[0] + vun[1]*vun[1] + vun[2]*vun[2] ) );
#else
	const Numeric & dthdvx( v[0] ), & dthdvy( v[1] ), & dthdvz( v[2] );
	const Numeric vnorm( theta );
#endif

	// pre-compute some partials and intermediates that will be used for all of
	// the derivatives
	const Numeric & dEdpx( res[0] ), dEdpy( res[1] );
	const Numeric P22( P[2]*P[2] );
	const Numeric pp00( pp[0]*pp[0] ), pp01( pp[0]*pp[1] ), pp11( pp[1]*pp[1] );

	const Numeric dpixxdppx( ( dstn + tmp1*pp00 ) );
	const Numeric dpixxdppy( ( tmp1*pp01 ) );
	const Numeric dpixydppx( ( tmp1*pp01 ) );
	const Numeric dpixydppy( ( dstn + tmp1*pp11 ) );

	const Numeric dPxdvpx( ( vdp + v[0]*q[0] ) * (1-costh) );
	const Numeric dPxdvpy(  q[2]*sinth + v[0]*q[1]*(1-costh) );
	const Numeric dPxdvpz( -q[1]*sinth + v[0]*q[2]*(1-costh) );
	const Numeric dPxdth(  -q[0]*sinth + vxp[0]*costh + v[0]*vdp*sinth );

	const Numeric dPydvpx( -q[2]*sinth + v[1]*q[0]*(1-costh) );
	const Numeric dPydvpy( ( vdp + v[1]*q[1] )*(1-costh) );
	const Numeric dPydvpz(  q[0]*sinth + v[1]*q[2]*(1-costh) );
	const Numeric dPydth(  -q[1]*sinth + vxp[1]*costh + v[1]*vdp*sinth );

	const Numeric dPzdvpx(  q[1]*sinth + v[2]*q[0]*(1-costh) );
	const Numeric dPzdvpy( -q[0]*sinth + v[2]*q[1]*(1-costh) );
	const Numeric dPzdvpz( ( vdp + v[2]*q[2] )*(1-costh) );
	const Numeric dPzdth(  -q[2]*sinth + vxp[2]*costh + v[2]*vdp*sinth );


	// dE / dvx
	const Numeric dvpxdvx( ( v[1]*v[1] + v[2]*v[2] ) / vnorm );
	const Numeric dvpydvx( -v[0]*v[1] / vnorm );
	const Numeric dvpzdvx( -v[0]*v[2] / vnorm );

	const Numeric dPxdvx( dPxdvpx*dvpxdvx + dPxdvpy*dvpydvx +
			dPxdvpz*dvpzdvx + dPxdth*dthdvx );
	const Numeric dPydvx( dPydvpx*dvpxdvx + dPydvpy*dvpydvx +
			dPydvpz*dvpzdvx + dPydth*dthdvx );
	const Numeric dPzdvx( dPzdvpx*dvpxdvx + dPzdvpy*dvpydvx +
			dPzdvpz*dvpzdvx + dPzdth*dthdvx );

	const Numeric dppxdvx( ( P[0]*dPzdvx - P[2]*dPxdvx ) / P22 );
	const Numeric dppydvx( ( P[1]*dPzdvx - P[2]*dPydvx ) / P22 );
	const Numeric drxdvx( dEdpx*( dpixxdppx*dppxdvx + dpixxdppy*dppydvx ) );
	const Numeric drydvx( dEdpy*( dpixydppx*dppxdvx + dpixydppy*dppydvx ) );
	const Numeric dEdvx( f * ( drxdvx + drydvx ) );

	// dE / dvy
	const Numeric dvpxdvy( -v[0]*v[1] / vnorm );
	const Numeric dvpydvy( ( v[0]*v[0] + v[2]*v[2] ) / vnorm );
	const Numeric dvpzdvy( -v[1]*v[2] / vnorm );

	const Numeric dPxdvy( dPxdvpx*dvpxdvy + dPxdvpy*dvpydvy +
			dPxdvpz*dvpzdvy + dPxdth*dthdvy );
	const Numeric dPydvy( dPydvpx*dvpxdvy + dPydvpy*dvpydvy +
			dPydvpz*dvpzdvy + dPydth*dthdvy );
	const Numeric dPzdvy( dPzdvpx*dvpxdvy + dPzdvpy*dvpydvy +
			dPzdvpz*dvpzdvy + dPzdth*dthdvy );

	const Numeric dppxdvy( ( P[0]*dPzdvy - P[2]*dPxdvy ) / P22 );
	const Numeric dppydvy( ( P[1]*dPzdvy - P[2]*dPydvy ) / P22 );
	const Numeric drxdvy( dEdpx*( dpixxdppx*dppxdvy + dpixxdppy*dppydvy ) );
	const Numeric drydvy( dEdpy*( dpixydppx*dppxdvy + dpixydppy*dppydvy ) );
	const Numeric dEdvy( f * ( drxdvy + drydvy ) );

	// dE / dvz
	const Numeric dvpxdvz( -v[0]*v[2] / vnorm );
	const Numeric dvpydvz( -v[1]*v[2] / vnorm );
	const Numeric dvpzdvz( ( v[0]*v[0] + v[1]*v[1] ) / vnorm );

	const Numeric dPxdvz( dPxdvpx*dvpxdvz + dPxdvpy*dvpydvz +
			dPxdvpz*dvpzdvz + dPxdth*dthdvz );
	const Numeric dPydvz( dPydvpx*dvpxdvz + dPydvpy*dvpydvz +
			dPydvpz*dvpzdvz + dPydth*dthdvz );
	const Numeric dPzdvz( dPzdvpx*dvpxdvz + dPzdvpy*dvpydvz +
			dPzdvpz*dvpzdvz + dPzdth*dthdvz );

	const Numeric dppxdvz( ( P[0]*dPzdvz - P[2]*dPxdvz ) / P22 );
	const Numeric dppydvz( ( P[1]*dPzdvz - P[2]*dPydvz ) / P22 );
	const Numeric drxdvz( dEdpx*( dpixxdppx*dppxdvz + dpixxdppy*dppydvz ) );
	const Numeric drydvz( dEdpy*( dpixydppx*dppxdvz + dpixydppy*dppydvz ) );
	const Numeric dEdvz( f * ( drxdvz + drydvz ) );

	// dE / dtheta
#ifdef SEPARATE_THETA_VAR
	const Numeric dcosthdth( -sinth );
	const Numeric dsinthdth( costh );

	const Numeric dPxdtheta( q[0]*dcosthdth + vxp[0]*dsinthdth - v[0]*vdp*dcosthdth );
	const Numeric dPydtheta( q[1]*dcosthdth + vxp[1]*dsinthdth - v[1]*vdp*dcosthdth );
	const Numeric dPzdtheta( q[2]*dcosthdth + vxp[2]*dsinthdth - v[2]*vdp*dcosthdth );

	const Numeric dppxdth( ( P[0]*dPzdtheta - P[2]*dPxdtheta ) / P22 );
	const Numeric dppydth( ( P[1]*dPzdtheta - P[2]*dPydtheta ) / P22 );
	const Numeric drxdth( dEdpx*( dpixxdppx*dppxdth + dpixxdppy*dppydth ) );
	const Numeric drydth( dEdpy*( dpixydppx*dppxdth + dpixydppy*dppydth ) );
	const Numeric dEdtheta( f * ( drxdth + drydth ) );
#endif

	// dE / dqx
	const Numeric dPxdqx( costh * ( 1.0 - v[0]*v[0] ) + v[0]*v[0] );
	const Numeric dPydqx( v[2]*sinth + v[0]*v[1]*( 1.0 - costh ) );
	const Numeric dPzdqx( -v[1]*sinth + v[0]*v[2]*( 1.0 - costh ) );

	const Numeric dppxdqx( ( P[0]*dPzdqx - P[2]*dPxdqx ) / P22 );
	const Numeric dppydqx( ( P[1]*dPzdqx - P[2]*dPydqx ) / P22 );
	const Numeric drxdqx( dEdpx*( dpixxdppx*dppxdqx + dpixxdppy*dppydqx ) );
	const Numeric drydqx( dEdpy*( dpixydppx*dppxdqx + dpixydppy*dppydqx ) );
	const Numeric dEdqx( f * ( drxdqx + drydqx ) );

	// dE / dqy
	const Numeric dPxdqy( -v[2]*sinth + v[0]*v[1]*( 1.0 - costh ) );
	const Numeric dPydqy( costh * ( 1.0 - v[1]*v[1] ) + v[1]*v[1] );
	const Numeric dPzdqy( v[0]*sinth + v[1]*v[2]*( 1.0 - costh ) );

	const Numeric dppxdqy( ( P[0]*dPzdqy - P[2]*dPxdqy ) / P22 );
	const Numeric dppydqy( ( P[1]*dPzdqy - P[2]*dPydqy ) / P22 );
	const Numeric drxdqy( dEdpx*( dpixxdppx*dppxdqy + dpixxdppy*dppydqy ) );
	const Numeric drydqy( dEdpy*( dpixydppx*dppxdqy + dpixydppy*dppydqy ) );
	const Numeric dEdqy( f * ( drxdqy + drydqy ) );

	// dE / dqz
	const Numeric dPxdqz( v[1]*sinth + v[0]*v[2]*( 1.0 - costh ) );
	const Numeric dPydqz( -v[0]*sinth + v[1]*v[2]*( 1.0 - costh ) );
	const Numeric dPzdqz( costh * ( 1.0 - v[2]*v[2] ) + v[2]*v[2] );

	const Numeric dppxdqz( ( P[0]*dPzdqz - P[2]*dPxdqz ) / P22 );
	const Numeric dppydqz( ( P[1]*dPzdqz - P[2]*dPydqz ) / P22 );
	const Numeric drxdqz( dEdpx*( dpixxdppx*dppxdqz + dpixxdppy*dppydqz ) );
	const Numeric drydqz( dEdpy*( dpixydppx*dppxdqz + dpixydppy*dppydqz ) );
	const Numeric dEdqz( f * ( drxdqz + drydqz ) );

	// dE / dtx
//	const Numeric dEdtx( ( tmp1*( dEdpx*pp00 + dEdpy*pp01 ) + dEdpx*dstn ) * -f / P[2] );
	const Numeric dEdtx( ( dEdpx*dpixxdppx + dEdpy*dpixydppx ) * -f / P[2] );

	// dE / dty
//	const Numeric dEdty( ( tmp1*( dEdpy*pp11 + dEdpx*pp01 ) + dEdpy*dstn ) * -f / P[2] );
	const Numeric dEdty( ( dEdpx*dpixxdppy + dEdpy*dpixydppy ) * -f / P[2] );

	// dE / dtz
	const Numeric dpxdtz( dpixxdppx*P[0] + dpixxdppy*P[1] );
	const Numeric dpydtz( dpixydppx*P[0] + dpixydppy*P[1] );
	const Numeric dEdtz( ( dEdpx*dpxdtz + dEdpy*dpydtz ) * f / P22 );

	// dE / df
	const Numeric dpxdf( dstn * pp[0] ), dpydf( dstn * pp[1] );
	const Numeric dEdf( dEdpx*dpxdf + dEdpy*dpydf );

	// dE / dk1
	const Numeric dpxdk1( f * r2 * pp[0] ), dpydk1( f * r2 * pp[1] );
	const Numeric dEdk1( dEdpx*dpxdk1 + dEdpy*dpydk1 );

	// dE / dk2
	const Numeric dpxdk2( f * r2*r2 * pp[0] ), dpydk2( f * r2*r2 * pp[1] );
	const Numeric dEdk2( dEdpx*dpxdk2 + dEdpy*dpydk2 );

	grad[CAM_ROT_X] = dEdvx;
	grad[CAM_ROT_Y] = dEdvy;
	grad[CAM_ROT_Z] = dEdvz;
#ifdef SEPARATE_THETA_VAR
	grad[CAM_THETA] = dEdtheta;
#endif

	grad[CAM_TRANS_X] = dEdtx;
	grad[CAM_TRANS_Y] = dEdty;
	grad[CAM_TRANS_Z] = dEdtz;

	grad[CAM_FOCAL] = dEdf;
	grad[CAM_RDL_K1] = dEdk1;
	grad[CAM_RDL_K2] = dEdk2;

	grad[POINT_X] = dEdqx;
	grad[POINT_Y] = dEdqy;
	grad[POINT_Z] = dEdqz;

	return feval;
}



void BundleAdjustmentFactor::print( std::ostream & os ) const {
	os << "factor " << id << ", camera " << m_cameraID << ", point " <<
			m_pointID << ", obs: (" << m_obsX << ", " << m_obsY << ")" << std::endl;
	const VariablePtrVec & vars( bundleAdjVars );
	os << "\t[" << vars[CAM_ROT_X]->getID() << ", " << vars[CAM_ROT_Y]->getID()
			<< ", " << vars[CAM_ROT_Z]->getID() << "], [";
	os << vars[CAM_TRANS_X]->getID() << ", " << vars[CAM_TRANS_Y]->getID()
			<< ", " << vars[CAM_TRANS_Z]->getID() << "], ";
	os << vars[CAM_FOCAL]->getID() << ", (" << vars[CAM_RDL_K1]->getID() << ", "
			<< vars[CAM_RDL_K2]->getID() << ")" << std::endl;
}

} // namespace rdis
