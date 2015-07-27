/*
 * BundleAdjustmentFactor.h
 *
 *  Created on: Apr 6, 2013
 *      Author: afriesen
 */

#ifndef RDIS_BUNDLEADJUSTMENTFACTOR_H_
#define RDIS_BUNDLEADJUSTMENTFACTOR_H_

#include "common.h"
#include "BundleAdjustmentCommon.h"
#include "Factor.h"

namespace rdis {

class BundleAdjustmentFactor : public Factor {

	typedef Factor super;

	typedef std::vector< NumericInterval > IntervalVec;
	typedef std::vector< Numeric > NumericVec;

public:
	BundleAdjustmentFactor( FactorID factorID, CameraID cameraID,
			PointID pointID, Numeric obsX, Numeric obsY );
	virtual ~BundleAdjustmentFactor();

	// add the specified variable to the set of variables in this factor
	void addVariable( Variable * vp, BundleAdjust::VariableType index );
	inline void addVariable( Variable * /*vp*/ ) { assert( false ); }

public:
	virtual Numeric evalNoCache() const;

protected:
	virtual NumericInterval computeFactorBounds(
			const VariableIDVec & vidsToIgnore = VariableIDVec() ) const;

public:
	virtual void computeGradient( PartialGradient & g, bool doGradCheck ) const;


public:
	CameraID getCameraID() const { return m_cameraID; }

	PointID getPointID() const { return m_pointID; }

	void getObservation( Numeric & obsX, Numeric & obsY ) const {
		obsX = m_obsX;
		obsY = m_obsY;
	}

	void print( std::ostream & os ) const;

protected:
	// evaluate the factor
	Numeric evalFactor() const;

	void getVarVals( IntervalVec & vals, bool useFullDomains,
			bool ignoreSpecificVarAssignments,
			const VariableIDVec & vidsToIgnore, bool useAssigned ) const;

	void getVarVals( NumericVec & vals ) const;

protected:

	NumericInterval evalFactor( const IntervalVec & vals ) const;
	Numeric evalFactor( const NumericVec & vals ) const;

//	void evalResiduals( const IntervalVec & vals, NumericInterval res[2] ) const;
	void evalPixelVals( const NumericVec & vals, Numeric v[3], Numeric & theta,
			Numeric P[3], Numeric pp[2], Numeric & r2, Numeric & dstn,
			Numeric pix[2] ) const;

	void angleAxisRotatePoint( const IntervalVec & vals, NumericInterval p[3] ) const;
	void angleAxisRotateTranslatePoint( const NumericVec & vals, Numeric v[3],
			Numeric & theta, Numeric vxp[3], Numeric & vdp, Numeric P[3] ) const;

	inline void perspectiveDivide( const Numeric P[3], Numeric pp[2] ) const {
		pp[0] = -P[0] / P[2];
		pp[1] = -P[1] / P[2];
	}

	inline void distort( const NumericVec & vals, const Numeric pp[2],
			Numeric & r2, Numeric & dstn, Numeric pix[2] ) const {
		using namespace BundleAdjust;

		// r^2 = px^2 + py^2, distortion = ( 1 + k1*r^2 + k2*r^4 )
		r2 = pp[0]*pp[0] + pp[1]*pp[1];
		dstn = 1 + r2*( vals[CAM_RDL_K1] + vals[CAM_RDL_K2]*r2 );

		pix[0] = vals[CAM_FOCAL] * dstn * pp[0];
		pix[1] = vals[CAM_FOCAL] * dstn * pp[1];
	}

	inline Numeric error( const Numeric pix[2], Numeric res[2] ) const {
		// error = (px - obsx)^2 + (py - obsy)^2, so compute the residuals
		res[0] = ( pix[0] - m_obsX );
		res[1] = ( pix[1] - m_obsY);

		const Numeric error( ( res[0]*res[0] + res[1]*res[1] ) / 2.0 );
//	std::cout << "factor " << id << ": " << error << " (" << residx << ", " <<
//			residy << ")" << std::endl;

		return error;
	}

//	void angleSphericalAxisRotatePoint( const IntervalVec & vals, NumericInterval p[3] ) const;
//	void angleAxisNormalizedRotatePoint( const IntervalVec & vals, NumericInterval p[3] ) const;


protected:
	Numeric computeGradient( const NumericVec & vals, NumericVec & grad ) const;

protected:
	// the ID of the camera that this factor is projecting the point onto
	const CameraID m_cameraID;

	// the ID of the point that this factor is projecting
	const PointID m_pointID;

	// the observed pixel values used in this factor
	const Numeric m_obsX;
	const Numeric m_obsY;

	// vector containing the variables in this factor, indexed by their
	// respective BundleAdjVarIndex
	VariablePtrVec bundleAdjVars;

protected:
	// keep this vector around so we don't have to re-allocate each time a call
	// to evalFactor is made -- never assume that this is set appropriately,
	// always reinitialize it (just preallocated mem. used wherever its needed)
	mutable IntervalVec m_tempVarIntervals;
	mutable NumericVec m_tempVarVals, m_tempGradient;
};

} // namespace rdis

#endif // RDIS_BUNDLEADJUSTMENTFACTOR_H_
