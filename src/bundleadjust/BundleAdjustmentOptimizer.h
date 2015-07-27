/*
 * BundleAdjustmentOptimizer.h
 *
 *  Created on: May 13, 2013
 *      Author: afriesen
 */

#ifndef BUNDLEADJUSTMENTOPTIMIZER_H_
#define BUNDLEADJUSTMENTOPTIMIZER_H_

#include "common.h"
#include "BundleAdjustmentCommon.h"
#include "bundleadjust/BundleAdjustmentFunction.h"
#include "RDISOptimizer.h"

#include <vector>

namespace rdis {


typedef std::vector< CameraID > CameraIDVec;
typedef std::vector< PointID > PointIDVec;

typedef std::vector< CameraIDVec > CamIDVecVec;
typedef std::vector< PointIDVec > PointIDVecVec;

//// class to optimize BundleAdjustmentFunctions -- main RDISOptimizer works well, but
//// this has specific heuristics, etc, for solving Bundle Adjustment
//class BundleAdjustmentOptimizer
//	: public RDISOptimizer {
//
//typedef RDISOptimizer super;
//
//public:
//	BundleAdjustmentOptimizer( BundleAdjustmentFunction & baf,
//			Numeric epsilon = 0.1, Numeric globalLipschitz = 20,
//			size_t cacheSize = 10000 );
//	virtual ~BundleAdjustmentOptimizer();
//
//
//protected:
//
//	virtual bool assignAllVariables();
//
//	void getCamOrPointToAssign( CameraID & cid, PointID & pid );
//	bool assignCamera( const CameraID & cid );
//	bool assignPoint( const PointID & pid );
//
//	virtual bool markDescendantsUnassigned( const VariableID & vid );
//
//protected:
//	virtual void init( size_t cacheSize );
//	virtual void prepareToOptimize( bool initVDNbrs = true );
//
//protected:
//
//	BundleAdjustmentFunction & m_BAfunc;
//
//	// lists of the unassigned cameras/points
//	CameraIDVec unassignedCams;
//	PointIDVec unassignedPts;
//
//	// list of the connected points (to each camera) and the connected cameras
//	// (to each point)
//	PointIDVecVec connectedPoints;
//	CamIDVecVec connectedCameras;
//
//	// a count of the number of connected (and unassigned) points, per camera
//	// (indexed by camera ID)
//	CountVec numNeighbourPoints;
//
//	// a count of the number of connected (and unassigned) cameras, per point
//	// (indexed by point ID)
//	CountVec numNeighbourCameras;
//
//};

} // namespace rdis
#endif // BUNDLEADJUSTMENTOPTIMIZER_H_
