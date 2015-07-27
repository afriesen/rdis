/*
 * BundleAdjustmentFunction.h
 *
 *  Created on: Apr 10, 2013
 *      Author: afriesen
 */

#ifndef BUNDLEADJUSTMENTFUNCTION_H_
#define BUNDLEADJUSTMENTFUNCTION_H_

#include "common.h"
#include "BundleAdjustmentCommon.h"
#include "OptimizableFunction.h"

namespace rdis {

class BundleAdjustmentFunction
		: public rdis::OptimizableFunction {

	typedef OptimizableFunction super;

	typedef std::vector< string > SplitVector;

protected:
	MinSumSemiring< Numeric > m_semiringNum;
	MinSumSemiring< NumericInterval > m_semiringBound;
public:
	virtual const Semiring< Numeric > & semiring() const { return m_semiringNum; }
	virtual const Semiring< NumericInterval > & semiringB() const { return m_semiringBound; }

public:
	BundleAdjustmentFunction();
	BundleAdjustmentFunction( const string & filename );
	virtual ~BundleAdjustmentFunction();

public:
	// load/save this bundle adjustment problem from/to a file (can specify the
	// number of cameras and points to load)
	virtual bool load( const string & file, VariableCount numcams,
			VariableCount numpoints );
	virtual bool save( const string &file, State * state = NULL );

	virtual bool load( const string & file ) {
		return load( file, 0, 0 );
	}

//	// print a representation of this to an output stream
//	virtual void print( std::ostream & os ) const;

public:
	inline virtual VariableCount getNumBlocks() const {
		return ( m_numCameras + m_numPoints );
	}

	inline virtual void getBlockRangeByBlkId( VariableCount blockid,
			VariableID & vidlower, VariableID & vidupper ) const {
		if ( blockid < m_numCameras ) {
			vidlower = getCamVID( blockid, 0 );
			vidupper = getCamVID( blockid, m_numCameraParams - 1 );
		} else {
			PointID pid = blockid - m_numCameras;
			vidlower = getPointVID( pid, 0 );
			vidupper = getPointVID( pid, m_numPointParams - 1 );
		}
	}

	inline virtual void getBlockRangeByVid( VariableID vid,
			VariableID & vidlower, VariableID & vidupper ) const {
		getBlockRangeByBlkId( getBlockID( vid ), vidlower, vidupper );
	}

	inline virtual VariableCount getBlockID( VariableID vid ) const {
		if ( vid < m_numCameras * m_numCameraParams ) {
			return getCameraID( vid );
		} else {
			return getPointID( vid ) + m_numCameras;
		}
	}

public:
	inline long long int getNumCameras() const { return m_numCameras; }
	inline long long int getNumPoints() const { return m_numPoints; }

	inline int getNumCameraParams() const { return m_numCameraParams; }
	inline int getNumPointParams() const { return m_numPointParams; }

	// get the variable ID  corresponding to this parameter of this camera
	inline VariableID getCamVID( CameraID cid, unsigned int camParam ) const {
		return cid * m_numCameraParams + camParam;
	}

	// get the variable ID corresponding to this parameter of this point
	inline VariableID getPointVID( PointID pid, unsigned int pointParam ) const {
		return m_numCameras*m_numCameraParams +
				pid*m_numPointParams + pointParam;
	}

	// get the type (BundleAdjVarType) of the specified variable
	BundleAdjust::VariableType getVarType( VariableID vid ) const {
		const long int numCamVars( m_numCameras * m_numCameraParams );
		long int type( vid );
		if ( type < numCamVars ) {
			type = type % (long int) m_numCameraParams;
		} else {
			type = BundleAdjust::POINT_X +
					( type - numCamVars ) % (long int) m_numPointParams;
		}
		return (BundleAdjust::VariableType) type;
	}

	// get the camera ID for the camera that this variable is in
	CameraID getCameraID( const VariableID & vid ) const {
		const long int numCamVars( m_numCameras * m_numCameraParams );
		assert( vid < numCamVars );
		return ( vid / m_numCameraParams );
	}

	// get the point ID for the point that this variable is in
	PointID getPointID( const VariableID & vid ) const {
		const long int numCamVars( m_numCameras * m_numCameraParams );
		assert( vid >= numCamVars );
		return ( ( vid - numCamVars ) / m_numPointParams );
	}

public:
	void graph( const std::string & fname = "BA_function.gv" );

	const State & getInitialState() const { return xinit; }

protected:
	void readHeader( BOOSTNS::filesystem::fstream & ifs, long long int & numCameras,
			long long int & numPoints, long long int & numObservations );

	void readObservation( BOOSTNS::filesystem::fstream & ifs, CameraID & cid,
			PointID & pid, Numeric & obsX, Numeric & obsY );

	void readCameraParams( BOOSTNS::filesystem::fstream & ifs, Numeric * camera );
	void readPoint( BOOSTNS::filesystem::fstream & ifs, Numeric * point );

	// set the domain of this variable based on its initial value
	void setDomain( const VariableID vid, const Numeric initialVal );

protected:
	  // Move the "center" of the reconstruction to the origin, where the
	  // center is determined by computing the marginal median of the
	  // points. The reconstruction is then scaled so that the median
	  // absolute deviation of the points measured from the origin is
	  // 100.0.
	  //
	  // The reprojection error of the problem remains the same.
	  // (from Sameer's Ceres code)
	void normalizeVariables();

	// TODO: use a proper vector library...
	typedef Numeric Vector3D[3];

	Numeric median( std::vector< Numeric > & data ) const;
	Numeric sum( Vector3D v ) const { return v[0] + v[1] + v[2]; }

	inline void getPoint( PointID /*pid*/, Vector3D /*pt*/ ) const {
		assert( false );
//		for ( int i = 0; i < 3; ++i ) {
//			pt[i] = variables[getPointVID( pid, i )]->getDomain().interval().lower();
//		}
	}

protected:
	const int m_numCameraParams;
	const int m_numPointParams;

	long long int m_numCameras;
	long long int m_numPoints;
	long long int m_numObservations;

	// initial state as loaded from the file
	State xinit;

};

//typedef BOOSTNS::shared_ptr< BundleAdjustmentFunction >
//		BundleAdjustmentFunctionSP;

} // namespace rdis

#endif // BUNDLEADJUSTMENTFUNCTION_H_
