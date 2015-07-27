/*
 * BundleAdjustmentFunction.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: afriesen
 */

#include "bundleadjust/BundleAdjustmentFunction.h"

#include "bundleadjust/BundleAdjustmentFactor.h"

#include BOOSTPATH/format.hpp>
#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/lexical_cast.hpp>
#include BOOSTPATH/algorithm/string.hpp>
#include BOOSTPATH/math/constants/constants.hpp>

#include BOOSTPATH/graph/graphviz.hpp>

#include <algorithm>

namespace fs = BOOSTNS::filesystem;


namespace rdis {

BundleAdjustmentFunction::BundleAdjustmentFunction()
	: super( true, true, false, true )
	, m_numCameraParams( ( (int) BundleAdjust::NUMPARAMS ) - 3 )
	, m_numPointParams( 3 )
	, m_numCameras( 0 )
	, m_numPoints( 0 )
	, m_numObservations( 0 )
{}

BundleAdjustmentFunction::BundleAdjustmentFunction( const string & filename )
		: super( true, true, false, true )
		, m_numCameraParams( ( (int) BundleAdjust::NUMPARAMS ) - 3 )
		, m_numPointParams( 3 )
		, m_numCameras( 0 )
		, m_numPoints( 0 )
		, m_numObservations( 0 ) {
	bool loaded( load( filename ) );
	assert( loaded );
}

BundleAdjustmentFunction::~BundleAdjustmentFunction() {}


bool BundleAdjustmentFunction::load( const string & file_in,
		VariableCount numcams, VariableCount numpoints ) {
	// open the specified file
	fs::fstream ifs( file_in, std::ios_base::in );
	if ( !ifs.is_open() ) {
		std::cerr << "BundleAdjustmentFunction::Load: Failed to open file for reading " <<
				file_in << std::endl;
		return false;
	}

	VariableID unique_var_id = 0;
	string line;
//	bool result = true;

	using namespace BOOSTNS;
	typedef BundleAdjust::VariableType BAVarType;

//	namespace bmc = BOOSTNS::math::constants;
//	const Numeric maxPosVal( 1e4 ), posSF( 1e3 );
//	const NumericInterval maxThetaDom( -bmc::pi<Numeric>(), bmc::pi<Numeric>() );
//#ifndef SEPARATE_THETA_VAR
//	const NumericInterval maxVDom( maxThetaDom );
//#else
//	const NumericInterval maxVDom( -1, 1 );
//#endif
//	const NumericInterval maxPosDom( -maxPosVal, maxPosVal );
//
//	const NumericInterval sampleDomains[BundleAdjust::NUMPARAMS] = {
//			maxVDom, maxVDom, maxVDom, 			// camera rotation params
//#ifdef SEPARATE_THETA_VAR
//			maxThetaDom,
//#endif // SEPARATE_THETA_VAR
//			maxPosDom, maxPosDom, maxPosDom, 	// camera translation params
//			maxPosDom,					 		// camera focal length
//			NumericInterval( -1e-3, 1e-2 ),		// camera distortion param 1
//			NumericInterval( -1e-5, 1e-4 ), 	// camera distortion param 2
//			maxPosDom, maxPosDom, maxPosDom		// point position variables
//	};
//
//	const NumericInterval fullDomains[BundleAdjust::NUMPARAMS] = {
//			maxVDom*10, maxVDom*10, maxVDom*10, // camera rotation params
//#ifdef SEPARATE_THETA_VAR
//			maxThetaDom*10,
//#endif // SEPARATE_THETA_VAR
//			maxPosDom*posSF, maxPosDom*posSF, maxPosDom*posSF, // camera translation params
//			NumericInterval( 0, maxPosVal*posSF ), // camera focal length
//			NumericInterval( -1e-2, 1e-2 ),		// camera distortion param 1
//			NumericInterval( -1e-2, 1e-2 ), 	// camera distortion param 2
//			maxPosDom, maxPosDom, maxPosDom		// point position variables
//	};
//
//	const bool sampleDomIsRelative[BundleAdjust::NUMPARAMS] = {
//			false, false, false,		// camera rotation params
//#ifdef SEPARATE_THETA_VAR
//			false,
//#endif // SEPARATE_THETA_VAR
//			true, true, true,			// camera translation params
//			true, false, false,			// camera focal length and distortions
//			true, true, true			// point position variables
//	};


	try {
		long long int numCameras( 0 ), numPoints( 0 );
		CameraID camID( 0 );
		PointID pointID( 0 );
		Numeric obsX( 0 ), obsY( 0 );

		readHeader( ifs, numCameras, numPoints, m_numObservations );

		m_numCameras = ( numcams <= 0 ? numCameras : numcams );
		m_numPoints  = ( numpoints <= 0 ? numPoints : numpoints );
//		m_numCameras = numCameras;
//		m_numPoints  = numPoints;
		assert( numCameras >= m_numCameras );
		assert( numPoints >= m_numPoints );

		// create all of the camera variables
		format camfmt( "C%1%.%2%" );
		for ( CameraID cid = 0; cid < numCameras; ++cid ) {
			if ( cid >= m_numCameras ) break;
			for ( int cpi = 0; cpi < m_numCameraParams; ++cpi ) {
				assert( ( getVarType( unique_var_id ) ) == (BundleAdjust::VariableType) cpi );
				addVariable( str( camfmt % cid % cpi ), "0:0", unique_var_id );
			}
		}

		// create all of the point variables
		format ptfmt( "x%1%.%2%" );
		for ( PointID pid = 0; pid < numPoints; ++pid ) {
			if ( pid >= m_numPoints ) break;
			for ( int ppi = 0; ppi < m_numPointParams; ++ppi ) {
				assert( unique_var_id == getPointVID( pid, ppi ) );
				assert( ( getVarType( unique_var_id ) ) ==
						BundleAdjust::POINT_X + (BundleAdjust::VariableType) ppi );

				addVariable( str( ptfmt % pid % (char) ( 'x' + ppi ) ),
						"0:0", unique_var_id );
			}
		}

		// create all of the factors
		for ( long long int i = 0; i < m_numObservations; ++i ) {
			readObservation( ifs, camID, pointID, obsX, obsY );

			if ( camID >= m_numCameras || pointID >= m_numPoints ) {
				continue;
			}

//			std::cout << "factor " << factors.size() << ": cam " << camID <<
//					", point " << pointID << " -- obs (" << obsX << ", " << obsY
//					<< ")" << std::endl;

			// create the corresponding factor
			BundleAdjustmentFactor * f = new BundleAdjustmentFactor(
					factors.size(), camID, pointID, obsX, obsY );
			factors.push_back( f );

			// add all of the camera variables
			for ( int cpi = 0; cpi < m_numCameraParams; ++cpi ) {
				f->addVariable( variables[getCamVID( camID, cpi )],
						(BAVarType) cpi );
			}

			// add all of the point variables
			for ( int ppi = 0; ppi < m_numPointParams; ++ppi ) {
				f->addVariable( variables[getPointVID( pointID, ppi )],
						(BAVarType) ( (int) BundleAdjust::POINT_X + ppi ) );
			}
		}

		xinit.resize( variables.size(), NumericMAX );

//		format domainfmt( "%1%:%1%" );
		Numeric cameraParams[m_numCameraParams];
		Numeric point[m_numPointParams];

		for ( CameraID cid = 0; cid < numCameras; ++cid ) {
			readCameraParams( ifs, cameraParams );

			if ( cid >= m_numCameras ) continue;

			// set the initial values for the camera parameters
			for ( int cpi = 0; cpi < m_numCameraParams; ++cpi ) {

				VariableID vid( getCamVID( cid, cpi ) );
////				variables[vid]->setDomain( cameraParams[cpi] );
//////				variables[vid]->setMaxDomain( maxDomains[cpi] );
//				NumericInterval itvl = cameraParams[cpi] + maxDomains[cpi];
//				variables[vid]->setDomain( VariableDomain( maxDomains[cpi].lower(), maxDomains[cpi].upper() ) );

				xinit[vid] = cameraParams[cpi];
				setDomain( vid, cameraParams[cpi] );

//				std::cout << "cam " << cid << "." << cpi << ": " <<
//						cameraParams[cpi] << std::endl;
			}
		}

		for ( PointID pid = 0; pid < numPoints; ++pid ) {
			readPoint( ifs, point );

			if ( pid >= m_numPoints ) continue;

			// set the initial values for the point parameters
			for ( int ppi = 0; ppi < m_numPointParams; ++ppi ) {
				VariableID vid( getPointVID( pid, ppi ) );
////				variables[vid]->setDomain( point[ppi] );
//////				variables[vid]->setMaxDomain(
//////						maxDomains[(int) BundleAdjust::POINT_X + ppi] );
//				NumericInterval itvl = point[ppi] +
//						maxDomains[(int) BundleAdjust::POINT_X + ppi];
//				variables[vid]->setDomain(
//						VariableDomain( itvl.lower(), itvl.upper() ) );

				xinit[vid] = point[ppi];
				setDomain( vid, point[ppi] );

//				std::cout << "point " << pid << "." << ppi << ": " <<
//						point[ppi] << std::endl;
			}
		}

	} catch ( std::exception & e ) {
		// catch anything
		std::cerr << "PolynomialFunction::load(" << file_in << ") caught exception: "
				<< std::endl << e.what() << std::endl;
//		result = false;
	}

	// close file
	if ( ifs.is_open() ) ifs.close();

	init();

//	for ( size_t i = 0; i < variables.size(); ++i ) {
//		std::cout << "var domain: " << variables[i]->getDomain().interval << std::endl;
//	}

	return ( variables.size() > 0 && factors.size() > 0 );
}


bool BundleAdjustmentFunction::save( const string & file_out,
		State * state ) {
	fs::fstream ofs( file_out, std::ios_base::out );
	if ( !ofs.is_open() ) {
		std::cerr << "BundleAdjustmentFunction::Load: Failed to open file for writing "
				<< file_out << std::endl;
		return false;
	}

	ofs << m_numCameras << " " << m_numPoints << " " << factors.size() << std::endl;

	// write out each <camID, pointID, observation> tuple
	for ( Factor * f : factors ) {
		const BundleAdjustmentFactor * baf(
				dynamic_cast< BundleAdjustmentFactor * >( f ) );
		if ( baf == NULL ) continue;
		Numeric obsx, obsy;
		baf->getObservation( obsx, obsy );
		ofs << baf->getCameraID() << " " << baf->getPointID() << " " << obsx
				<< " " << obsy << std::endl;
	}

	using BOOSTNS::numeric::median;

	// write out each camera variable value
	for ( CameraID cid = 0; cid < m_numCameras; ++cid ) {
		for ( int cpi = 0; cpi < m_numCameraParams; ++cpi ) {
			VariableID vid( getCamVID( cid, cpi ) );
			Numeric val = ( state != NULL ) ? (*state)[vid] :
						  variables[vid]->getDomain().median();

#ifdef SEPARATE_THETA_VAR
			switch( getVarType( vid ) ) {
			case BundleAdjust::CAM_ROT_X:
			case BundleAdjust::CAM_ROT_Y:
			case BundleAdjust::CAM_ROT_Z: {
				// write these out in Rodrigues' representation
				VariableID thetavid( getCamVID( cid, BundleAdjust::CAM_THETA ) );
				Numeric theta( state != NULL ? (*state)[thetavid] :
						variables[thetavid]->getDomain().median() );
				val *= theta;
				break;
			}
			case BundleAdjust::CAM_THETA:
				// don't write out THETA
				continue;
			default:
				// do nothing
				break;
			}
#endif // SEPARATE_THETA_VAR

			ofs << val << std::endl;
		}
	}

	// write out each point dimension value
	for ( PointID pid = 0; pid < m_numPoints; ++pid ) {
		for ( int ppi = 0; ppi < m_numPointParams; ++ppi ) {
			VariableID vid( getPointVID( pid, ppi ) );
			Numeric val = ( state != NULL ) ? (*state)[vid] :
						  variables[vid]->getDomain().median();
			ofs << val << std::endl;
		}
	}

	return true;
}


void BundleAdjustmentFunction::readHeader( fs::fstream & ifs, long long int & numCams,
		long long int & numPoints, long long int & numObs ) {

	SplitVector sv;
	bool result( readAndSplit( ifs, sv ) );

	assert( result );
	assert( sv.size() == 3 );

	numCams = BOOSTNS::lexical_cast< long long int >( sv[0] );
	numPoints = BOOSTNS::lexical_cast< long long int >( sv[1] );
	numObs = BOOSTNS::lexical_cast< long long int >( sv[2] );
}


void BundleAdjustmentFunction::readObservation( fs::fstream & ifs,
		CameraID & cid, PointID & pid, Numeric & obsX, Numeric & obsY ) {

	SplitVector sv;
	bool result( readAndSplit( ifs, sv ) );

	assert( result );
	assert( sv.size() == 4 );

	cid = BOOSTNS::lexical_cast< CameraID >( sv[0] );
	pid = BOOSTNS::lexical_cast< PointID >( sv[1] );
	obsX = BOOSTNS::lexical_cast< Numeric >( sv[2] );
	obsY = BOOSTNS::lexical_cast< Numeric >( sv[3] );
}


void BundleAdjustmentFunction::readCameraParams( fs::fstream & ifs,
		Numeric * camera ) {

	using namespace BundleAdjust;
	SplitVector sv;

	for ( int cpi = 0; cpi < m_numCameraParams; ++cpi ) {
#ifdef SEPARATE_THETA_VAR
		if ( cpi == (int) CAM_THETA ) {
			Numeric theta2( dot( &camera[CAM_ROT_X], &camera[CAM_ROT_X] ) );

			assert( theta2 == ( camera[CAM_ROT_X]*camera[CAM_ROT_X] +
					camera[CAM_ROT_Y]*camera[CAM_ROT_Y] +
					camera[CAM_ROT_Z]*camera[CAM_ROT_Z] ) );

			Numeric theta( std::sqrt( theta2 ) );
			camera[CAM_ROT_X] /= theta;
			camera[CAM_ROT_Y] /= theta;
			camera[CAM_ROT_Z] /= theta;
			camera[CAM_THETA] = theta;
			continue;
		}
#endif // SEPARATE_THETA_VAR

		bool result( readAndSplit( ifs, sv ) );
		if ( !result || sv.size() != 1 ) {
			std::cout << "readCameraParams failed: "; printContainer( sv );
			std::cout << std::endl;
		}
		assert( result );
		assert( sv.size() == 1 );
		camera[cpi] = BOOSTNS::lexical_cast< Numeric >( sv[0] );
	}
}


void BundleAdjustmentFunction::readPoint( fs::fstream & ifs, Numeric * point ) {

	SplitVector sv;
	for ( int ppi = 0; ppi < m_numPointParams; ++ppi ) {
		bool result( readAndSplit( ifs, sv ) );
		assert( result );
		assert( sv.size() == 1 );
		point[ppi] = BOOSTNS::lexical_cast< Numeric >( sv[0] );
	}
}


void BundleAdjustmentFunction::setDomain( const VariableID vid,
		const Numeric initialVal ) {

	BundleAdjust::VariableType type = getVarType( vid );
	Variable * v = variables[vid];

	namespace bmc = BOOSTNS::math::constants;

	NumericInterval dom( NumericMAX ), sit( NumericMAX );
	const Numeric dsf = 1000.0;

#ifdef SEPARATE_THETA_VAR
	const Numeric rotscale = 1;
#else
	const Numeric rotscale = bmc::pi< Numeric >();
#endif // SEPARATE_THETA_VAR

	switch( type ) {
	case BundleAdjust::CAM_ROT_X:
	case BundleAdjust::CAM_ROT_Y:
	case BundleAdjust::CAM_ROT_Z:
		sit = NumericInterval( -1, 1 ) * rotscale;
		dom = sit * dsf;
		break;

	case BundleAdjust::CAM_TRANS_X:
	case BundleAdjust::CAM_TRANS_Y:
	case BundleAdjust::CAM_TRANS_Z:
	case BundleAdjust::POINT_X:
	case BundleAdjust::POINT_Y:
	case BundleAdjust::POINT_Z:
		sit = initialVal + NumericInterval( -1, 1 )*1e2;
		dom = sit * dsf;
		break;

#ifdef SEPARATE_THETA_VAR
	case BundleAdjust::CAM_THETA:
		sit = NumericInterval( -1, 1 )*bmc::pi<Numeric>();
		dom = sit * dsf;
		break;
#endif // SEPARATE_THETA_VAR

	case BundleAdjust::CAM_FOCAL:
	{
		sit = initialVal + NumericInterval( -1, 1 )*1e2;
		Numeric lower = std::min( sit.lower(), sit.lower()*dsf );
		lower = std::max( lower, 0.0 );
		dom.assign( lower, sit.upper()*dsf );
		break;
	}

	case BundleAdjust::CAM_RDL_K1:
		sit = initialVal + NumericInterval( -1, 1 )*1e-4;
		sit.assign( sit.lower(), sit.upper() );
		dom.assign( -1e-1, 1e-1 );
		break;
	case BundleAdjust::CAM_RDL_K2:
		sit = initialVal + NumericInterval( -1, 1 )*1e-6;
		sit.assign( sit.lower(), sit.upper() );
		dom.assign( -1e-3, 1e-3 );
		break;

	default:
		assert( false );
	}

	// ensure overlap between the two intervals
	dom = BOOSTNS::numeric::hull( dom, sit );

	v->setDomain( VariableDomain( dom.lower(), dom.upper()));
	v->setSamplingInterval( sit );
//	std::cout << "sampling interval: " << sit << " (init: " << initialVal << std::endl;

	assert( !std::isnan( sit.lower() ) && !std::isnan( sit.upper() ) );
	assert( !std::isnan( dom.lower() ) && !std::isnan( dom.upper() ) );
}


void BundleAdjustmentFunction::normalizeVariables() {

	assert( false );

	std::vector< Numeric > tmp( m_numPoints );
	Numeric med[3];

	// compute the marginal median of the geometry
	for ( int i = 0; i < m_numPointParams; ++i ) {
		for ( PointID j = 0; j < m_numPoints; ++j ) {
			assert( BOOSTNS::numeric::singleton(
					variables[getPointVID( j, i )]->getDomain().interval() ) );
			tmp[j] = variables[getPointVID( j, i )]->getDomain().min();
		}
		med[i] = median( tmp );
	}

	for ( PointID i = 0; i < m_numPoints; ++i ) {
		Numeric point[3];
		getPoint( i, point );
		tmp[i] = 0;
		for ( int j = 0; j < m_numPointParams; ++j ) {
			tmp[i] += ( point[j] - med[j] );
		}
	}


	const double median_abs_deviation = median( tmp );

//	  // scale so that the median absolute deviation of the resulting
//	  // reconstruction is 100
	const double scale = 100.0 / median_abs_deviation;


	std::cerr << "median: (" << med[0] << ", " << med[1] << ", " <<
			med[2] << ")" << std::endl;
	std::cerr << "median absolute deviation: " << median_abs_deviation << std::endl;
	std::cerr << "scale: " << scale;

	// x = scale * (x - median)
	for ( PointID pid = 0; pid < m_numPoints; ++pid ) {
		Numeric point[3];
		getPoint( pid, point );
		for ( int j = 0; j < m_numPointParams; ++j ) {
			variables[getPointVID( pid, j )]->setDomain(
					scale * ( point[j] - med[j] ));
		}
	}

	// TODO: finish normalization function
	assert( false );

//	  double* cameras = mutable_cameras();
//	  double angle_axis[3];
//	  double center[3];
//	  for (int i = 0; i < num_cameras_; ++i) {
//	    double* camera = cameras + camera_block_size() * i;
//	    CameraToAngleAxisAndCenter(camera, angle_axis, center);
//	    // center = scale * (center - median)
//	    VectorRef(center, 3) = scale * (VectorRef(center, 3) - median);
//	    AngleAxisAndCenterToCamera(angle_axis, center, camera);
//	  }
//	}
}

Numeric BundleAdjustmentFunction::median( std::vector< Numeric > & data ) const {
	std::vector< Numeric >::iterator mid( data.begin() + data.size() / 2 );
	std::nth_element( data.begin(), mid, data.end() );
	return *mid;
}


void BundleAdjustmentFunction::graph( const std::string & fname ) {
	typedef BOOSTNS::adjacency_list< BOOSTNS::setS, BOOSTNS::vecS, BOOSTNS::undirectedS > Graph;
	fs::fstream ofs( fname, std::ios_base::out );
	if ( ofs.is_open() ) {
		Graph g( m_numCameras + m_numPoints );
		BundleAdjustmentFactor * baf( NULL );

		for ( size_t i = 0; i < factors.size(); ++i ) {
			baf = dynamic_cast< BundleAdjustmentFactor * >( factors[i] );
			if ( baf == NULL ) continue;
			BOOSTNS::add_edge( baf->getCameraID(),
					m_numCameras + baf->getPointID(), g );
		}

		BOOSTNS::write_graphviz( ofs, g );
		ofs.close();
	}
}


} // namespace rdis
