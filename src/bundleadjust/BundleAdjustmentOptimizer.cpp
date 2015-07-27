/*
 * BundleAdjustmentOptimizer.cpp
 *
 *  Created on: May 13, 2013
 *      Author: afriesen
 */

#include "BundleAdjustmentOptimizer.h"
#include "bundleadjust/BundleAdjustmentFactor.h"

namespace rdis {

//BundleAdjustmentOptimizer::BundleAdjustmentOptimizer(
//		BundleAdjustmentFunction & baf, Numeric epsilon,
//		Numeric globalLipschitz, size_t cacheSize )
//	: super( baf, epsilon, globalLipschitz, cacheSize )
//	, m_BAfunc( baf )
//{
////	init( cacheSize ); // super() call above calls this already
//}
//
//
//BundleAdjustmentOptimizer::~BundleAdjustmentOptimizer() {}
//
//
//
//bool BundleAdjustmentOptimizer::assignAllVariables() {
//	CameraID cid( -1 );
//	PointID pid( -1 );
//	bool result( true );
////	while ( !vardata.getUnassignedVars().empty() ) {
//	while ( !unassignedCams.empty() || !unassignedPts.empty() ) {
//		getCamOrPointToAssign( cid, pid );
//
//		if ( cid >= 0 ) result = assignCamera( cid );
//		else if ( pid >= 0 ) result = assignPoint( pid );
//		else result = false;
//
//		if ( !result ) return false;
//	}
//
//	if ( !vardata.getUnassignedVars().empty() ) {
//		return super::assignAllVariables();
//	}
//
//	return true;
//}
//
//
//void BundleAdjustmentOptimizer::getCamOrPointToAssign( CameraID & cid,
//		PointID & pid ) {
//
//	CameraIDVec::iterator icam( unassignedCams.end() );
//	PointIDVec::iterator ipnt( unassignedPts.end() );
//	Numeric temphv( 0 ), camMaxHVal( -std::numeric_limits< Numeric >::max() ),
//			pntMaxHVal( -std::numeric_limits< Numeric >::max() );
//	cid = pid = -1;
//
//	assert( !unassignedCams.empty() || !unassignedPts.empty() );
//
//	// find the best camera to assign
//	for ( CameraIDVec::iterator it( unassignedCams.begin() ),
//			end( unassignedCams.end() ); it != end; ++it ) {
////		temphv = getCamHeuristicVal( *it );
//		temphv = numNeighbourPoints[*it];
//		assert( !std::isnan( temphv ) );
//		if ( temphv > camMaxHVal ) {
//			camMaxHVal = temphv;
//			icam = it;
//		}
//	}
//
//	// find the best point to assign
//	for ( PointIDVec::iterator it( unassignedPts.begin() ),
//			end( unassignedPts.end() ); it != end; ++it ) {
////		temphv = getPointHeuristicVal( *it );
//		temphv = numNeighbourCameras[*it];
//		assert( !std::isnan( temphv ) );
//		if ( temphv > pntMaxHVal ) {
//			pntMaxHVal = temphv;
//			ipnt = it;
//		}
//	}
//
//	// scale to the number of variables that each is actually connected to
//	camMaxHVal *= (Numeric) m_BAfunc.getNumPointParams();
//	pntMaxHVal *= (Numeric) m_BAfunc.getNumCameraParams();
//
//	if ( camMaxHVal > pntMaxHVal ) {
//		cid = *icam;
//		unassignedCams.erase( icam );
//	} else {
//		pid = *ipnt;
//		unassignedPts.erase( ipnt );
//	}
//}
//
//
//bool BundleAdjustmentOptimizer::assignCamera( const CameraID & cid ) {
//	const int ncparams( m_BAfunc.getNumCameraParams() );
//	VariableID vid( -1 );
//	Numeric value( 0 );
//
//	std::cout << "assigning camera " << cid << ": { ";
//	for ( int cpi = 0; cpi < ncparams; ++cpi ) {
//		vid = m_BAfunc.getCamVID( cid, cpi );
//		std::cout << (cpi == 0 ? "" : ", " ) << vid;
//		assert( vardata.isUnassigned( vid ) );
//
//		if ( !getValueFromDomain( vid, value ) ) return false;
//		assign( vid, value );
//	}
//	std::cout << " }" << std::endl;
//
//	assert( false );
//
////	VariableIDVec & unassigned( vardata.getUnassignedVars() );
////	std::sort( unassigned.begin(), unassigned.end() );
////	vid = m_BAfunc.getCamVID( cid, 0 );
////	VariableIDVec::iterator it( std::lower_bound( unassigned.begin(),
////			unassigned.end(), vid ) );
////
////	assert( *it == vid );
////	assert( *(it + ncparams - 1 ) == vid + ncparams - 1 );
////	unassigned.erase( it, it + ncparams );
//
//	return true;
//}
//
//
//bool BundleAdjustmentOptimizer::assignPoint( const PointID & pid ) {
//	const int npparams( m_BAfunc.getNumPointParams() );
//	VariableID vid( -1 );
//	Numeric value( 0 );
//
//	std::cout << "assigning point " << pid << ": { ";
//	for ( int ppi = 0; ppi < npparams; ++ppi ) {
//		vid = m_BAfunc.getPointVID( pid, ppi );
//		std::cout << (ppi == 0 ? "" : ", " ) << vid;
//		assert( vardata.isUnassigned( vid ) );
//
//		if ( !getValueFromDomain( vid, value ) ) return false;
//		assign( vid, value );
//	}
//	std::cout << " }"<< std::endl;
//
//	assert( false );
////	VariableIDVec & unassigned( vardata.getUnassignedVars() );
////	std::sort( unassigned.begin(), unassigned.end() );
////	vid = m_BAfunc.getPointVID( pid, 0 );
////	VariableIDVec::iterator it( std::lower_bound( unassigned.begin(),
////			unassigned.end(), vid ) );
////
////	assert( *it == vid );
////	assert( *(it + npparams - 1) == vid + npparams - 1 );
////	unassigned.erase( it, it + npparams );
//
//	return true;
//}
//
//
//bool BundleAdjustmentOptimizer::markDescendantsUnassigned(
//		const VariableID & vid ) {
//
//	const VariableIDVec & c( vardata.get( vid ).children );
//	if ( c.empty() ) return false;
//
//	for ( VariableIDVec::const_iterator it( c.begin() ), end( c.end() );
//			it != end; ++it ) {
//		vardata.addUnassignedVar( *it );
//		m_components.unmerge( *it );
////		std::cout << "marking " << *it << " as unassigned" << std::endl;
//
//		BundleAdjust::VariableType type( m_BAfunc.getVarType( vid ) );
//
//		if ( type == BundleAdjust::CAM_ROT_X ) {
//			CameraID cid( m_BAfunc.getCameraID( vid ) );
//			assert( std::find( unassignedCams.begin(), unassignedCams.end(),
//					cid ) == unassignedCams.end() );
//			unassignedCams.push_back( cid );
//
//			// tell each connected point that it has another neighbour camera
//			for ( PointIDVec::iterator it( connectedPoints[cid].begin() ),
//					end( connectedPoints[cid].end() ); it != end; ++it ) {
//				++numNeighbourCameras[*it];
//			}
//
//		} else if ( type == BundleAdjust::POINT_X ) {
//			PointID pid( m_BAfunc.getPointID( vid ) );
//			assert( std::find( unassignedPts.begin(), unassignedPts.end(),
//					pid ) == unassignedPts.end() );
//			unassignedPts.push_back( pid );
//
//			// tell each connected camera that it has another neighbour point
//			for ( CameraIDVec::iterator it( connectedCameras[pid].begin() ),
//					end( connectedCameras[pid].end() ); it != end; ++it ) {
//				++numNeighbourPoints[*it];
//			}
//		}
//
//		markDescendantsUnassigned( *it );
//	}
//
//	return true;
//}
//
//
//void BundleAdjustmentOptimizer::init( size_t cacheSize ) {
//	super::init( cacheSize );
//
//	unassignedCams.reserve( m_BAfunc.getNumCameras() );
//	unassignedPts.reserve( m_BAfunc.getNumPoints() );
//
//	connectedPoints.resize( m_BAfunc.getNumCameras() );
//	connectedCameras.resize( m_BAfunc.getNumPoints() );
//
//	numNeighbourPoints.resize( m_BAfunc.getNumCameras(), 0 );
//	numNeighbourCameras.resize( m_BAfunc.getNumPoints(), 0 );
//
//	CameraID cid( -1 );
//	PointID pid( -1 );
//	for ( FactorID fid = 0; fid < factors.size(); ++fid ) {
//		BundleAdjustmentFactor & baf( *( (BundleAdjustmentFactor *) factors[fid] ) );
//		cid = baf.getCameraID();
//		pid = baf.getPointID();
//
//		connectedPoints[cid].push_back( pid );
//		connectedCameras[pid].push_back( cid );
//
//		++numNeighbourPoints[cid];
//		++numNeighbourCameras[pid];
//	}
//}
//
//
//void BundleAdjustmentOptimizer::prepareToOptimize( bool initVDNbrs ) {
//	super::prepareToOptimize( initVDNbrs );
//
//	unassignedCams.clear();
//	unassignedPts.clear();
//
//	for ( CameraID cid = 0; cid < m_BAfunc.getNumCameras(); ++cid ) {
//		unassignedCams.push_back( cid );
//	}
//	for ( PointID pid = 0; pid < m_BAfunc.getNumPoints(); ++pid ) {
//		unassignedPts.push_back( pid );
//	}
//}


} // namespace rdis
