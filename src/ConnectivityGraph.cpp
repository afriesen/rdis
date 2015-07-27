/*
 * ConnectivityGraph.cpp
 *
 *  Created on: Jun 14, 2013
 *      Author: afriesen
 */

#include "ConnectivityGraph.h"

#include BOOSTPATH/bind.hpp>
#include BOOSTPATH/filesystem/fstream.hpp>

namespace rdis {


void ConnectivityGraph::Edge::set( VertexID vxid1_, VertexID vxid2_,
								   ETForestLevel l_ ) {
	vxid1 = vxid1_;
	vxid2 = vxid2_;
	level = l_;
	isTree = false;

	arcs.clear();
	arcs.reserve( 1 );
}

void ConnectivityGraph::Edge::clear() {
	level = -1;

	// note that we're assuming the tree gets destroyed elsewhere
	arcs.clear();
}


ConnectivityGraph::Vertex::Vertex( Variable * vp_, VariableCount idOffset_,
		VariableCount nlevels )
	: vp( vp_ )
	, fp( NULL )
	, idOffset( idOffset_ )
{
	initialize( nlevels );
}

ConnectivityGraph::Vertex::Vertex( Factor * fp_, VariableCount idOffset_,
		VariableCount nlevels )
	: vp( NULL )
	, fp( fp_ )
	, idOffset( idOffset_ )
{
	initialize( nlevels );
}

void ConnectivityGraph::Vertex::initialize( VariableCount nlevels ) {
	ETvertices.reserve( nlevels );
	ETVertexID vxid = getETVertexID();

	// create an Euler tour vertex for each level and put each of these in their
	// own tour
	for ( VariableCount i = 0; i < nlevels; ++i ) {
		ETvertices.emplace_back( vxid, i );
		EulerTourTree::createTour( &ETvertices[i] );
	}

	adjTreeEdges.resize( nlevels );
	adjNontreeEdges.resize( nlevels );
}



ConnectivityGraph::ConnectivityGraph( VariableCount numvars )
	: numvariables( numvars )
	, numfactors( -1 )
	, numvertices( -1 )
	, numlevels( -1 )
{}

ConnectivityGraph::~ConnectivityGraph() {
	clear();
}


void ConnectivityGraph::init( const OptimizableFunction & func,
		const VariablePtrVec & vars, const FactorPtrVec & factors ) {

	assert( (VariableCount) vars.size() == numvariables );

	numfactors = factors.size();
	numvertices = numvariables + numfactors;

	numlevels = (VariableCount) ( std::log( (double) numvertices ) /
				std::log( 2.0 ) );

	const Clock::time_point starttime = Clock::now();

	std::cout << "connectivity graph init: " << numvariables << " vars, " <<
			numfactors << " factors, " << numvertices << " vertices, " <<
			numlevels << " levels" << std::endl;

	vertices.reserve( numvertices );

	// create all the forests (empty)
	forests.resize( numlevels );

//	edgepool.setMaxSize( numfactors * numvertices * numvertices );

#ifdef DEBUG
	assignedVars.clear();
	assignedFactors.clear();
	assignedVars.resize( numvertices );
	assignedFactors.resize( numfactors );

//	existingEdges.clear();
//	existingEdges.resize( numvertices,
//			std::vector< BOOSTNS::dynamic_bitset<> >( numvertices,
//					BOOSTNS::dynamic_bitset<>( numfactors )  ) );
#endif // DEBUG

	// create a vertex for each variable and put all the loop nodes in forests
	for ( Variable * vp : vars ) {
		vertices.emplace_back( vp, 0, numlevels );
		Vertex & vx = vertices.back();

		for ( ETForestLevel l = 0; l < numlevels; ++l ) {
			forests[l].emplace( vx.getETVertexID(), vx.ETvertices[l].loopNode );
		}
	}

	// create a vertex for each factor and put all the loop nodes in forests
	for ( Factor * fp : factors ) {
		vertices.emplace_back( fp, numvariables, numlevels );
		Vertex & vx = vertices.back();

		for ( ETForestLevel l = 0; l < numlevels; ++l ) {
			forests[l].emplace( vx.getETVertexID(), vx.ETvertices[l].loopNode );
		}
	}

	size_t nedges = 0;

	// add all of the edges -- this is a factor graph, so connect all factor
	// vertices to all variables contained in that factor
	for ( const Factor * f : factors ) {
		if ( f->isAssigned() ) continue;
		const VariablePtrVec & fvars( f->getVariables() );
		for ( const Variable * v : fvars ) {
			if ( v->isAssigned() ) continue;
			if ( f->containsVar( v->getID() ) ) {
//				std::cout << "inserting edge from vid: " << v->getID() <<
//						" to fid " << f->getID() << std::endl;
				insertEdge( vidToVertexID( v->getID() ),
						fidToVertexID( f->getID() ) );
				++nedges;
			}
		}
	}

	BOOSTNS::dynamic_bitset<> addedBlockEdges( func.getNumBlocks(), false );
	VariableID vidl, vidu;
	VariableCount bid;

	// add edges between all the variables in each block (so they never get
	// separated)
	for ( const Variable * vp : vars ) {
		if ( vp->isAssigned() ) continue;
		bid = func.getBlockID( vp->getID() );
		if ( addedBlockEdges[bid] ) continue;

		func.getBlockRangeByBlkId( bid, vidl, vidu );

		for ( VariableID vid1 = vidl; vid1 < vidu; ++vid1 ) {
			for ( VariableID vid2 = vid1+1; vid2 <= vidu; ++vid2 ) {
//				std::cout << "inserting edge from vid: " << vid1 <<
//						" to vid " << vid2 << " (block " << bid << ")"  << std::endl;
				insertEdge( vidToVertexID( vid1 ), vidToVertexID( vid2 ) );
				++nedges;
			}
		}

		addedBlockEdges[bid] = true;
	}

//	drawgraph();

	std::cout << "connectivity graph initialization completed in " <<
			( (Duration) ( Clock::now() - starttime ) ).count() <<
			" seconds (" << nedges << " edges)" << std::endl;
}


void ConnectivityGraph::clear() {
	for ( ETForestLevel i = 0; i < numlevels; ++i ) {
		for ( ETForest & f : forests ) {
			for ( ETForest::value_type & root : f ) {
				EulerTourTree::destroyTour( root.second->findRoot() );
			}
			f.clear();
		}
	}

	forests.clear();

	for ( Vertex & v : vertices ) {
		v.adjNontreeEdges.clear();
		v.adjTreeEdges.clear();
	}

	vertices.clear();
}


void ConnectivityGraph::onVariableAssigned( const Variable * v ) {
	// for each factor f that v is in and is not assigned
	// 		deleteEdge( v, f )

//	std::cout << "var " << v->getID() << " assigned (CG)" << std::endl;

#ifdef DEBUG
	assert( !assignedVars[v->getID()] );
	assignedVars.set( v->getID() );
#endif // DEBUG

	// note: we do assign/unassign in reverse order of each other, so we don't
	// get an inefficient chaining of edge replacements when all we really want
	// is a disconnect

	for ( auto fit = v->getFactors().begin(), fend = v->getFactors().end();
			fit != fend; ++fit ) {
		const Factor * f = *fit;
		if ( f->isAssigned() ) continue;
		if ( !f->containsVar( v->getID() ) ) continue;
//		std::cout << "delete edge v " << v->getID() << " to f " << f->getID() << std::endl;
		deleteEdge( vidToVertexID( v->getID() ), fidToVertexID( f->getID() ) );
	}
}


void ConnectivityGraph::onVariableUnassigned( const Variable * v ) {
	// for each factor f that v is in and is not assigned
	// 		addEdge( v, f )

//	std::cout << "var " << v->getID() << " unassigned (CG)" << std::endl;

#ifdef DEBUG
	assert( assignedVars[v->getID()] );
	assignedVars.reset( v->getID() );
#endif // DEBUG

	for ( auto fit = v->getFactors().rbegin(), fend = v->getFactors().rend();
			fit != fend; ++fit ) {
		const Factor * f = *fit;
		if ( f->isAssigned() ) continue;
		if ( !f->containsVar( v->getID() ) ) continue;
		insertEdge( vidToVertexID( v->getID() ), fidToVertexID( f->getID() ) );
	}
}


void ConnectivityGraph::onFactorAssigned( const Factor * f ) {
	// for each v in this factor that is not assigned
	//		deleteEdge( v, fid );

//	std::cout << "factor " << f->getID() << " assigned (CG) (" <<
// 		f->getVidAssigned() << ")" << std::endl;

#ifdef DEBUG
	assert( !assignedFactors[f->getID()] );
	assignedFactors.set( f->getID() );
#endif // DEBUG

	disconnectVars( f->getVariables(), f );
}

void ConnectivityGraph::onFactorUnassigned( const Factor * f ) {
	// for each v in this factor that is not assigned
	// 		addEdge( v, fid );

//	std::cout << "factor " << f->getID() << " unassigned (CG)" << std::endl;

#ifdef DEBUG
	assert( assignedFactors[f->getID()] );
	assignedFactors.reset( f->getID() );
#endif // DEBUG

	connectVars( f->getVariables(), f );
}


void ConnectivityGraph::connectVars( const VariablePtrVec & vars,
		const Factor * f ) {

	for ( auto vit = vars.begin(), vend = vars.end(); vit != vend; ++vit ) {
		const Variable & v = **vit;
		if ( !v.isAssigned() ) {
			insertEdge( vidToVertexID( v.getID() ),
					fidToVertexID( f->getID() ) );
		}
	}

}


void ConnectivityGraph::disconnectVars( const VariablePtrVec & vars,
		const Factor * f ) {

	for ( auto vit = vars.rbegin(), vend = vars.rend(); vit != vend; ++vit ) {
		const Variable & v = **vit;
		if ( !v.isAssigned() ) {
			deleteEdge( vidToVertexID( v.getID() ),
					fidToVertexID( f->getID() ) );
		}
	}
}


//void ConnectivityGraph::insertEdge( VariableID vid, FactorID fid ) {
void ConnectivityGraph::insertEdge( VertexID vxid1, VertexID vxid2 ) {

	assert( vxid1 != vxid2 );

	Vertex & u( vertices[vxid1] );
	Vertex & v( vertices[vxid2] );

	ETNode * ru = u.ETvertices[0].loopNode->findRoot();
	ETNode * rv = v.ETvertices[0].loopNode->findRoot();

	EdgeIP e = edgepool.get();
	e->set( /*vid, fid,*/ vxid1, vxid2, /*0,*/ 0 );

//	std::cout << "inserting/incrementing edge: " << *e << std::endl;


//#ifdef DEBUG
//	assert( !existingEdges[vxid][fxid][fid] );
//	assert( !existingEdges[fxid][vxid][fid] );
//	existingEdges[vxid][fxid].set( fid );
//	existingEdges[fxid][vxid].set( fid );
//#endif // DEBUG

//	if ( e->count == 0 ) {
//		e->set( vid, fid, vxid, fxid, 1, 0 ); // ensure level is 0

//		std::cout << "adding edge to graph: " << *e << std::endl;

		// check if they're already connected
		if ( ru == rv ) {
			// add edge to non-tree adjacency lists and update weights, which
			// correspond to edge counts
			addNontreeEdge( u, v, e, 0 );

		} else {
			e->isTree = true;

			// find and remove the two previous roots from the level 0 forest
			removeRoot( ru, 0 );
			removeRoot( rv, 0 );

			ETNode * newroot = addTreeEdge( u, v, e, 0 );

			// put the new root into the level 0 forest -- note that we use the
			// root's vertex ID as the index
			insertRoot( newroot, 0 );
		}
//	} else {
//		++e->count;
//	}
}


//void ConnectivityGraph::deleteEdge( VariableID vid, FactorID fid ) {
void ConnectivityGraph::deleteEdge( VertexID vxid1, VertexID vxid2 ) {

	assert( vxid1 != vxid2 );

//	std::cout << "deleting / decrementing edge (" << vid << ", " << fid << " / " << fxid << ")" << std::endl;

////#ifdef DEBUG
////	assert( existingEdges[e->vxid][e->fxid][fid] );
////	assert( existingEdges[e->fxid][e->vxid][fid] );
////	existingEdges[e->vxid][e->fxid].reset( fid );
////	existingEdges[e->fxid][e->vxid].reset( fid );
////#endif // DEBUG
//
//
//	if ( --e->count > 0 ) return;

	Vertex & u( vertices[vxid1] );
	Vertex & v( vertices[vxid2] );

	// find (and remove) the edge from the adjacency lists
//	EdgeIP ee = getFromAdjLists( u, v, e->level, e->isTree );
	EdgeIP ee = getFromAdjLists( u, v );

	if ( ee == NULL ) {
		std::cout << "Can't delete edge " <<
				Edge( /*vid, fid,*/ vxid1, vxid2, -1 ) << " (not found) -- " <<
				vxid1 << " / " << vxid2 << std::endl;
		assert( ee != NULL );
		return;
	}

	assert( ee->isequal( vxid1, vxid2 ) /*&& *e == *ee*/ );

	deleteEdge( ee );
}


void ConnectivityGraph::deleteEdge( EdgeIP e ) {
//	std::cout << "deleting edge " << *e <<
//			( e->isTree ? " (tree)" : " (nontree)" ) << std::endl;

//	assert( e->count == 0 );

	Vertex & u( vertices[e->vxid1] );
	Vertex & v( vertices[e->vxid2] );
	assert( EulerTourTree::connected( &u.ETvertices[0], &v.ETvertices[0] ) );

	if ( e->isTree ) {
		// cut all forests F_i, for 0 <= i <= level(e)
		for ( ETForestLevel i = e->level; i >= 0; --i ) {
			ETNode * oldroot = e->arcs[i].first->findRoot();
			assert( oldroot == e->arcs[i].second->findRoot() );

			// remove the old root from the forest
			removeRoot( oldroot, i );

			ETNode * r1 = NULL, * r2 = NULL;
			EulerTourTree::cut( e->arcs[i].first, e->arcs[i].second, &r1, &r2 );

#ifdef DEBUG
			ETNode * ru = u.ETvertices[i].loopNode->findRoot();
			ETNode * rv = v.ETvertices[i].loopNode->findRoot();
			assert( ru != rv );
			assert( ru == r1 || ru == r2 );
			assert( rv == r1 || rv == r2 );
#endif

			// put the new roots in the forest
			insertRoot( r1, i );
			insertRoot( r2, i );
		}

		// attempt to find a replacement edge
		for ( ETForestLevel i = e->level; i >= 0; --i ) {
			if ( replaceEdge( *e, i ) ) break;
		}
	}

//	delete e;
}


bool ConnectivityGraph::replaceEdge( const Edge & e, ETForestLevel lvl ) {

	Vertex * u = &vertices[e.vxid1];
	Vertex * v = &vertices[e.vxid2];

	ETNode * rTu = u->ETvertices[lvl].loopNode->findRoot();
	ETNode * rTv = v->ETvertices[lvl].loopNode->findRoot();

	assert( u != v );
	assert( rTu != rTv );

	assert( forests[lvl].find( rTu->vx1->vxid ) != forests[lvl].end() );
	assert( forests[lvl].find( rTv->vx1->vxid ) != forests[lvl].end() );

	// Tu = the smallest of the two subtrees resulting from the deletion of e
	if ( rTv->stsize < rTu->stsize ) {
		std::swap( u, v );
		std::swap( rTu, rTv );
	}

	// ------ outline of algorithm -------
	// raise all edges in T (i.e., tree edges) at level i to i+1
	// --> move them up a level and add them all to F_{i+1} (but don't change F_i)
	// for all level i non-tree edges incident to T_i (skip subtrees w/o level i edges)
	//		if edge f doesn't connect Tv and Tw
	//			increase its level to i+1
	//		else
	//			insert f as replacement edge for e (for all F_i, i <= l(e))
	//			return true
	//
	// return false

	// walk over all level i edges incident to T -- tree and non-tree

	ETNode * n = rTu;
	assert( n != NULL );

	ETSplayTree::SplayNodeStack stack;
	stack.reserve( 5 );

	// pre-order visit all loop nodes with tree edges in the current ET tree and
	// raise these to the next level
	while ( n != NULL ) {
		if ( n->isLoop() && n->stweight[ETWeight::Tree] > 0 ) {
			raiseTreeEdges( n, lvl );
		}

		n = (ETNode *) n->nextPreOrder( stack, true, ETWeight::Tree );
	}

	assert( stack.empty() );
	stack.clear();
	n = rTu;

	// check all non-tree edges on this level for a replacement edge
	while ( n != NULL ) {
		if ( n->isLoop() && n->stweight[ETWeight::Nontree] > 0 ) {
			if ( tryToReplace( *u, *v, lvl, rTu, rTv, n ) ) return true;
		}

		n = (ETNode *) n->nextPreOrder( stack, true, ETWeight::Nontree );
	}

//	std::cout << "unable to replace " << e << " at level " << lvl << ", giving: \n\t";
//	EulerTourTree::print( rTu ); std::cout << "\n\t";
//	EulerTourTree::print( rTv ); std::cout << std::endl;

	return false;
}


void ConnectivityGraph::addNontreeEdge( Vertex & u, Vertex & v, EdgeIP e,
		ETForestLevel l ) {

	assert( l < numlevels );
	assert( !e->isTree );

	// put edge into the adjacency lists (for quick iteration)
	u.adjNontreeEdges[l].push_back( e );
	v.adjNontreeEdges[l].push_back( e );

	// update the weight of the relevant node (weight == number of children)
	u.ETvertices[l].loopNode->addWeight( 1, ETWeight::Nontree );
	v.ETvertices[l].loopNode->addWeight( 1, ETWeight::Nontree );

//	std::cout << "added non-tree edge " << *e << " to level "<< l << std::endl;
}


ETNode * ConnectivityGraph::addTreeEdge( Vertex & u, Vertex & v, EdgeIP e,
		ETForestLevel l ) {
	assert( l < numlevels );
	assert( e->level >= l );
	assert( e->isTree );

	if ( e->level == l ) {
		// put edge into the adjacency lists (for quick iteration)
		u.adjTreeEdges[l].push_back( e );
		v.adjTreeEdges[l].push_back( e );

		u.ETvertices[l].loopNode->addWeight( 1, ETWeight::Tree );
		v.ETvertices[l].loopNode->addWeight( 1, ETWeight::Tree );
	}

	// connect the two vertices in the level 0 forest (F_0)
	ETNode * neuv = NULL, * nevu = NULL;
	ETNode * newroot = EulerTourTree::link( &u.ETvertices[l],
			&v.ETvertices[l], &neuv, &nevu );

//	std::cout << "new ET tree after adding edge " << *e << ": ";
//	EulerTourTree::print( newroot ); std::cout << std::endl;

	// store the pointers in the edge
	if ( l >= (ETForestLevel) e->arcs.size() ) {
		e->arcs.resize( l+1 );
	}
	e->arcs[l].first = neuv;
	e->arcs[l].second = nevu;

	return newroot;
}


void ConnectivityGraph::raiseTreeEdges( ETNode * n, ETForestLevel l ) {
	assert( n->isLoop() );
	Vertex & u( vertices[n->vx1->vxid] );

	auto & list = u.adjTreeEdges[l];

	while ( !list.empty() ) {
		// get and remove the last edge (efficient in a vector)
		EdgeIP e = list.back();

		assert( e->isTree );
		assert( e->level == l );

//		removeTreeEdge( e.get(), l );

		list.pop_back();

		Vertex & v( e->vxid1 == n->vx1->vxid ? vertices[e->vxid2] :
				vertices[e->vxid1] );

		assert( EulerTourTree::connected( &u.ETvertices[l], &v.ETvertices[l] ) );

		// remove e from v's adjacency lists
		EdgeIP ev = getFromList( e->vxid1, e->vxid2, v.adjTreeEdges[l], l );
		assert( ev == e );

		u.ETvertices[l].loopNode->addWeight( -1, ETWeight::Tree );
		v.ETvertices[l].loopNode->addWeight( -1, ETWeight::Tree );

		// get the previous roots on level l+1
		ETNode * ru = u.ETvertices[l+1].loopNode->findRoot();
		ETNode * rv = v.ETvertices[l+1].loopNode->findRoot();

		assert( ru != rv );

		// remove the previous roots
		removeRoot( ru, l+1 );
		removeRoot( rv, l+1 );

		// move the edge up a level, add it to this new level (and the forest)
		e->level = l + 1;
		ETNode * newroot = addTreeEdge( u, v, e, l+1 );

		// remember the new root
		insertRoot( newroot, l+1 );
	}
}


void ConnectivityGraph::raiseNontreeEdge( Edge * e, ETForestLevel l ) {
	Vertex & u( vertices[e->vxid1] );
	Vertex & v( vertices[e->vxid2] );

	assert( e->level == l );
	assert( !e->isTree );
	assert( EulerTourTree::connected( &u.ETvertices[l], &v.ETvertices[l] ) );

	removeNontreeEdge( e, l );
	e->level = l+1;
	addNontreeEdge( u, v, e, l+1 );
}


void ConnectivityGraph::removeNontreeEdge( Edge * e, ETForestLevel l ) {
	Vertex & u( vertices[e->vxid1] );
	Vertex & v( vertices[e->vxid2] );

	// remove e from adjacency lists
	EdgeIP eu = getFromList( e->vxid1, e->vxid2, u.adjNontreeEdges[l], l );
	EdgeIP ev = getFromList( e->vxid1, e->vxid2, v.adjNontreeEdges[l], l );
	assert( eu != NULL && eu == ev );

	u.ETvertices[l].loopNode->addWeight( -1, ETWeight::Nontree );
	v.ETvertices[l].loopNode->addWeight( -1, ETWeight::Nontree );
}


bool ConnectivityGraph::tryToReplace( const Vertex & u, const Vertex & v,
		ETForestLevel lvl, ETNode * rTu, ETNode * rTv, ETNode * n ) {

	assert( n->isLoop() );
	assert( n->vx1->level == lvl );
	assert( rTu->stsize <= rTv->stsize );

	Vertex & x( vertices[n->vx1->vxid] );

	assert( EulerTourTree::connected( &x.ETvertices[lvl], &u.ETvertices[lvl] ) );

	// iterate over all non-tree edges ( f = (x,y) ) incident to n (which is a
	// node in u's Euler tour tree)

	for ( auto reit = x.adjNontreeEdges[lvl].rbegin(),
			rend = x.adjNontreeEdges[lvl].rend(); reit != rend; ) {

		EdgeIP f = *reit;

		assert( f->level == lvl );
		assert( (ETForestLevel) f->arcs.size() <= lvl ||
				f->arcs[lvl].first == NULL );
		assert( f->vxid1 == x.getETVertexID() || f->vxid2 == x.getETVertexID() );

		VertexID yvid = ( f->vxid1 == x.getETVertexID() ? f->vxid2 : f->vxid1 );
		assert( yvid != n->vx1->vxid );
		ETVertex & yET( vertices[yvid].ETvertices[lvl] );

		assert( EulerTourTree::connected( &v.ETvertices[lvl], &yET ) ==
				( rTv == yET.loopNode->findRoot() ) );

		// must do this before removing this edge from the list we're
		// iterating over
		++reit;

		// check if we found a replacement or not (i.e., this edge is (x,y) and
		// y is connected to v, thus connecting u to v through x)
		if ( rTv == yET.loopNode->findRoot() ) {
			assert( !EulerTourTree::connected( &u.ETvertices[0],
					&v.ETvertices[0] ) );

			// found a replacement, so insert it into all lower trees
			setAsReplacement( f, lvl );

			assert( EulerTourTree::connected( &u.ETvertices[0],
					&v.ETvertices[0] ) );
			return true;

		} else {
			// didn't find a replacement, so increase the level of this edge
			raiseNontreeEdge( f.get(), lvl );
		}
	}

	return false;
}


void ConnectivityGraph::setAsReplacement( EdgeIP e, ETForestLevel l ) {
	assert( e->level == l );
	assert( !e->isTree );

	Vertex & u( vertices[e->vxid1] );
	Vertex & v( vertices[e->vxid2] );

	// remove from non-tree list and update weights
	removeNontreeEdge( e.get(), l );

	e->isTree = true;

//	std::cout << "adding edge " << *e << " as replacement, giving " << std::endl;

	for ( int i = 0; i <= l; ++i ) {
	 	// get the previous roots on level l
		ETNode * ru = u.ETvertices[i].loopNode->findRoot();
		ETNode * rv = v.ETvertices[i].loopNode->findRoot();

		assert( ru != rv );

		// remove the previous roots
		removeRoot( ru, i );
		removeRoot( rv, i );

		// add the edge to this level (and the forest)
		ETNode * newroot = addTreeEdge( u, v, e, i );

		// remember the new root
		insertRoot( newroot, i );
	}

//	EulerTourTree::print( u.ETvertices[0].loopNode->findRoot() );
//	std::cout << std::endl;
}


void ConnectivityGraph::removeRoot( ETNode * r, ETForestLevel l ) {
	assert( l < (ETForestLevel) forests.size() );
	ETForest::iterator itr = forests[l].find( r->vx1->vxid );
	assert( itr != forests[l].end() );
	forests[l].erase( itr );
}


void ConnectivityGraph::insertRoot( ETNode * r, ETForestLevel l ) {
	assert( forests[l].find( r->vx1->vxid ) == forests[l].end() );
	forests[l][r->vx1->vxid] = r;
}


ConnectivityGraph::EdgeIP ConnectivityGraph::getFromAdjLists( Vertex & u,
		Vertex & v, ETForestLevel l, bool isTree ) {

	VertexID uvid( u.getETVertexID() );
	VertexID vvid( v.getETVertexID() );

	// look for the specified edge in the set of tree edges, organized by level
	// (note that this erases the edge from the adjacency lists)

	EdgeIP e = NULL;
	EdgeIP eu = NULL, ev = NULL;

	AdjacencyList * l1 = ( isTree ? &u.adjTreeEdges[l] : &u.adjNontreeEdges[l] );
	AdjacencyList * l2 = ( isTree ? &v.adjTreeEdges[l] : &v.adjNontreeEdges[l] );

	if ( l2->size() < l1->size() ) std::swap( l1, l2 );

	eu = getFromList( uvid, vvid, *l1, l );
	ev = getFromList( uvid, vvid, *l2, l );

//	if ( eu == NULL || ev == NULL ) e = NULL;

	if ( eu == NULL && ev == NULL ) {
		std::cout << "getFromAdjList: " << uvid << ", " << vvid <<
				" for vertices "<< std::endl;
		std::cout << "list l1: "; printContainer( *l1 ); std::cout << "\n";
		std::cout << "list l2: "; printContainer( *l2 ); std::cout << "\n";
	}
	assert( eu != NULL && ev != NULL );

	if ( eu != ev ) {
		std::cout << "nt w: " << v.ETvertices[l].loopNode->weight[ETWeight::Nontree] <<
				", t w: " << v.ETvertices[l].loopNode->weight[ETWeight::Nontree] << std::endl;
		std::cout << *ev << ", l = " << l << std::endl;
	}

	assert( eu == ev );
	assert( eu->level == l );
//	assert( eu->count == 0 );

	if ( eu != NULL ) {
		e = eu;

		// decrease the appropriate weight now that we've removed it
		ETWeight::Type t = ( e->isTree ? ETWeight::Tree : ETWeight::Nontree );
		u.ETvertices[e->level].loopNode->addWeight( -1, t );
		v.ETvertices[e->level].loopNode->addWeight( -1, t );
	}

	return e;
}


ConnectivityGraph::EdgeIP ConnectivityGraph::getFromAdjLists(
		Vertex & u, Vertex & v ) {

	VertexID uvid( u.getETVertexID() );
	VertexID vvid( v.getETVertexID() );

	EdgeIP e, eu, ev;

	bool istree = false;
	ETForestLevel l = 0;

//	std::cout << "trying to getFromAdjLists edge " << uvid << " -- " << vvid << std::endl;

	AdjacencyList * l1 = NULL, * l2 = NULL;

	// look in all of the lists for this edge
	for ( l = numlevels-1; l >= 0 && !eu && !ev; --l ) {
		for ( int isT = 0; isT <= 1; ++isT ) {
			istree = (bool) isT;
			l1 = ( istree ? &u.adjTreeEdges[l] : &u.adjNontreeEdges[l] );
			l2 = ( istree ? &v.adjTreeEdges[l] : &v.adjNontreeEdges[l] );

			// find (and remove) the edge from the smaller list
			if ( !l1->empty() && !l2->empty() ) {
				if ( l1->size() <= l2->size() ) {
					eu = getFromList( uvid, vvid, *l1, l );
					if ( eu != NULL ) {
						ev = getFromList( uvid, vvid, *l2, l );
						break;
					}
				} else {
					ev = getFromList( uvid, vvid, *l2, l );
					if ( ev != NULL ) {
						eu = getFromList( uvid, vvid, *l1, l );
						break;
					}
				}
			}
		}
	}

	if ( eu == NULL || ev == NULL ) {
		std::cout << "eu NULL? " << (eu == NULL) << ", ev NULL? " << (ev==NULL)
				<< std::endl;
		if ( eu != NULL ) std::cout << "eu " << *eu << std::endl;
		if ( ev != NULL ) std::cout << "ev " << *ev << std::endl;
		std::cout << "istree " << istree << ", lvl " << l << std::endl;
		std::cout << "l1: "; printContainer( *l1 ); std::cout << std::endl;
		std::cout << "l2: "; printContainer( *l2 ); std::cout << std::endl;
	}

	// make sure we found it, otherwise something went wrong...
	assert( eu != NULL && ev != NULL );
	assert( eu == ev );

	if ( eu != NULL ) {
		e = eu;

		// decrease the appropriate weight now that we've removed it
		ETWeight::Type t = ( e->isTree ? ETWeight::Tree : ETWeight::Nontree );
		u.ETvertices[e->level].loopNode->addWeight( -1, t );
		v.ETvertices[e->level].loopNode->addWeight( -1, t );
	}

	return e;
}


ConnectivityGraph::EdgeIP ConnectivityGraph::getFromList( VertexID uvxid,
		VertexID vvxid, AdjacencyList & edges, ETForestLevel lvl ) const {

	EdgeIP e = NULL;

	for ( auto rit = edges.end(), rend = edges.begin(); rit != rend; ) {
		--rit;
		if ( (*rit)->level != lvl ) {
			std::cout << "CG ERR: getFromList " << uvxid << "," << vvxid << ": "
					<< **rit << ", lvl " << lvl << std::endl;
		}
		assert( (*rit)->level == lvl );

		if ( (*rit)->isequal( uvxid, vvxid ) ) {
			e = *rit;
			edges.erase( rit );
			break;
		}
	}

	return e;
}


void ConnectivityGraph::drawgraph( ETForestLevel level ) const {
	BOOSTNS::filesystem::fstream ofs( "splaytree.gv", std::ios_base::out );
	ofs << "graph G {" << std::endl;
	size_t curid = 0;
	for ( auto it = forests[level].begin(); it != forests[level].end(); ++it ) {
		it->second->findRoot()->graphsubtree( ofs, curid );
	}
	ofs << "}" << std::endl;
	ofs.close();
}


VariableCount ConnectivityGraph::getNumEdges( VertexID vxid ) {
	return vertices[vxid].adjNontreeEdges[0].size() +
			vertices[vxid].adjTreeEdges[0].size();
}


ConnectivityGraph::IntraCCompIterator::IntraCCompIterator( ETNode * n )
		: cur( (ETNode *) n->findRoot()->subtree_minimum() ) {
	if ( !cur->isLoop() ) step();
	assert( cur->isLoop() );
}


void ConnectivityGraph::IntraCCompIterator::step() {
	do {
		cur = (ETNode *) cur->next();
	} while ( cur != NULL && !cur->isLoop() );
}


std::ostream & operator<<( std::ostream & os, const ConnectivityGraph::Edge & e ) {
	std::cout << "{(" << /*e.vid << "," << e.fid << "--" <<*/ e.vxid1 << "," <<
			e.vxid2 << ")," /*<< e.count << "|"*/ << e.level <<
			( e.isTree ? ",t" : "" ) << "}";
	return os;
}

} // namespace rdis
