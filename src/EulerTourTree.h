/*
 * EulerTourTree.h
 *
 *  Created on: Jun 14, 2013
 *      Author: afriesen
 */

#ifndef RDIS_EULERTOURTREE_H_
#define RDIS_EULERTOURTREE_H_

#include "common.h"
#include "SplayTree.h"
#include "IntrusivePtrPool.h"

#include BOOSTPATH/container/vector.hpp>


namespace rdis {


// NOTE: We are using Tarjan's specification of Euler tours from "Dynamic trees
// as search trees via Euler tours, applied to the network simplex algorithm.",
// Tarjan, 1996 -- here, the Euler tour is explicitly represented as a sequence
// of arcs (directed edges) over a tree (a spanning tree in this case) with one
// "loop" arc per vertex visited by the tour which contains subtree size and
// weight information.


typedef VariableCount ETForestLevel;

struct ETComp {
	bool operator()() const {
		assert( false );
		return false;
	}
};


// typedefs and forward declarations
typedef long long int ETSplayKey;
typedef long long int ETSplayWeight;
typedef SplayTree< ETSplayKey, ETSplayWeight, ETComp, 2 > ETSplayTree;

typedef long long int ETVertexID;

namespace ETWeight {
enum Type {
	Tree = 0,
	Nontree = 1
};
}

struct ETVertex;
struct ETNode;

//typedef BOOSTNS::intrusive_ptr< ETNode > ETNodeIP;
//typedef BOOSTNS::intrusive_ptr< const ETNode > ETNodeCIP;

typedef IntrusivePtrPool< ETNode > ETNodePool;


// ET nodes represent the arcs (directed edges) and loops (the visit of the
// corresponding vertex in the tour) in the Euler tour of the spanning tree
struct ETNode
		: public ETSplayTree::SplayNode
		, public IntrusivePtrPoolObj< ETNode, ETNodePool > {

	typedef ETSplayTree::SplayNode super;

	ETNode( ETVertex * v1, ETVertex * v2 )
			: super()
			, vx1( v1 )
			, vx2( v2 )
	{
		assert( vx1 != NULL && vx2 != NULL );
//		std::cout << "created ETNode: "; this->print(); std::cout << std::endl;
	}

	virtual ~ETNode() {
//		std::cout << "deleting ETNode: "; this->print(); std::cout << std::endl;
		assert( isDisconnected() );
	}

	void set( ETVertex * v1, ETVertex * v2 ) {
		vx1 = v1;
		vx2 = v2;
	}

public:

	// return true if this is a loop node
	bool isLoop() const { return ( vx1 == vx2 ); }

	bool operator==( const ETNode & rhs ) const {
		return ( vx1 == rhs.vx1 && vx2 == rhs.vx2 );
	}

	// return the root of this node's tree
	virtual const ETNode * findRoot() const {
		return (ETNode *) super::findRoot();
	}

	virtual ETNode * findRoot() {
		return (ETNode *) super::findRoot();
	}
	// TODO: get rid of these overrides

public:
	virtual void print( std::ostream & os = std::cout ) const;

public:
	virtual void clear();

public:
	// pointers to the two vertices that are at the ends of the arc (directed
	// edge) that this node represents in the Euler tour
	ETVertex * vx1;
	ETVertex * vx2;
};

std::ostream & operator<<( std::ostream & os, const ETNode & n );


// this struct represents the unique vertices in the graph that the Euler tour
// visits -- there will (in general) be multiple ETNodes pointing to each ETVertex
struct ETVertex {

	ETVertex( ETVertexID vxid_, ETForestLevel level_ )
		: vxid( vxid_ )
		, level( level_ )
		, loopNode( NULL )
	{}

	virtual ~ETVertex() {}


	// the ID of this vertex
	ETVertexID vxid;

	// the level of the forest in which this vertex exists
	ETForestLevel level;

	// the tree node representing the active occurrence of this vertex
	ETNode * loopNode;
};


typedef BOOSTNS::container::vector< ETVertex > ETVertexVec;


// extends a balanced n-ary search tree (in this case a binary splay tree) to
// provide O(log(n)) manipulation of a tree's Euler tour by putting the tour in
// a search tree (with somewhat extended functionality)
class EulerTourTree {
public:

	// called to create a tour containing only this vertex
	static void createTour( ETVertex * v );

	// delete the entire Euler tour rooted at this node
	static void destroyTour( ETNode * root );


	// remove the edge {u, v}, to which neuv and nevu are the corresponding arcs
	// in the Euler tour -- set r1 and r2 to point to the roots of the resulting
	// trees
	static void cut( ETNode * neuv, ETNode * nevu, ETNode ** r1, ETNode ** r2 );

	// create an edge connecting u and v, setting neuv and nevu to point to the
	// two arcs representing that edge -- return the root of the tree that is
	// the combination of the existing two trees
	static ETNode * link( ETVertex * u, ETVertex * v,
			ETNode ** neuv, ETNode ** nevu );

	// return true if the specified nodes are connected (i.e. in the same ET tree)
	static bool connected( const ETVertex * u, const ETVertex * v );
	static bool connected( const ETNode * u, const ETNode * v );

	// return the root of this node's Euler tour -- the first (or last) node in
	// the sequence representing the tour
	static ETNode * findTourRoot( const ETNode * n ) {
		ETNode * r( (ETNode *) n->findRoot()->subtree_minimum() );
		assert( r );
		return r;
	}

	static void print( const ETNode * n );


private:
	// join the two trees represented by nodes p and q (duplicates functionality
	// in splay tree, but is more convenient here) and return the root of the
	// resulting tree
	static ETNode * join( ETNode * p, ETNode * q );

};

} // namespace rdis

#endif // RDIS_EULERTOURTREE_H_
