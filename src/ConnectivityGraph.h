/*
 * ConnectivityGraph.h
 *
 *  Created on: Jun 14, 2013
 *      Author: afriesen
 */

#ifndef RDIS_CONNECTIVITYGRAPH_H_
#define RDIS_CONNECTIVITYGRAPH_H_

#include "common.h"
#include "IntrusivePtrPool.h"
#include "Variable.h"
#include "Factor.h"
#include "EulerTourTree.h"

#include "OptimizableFunction.h"

#include BOOSTPATH/multi_array.hpp>
#include BOOSTPATH/dynamic_bitset.hpp>
#include BOOSTPATH/container/list.hpp>
#include BOOSTPATH/container/vector.hpp>
#include BOOSTPATH/container/flat_map.hpp>

#include <list>

namespace rdis {

class ConnectivityGraph {
public:
	typedef long long int VertexID;

protected:
	typedef std::pair< ETNode *, ETNode * > ETArcPair;
	typedef BOOSTNS::container::vector< ETArcPair > ETArcPairs;

	typedef BOOSTNS::container::flat_map< ETVertexID, ETNode * > ETForest;
	typedef BOOSTNS::container::vector< ETForest > ETForestVec;


	struct Edge;
	typedef BOOSTNS::intrusive_ptr< Edge > EdgeIP;
	typedef BOOSTNS::intrusive_ptr< const Edge > EdgeCIP;
	typedef IntrusivePtrPool< Edge > EdgePool;

	typedef BOOSTNS::container::vector< EdgeIP > AdjacencyList;
//	typedef BOOSTNS::container::vector< EdgeIP > AdjacencyVec;
	typedef BOOSTNS::container::vector< AdjacencyList > AdjacencyListVec; // index by level
//	typedef BOOSTNS::multi_array< EdgeIP, 2 > AdjacencyMatrix;


	// edges between variables (or, equivalently, between their corresponding
	// vertices)
	struct Edge
			: public IntrusivePtrPoolObj< Edge, EdgePool > {

		Edge()
			: vxid1( -1 )
			, vxid2( -1 )
			, level( -1 )
			, isTree( false )
		{}

		Edge( VertexID vxid1_, VertexID vxid2_, ETForestLevel level_ )
				: vxid1( vxid1_ )
				, vxid2( vxid2_ )
				, level( level_ )
//				, count( 1 )
				, isTree( false )
		{
			arcs.reserve( 1 );
		}

		~Edge() {
			// note: arcs should already have been deleted when the tree was
			// cleaned up, so don't re-delete here (that memory may have already
			// been reallocated so can't do checks here either)
		}

	public:
		void set( VertexID vxid1_, VertexID vxid2_, ETForestLevel l_ );
		void clear();

	public:
		bool isequal( VertexID uu, VertexID vv ) const {
			return ( ( vxid1 == uu && vxid2 == vv ) ||
					( vxid1 == vv && vxid2 == uu ) );
		}

		bool operator==( const Edge & rhs ) const {
			assert( level == rhs.level ); // shouldn't compare across levels?
			bool res = isequal( rhs.vxid1, rhs.vxid2 );
//			assert( !res || ( rhs.vid == vid && rhs.fid == fid ) );
			return res;
		}

	public:
		// the corresponding vertex IDs for this edge
		VertexID vxid1;
		VertexID vxid2;

		// the level that this edge is at
		ETForestLevel level;

		// true if this is a tree edge, false otherwise
		bool isTree;

		// list of pairs of arcs that correspond to this edge (indexed by level,
		// for all levels <= this->level)
		ETArcPairs arcs;
	};


	// a vertex in the connectivity graph is a representation of a variable or
	// factor and contains a corresponding EulerTour vertex for each EulerTour
	// forest (and there is an ET forest for each level, of which there are
	// log(|V|) levels, where |V| is the number of vertices)
	struct Vertex {
		Vertex( Variable * vp_, VariableCount idOffset_, VariableCount nlevels );
		Vertex( Factor * fp_, VariableCount idOffset_, VariableCount nlevels );
		~Vertex() {}

	private:
		void initialize( VariableCount nlevels );

	public:
		ETVertexID getETVertexID() const {
			assert( vp != NULL || fp != NULL );
			return ( ( vp != NULL ? vp->getID() : fp->getID() ) + idOffset );
		}

	public:
		// pointer to the variable/factor that this vertex represents
		const Variable * vp;
		const Factor * fp;

		// offset for the ETvertex ID of this vertex
		VariableCount idOffset;

		// the copies of this vertex -- one for each level
		ETVertexVec ETvertices;

		// adjacency list of edges incident to this vertex that are in an ET
		// tree (one list per level)
		AdjacencyListVec adjTreeEdges;

		// adjacency list containing the incident edges that are not in an ET
		// tree (one list per level)
		AdjacencyListVec adjNontreeEdges;
	};


	typedef BOOSTNS::container::vector< Vertex > VertexVec;

	friend std::ostream & operator<<( std::ostream & os,
			const ConnectivityGraph::Edge & e );

public:
	class InterCCompIterator;
	class IntraCCompIterator;


public:
	ConnectivityGraph( VariableCount numvars );
	~ConnectivityGraph();

	void init( const OptimizableFunction & func,
			const VariablePtrVec & vars, const FactorPtrVec & factors );

	void clear();

public:
	// called when a variable is assigned/unassigned
	void onVariableAssigned( const Variable * v );
	void onVariableUnassigned( const Variable * v );

	// called when a factor is assigned/unassigned
	void onFactorAssigned( const Factor * f );
	void onFactorUnassigned( const Factor * f );

	// connect/disconnect the specified set of vars through this factor id
	void connectVars( const VariablePtrVec & vars, const Factor * f );
	void disconnectVars( const VariablePtrVec & vars, const Factor * f );

private:
	void insertEdge( VertexID vxid1, VertexID vxid2 );

	void deleteEdge( VertexID vxid1, VertexID vxid2 );
	void deleteEdge( EdgeIP e );

	// conversion helper functions
	VertexID vidToVertexID( VariableID vid ) const { return vid; }
	VertexID fidToVertexID( FactorID fid ) const {
		return ( numvariables + fid );
	}

public:
	bool isVariableID( VertexID vxid ) const { return ( vxid < numvariables ); }
	bool isFactorID( VertexID vxid ) const { return ( vxid >= numvariables ); }

	VariableID vxidToVID( VertexID vxid ) const { return vxid; }
	FactorID vxidToFID( VertexID vxid ) const { return ( vxid - numvariables ); }

private:
	bool replaceEdge( const Edge & e, ETForestLevel lvl );

private:

	void addNontreeEdge( Vertex & u, Vertex & v, EdgeIP e, ETForestLevel l );
	ETNode * addTreeEdge( Vertex & u, Vertex & v, EdgeIP e, ETForestLevel lvl );

	void raiseTreeEdges( ETNode * n, ETForestLevel l );
//	void removeTreeEdge( Edge * e, ETForestLevel l );

	void raiseNontreeEdge( Edge * e, ETForestLevel l );
	void removeNontreeEdge( Edge * e, ETForestLevel l );

	bool tryToReplace( const Vertex & u, const Vertex & v, ETForestLevel lvl,
			ETNode * rTu, ETNode * rTv, ETNode * n );

	void setAsReplacement( EdgeIP e, ETForestLevel l );

private:
	// find and remove the specified root from the level l forest
	void removeRoot( ETNode * r, ETForestLevel l );

	// insert the specified root into the level l forest -- note that we use the
	// root's vertex ID as the index
	void insertRoot( ETNode * r, ETForestLevel l );

	// find the edge between these vertices and remove it from their adjacency
	// lists, given that we know the level and whether or not it's a tree edge
	EdgeIP getFromAdjLists( Vertex & u, Vertex & v, ETForestLevel l,
			bool isTree );

	// find the edge between these vertices and remove it from their adjacency
	// lists, not knowing the level or whether it's in the tree
	EdgeIP getFromAdjLists( Vertex & u, Vertex & v );

	// find the corresponding edge in the edge list, remove it and return it
	EdgeIP getFromList( VertexID uvxid, VertexID vvxid,
			AdjacencyList & edges, ETForestLevel lvl ) const;

	void drawgraph( ETForestLevel level = 0 ) const;

public:
	VariableCount getNumComponents() const {
		return forests[0].size();
	}

	InterCCompIterator getConnCompIterator() {
		return InterCCompIterator( forests[0] );
	}

	IntraCCompIterator getIntraCCompIterator( VertexID vxid ) {
		return IntraCCompIterator( vertices[vxid].ETvertices[0].loopNode );
	}

	VertexID getComponentLabel( VertexID vxid ) {
		return vertices[vxid].ETvertices[0].loopNode->findRoot()->vx1->vxid;
	}

	// return the number of edges connected to the specified variable
	VariableCount getNumEdges( VertexID vxid );


private:
	// the number of variables and factors that in this graph (summing these
	// gives the number of vertices)
	const VariableCount numvariables;
	VariableCount numfactors;

	// the number of vertices in this graph
	VariableCount numvertices;

	// the number of levels used to maintain the connectivity graph
	VariableCount numlevels;

	// a pool of edges (so we're not constantly creating/deleting these)
	// -- KEEP AT TOP
	EdgePool edgepool;

	// the set of vertices that exist in the graph
	VertexVec vertices;

	// the set of forests (one for each level)
	ETForestVec forests;

#ifdef DEBUG
	BOOSTNS::dynamic_bitset<> assignedVars;
	BOOSTNS::dynamic_bitset<> assignedFactors;
//	std::vector< std::vector< BOOSTNS::dynamic_bitset<> > > existingEdges;
#endif // DEBUG


public:
	// iterator to iterate across connected components
	class InterCCompIterator {
	public:
		InterCCompIterator( const ETForest & forest )
			: fiter( forest.begin() )
			, fend( forest.end() ) {}
		~InterCCompIterator() {}

		IntraCCompIterator operator*() {
			return IntraCCompIterator( fiter->second );
		}

		InterCCompIterator & operator++() { // prefix
			++fiter;
			return *this;
		}

		InterCCompIterator operator++( int /*unused*/ ) { // postfix
			InterCCompIterator result = *this;
			++( *this );
			return result;
		}

		bool operator==( const InterCCompIterator & rhs ) {
			return ( fiter == rhs.fiter );
		}
		bool operator!=( const InterCCompIterator & rhs ) {
			return ( fiter != rhs.fiter );
		}

		bool end() { return ( fiter == fend ); }

	private:
		ETForest::const_iterator fiter;
		ETForest::const_iterator fend;
	};

	// iterator to iterate within connected components
	class IntraCCompIterator {
	public:
		IntraCCompIterator( ETNode * n );
		~IntraCCompIterator() {}

		inline VertexID operator*() {
			assert( cur->isLoop() );
			return cur->vx1->vxid;
		}

		inline IntraCCompIterator & operator++() { // prefix
			step();
			return *this;
		}

		IntraCCompIterator operator++( int ) { // postfix
			IntraCCompIterator result = *this;
			step();
			return result;
		}

		inline bool operator==( const IntraCCompIterator & rhs ) {
			return ( cur == rhs.cur );
		}
		inline bool operator!=( const IntraCCompIterator & rhs ) {
			return ( cur != rhs.cur );
		}

		inline bool end() const { return ( cur == NULL ); }

	private:
		void step();

	private:
		ETNode * cur;
	};
};

std::ostream & operator<<( std::ostream & os, const ConnectivityGraph::Edge & e );

typedef BOOSTNS::container::flat_set< ConnectivityGraph::VertexID > CGVertexIDSet;

} // namespace rdis

#endif // RDIS_CONNECTIVITYGRAPH_H_
