/*
 * EulerTourTree.cpp
 *
 *  Created on: Jun 14, 2013
 *      Author: afriesen
 */

#include "EulerTourTree.h"

namespace rdis {


void ETNode::print( std::ostream & os ) const {
	os << "(";

	if ( isLoop() ) os << vx1->vxid;
	else os << vx1->vxid << "," << vx2->vxid;

	os << ")";
}


void ETNode::clear() {
	assert( isDisconnected() );
	vx1 = NULL;
	vx2 = NULL;
}


std::ostream & operator<<( std::ostream & os, const ETNode & n ) {
	n.print( os );
	return os;
}


void EulerTourTree::createTour( ETVertex * v ) {
	assert( v != NULL && v->loopNode == NULL );
//	v->loopNode = EulerTourTree::nodepool.get( v, v );
	v->loopNode = new ETNode( v, v );
}


void EulerTourTree::destroyTour( ETNode * root ) {
	ETSplayTree::destroy( root );
}


void EulerTourTree::cut( ETNode * n1, ETNode * n2,
		ETNode ** r1, ETNode ** r2 ) {
	assert( n1->findRoot() == n2->findRoot() );
	assert( !n1->isLoop() && !n2->isLoop() );
	assert( n1->vx1 == n2->vx2 && n1->vx2 == n2->vx1 );

	// cut out the two edges in the tour (represented by n1 and n2)

	n1->splay();
	assert( n1->parent == NULL );
//	std::cout << "splitting 1 "; n1->print(); std::cout << " from ";
//	if ( n1->left ) n1->left->print();
//	else std::cout << "null";
//	std::cout << " and ";
//	if ( n1->right ) n1->right->print();
//	else std::cout << "null";
//	std::cout << std::endl;
	ETNode * l1 = (ETNode *) n1->splitLeft();
	ETNode * l2 = (ETNode *) n1->splitRight();

	// get the actual previous/next node (not just the root of the pre-split subtree)
	if ( l1 != NULL ) l1 = (ETNode *) l1->subtree_maximum();
	if ( l2 != NULL ) l2 = (ETNode *) l2->subtree_minimum();

	n2->splay();
	assert( n2->parent == NULL );
//	std::cout << "splitting 2 "; n2->print(); std::cout << " from ";
//	if ( n2->left ) n2->left->print();
//	else std::cout << "null";
//	std::cout << " and ";
//	if ( n2->right ) n2->right->print();
//	else std::cout << "null";
//	std::cout << std::endl;
	ETNode * l3 = (ETNode *) n2->splitLeft();
	ETNode * l4 = (ETNode *) n2->splitRight();

	// get the actual previous/next node (not just the root of the pre-split subtree)
	if ( l3 != NULL ) l3 = (ETNode *) l3->subtree_maximum();
	if ( l4 != NULL ) l4 = (ETNode *) l4->subtree_minimum();

	// check whether n1 and n2 were in the correct order or reversed, such that
	// their order in the tour was actually ..., n2, ..., n1, ...
	if ( !l2 || !l3 || l2->findRoot() != l3->findRoot() ) {
		assert( l1->findRoot() == l4->findRoot() );
		std::swap( l1, l3 );
		std::swap( l2, l4 );
	}

	// join l1 and l4 (l2 and l3 should already be the same subtree)
	l1 = EulerTourTree::join( l1, l4 );

	// delete these nodes that are no longer needed
	delete n1;
	delete n2;

	// return the roots of the two trees
	*r1 = l1->findRoot();
	*r2 = l3->findRoot();
}


ETNode * EulerTourTree::link( ETVertex * u, ETVertex * v,
		ETNode ** neuv, ETNode ** nevu ) {
	ETNode * nluu = u->loopNode;
	ETNode * nlvv = v->loopNode;

	assert( nluu != NULL && nlvv != NULL );
	assert( !EulerTourTree::connected( nluu, nlvv ) );
	assert( *neuv == NULL && *nevu == NULL );

	// split each tour (represented as a list embedded in a BST), after each
	// loop node, to get L1 -> (L11, L12) and L2 -> (L21, L22)

	nluu->splay();
	ETNode * l11 = nluu;
	ETNode * l12 = (ETNode *) nluu->splitRight();

	nlvv->splay();
	ETNode * l21 = nlvv;
	ETNode * l22 = (ETNode *) nlvv->splitRight();

	// create the new nodes for the two arcs representing the new edge
	*neuv = new ETNode( u, v );
	*nevu = new ETNode( v, u );

	// concat as: ( l12, l11, neuv, l22, l21, nevu ) -- takes 5 joins

	ETNode * tmp1 = EulerTourTree::join( l11, *neuv );
	ETNode * tmp2 = EulerTourTree::join( l21, *nevu );

	ETNode * rl = EulerTourTree::join( l12, tmp1 );
	ETNode * rr = EulerTourTree::join( l22, tmp2 );

	return EulerTourTree::join( rl, rr );
}


bool EulerTourTree::connected( const ETVertex * u, const ETVertex * v ) {
	return EulerTourTree::connected( u->loopNode, v->loopNode );
}


bool EulerTourTree::connected( const ETNode * u, const ETNode * v ) {
	return ( u->findRoot() == v->findRoot() );
}


void EulerTourTree::print( const ETNode * n ) {

	VariableCount count( 0 );
	n = (ETNode *) n->findRoot()->subtree_minimum();

	assert( n != NULL );

	while ( n != NULL ) {

		std::cout << ( count > 0 ? ", " : "" );
		n->print( std::cout );

		assert( n->checkSizeAndWeight() );

		++count;
		n = (ETNode *) n->next();
	}

//	std::cout << std::endl;
}


ETNode * EulerTourTree::join( ETNode * p, ETNode * q ) {
	if ( p == NULL ) return q->findRoot();
	if ( q == NULL ) return p->findRoot();
	p = (ETNode *) p->findRoot()->subtree_maximum()->splay();
	assert( p->right == NULL );
	p->joinRight( q->findRoot() );
	return p;
}


} // namespace rdis
