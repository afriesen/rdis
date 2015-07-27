/*
 * SplayTree.h
 *
 *  Created on: Jun 14, 2013
 *      Author: afriesen
 */

#ifndef RDIS_SPLAYTREE_H_
#define RDIS_SPLAYTREE_H_

#include BOOSTPATH/container/vector.hpp>

#include <iostream>
#include <functional>

namespace rdis {

// the splaytree class (and those derived from it) should be lightweight
// wrappers around the root node
template< typename T, typename Weight = double, typename Comp = std::less<T>,
	int NumWeights = 1 >
class SplayTree {

public:
	struct SplayNode;
	typedef BOOSTNS::container::vector< SplayNode * > SplayNodeStack;

	struct SplayNode {
	public:
		SplayNode( const T & init = T(), const Weight & weight_ = 0.0 )
				: parent( NULL )
				, left( NULL )
				, right( NULL )
				, key( init )
				, stsize( 1 )
		{
			for ( int i = 0; i < NumWeights; ++i ) {
				weight[i] = weight_;
				stweight[i] = weight[i];
			}
		}

		virtual ~SplayNode() {}


	public:
		// split to the left of this node and return the new tree (which was
		// the left subtree of this node)
		virtual SplayNode * splitLeft() {
			SplayNode * z( left );
			left = NULL;
			if ( z ) {
				z->parent = NULL;
				stsize -= z->stsize;
				for ( int i = 0; i < NumWeights; ++i ) {
					stweight[i] -= z->stweight[i];
				}
			}
			return z;
		}

		// split to the right of this node and return the new tree (which was
		// the right subtree of this node)
		virtual SplayNode * splitRight() {
			SplayNode * z( right );
			right = NULL;
			if ( z ) {
				z->parent = NULL;
				stsize -= z->stsize;
				for ( int i = 0; i < NumWeights; ++i ) {
					stweight[i] -= z->stweight[i];
				}
			}
			return z;
		}


	public:
		virtual void joinLeft( SplayNode * leftRoot ) {
			assert( !left && !leftRoot->parent );
			left = leftRoot;
			if ( left ) {
				stsize += left->stsize;
				for ( int i = 0; i < NumWeights; ++i ) {
					stweight[i] += left->stweight[i];
				}
				left->parent = this;
			}
		}

		virtual void joinRight( SplayNode * rightRoot ) {
			assert( !right && !rightRoot->parent );
			right = rightRoot;
			if ( right ) {
				stsize += right->stsize;
				for ( int i = 0; i < NumWeights; ++i ) {
					stweight[i] += right->stweight[i];
				}
				right->parent = this;
			}
		}


	public:
		// return the next node following this one (or null if there aren't any)
		// in an in-order traversal of the tree
		virtual SplayNode * next() const {
			SplayNode * nxt( NULL );

			if ( right ) {
				nxt = right->subtree_minimum();
			} else if ( parent ) {
				if ( parent->left == this ) {
					nxt = parent;
				} else {
					assert( parent->right == this );
					for ( SplayNode * n = parent; n->parent; n = n->parent ) {
						if ( n == n->parent->left ) {
							nxt = n->parent;
							break;
						}
					}
				}
			}

			return nxt;
		}

		// return the node that precedes the specified node (or null) in an
		// in-order traversal of the tree
		virtual SplayNode * prev() const {

			SplayNode * prv( NULL );

			if ( left ) {
				prv = left->subtree_maximum();
			} else if ( parent ) {
				if ( parent->right == this ) {
					prv = parent;
				} else {
					assert( parent->left == this );
					for ( SplayNode * n = parent; n->parent; n = n->parent ) {
						if ( n == n->parent->right ) {
							prv = n->parent;
							break;
						}
					}
//					for ( rnbnode aux = par; aux->par; aux = aux->par ) {
//						if ( aux == aux->par->child[rnbright] ) return aux->par;
//					}
				}
			}

			return prv;
		}

	public:
		// return the next node in a pre-order traversal of this node's tree,
		// using the stack to maintain state across calls to this function (can
		// also use to traverse only those nodes with positive weights)
		virtual SplayNode * nextPreOrder( SplayNodeStack & stack,
				bool onlyWeighted = false, int weightIdx = 0 ) const {

			SplayNode * nxt( NULL );

			// remember the right node for later
			if ( right != NULL && ( !onlyWeighted ||
					right->stweight[weightIdx] > 0 ) ) {
				stack.push_back( right );
			}

			// try to move to the left node now, or else get one from the stack
			if ( left != NULL &&
					( !onlyWeighted || left->stweight[weightIdx] > 0 ) ) {
				nxt = left;

			} else if ( !stack.empty() ) {
				nxt = stack.back();
				stack.pop_back();
			}

			assert( nxt != NULL || stack.empty() );

			return nxt;
		}

	public:
		virtual void graphsubtree( std::ostream & os, size_t & curid ) const {
			size_t myid = curid;
			os << myid << " ; " << std::endl;

			if ( this->left ) {
				++curid;
				os << myid << " -- " << curid << " ; " << std::endl;
				left->graphsubtree( os, curid );
			}

			if ( this->right ) {
				++curid;
				os << myid << " -- " << curid << " ; " << std::endl;
				right->graphsubtree( os, curid );
			}
		}

		virtual void printsubtree( std::ostream & os ) const {
			if ( this->left == NULL && this->right == NULL && this->parent == NULL ) {
				os << key << std::endl;
			}

			if ( this->left ) {
				os << key << " --> " << this->left->key << " (L)" << std::endl;
				left->printsubtree( os );
			}

			if ( this->right ) {
				os << key << " --> " << this->right->key << " (R)" << std::endl;
				right->printsubtree( os );
			}
		}

		// return the root of this node's tree
		virtual const SplayNode * findRoot() const {
			const SplayNode * z = this;

			while ( z->parent ) { z = z->parent; }
			return z;
		}

		virtual SplayNode * findRoot() {
			SplayNode * z = this;
			while ( z->parent ) { z = z->parent; }
			return z;
		}

		// descend the left branch until we find the minimum node
		virtual const SplayNode * subtree_minimum() const {
			const SplayNode * u( this );
			while ( u->left ) { u = u->left; }
			return u;
		}

		virtual SplayNode * subtree_minimum() {
			SplayNode * u( this );
			while ( u->left ) { u = u->left; }
			return u;
		}

		// descend the right branch until we find the maximum node
		virtual const SplayNode * subtree_maximum() const {
			const SplayNode * u( this );
			while ( u->right ) { u = u->right; }
			return u;
		}

		virtual SplayNode * subtree_maximum() {
			SplayNode * u( this );
			while ( u->right ) { u = u->right; }
			return u;
		}

	public:
		virtual void changeWeight( Weight wnew, int widx ) {
			Weight diff( wnew - weight[widx] );
			weight[widx] = wnew;
			stweight[widx] += diff;

			assert( weight[widx] >= 0 && stweight[widx] >= 0 );
			assert( checkSizeAndWeight() );

			// walk up the tree and adjust all subtree weights accordingly
			SplayNode * n( this );
			while ( ( n = n->parent ) ) {
				n->stweight[widx] += diff;
//				assert( n->weight[widx] >= 0 && n->stweight[widx] >= 0 );
				n->checkSizeAndWeight();
			}
		}

		virtual void addWeight( Weight wdelta, int widx ) {
			changeWeight( weight[widx] + wdelta, widx );
		}

		inline virtual bool checkSizeAndWeight() const {
#ifdef DEBUG
			bool sz( stsize == ( 1 + ( left ? left->stsize : 0 ) +
					( right ? right->stsize : 0 ) ) );

			bool wt( true );
			for ( int i = 0; i < NumWeights; ++i ) {
				wt = wt && ( stweight[i] == ( weight[i] +
					( left ? left->stweight[i] : 0 ) +
					( right ? right->stweight[i] : 0 ) ) );

				assert( weight[i] >= 0 && weight[i] <= 1000000 );
			}

			return ( sz && wt );
#else
			return true;
#endif // DEBUG
		}


	public:
		// splay this node to the top of the tree that it's currently in and
		// return this node (for convenience when stringing together node ops.)
		// NOTE: This is not a very clean implementation, as splaying a node
		// makes it the new root of its tree, but this can't actually update an
		// instance of SplayTree if it exists, so that will no longer point to
		// the correct root if a node is splayed without that being updated, so
		// use this carefully.
		virtual SplayNode * splay() {
			SplayNode * x = this;
			while (x->parent) {
				if (!x->parent->parent) {
					if (x->parent->left == x) x->parent->right_rotate();
					else x->parent->left_rotate();

				} else if (x->parent->left == x
						&& x->parent->parent->left == x->parent) {
					x->parent->parent->right_rotate();
					x->parent->right_rotate();

				} else if (x->parent->right == x
						&& x->parent->parent->right == x->parent) {
					x->parent->parent->left_rotate();
					x->parent->left_rotate();

				} else if (x->parent->left == x
						&& x->parent->parent->right == x->parent) {
					x->parent->right_rotate();
					x->parent->left_rotate();

				} else {
					x->parent->left_rotate();
					x->parent->right_rotate();
				}
			}

			assert( checkSizeAndWeight() );
			assert( this->parent == NULL );
//			assert( root == x );
			return this;
		}

	protected:
		virtual void left_rotate() {
			SplayNode *x = this;
			SplayNode *y = x->right;

			x->right = y->left;
			if (y->left) y->left->parent = x;

			y->parent = x->parent;
//			if (!x->parent) root = y;
			if ( x->parent ) {
				if (x == x->parent->left) x->parent->left = y;
				else x->parent->right = y;
			}

			// update sizes before setting y->left
			x->stsize -= y->stsize - ( y->left ? y->left->stsize : 0 );
			y->stsize += 1 + ( x->left ? x->left->stsize : 0 );

			y->left = x;
			x->parent = y;

			// update weights
			for ( int i = 0; i < NumWeights; ++i ) {
				x->stweight[i] = x->weight[i] +
						( x->left ? x->left->stweight[i] : 0 ) +
						( x->right ? x->right->stweight[i] : 0 );
				y->stweight[i] = y->weight[i] +
						( y->left ? y->left->stweight[i] : 0 ) +
						( y->right ? y->right->stweight[i] : 0 );
			}

			assert( !x || x->checkSizeAndWeight() );
			assert( !y || y->checkSizeAndWeight() );
			assert( !y->parent || y->parent->checkSizeAndWeight() );
		}

		virtual void right_rotate() {
			SplayNode *x = this;
			SplayNode *y = x->left;
			x->left = y->right;
			if (y->right) y->right->parent = x;
			y->parent = x->parent;

//			if (!x->parent) root = y;
			if ( x->parent ) {
				if (x == x->parent->left) x->parent->left = y;
				else x->parent->right = y;
			}

			// update the sizes before setting y->right
			x->stsize -= y->stsize - ( y->right ? y->right->stsize : 0 );
			y->stsize += 1 + ( x->right ? x->right->stsize : 0 );

			y->right = x;
			x->parent = y;

			// update weights after rotation
			for ( int i = 0; i < NumWeights; ++i ) {
				x->stweight[i] = x->weight[i] +
						( x->left ? x->left->stweight[i] : 0 ) +
						( x->right ? x->right->stweight[i] : 0 );
				y->stweight[i] = y->weight[i] +
						( y->left ? y->left->stweight[i] : 0 ) +
						( y->right ? y->right->stweight[i] : 0 );
			}

			assert( !x || x->checkSizeAndWeight() );
			assert( !y || y->checkSizeAndWeight() );
			assert( !y->parent || y->parent->checkSizeAndWeight() );
		}

	public:
		virtual void disconnect() {
			// NOTE: This doesn't change the nodes connected to this node!
			parent = left = right = NULL;
		}

		virtual bool isDisconnected() const {
			return ( parent == NULL && left == NULL && right == NULL );
		}

	public:
		virtual void print( std::ostream & os = std::cout ) const {
			os << key;
		}

	public:
		// pointers to the attached nodes in the tree
		SplayNode *parent;
		SplayNode *left, *right;

		// the key contained at this node
		T key;

		// the weight(s) of this node
		Weight weight[NumWeights];

		// the weight(s) of the subtree rooted at this node (includes this node)
		Weight stweight[NumWeights];

		// the size of the subtree rooted at this node (includes this node)
		std::size_t stsize;
	}; // struct SplayNode


protected:
	typedef SplayTree< T, Weight, Comp > SplayTree_t;

	typedef void (SplayTree::*SplayTreeMemFn)( SplayNode * u );

public:
	SplayTree()
		: root( NULL )
	{}

	virtual ~SplayTree() {}

	// delete all of the nodes in this tree
	static void destroy( SplayNode * r ) {
//		SplayTree_t::postOrderTraverse( root, &SplayTree::deleteNode );
		SplayTree::postOrderDelete( r );
	}


	virtual void insert( const T & key, const Weight & w = 0.0 ) {
		// move the closest node to the root, split the tree, and then insert
		// the new node as the new root
		SplayNode * z = new SplayNode( key, w );

		SplayNode * rn = NULL;
		if ( root != NULL ) rn = split( key );

		z->left = root;
		if ( root ) root->parent = z;

		z->right = rn;
		if ( rn ) rn->parent = z;

		root = z;

		if ( z->left ) {
			for ( int i = 0; i < NumWeights; ++i ) {
				z->stweight[i] += z->left->stweight[i];
			}
			z->stsize += z->left->stsize;
		}

		if ( z->right ) {
			for ( int i = 0; i < NumWeights; ++i ) {
				z->stweight[i] += z->right->stweight[i];
			}
			z->stsize += z->right->stsize;
		}

		assert( !root || root->checkSizeAndWeight() );
	}


	// return the node that matches the specified key
	virtual SplayNode * find( const T & key ) {
		return access( key );
	}


	// splay the node we're erasing to the top, remove it, and then connect the
	// resulting subtrees such that the min element of the right subtree is the
	// new root
	virtual void erase( const T & key ) {

		// find the node to delete
		SplayNode * z = access( key );
		if ( z ) {
			assert( !z->parent );

			// split away the left subtree
			SplayNode * lr = z->splitLeft();
			assert( z != lr );

			// split away the right subtree
			SplayNode * rr = z->splitRight();
			assert( z != rr );

			// re-connect the two trees
			assert( root == z );
			root = lr;
			join( rr );

			// delete the node
			assert( !comp( z->key, key ) && !comp( key, z->key ) );
			delete z;
		}
	}

	virtual const SplayNode * const getroot() const { return root; }

	virtual const T& minimum() const {
		return root->subtree_minimum()->key;
	}

	virtual const T& maximum() const {
		return root->subtree_maximum()->key;
	}

	virtual bool empty() const {
		return ( root == NULL );
	}

	virtual unsigned long size() const {
		return ( !root ? 0 : root->stsize );
	}

protected:

	// return a pointer to the node with the specified key (or null) -- splay
	// the found node, or the last accessed node if it wasn't found
	virtual SplayNode * access( const T & key ) {
		bool found( false );
		SplayNode * z = root;

		while ( z ) {
			if ( comp(z->key, key) ) {
				if ( z->right ) z = z->right;
				else break;

			} else if ( comp(key, z->key) ) {
				if ( z->left ) z = z->left;
				else break;

			} else {
				found = true;
				break;
			}
		}

		z->splay();
		root = z;
		assert( !found || ( !comp( z->key, key ) && !comp( key, z->key ) ) );

		return ( found ? z : NULL );
	}


	// split this tree into t1,t2 at the specified key, returning t2 (where t2
	// has keys greater than those in t1)
	virtual SplayNode * split( const T & key ) {
		SplayNode * n( NULL );
		SplayNode * z = access( key );
		if ( z ) {
			assert( z == root );
			assert( root == NULL || root->right == NULL || comp(key, root->right->key) );
			n = root;
			root = root->splitLeft();
		} else if ( comp(key, root->key) ) {
			n = root;
			root = root->splitLeft();
		} else {
			n = root->splitRight();
		}
		return n;
	}


	// join <rightTree> to this tree by splaying the max node and then inserting
	// <rightTree> as the right subtree of the splayed node
	virtual void join( SplayNode * rightRoot ) {
		SplayNode * z = root->subtree_maximum();
		z->splay();
		root = z;

		assert( !rightRoot->parent );
		assert( !z->right );
		z->joinRight( rightRoot );
	}


	static void postOrderDelete( SplayNode * root ) {
		if ( !root ) return;
		SplayNode * cur = root;
		SplayNode * prev = NULL;

		while ( cur ) {
			if ( prev == cur->parent || prev == NULL ) { // traversing down the tree
				prev = cur;
				if ( cur->left )
					cur = cur->left;
				else if ( cur->right )
					cur = cur->right;
				else
					cur = cur->parent; // reached the bottom

			} else if ( prev == cur->left ) { // traversing up from the left
				// perform the visit of the node we are coming from (left)
//				( this->*(onVisit) )( prev );
				prev->disconnect();
				delete prev;

				prev = cur;
				if ( cur->right ) cur = cur->right;
				else cur = cur->parent;

			} else if ( prev == cur->right ) { // traversing up from the right
				// perform the visit of the node we are coming from (right)
//				( this->*(onVisit) )( prev );
				prev->disconnect();
				delete prev;

				prev = cur;
				cur = cur->parent;
			} else {
				assert( false );
			}
		}

		// finally, visit the root
		assert( prev == root );
		prev->disconnect();
//		( this->*(onVisit) )( prev );
		delete prev;
	}


protected:

	// the root node of this tree
	SplayNode * root;

	// comparator used to order the nodes of this tree
	Comp comp;
};

} // namespace rdis

#endif // RDIS_SPLAYTREE_H_
