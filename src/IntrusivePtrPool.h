/*
 * IntrusivePtrPool.h
 *
 *  Created on: Aug 22, 2013
 *      Author: afriesen
 */

#ifndef INTRUSIVEPTRPOOL_H_
#define INTRUSIVEPTRPOOL_H_

#include "common.h"

#include BOOSTPATH/intrusive_ptr.hpp>

namespace rdis {

template< class Object, typename... ConstructorArgs >
class IntrusivePtrPool;

template< class Object >
void intrusive_ptr_add_ref( const Object * o );

template< class Object >
void intrusive_ptr_release( const Object * o );


// NOTE: IMPROPER ORDER OF RELEASE OF RESOURCES CAN RESULT IN MEMORY LEAKS AND
// CRASHES -- THIS IS NOT NECESSARILY A SAFE CLASS

template< class Object, typename PoolType >
class IntrusivePtrPoolObj {

public:
	IntrusivePtrPoolObj()
		: pool( NULL )
		, m_refCount( 0 )
	{}

	virtual ~IntrusivePtrPoolObj() {}

public:
	// clear member data so this object can be put in the pool and reused
	virtual void clear() = 0;


public:
	size_t refCount() const { return m_refCount; }
	bool refUnique() const { return m_refCount == 1; }

	bool hasPool() const { return ( pool != NULL ); }

protected:
	// TODO: These won't be safe with threading
	void ref() const { ++m_refCount; }
	void unref() const { --m_refCount; }

	template< class Object1 >
	friend void intrusive_ptr_add_ref( const Object1 * o );
	template< class Object1 >
	friend void intrusive_ptr_release( const Object1 * o );

	friend PoolType;


protected:
	// pointer to the pool of subdomains that this is from (and to return this
	// to on its "destruction")
	PoolType * pool;

	// used for boost intrusive_ptr
	mutable size_t m_refCount;
};



template< class Object, typename... ConstructorArgs >
class IntrusivePtrPool {
protected:
	typedef BOOSTNS::intrusive_ptr< Object > ObjectIP;
	typedef BOOSTNS::intrusive_ptr< const Object > ObjectCIP;

public:
	// a size of 0 means unlimited (except by memory)
	IntrusivePtrPool( long long int maxUnusedObjects = -1 )
		: maxsize( maxUnusedObjects ) {
		if ( maxsize > 0 ) pool.reserve( std::log( maxsize ) );
	}

	virtual ~IntrusivePtrPool() {
		while ( !pool.empty() ) {
			delete pool.back();
			pool.pop_back();
		}
	}

	virtual void setMaxSize( size_t maxUnusedObjects ) {
		maxsize = maxUnusedObjects;
		if ( maxsize > 0 ) pool.reserve( std::log( maxsize ) );
		while ( (long long int) pool.size() > maxsize ) {
			delete pool.back();
			pool.pop_back();
		}
	}


public:
	friend Object;


public:
	// get an object from this pool (pool gives up ownership)
	virtual ObjectIP get( ConstructorArgs... args ) {
		Object * o( NULL );

		if ( !pool.empty() ) {
			o = pool.back();
			pool.pop_back();
		} else {
			o = new Object( args... );
			o->pool = this;
		}

		return ObjectIP( o );
	}

	virtual ObjectCIP getc( ConstructorArgs... args ) {
		Object * o( NULL );

		if ( !pool.empty() ) {
			o = pool.back();
			pool.pop_back();
		} else {
			o = new Object( args... );
			o->pool = this;
		}

		return ObjectCIP( o );
	}

public:
	virtual size_t size() const { return pool.size(); }

protected:
	// release an object to this pool (pool takes ownership)
	virtual void release( Object * o ) {
		assert( o->refCount() == 0 );

		if ( maxsize < 0 || (long long int) pool.size() < maxsize ) {
			o->clear();
			pool.push_back( o );
		} else {
			delete o;
		}

	//	if ( maxsize == 0 && pool.size() > ( lastsize + 100 ) ) {
	//		std::cout << "pool size: " << pool.size() << std::endl;
	//		lastsize = pool.size();
	//	}
	}

protected:
	template< class Object1 >
	friend void intrusive_ptr_add_ref( const Object1 * o );
	template< class Object1 >
	friend void intrusive_ptr_release( const Object1 * o );

protected:
	long long int maxsize;

	std::vector< Object * > pool;

//	size_t lastsize;
};


template< class Object >
void intrusive_ptr_add_ref( const Object * o ) {
	o->ref();
}

template< class Object >
void intrusive_ptr_release( const Object * o ) {
	o->unref();
	if ( o->refCount() == 0 ) {
		// if there's a pool, release to it, otherwise delete
		if ( o->pool != nullptr ) {
			o->pool->release( const_cast< Object * >( o ) );
//			std::cout << "releasing IPPObj to pool : count is now " << o->pool->pool.size() << std::endl;
		} else {
			delete o;
		}
	}
}

} // namespace rdis

#endif // INTRUSIVEPTRPOOL_H_
