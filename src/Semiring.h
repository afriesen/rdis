/*
 * SumProduct.h
 *
 *  Created on: Apr 17, 2013
 *      Author: afriesen
 */

#ifndef RDIS_SEMIRING_H_
#define RDIS_SEMIRING_H_

#include <limits>
#include <algorithm>


//-----------------------------------------------------------------------------
// Use this by specifying the appropriate semiring struct as the return of the
// OptimizableFunction::getSemiring() function. Choices are:
//  - MinSum
//  - MinProduct
//  - MaxSum
//
//// ------- OUTDATED --------------
//// Use this by using the appropriate semiring namespace. Choices are:
////  - SumProduct
////	- MaxSum
////	- MaxProduct -- not supported
////	- MinSum
////	- MinProduct -- not supported
//// -------------------------------
//
// Then just call Semiring::Sum(x, y), Semiring::Product(x,y),
// Semiring::InverseProduct(x,y), Semiring::SumIdentity(), or
// Semiring::ProductIdentity().
//-----------------------------------------------------------------------------


namespace rdis {

// the different possible semiring operations
namespace SemiringOps {

// the available semiring operations
template< class T >
inline T Sum( const T & x, const T & y ) { return x + y; }
template< class T >
inline T Product( const T & x, const T & y ) { return x*y; }
template< class T >
inline T Max( const T & x, const T & y ) { return std::max( x, y ); }
template< class T >
inline T Min( const T & x, const T & y ) { return std::min( x, y ); }

// the available inverse operations
template< class T >
inline T Subtract( const T & x, const T & y ) { return x - y; }
template< class T >
inline T Divide( const T & x, const T & y ) { return x / y; }

// the identity elements of each type of operation
template< class T >
inline T SumIdentity() { return 0; }
template< class T >
inline T ProductIdentity() { return 1; }
template< class T >
inline T MaxIdentity() { return -std::numeric_limits< T >::max(); } // need C++11 to get ::lowest()
template< class T >
inline T MinIdentity() { return std::numeric_limits< T >::max(); }

} // namespace SemiringOps


// define the different semirings here

template< class T >
struct Semiring {
	virtual ~Semiring() {};
	virtual T SumIdentity() const = 0;
	virtual T ProductIdentity() const = 0;
	virtual T Sum( const T & x, const T & y ) const = 0;
	virtual T Product( const T & x, const T & y ) const = 0;
	virtual T InverseProduct( const T & x, const T & y ) const = 0;
	virtual Semiring * clone() const = 0;
};

template< class T >
struct SumProductSemiring
		: public Semiring< T > {
	virtual ~SumProductSemiring() {};
	inline T SumIdentity() const { return SemiringOps::SumIdentity< T >(); }
	inline T ProductIdentity() const { return SemiringOps::ProductIdentity< T >(); }
	inline T Sum( const T & x, const T & y ) const { return SemiringOps::Sum( x, y ); }
	inline T Product( const T & x, const T & y ) const { return SemiringOps::Product( x, y ); }
	inline T InverseProduct( const T & x, const T & y ) const { return SemiringOps::Divide( x, y ); }
	inline Semiring< T > * clone() const { return new SumProductSemiring< T >(); }
};

template< class T >
struct MinSumSemiring
		: public Semiring< T > {
	virtual ~MinSumSemiring() {};
	inline T SumIdentity() const { return SemiringOps::MinIdentity< T >(); }
	inline T ProductIdentity() const { return SemiringOps::SumIdentity< T >(); }
	inline T Sum( const T & x, const T & y ) const { return SemiringOps::Min( x, y ); }
	inline T Product( const T & x, const T & y ) const { return SemiringOps::Sum( x, y ); }
	inline T InverseProduct( const T & x, const T & y ) const { return SemiringOps::Subtract( x, y ); }
	inline Semiring< T > * clone() const { return new MinSumSemiring< T >(); }
};

template< class T >
struct MinProductSemiring
		: public Semiring< T > {
	virtual ~MinProductSemiring() {};
	inline T SumIdentity() const { return SemiringOps::MinIdentity< T >(); }
	inline T ProductIdentity() const { return SemiringOps::ProductIdentity< T >(); }
	inline T Sum( const T & x, const T & y ) const { return SemiringOps::Min( x, y ); }
	inline T Product( const T & x, const T & y ) const { return SemiringOps::Product( x, y ); }
	inline T InverseProduct( const T & x, const T & y ) const { return SemiringOps::Divide( x, y ); }
	inline Semiring< T > * clone() const { return new MinProductSemiring< T >(); }
};

template< class T >
struct MaxSumSemiring
		: public Semiring< T > {
	virtual ~MaxSumSemiring() {};
	inline T SumIdentity() const { return SemiringOps::MaxIdentity< T >(); }
	inline T ProductIdentity() const { return SemiringOps::SumIdentity< T >(); }
	inline T Sum( const T & x, const T & y ) const { return SemiringOps::Max( x, y ); }
	inline T Product( const T & x, const T & y ) const { return SemiringOps::Sum( x, y ); }
	inline T InverseProduct( const T & x, const T & y ) const { return SemiringOps::Subtract( x, y ); }
	inline Semiring< T > * clone() const { return new MaxSumSemiring< T >(); }
};

template< class T >
struct MaxProductSemiring
		: public Semiring< T > {
	virtual ~MaxProductSemiring() {};
	inline T SumIdentity() const { return SemiringOps::MaxIdentity< T >(); }
	inline T ProductIdentity() const { return SemiringOps::ProductIdentity< T >(); }
	inline T Sum( const T & x, const T & y ) const { return SemiringOps::Max( x, y ); }
	inline T Product( const T & x, const T & y ) const { return SemiringOps::Product( x, y ); }
	inline T InverseProduct( const T & x, const T & y ) const { return SemiringOps::Divide( x, y ); }
	inline Semiring< T > * clone() const { return new MaxProductSemiring< T >(); }
};

//namespace SumProduct {
//namespace Semiring {
//
//template< class T >
//inline T SumIdentity() { return SemiringOps::SumIdentity< T >(); }
//template< class T >
//inline T ProductIdentity() { return SemiringOps::ProductIdentity< T >(); }
//template< class T >
//inline T Sum( const T & x, const T & y ) { return SemiringOps::Sum( x, y ); }
//template< class T >
//inline T Product( const T & x, const T & y ) { return SemiringOps::Product( x, y ); }
//template< class T >
//inline T InverseProduct( const T & x, const T & y ) { return SemiringOps::Divide( x, y ); }
//
//} // namespace Semiring
//} // namespace SumProduct
//
//
//namespace MaxSum {
//namespace Semiring {
//
//template< class T >
//inline T SumIdentity() { return SemiringOps::MaxIdentity< T >(); }
//template< class T >
//inline T ProductIdentity() { return SemiringOps::SumIdentity< T >(); }
//template< class T >
//inline T Sum( const T & x, const T & y ) { return SemiringOps::Max( x, y ); }
//template< class T >
//inline T Product( const T & x, const T & y ) { return SemiringOps::Sum( x, y ); }
//template< class T >
//inline T InverseProduct( const T & x, const T & y ) { return SemiringOps::Subtract( x, y ); }
//
//} // namespace Semiring
//} // namespace MaxSum
//
//
////namespace MaxProduct {
////namespace Semiring {
////
////template< class T >
////T SRSumIdentity() { return SemiringOps::MaxIdentity< T >(); }
////template< class T >
////T SRProductIdentity() { return SemiringOps::ProductIdentity< T >(); }
////template< class T >
////T SRSum( const T & x, const T & y ) { return SemiringOps::Max( x, y ); }
////template< class T >
////T SRProduct( const T & x, const T & y ) { return SemiringOps::Product( x, y ); }
////
////} // namespace Semiring
////} // namespace MaxProduct
//
//
//namespace MinSum {
//namespace Semiring {
//
//template< class T >
//T SumIdentity() { return SemiringOps::MinIdentity< T >(); }
//template< class T >
//T ProductIdentity() { return SemiringOps::SumIdentity< T >(); }
//template< class T >
//T Sum( const T & x, const T & y ) { return SemiringOps::Min( x, y ); }
//template< class T >
//T Product( const T & x, const T & y ) { return SemiringOps::Sum( x, y ); }
//template< class T >
//T InverseProduct( const T & x, const T & y ) { return SemiringOps::Subtract( x, y ); }
//
//} // namespace Semiring
//} // namespace MaxSum
//
//
////namespace MinProduct {
////namespace Semiring {
////
////template< class T >
////T SumIdentity() { return SemiringOps::MinIdentity< T >(); }
////template< class T >
////T ProductIdentity() { return SemiringOps::ProductIdentity< T >(); }
////template< class T >
////T Sum( const T & x, const T & y ) { return SemiringOps::Min( x, y ); }
////template< class T >
////T Product( const T & x, const T & y ) { return SemiringOps::Product( x, y ); }
////
////} // namespace Semiring
////} // namespace MaxProduct

} // namespace rdis

#endif // RDIS_SEMIRING_H_
