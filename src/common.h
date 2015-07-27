// file containing common and shared definitions

#ifndef RDIS_COMMON_H_
#define RDIS_COMMON_H_

// *necessary* common includes only

#include "Semiring.h"

#include BOOSTPATH/numeric/interval.hpp>
#include BOOSTPATH/chrono/system_clocks.hpp>

#include <assert.h>
#include <string>
#include <vector>


using std::size_t;
using std::string;


namespace rdis {

// common typedefs
typedef double Numeric;

typedef BOOSTNS::chrono::steady_clock Clock;
typedef BOOSTNS::chrono::duration< Numeric > Duration; // in seconds, floating point

typedef long long int VariableCount;
typedef long long int VariableID;
typedef long long int FactorID;

typedef std::vector< size_t > SizeVec;
typedef std::vector< Numeric > NumericVec;
typedef std::vector< NumericVec > NumericVecVec;
typedef std::vector< VariableID > VariableIDVec;

// useful constants
extern const Numeric NumericMAX;
extern const Numeric NumericMIN;

// faster interval arithmetic (more errors though), without rounding/checking
namespace IntervalDefines {

using namespace BOOSTNS::numeric;
typedef interval_lib::rounded_transc_exact< Numeric >
	ExactNumericRounding;
typedef interval_lib::save_state_nothing< ExactNumericRounding >
	SaveNothingExactNumericRounding;

typedef interval_lib::checking_base< Numeric > CheckingBaseNumeric;

typedef interval_lib::policies< SaveNothingExactNumericRounding,
		CheckingBaseNumeric > PolicyNoRoundingNoChecking;
} // namespace IntervalDefines

// without rounding or checking (higher performance)
typedef BOOSTNS::numeric::interval< Numeric,
		IntervalDefines::PolicyNoRoundingNoChecking > NumericInterval;
//// with rounding + checking
//typedef BOOSTNS::numeric::interval< Numeric > NumericInterval;



inline bool approxeq( const Numeric & x, const Numeric & y, Numeric tol = 1e-10 ) {
	return ( std::abs( x - y ) < tol );
}

inline bool approxleq( const Numeric & x, const Numeric & y, Numeric tol = 1e-10 ) {
	return ( ( y - x ) >= -tol );
}

inline bool approxgeq( const Numeric x, const Numeric y, Numeric tol = 1e-10 ) {
	return ( x >= ( y - tol ) );
}

inline bool approxin( const NumericInterval & x, const NumericInterval & y, Numeric tol = 1e-10 ) {
	return ( approxleq( y.lower(), x.lower(), tol ) &&
			approxleq( x.upper(), y.upper(), tol ) );
}


inline bool approxeq( const NumericVec & x, const NumericVec & y, Numeric tol = 1e-10 ) {
	assert( x.size() == y.size() );
	for ( size_t i = 0; i < x.size(); ++i ) {
		if ( !approxeq( x[i], y[i], tol ) ) return false;
	}
	return true;
}


template< class T >
std::ostream & operator<<( std::ostream & out, const std::vector< T > & v ) {
	for ( size_t i = 0; i < v.size(); ++i ) {
		out << ( i == 0 ? "" : ", " ) << v[i];
	}
	return out;
}


namespace Heuristic {
enum Type {
	MOMS3 = 0,
	MOMS4,
	MOMS5,
	BottomUp,
	VarPriorities,

	NumHeuristics
};
} // namespace Heuristic


namespace instrumentation {

// global variables containing the number of operations performed during an
// integration (used for instrumentation)
extern size_t g_numFactEvals;
extern size_t g_numVarOps;
extern size_t g_numProducts;
extern size_t g_numSums;
//extern size_t g_numExp;
//extern size_t g_numMax;
extern Duration g_runTime;

extern size_t g_numDownhillSteps;
extern size_t g_numSteps;

extern size_t g_cacheInserts;
extern size_t g_cacheHits;
extern size_t g_cacheLookups;
extern size_t g_cacheMidInserts;

extern size_t g_counter1;
extern size_t g_counter2;

extern void clearInstrumentation();

} // namespace instrumentation


inline std::ostream & operator<<( std::ostream & os, NumericInterval ni) {
	os << "[" << ni.lower() << "," << ni.upper() << "]";
	return os;
}


template< class Container >
void printContainer( const Container & v ) {
	for ( auto it( v.begin() ), end( v.end() );
			it != end; ++it ) {
		std::cout << ( it == v.begin() ? "" : ", " ) << *it;
	}
}


template< class Container >
void printPtrContainer( const Container & c ) {
	for ( auto it( c.begin() ); it != c.end();
			++it ) {
		std::cout << ( it == c.begin() ? "" : ", " ) << **it;
	}
}

// Floating-point modulo -- from stackoverflow
// The result (the remainder) has same sign as the divisor.
// Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
//template<typename T>
inline Numeric Mod(Numeric x, Numeric y) {
//    static_assert(!std::numeric_limits<T>::is_exact , "Mod: floating-point type expected");

    if ( 0. == y ) return x;

    double m = x - y * floor(x/y);

    // handle boundary cases resulted from floating-point cut off:

    if (y > 0) {            // modulo range: [0..y)
        if ( m >= y ) return 0; // Mod(-1e-16             , 360.    ): m= 360.
        if ( m < 0 ) return ( y+m == y ? 0 : y+m );
    } else {                // modulo range: (y..0]
        if ( m <= y ) return 0; // Mod(1e-16              , -360.   ): m= -360.
        if ( m > 0 ) return ( y+m == y ? 0 : y+m );
    }

    return m;
}


//// The result (the remainder) has same sign as the divisor.
//// Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
//inline NumericInterval ModItvl( NumericInterval x, Numeric y) {
//    if ( 0. == y ) return x;
//
//	x = BOOSTNS::numeric::fmod( x, y );
//
//    double ml = x.lower() - y * floor( x.lower() / y );
//    double mu = x.upper() - y * floor( x.upper() / y );
//
//    // handle boundary cases resulted from floating-point cut off:
//
//    if (y > 0) {            // modulo range: [0..y)
//        if ( ml >= y ) ml = 0.0; // Mod(-1e-16             , 360.    ): m= 360.
//        if ( ml < 0 ) {
//        	if ( y+ml == y ) {
//        		ml = 0.0;
//        	} else {
//        		ml += y;
//        		mu += y;
//        	}
////        	return ( y+m == y ? 0 : y+m );
//        }
//
//        if ( mu >= y ) mu = 0.0;
//        if ( mu < 0 ) mu = ( y+mu == y ? 0 : y+mu );
//
//    } else {                // modulo range: (y..0]
//        if ( ml <= y ) ml = 0.0; // Mod(1e-16              , -360.   ): m= -360.
//        if ( ml > 0 ) {
//        	if ( y+ml == y ) {
//        		ml = 0.0;
//        	} else {
//        		ml += y;
//        		mu += y;
//        	}
////        	return ( y+m == y ? 0 : y+m );
//        }
//
//        if ( mu <= y ) mu = 0.0;
//        if ( mu > 0 ) mu = ( y+mu == y ? 0 : y+mu );
//    }
//
//    if ( mu < ml ) mu += y;
//
//    return NumericInterval( ml, mu );
//}

} // namespace rdis

#endif // RDIS_COMMON_H_





