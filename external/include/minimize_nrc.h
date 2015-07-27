/*
 * minimize_nrc.h
 *
 *  Created on: Jan 16, 2014
 *      Author: afriesen
 */

#ifndef RDIS_MINIMIZE_NRC_H_
#define RDIS_MINIMIZE_NRC_H_

#include <vector>
#include <math.h>

/*
 * this file contains code for minimization of functions taken from Numerical
 * Recipes in C, 3rd edition
 */

namespace rdis {
namespace nrc {

//template< class T >
//typedef < T > NRVector;

typedef bool Bool;
typedef int Int;
typedef Numeric Doub;

typedef std::vector< Doub > VecDoub, VecDoub_O, VecDoub_IO;
typedef const std::vector< Doub > VecDoub_I;

//////// forward declarations ///////
struct Bracketmethod;
struct Brent;
struct Dbrent;

template< class T >
struct Linemethod;

template< class T >
struct F1dim;

template< class T >
struct Dlinemethod;

template< class T >
struct Df1dim;


template< class T >
inline void SWAP( T & t1, T & t2 ) {
	std::swap( t1, t2 );
}

template< class T >
inline T SIGN( const T & t1, const T & t2 ) {
	return copysign( t1, t2 );
}

#ifndef MAX
template< class T >
inline const T & MAX( const T & t1, const T & t2 ) {
	return std::max( t1, t2 );
}
#endif




// Base class for one-dimensional minimization routines. Provides a routine to
// bracket a minimum and several utility functions.
struct Bracketmethod {
	Doub ax, bx, cx, fa, fb, fc;

// Given a function or functor func, and given distinct initial points ax and
// bx, this routine searches in the downhill direction (defined by the function
// as evaluated at the initial points) and returns new points ax, bx, cx that
// bracket a minimum of the function. Also returned are the function values at
// the three points, fa, fb, and fc.
	template< class T >
	void bracket( const Doub a, const Doub b, T &func ) {
		const Doub GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
		//Here GOLD is the default ratio by which successive intervals are
		// magnified and GLIMIT is the maximum magnification allowed for a
		// parabolic-fit step.
		ax = a;
		bx = b;
		Doub fu;
		fa = func( ax );
		fb = func( bx );

		// Switch roles of a and b so that we can go downhill in the
		// direction from a to b.
		if ( fb > fa ) {
			SWAP( ax, bx );
			SWAP( fb, fa );
		}

		// First guess for c.
		cx = bx + GOLD * ( bx - ax );
		fc = func( cx );

		while ( fb > fc ) {
			// Keep returning here until we bracket. Compute u by parabolic
			// extrapolation
			// from a; b; c. TINY is used to prevent any possible division by zero.
			Doub r = ( bx - ax ) * ( fb - fc );
			Doub q = ( bx - cx ) * ( fb - fa );
			Doub u = bx - ( ( bx - cx ) * q - ( bx - ax ) * r )
					/ ( 2.0 * SIGN( MAX( std::abs( q - r ), TINY ),	q - r ) );
			Doub ulim = bx + GLIMIT * ( cx - bx );

			//We won’t go farther than this. Test various possibilities:
			if ( ( bx - u ) * ( u - cx ) > 0.0 ) {
				//Parabolic u is between b and c: try it.
				fu = func( u );
				if ( fu < fc ) {  	//Got a minimum between b and c.
					ax = bx;
					bx = u;
					fa = fb;
					fb = fu;
					return;
				} else if ( fu > fb ) {	//Got a minimum between between a and u.
					cx = u;
					fc = fu;
					return;
				}

				u = cx + GOLD * ( cx - bx ); //Parabolic fit was no use. Use default magnification.
				fu = func( u );
			} else if ( ( cx - u ) * ( u - ulim ) > 0.0 ) {
				// Parabolic fit is between c and its allowed limit.
				fu = func( u );
				if ( fu < fc ) {
					shft3( bx, cx, u, u + GOLD * ( u - cx ) );
					shft3( fb, fc, fu, func( u ) );
				}
			} else if ( ( u - ulim ) * ( ulim - cx ) >= 0.0 ) {
				// Limit parabolic u to maximum allowed value.
				u = ulim;
				fu = func( u );
			} else {
				//Reject parabolic u, use default magnification.
				u = cx + GOLD * ( cx - bx );
				fu = func( u );
			}
			//Eliminate oldest point and continue.
			shft3( ax, bx, cx, u );
			shft3( fa, fb, fc, fu );
		}
	}

	// Utility function used in this structure or others derived from it.
	inline void shft2( Doub &a, Doub &b, const Doub c )
			{
		a = b;
		b = c;
	}
	inline void shft3( Doub &a, Doub &b, Doub &c, const Doub d ) {
		a = b;
		b = c;
		c = d;
	}
	inline void mov3( Doub &a, Doub &b, Doub &c, const Doub d, const Doub e,
			const Doub f ) {
		a = d;
		b = e;
		c = f;
	}
};

// Brent’s method to find a minimum.
struct Brent : Bracketmethod {
	Doub xmin, fmin;
	const Doub tol;
	Brent( const Doub toll = 3.0e-8 )
			: tol( toll ) {
	}

	// Given a function or functor f, and given a bracketing triplet of
	// abscissas ax, bx, cx (such that bx is between ax and cx, and f(bx) is
	// less than both f(ax) and f(cx)), this routine isolates the minimum to a
	// fractional precision of about tol using Brent’s method. The abscissa of
	// the minimum is returned as xmin, and the function value at the minimum is
	// returned as min, the returned function value.
	template< class T >
	Doub minimize( T &func ) {
		const Int ITMAX = 100;
		const Doub CGOLD = 0.3819660;
		const Doub ZEPS = std::numeric_limits< Doub >::epsilon() * 1.0e-3;
		// Here ITMAX is the maximum allowed number of iterations; CGOLD is the
		// golden ratio; and ZEPS is a small number that protects against for a
		// minimum that happens to be exactly zero.

		Doub a, b, d = 0.0, etemp, fu, fv, fw, fx;
		Doub p, q, r, tol1, tol2, u, v, w, x, xm;

		//This will be the distance moved on the step before last.
		Doub e = 0.0;

		// a and b must be in ascending order, but input abscissas need not be.
		a = ( ax < cx ? ax : cx );
		b = ( ax > cx ? ax : cx );

		// Initializations...
		x = w = v = bx;
		fw = fv = fx = func( x );

		// Main program loop.
		for ( Int iter = 0; iter < ITMAX; iter++ ) {
			xm = 0.5 * ( a + b );
			tol2 = 2.0 * ( tol1 = tol * std::abs( x ) + ZEPS );

			// Test for done here.
			if ( std::abs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) ) {
				fmin = fx;
				return xmin = x;
			}

			// Construct a trial parabolic fit.
			if ( std::abs( e ) > tol1 ) {
				r = ( x - w ) * ( fx - fv );
				q = ( x - v ) * ( fx - fw );
				p = ( x - v ) * q - ( x - w ) * r;
				q = 2.0 * ( q - r );
				if ( q > 0.0 ) p = -p;
				q = std::abs( q );
				etemp = e;
				e = d;

				// These conditions determine the acceptability of the parabolic
				// fit. Here we take the golden section step into the larger of
				// the two segments.
				if ( std::abs( p ) >= std::abs( 0.5 * q * etemp ) || p <= q * ( a - x )
						|| p >= q * ( b - x ) )
					d = CGOLD * ( e = ( x >= xm ? a - x : b - x ) );
				else {
					// Take the parabolic step.
					d = p / q;
					u = x + d;
					if ( u - a < tol2 || b - u < tol2 )
						d = SIGN( tol1, xm - x );
				}
			} else {
				d = CGOLD * ( e = ( x >= xm ? a - x : b - x ) );
			}

			u = ( std::abs( d ) >= tol1 ? x + d : x + SIGN( tol1, d ) );

			// This is the one function evaluation per iteration.
			fu = func( u );

			// Now decide what to do with our function evaluation.
			// Housekeeping follows:
			if ( fu <= fx ) {
				if ( u >= x )
					a = x;
				else
					b = x;
				shft3( v, w, x, u );
				shft3( fv, fw, fx, fu );
			} else {
				if ( u < x )
					a = u;
				else
					b = u;
				if ( fu <= fw || w == x ) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				} else if ( fu <= fv || v == x || v == w ) {
					v = u;
					fv = fu;
				}
			}
			// Done with housekeeping. Back for another iteration.
		}
		throw( "Too many iterations in brent" );
	}
};

// Brent’s method to find a minimum, modified to use derivatives.
struct Dbrent : Bracketmethod {
	Doub xmin, fmin;
	const Doub tol;

	Dbrent( const Doub toll = 3.0e-8 )
			: tol( toll ) {
	}

	// Given a functor funcd that computes a function and also its derivative
	// function df, and given a bracketing triplet of abscissas ax, bx, cx
	// [such that bx is between ax and cx, and f(bx) is less than both f(ax) and
	// f(cx)], this routine isolates the minimum to a fractional precision of
	// about tol using a modification of Brent’s method that uses derivatives.
	// The abscissa of the minimum is returned as xmin, and the minimum function
	// value is returned as min, the returned function value.
	template< class T >
	Doub minimize( T &funcd ) {
		const Int ITMAX = 100;
		const Doub ZEPS = std::numeric_limits< Doub >::epsilon() * 1.0e-3;

		// Will be used as flags for whether proposed steps are acceptable or not.
		Bool ok1, ok2;
		Doub a, b, d = 0.0, d1, d2, du, dv, dw, dx, e = 0.0;
		Doub fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

		// Comments following will point out only differences from the routine
		// in Brent. Read that routine first.
		a = ( ax < cx ? ax : cx );
		b = ( ax > cx ? ax : cx );
		x = w = v = bx;
		fw = fv = fx = funcd( x );
		dw = dv = dx = funcd.df( x );
		// All our housekeeping chores are doubled by the necessity of moving
		// around derivative values as well as function values.

		for ( Int iter = 0; iter < ITMAX; iter++ ) {
			xm = 0.5 * ( a + b );
			tol1 = tol * std::abs( x ) + ZEPS;
			tol2 = 2.0 * tol1;
			if ( std::abs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) ) {
				fmin = fx;
				return xmin = x;
			}

			if ( std::abs( e ) > tol1 ) {
				d1 = 2.0 * ( b - a );
				d2 = d1;
				if ( dw != dx ) d1 = ( w - x ) * dx / ( dx - dw ); // Secant method with one point.
				if ( dv != dx ) d2 = ( v - x ) * dx / ( dx - dv ); // And the other.
				// Which of these two estimates of d shall we take? We will
				// insist that they be within the bracket, and on the side
				// pointed to by the derivative at x:
				u1 = x + d1;
				u2 = x + d2;
				ok1 = ( a - u1 ) * ( u1 - b ) > 0.0 && dx * d1 <= 0.0;
				ok2 = ( a - u2 ) * ( u2 - b ) > 0.0 && dx * d2 <= 0.0;
				olde = e;
				e = d;
				if ( ok1 || ok2 ) {
					if ( ok1 && ok2 )
						d = ( std::abs( d1 ) < std::abs( d2 ) ? d1 : d2 );
					else if ( ok1 )
						d = d1;
					else
						d = d2;
					if ( std::abs( d ) <= std::abs( 0.5 * olde ) ) {
						u = x + d;
						if ( u - a < tol2 || b - u < tol2 )
							d = SIGN( tol1, xm - x );
					} else { //Bisect, not golden section.
						d = 0.5 * ( e = ( dx >= 0.0 ? a - x : b - x ) );
						// Decide which segment by the sign of the derivative.
					}
				} else {
					d = 0.5 * ( e = ( dx >= 0.0 ? a - x : b - x ) );
				}
			} else {
				d = 0.5 * ( e = ( dx >= 0.0 ? a - x : b - x ) );
			}

			if ( std::abs( d ) >= tol1 ) {
				u = x + d;
				fu = funcd( u );
			} else {
				u = x + SIGN( tol1, d );
				fu = funcd( u );
				if ( fu > fx ) {
					fmin = fx;
					return xmin = x;
					// If the minimum step in the downhill direction takes us
					// uphill, then we are done.
				}
			}

			du = funcd.df( u );

			// Now all the housekeeping, sigh.
			if ( fu <= fx ) {
				if ( u >= x )
					a = x;
				else
					b = x;
				mov3( v, fv, dv, w, fw, dw );
				mov3( w, fw, dw, x, fx, dx );
				mov3( x, fx, dx, u, fu, du );
			} else {
				if ( u < x )
					a = u;
				else
					b = u;
				if ( fu <= fw || w == x ) {
					mov3( v, fv, dv, w, fw, dw );
					mov3( w, fw, dw, u, fu, du );
				} else if ( fu < fv || v == x || v == w ) {
					mov3( v, fv, dv, u, fu, du );
				}
			}
		}
		throw( "Too many iterations in routine dbrent" );
	}
};



//Must accompany linmin in Dlinemethod.
template< class T >
struct Df1dim {

	const VecDoub &p;
	const VecDoub &xi;
	Int n;
	T &funcd;
	VecDoub xt;
	VecDoub dft;

	// Constructor takes as inputs an n-dimensional point p[0..n-1] and an
	// n-dimensional direc- tion xi[0..n-1] from linmin, as well as the functor
	// funcd.
	Df1dim( VecDoub_I &pp, VecDoub_I &xii, T &funcdd )
			: p( pp )
			, xi( xii )
			, n( pp.size() )
			, funcd( funcdd )
			, xt( n )
			, dft( n ) {
	}

	// Functor returning value of the given function along a one-dimensional line.
	Doub operator()( const Doub x ) {
		for ( Int j = 0; j < n; j++ )
			xt[j] = p[j] + x * xi[j];
		return funcd( xt );
	}

	// Returns the derivative along the line.
	Doub df( const Doub x ) {
		Doub df1 = 0.0;
		funcd.df( xt, dft );
		for ( Int j = 0; j < n; j++ )
			df1 += dft[j] * xi[j];
		return df1;
		// Dbrent always evaluates the derivative at the same value as the function,
		// so xt is un- changed.
	}



//	void clip( Doub & ax, Doub & bx, Doub & cx ) {
//		// order the points
//		if ( cx < bx ) SWAP( bx, cx );
//		if ( bx < ax ) SWAP( ax, bx );
//		if ( cx < bx ) SWAP( bx, cx );
//		assert( ax <= bx && bx <= cx );
//
//
//
//		VecDoub p( n ), q( n );
//		for ( Int j = 0; j < n; ++j ) {
//			p[j] =
//
//		}
//	}
};


// Base class for line-minimization algorithms using derivative information.
// Provides the line- minimization routine linmin.
template< class T >
struct Dlinemethod {

	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;

// Constructor argument is the user-supplied function or functor to be minimized
	Dlinemethod( T &funcc )
			: func( funcc )
			, n( -1 ) {
	}

	// Line-minimization routine. Given an n-dimensional point p[0..n-1] and an
	// n-dimensional direction xi[0..n-1], moves and resets p to where the
	// function or functor func(p) takes on a minimum along the direction xi
	// from p, and replaces xi by the actual vector displacement that p was
	// moved. Also returns the value of func at the returned location p. All of
	// this is actually accomplished by calling the routines bracket and
	// minimize of Dbrent.
	Doub linmin()
	{
		Doub ax, xx, xmin;
		n = p.size();
		Df1dim< T > df1dim( p, xi, func );
		ax = 0.0;
		xx = 1.0;
		Dbrent dbrent;
		dbrent.bracket( ax, xx, df1dim );
//		std::cout << "dbrent.bracket finished: " << ax << " / " << xx << std::endl;
//		std::cout << "ax : bx : cx -- " << dbrent.ax << " : " << dbrent.bx <<
//				" : " << dbrent.cx << std::endl;

//		df1dim.clip(  )

		xmin = dbrent.minimize( df1dim );
		for ( Int j = 0; j < n; j++ ) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return dbrent.fmin;
	}
};


// Base class for line-minimization algorithms. Provides the line-minimization
// routine linmin.
template< class T >
struct Linemethod {
	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;

	// Constructor argument is the user-supplied function or functor to be
	// minimized.
	Linemethod( T & funcc )
			: func( funcc ) {
	}

	// Line-minimization routine. Given an n-dimensional point p[0..n-1] and an
	// n-dimensional direction xi[0..n-1], moves and resets p to where the
	// function or functor func(p) takes on a minimum along the direction xi
	// from p, and replaces xi by the actual vector displacement that p was
	// moved. Also returns the value of func at the returned location p. This is
	// actually all accomplished by calling the routines bracket and minimize of
	// Brent.
	Doub linmin() {
		Doub ax, xx, xmin;
		n = p.size();
		F1dim< T > f1dim( p, xi, func );

		// Initial guess for brackets.
		ax = 0.0;
		xx = 1.0;
		Brent brent;
		brent.bracket( ax, xx, f1dim );
		xmin = brent.minimize( f1dim );

		// Construct the vector results to return
		for ( Int j = 0; j < n; j++ ) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return brent.fmin;
	}
};

template< class T >
struct F1dim {
	// Must accompany linmin in Linemethod.
	const VecDoub &p;
	const VecDoub &xi;
	Int n;
	T &func;
	VecDoub xt;

	// Constructor takes as inputs an n-dimensional point p[0..n-1] and an
	// n-dimensional di- rection xi[0..n-1] from linmin, as well as the function
	// or functor that takes a vector argument.
	F1dim( VecDoub_I &pp, VecDoub_I &xii, T &funcc )
			: p( pp ), xi( xii ), n( pp.size() ), func( funcc ), xt( n ) {
	}

	// Functor returning value of the given function along a one-dimensional line.
	Doub operator()( const Doub x ) {
		for ( Int j = 0; j < n; j++ )
			xt[j] = p[j] + x * xi[j];
		return func( xt );
	}
};


template< class T >
struct Frprmn : Dlinemethod< T > {

	// number of iterations taken to converge
	Int iter;

	// final function value
	Doub fret;

	// convergence tolerance
	const Doub ftol;

	// maximum number of iterations to take
	const Int maxiters;

	using Dlinemethod< T >::func;
	using Dlinemethod< T >::linmin;
	using Dlinemethod< T >::p;
	using Dlinemethod< T >::xi;

	//Constructor arguments are funcd, the function or functor to be minimized, and
	// an optional argument ftoll, the fractional tolerance in the function value
	// such that failure to decrease by more than this amount on one iteration
	// signals doneness.
	Frprmn( T & funcd, const Int maxiterss = 300, const Doub ftoll = 3.0e-8 )
			: Dlinemethod< T >( funcd )
			, iter( 0 )
			, fret( std::numeric_limits< Doub >::max() )
			, ftol( ftoll )
			, maxiters( maxiterss ) {
	}

	//Given a starting point pp[0..n-1], performs the minimization on a function
	// whose value and gradient are provided by a functor funcd (see text).
	VecDoub minimize( VecDoub_I & pp ) {
		const Int ITMAX = maxiters;
		const Doub EPS = 1.0e-18;
		const Doub GTOL = 1.0e-8;
		//Here ITMAX is the maximum allowed number of iterations; EPS is a small
		// number to rectify the special case of converging to exactly zero function
		// value; and GTOL is the convergence criterion for the zero gradient test.

		Doub gg, dgg;

		// Initializations.
		Int n = pp.size();
		p = pp;
		VecDoub g( n ), h( n );
		xi.resize( n );
		Doub fp = func( p );
		func.df( p, xi );

		for ( Int j = 0; j < n; j++ ) {
			g[j] = -xi[j];
			xi[j] = h[j] = g[j];
//			std::cout << "g[" << j << "] = " << g[j] << std::endl;
		}

		//Loop over iterations.
		for ( Int its = 0; its < ITMAX; its++ ) {
			iter = its;
			fret = linmin();
			//Next statement is one possible return:
			if ( 2.0 * std::abs( fret - fp ) <=
					ftol * ( std::abs( fret ) + std::abs( fp ) + EPS ) ) {
//				std::cout << "RETURNING CGD FTOL " << fret - fp << " (" << its << ")" << std::endl;
				return p;
			}
			fp = fret;
			func.df( p, xi );
			Doub test = 0.0;
			Doub den = MAX( std::abs( fp ), 1.0 );
			for ( Int j = 0; j < n; j++ ) {
				//Test for convergence on zero gradient.
				Doub temp = std::abs( xi[j] ) * MAX( std::abs( p[j] ), 1.0 ) / den;
				if ( temp > test ) test = temp;
			}

			if ( test < GTOL ) {
//				std::cout << "RETURNING CGD GTOL: " << test << std::endl;
				return p;
			}

			dgg = gg = 0.0;
			for ( Int j = 0; j < n; j++ ) {
				gg += g[j] * g[j];
				// dgg += xi[j] * xi[j];  // This statement for Fletcher-Reeves.
				dgg += ( xi[j] + g[j] ) * xi[j]; // This statement for Polak-Ribiere.
			}

			// Unlikely. If gradient is exactly zero, then we are already done.
			if ( gg == 0.0 ) {
//				std::cout << "RETURNING CGD GG==0.0: " << gg << std::endl;
				return p;
			}

			Doub gam = dgg / gg;
			for ( Int j = 0; j < n; j++ ) {
				g[j] = -xi[j];
				xi[j] = h[j] = g[j] + gam * h[j];
			}

//			std::cout << "gam: " << gam << std::endl;
		}

		throw( "Too many iterations in frprmn" );
	}
};


} // namespace nrc
} // namespace rdis
#endif // RDIS_MINIMIZE_NRC_H_
