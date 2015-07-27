/*
 * Optimizer.h
 *
 *  Created on: Feb 17, 2013
 *      Author: afriesen
 */

#ifndef RDIS_OPTIMIZER_H_
#define RDIS_OPTIMIZER_H_

#include "common.h"
#include "State.h"
#include "OptimizableFunction.h"

namespace rdis {

// utility function to print a vector
template< class Vector >
void printVector( const Vector & v ) {
    for ( size_t i = 0; i < v.size(); ++i ) {
        std::cout << ( i == 0 ? "" : ", " ) << v[i];
    }
}

// abstract base class for an Optimizer
class Optimizer {
public:
	Optimizer( OptimizableFunction & of );

	virtual ~Optimizer();

public:
	// optimize the function f(x) this Optimizer was created with
	// xinitial: initial state to start optimization at
	// finit: known evaluation of f(x) to use for BnB
	// timeout: time at which to abort the optimization
	// printInfo: print some output during & after the optimization
	// actualStartTime: used to keep printed timings consistent across
	// 					repeated calls to optimize())
	virtual Numeric optimize( const State * xinitial = NULL,
			Numeric finit = std::numeric_limits< Numeric >::max(),
			Clock::time_point timeout = Clock::time_point::max(),
			bool printInfo = true,
			Clock::time_point actualStartTime = Clock::time_point::max() ) = 0;

protected:
	virtual bool checkTimedOut();

public:

	// return true if the last optimization timed out, stopped, or was a bound
	virtual bool timedOut() { return m_timedOut; }
	virtual bool wasStopped() { return m_stopped; }
	virtual bool wasBound() { return m_wasBound; }

	// called to cleanly stop optimization (e.g., by a signal handler)
	virtual void stop( int sig ) {
		std::cout << "Optimizer received stop signal (" << sig << ")" << std::endl;
		m_stopped = true;
	}

	// return the value / state / error / time-taken from the last call to optimize()
	virtual Numeric getOptVal() const { return m_optVal; }
	virtual const State & getOptState() const { return m_optState; }
	virtual Numeric getOptError() const { return m_optError; }
	virtual Duration getOptTime() const { return m_optTime; }

    // print the statistics from the last call to optimize()
	virtual void printStats( std::ostream & os = std::cout, bool printInfo = true ) const = 0;
	friend std::ostream & operator<<( std::ostream & out, const Optimizer & p );

protected:
	// the function being optimized
	OptimizableFunction & m_func;

	// the optimal value, state, and error from the optimization
	Numeric m_optVal;
	State m_optState;
	Numeric m_optError;

	// the time taken to perform the most recent optimization
	Duration m_optTime;

	// set to true when integration times out (or is stopped)
	bool m_timedOut;
	bool m_stopped;
    bool m_wasBound;

	// how often to check the timer (too frequently is quite slow)
	int m_timerCheckFrequency;
	int m_timerCheckCount;

	// the start time and the timeout time of this optimization
	Clock::time_point m_starttime;
	Clock::time_point m_timeout;
};

} // namespace rdis

#endif // RDIS_OPTIMIZER_H_
