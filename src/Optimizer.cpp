/*
 * Optimizer.cpp
 *
 *  Created on: Feb 17, 2013
 *      Author: afriesen
 */

#include "common.h"
#include "Optimizer.h"

namespace rdis {

Optimizer::Optimizer( OptimizableFunction & of )
        : m_func( of )
        , m_optVal( 0.0 )
        , m_optError( 0.0 )
        , m_timedOut( false )
        , m_stopped( false )
        , m_timerCheckFrequency( 100 )
        , m_timerCheckCount( 0 )
        , m_wasBound( false )
{
}

Optimizer::~Optimizer() {}


bool Optimizer::checkTimedOut() {
    if ( ++m_timerCheckCount >= m_timerCheckFrequency ) {
        m_timerCheckCount = 0;
        // set this flag like this so we don't have a race condition with the
        // stop() function call
        if ( Clock::now() > m_timeout ) m_timedOut = true;
    }
    return ( m_timedOut || m_stopped );
}


std::ostream & operator<<( std::ostream & out, const Optimizer & o ) {
	o.printStats( out, true );
	return out;
}

} // namespace rdis
