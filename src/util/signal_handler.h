/*
 * signal_handler.h
 *
 *  Created on: Aug 20, 2013
 *      Author: afriesen
 */

#ifndef RDIS_SIGNAL_HANDLER_H_
#define RDIS_SIGNAL_HANDLER_H_

#include "Optimizer.h"

namespace rdis {

void registerHandler( bool clear = false );

void setSignalHandler( Optimizer * opt );

void clearSignalHandler();

extern "C" void signalHandler( int signal );

} // namespace rdis

#endif // RDIS_SIGNAL_HANDLER_H_
