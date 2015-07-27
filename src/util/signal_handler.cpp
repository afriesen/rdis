/*
 * signal_handler.cpp
 *
 *  Created on: Aug 20, 2013
 *      Author: afriesen
 */

#include "util/signal_handler.h"

#include <signal.h>


namespace rdis {

Optimizer * g_signalOptimizer = NULL;


void registerHandler( bool clear ) {
	struct sigaction action;
	if ( clear ) action.sa_handler = SIG_DFL;
	else action.sa_handler = signalHandler;
	sigemptyset( &action.sa_mask );
	action.sa_flags = 0;
	int result = sigaction( SIGINT, &action, NULL );
	int result2 = sigaction( SIGUSR1, &action, NULL );
	int result3 = sigaction( SIGTERM, &action, NULL );
	if ( result != 0 ) std::cout << "register handler returned code " << result << std::endl;
	if ( result2 != 0 ) std::cout << "register handler returned code " << result << std::endl;
	if ( result3 != 0 ) std::cout << "register handler returned code " << result << std::endl;
//	result = sigaction( SIGPIPE, &action, NULL );
//	sigaction( SIGTSTP, &action, NULL );
}

void setSignalHandler( Optimizer * opt ) {
	g_signalOptimizer = opt;
	registerHandler( false );
}

void clearSignalHandler() {
	g_signalOptimizer = NULL;
	registerHandler( true );
}

void signalHandler( int signal ) {
	std::cout << "stop signal handler called (null? " <<
			( ( g_signalOptimizer != NULL ) ? "no)" : "yes)" ) << std::endl;
	if ( g_signalOptimizer != NULL ) g_signalOptimizer->stop( signal );
	clearSignalHandler();
}

} // namespace rdis
