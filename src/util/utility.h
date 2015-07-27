/*
 * utility.h
 *
 *  Created on: Apr 15, 2013
 *      Author: afriesen
 */

#include "common.h"

#include BOOSTPATH/swap.hpp>
#include BOOSTPATH/random/mersenne_twister.hpp>
#include BOOSTPATH/random/uniform_int_distribution.hpp>

#include BOOSTPATH/program_options.hpp>

#ifndef RDIS_UTILITY_H_
#define RDIS_UTILITY_H_

namespace rdis {


bool parse_command_line( int argc, char * argv[],
		const BOOSTNS::program_options::options_description & desc,
		BOOSTNS::program_options::variables_map & options );


template < class Container >
void permute( Container & c, BOOSTNS::random::mt19937 & rng ) {

	// modern Fischer-Yates shuffle algorithm
	for ( int i = c.size() - 1; i > 0; --i ) {
		BOOSTNS::random::uniform_int_distribution<> unif(0, i);
		BOOSTNS::swap( c[i], c[unif( rng )] );
	}
}

} // namespace rdis

#endif // RDIS_UTILITY_H_
