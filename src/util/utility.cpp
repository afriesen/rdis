/*
 * utility.cpp
 *
 *  Created on: May 8, 2014
 *      Author: afriesen
 */

#include "util/utility.h"

namespace rdis {

bool parse_command_line( int argc, char * argv[],
		const BOOSTNS::program_options::options_description & desc,
		BOOSTNS::program_options::variables_map & options ) {

	namespace po = BOOSTNS::program_options;

	std::cout << "command line: " << std::endl << "\t";
	for ( int i = 0; i < argc; ++i ) {
		std::cout << argv[i] << " ";
	}
	std::cout << std::endl;

	po::store( po::parse_command_line( argc, argv, desc ), options );

	try {
		po::notify( options );
	} catch ( po::error & e) {
		  std::cerr << "\nCommand line error: " << e.what() << std::endl << std::endl;
		  std::cerr << desc << std::endl;
		  return false;
	}

	if ( options.count( "help" ) ) {
		std::cout << desc << std::endl;
		return false;
	}

	return true;
}

} // namespace rdis
