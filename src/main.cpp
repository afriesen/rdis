/*
 * main.cpp
 *
 *  Created on: Jan 10, 2013
 *      Author: afriesen
 */

#include "common.h"

#include "OptimizableFunctionGenerator.h"
#include "PolynomialFunction.h"
#include "optimizers/CGDSubspaceOptimizer.h"
#include "RDISOptimizer.h"

#include "util/utility.h"
#include "util/signal_handler.h"

namespace po = BOOSTNS::program_options;


using namespace rdis;

void run_debug_test( po::variables_map & options );
void optimize_poly( const string & fname, po::variables_map & options );
void addOptions( po::options_description & desc );


const Numeric g_epsilon( 0.05 );
typedef std::vector< NumericInterval > NumericIntervalVec;



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char *argv[] )
{
	const Clock::time_point start = Clock::now();

	std::cout << "called with: " << std::endl << "\t";
	for ( int i = 0; i < argc; ++i ) {
		std::cout << argv[i] << " ";
	}
	std::cout << std::endl;

	po::options_description desc( "testRDIS options" );

	try {
		addOptions( desc );

		po::variables_map options;
		po::store( po::parse_command_line( argc, argv, desc ), options );

		po::notify( options );

		if ( options.count( "help" ) ) {
			std::cout << desc << std::endl;

		} else {
			std::vector<string> files;

			if ( options.count( "file" ) ) {
				// get the set of files to optimize
				const std::vector<string> tmpf =
						options["file"].as< std::vector<string> >();
				files.assign( tmpf.begin(), tmpf.end() );
			}

			if ( files.empty() ) {
				// run the pre-defined debug tests
				run_debug_test( options );

			} else {
				// optimize each input polynomial
				for ( string file : files ) {
					optimize_poly( file, options );
				}
			}
		}

	} catch ( const std::exception & e ) {
		std::cerr << "Exception occurred in main: \n\t" << e.what() << std::endl;
	} catch ( const std::string & s ) {
		std::cerr << "Exception occurred in main: \n\t" << s << std::endl;
	} catch ( const char * c ) {
		std::cerr << "Exception occurred in main: \n\t" << c << std::endl;
	}

	const Duration elapsed = Clock::now() - start;
	std::cout << "total elapsed time: " << elapsed.count() << " seconds" << std::endl;

	return 0;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void run_debug_test( po::variables_map & options ) {
	std::cout << "RUN DEBUG TEST" << std::endl;
	State x;
	Numeric val( 0 );
	NumericInterval bounds;

	std::vector< OptimizableFunctionSP > functions;
	NumericVec expectedResults;

	functions.push_back( OptimizableFunctionGenerator::makeSimplePoly() );
	expectedResults.push_back( -12 );

	functions.push_back( OptimizableFunctionGenerator::makeNonDecompPoly() );
	expectedResults.push_back( -12 );

	functions.push_back( OptimizableFunctionGenerator::makeSetToConstFactor() );
	expectedResults.push_back( -7 );

	functions.push_back( OptimizableFunctionGenerator::makeSetToConstFactor2() );
	expectedResults.push_back( 5.01399 );

	functions.push_back( OptimizableFunctionGenerator::make2ComponentPoly() );
	expectedResults.push_back( 8 );

	functions.push_back( OptimizableFunctionGenerator::makeTreePoly() );
	expectedResults.push_back( 8 );

	functions.push_back( OptimizableFunctionGenerator::makeMinStateWrongPoly() );
	expectedResults.push_back( -3.0625 );

	functions.push_back( OptimizableFunctionGenerator::makeCrossPoly() );
	expectedResults.push_back( 16 );

	functions.push_back( OptimizableFunctionGenerator::makePowellsFunction() );
	expectedResults.push_back( 0 );

	functions.push_back( OptimizableFunctionGenerator::makeTreePoly3() );
	expectedResults.push_back( 0 );

	BOOSTNS::program_options::options_description desc;
	RDISOptimizer::addOptions( desc );


	bool failed = false;

	for ( size_t i = 0; i < functions.size(); ++i ) {
		OptimizableFunction & testPoly( *functions[i] );
		std::cout << testPoly << std::endl;

		CGDSubspaceOptimizer ssopt( testPoly );
		RDISOptimizer opt( testPoly, ssopt, options);
		opt.setEpsilon( g_epsilon );
		opt.setApproxFactors( true );
		opt.setUseCache( false );

		bounds = testPoly.computeBounds();
		val = opt.optimize();

		x = opt.getOptState();
		std::cout << "[" << i << "] result: " << val <<  " at " << x
				<< " -- bounds: [" << bounds.lower() << ", " << bounds.upper()
				<< "]" << std::endl;
		opt.printStats( std::cout, true );

		if ( !approxeq( expectedResults[i], val, 1e-5 ) ) {
			failed = true;
			std::cout << "RESULT FOR FUNCTION " << i << " WAS INCORRECT: expected "
				<< expectedResults[i] << ", actual: " << val << ", diff: "
				<< expectedResults[i] - val << std::endl;

		}
	}

	if ( failed ) {
		std::cout << "AT LEAST ONE OPTIMIZATION RESULT WAS INCORRECT" << std::endl;
	} else {
		std::cout << "DEBUG TESTS COMPLETED SUCCESSFULLY" << std::endl;
	}
}



void optimize_poly( const string & fname, po::variables_map & options ) {
	PolynomialFunction fpoly( fname );

	std::cout << "polynomial to optimize (subspace opt. is CGD): " << std::endl;
	std::cout << fpoly << std::endl;

	CGDSubspaceOptimizer ssopt( fpoly );
	RDISOptimizer opt( fpoly, ssopt, options );

	NumericInterval bounds = fpoly.computeBounds();
	Numeric val = opt.optimize();

	State x = opt.getOptState();
	std::cout << "optimization result: " << val <<  " at " << x
		<< " -- bounds: [" << bounds.lower() << ", " << bounds.upper()
		<< "]" << std::endl;
	opt.printStats( std::cout, true );
}


void addOptions( po::options_description & desc ) {

	desc.add_options()
			( "file,f", po::value< std::vector< string > >(),
				  "input file(s) containing polynomial problem(s) to optimize (can be multiple(" )
			;

	RDISOptimizer::addOptions( desc );
}

