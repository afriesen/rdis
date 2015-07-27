/*
 * PolynomialFunction.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: afriesen
 */

#include "PolynomialFunction.h"

//#include "SimpleSumFactor.h"
#include "NonlinearProductFactor.h"

#include BOOSTPATH/tokenizer.hpp>
#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/algorithm/string.hpp>
#include BOOSTPATH/filesystem/fstream.hpp>

#include <iostream>
#include <algorithm>
#include <boost/lexical_cast.hpp>

namespace fs = BOOSTNS::filesystem;


namespace rdis {

PolynomialFunction::PolynomialFunction( const Semiring< Numeric > & sr,
		const Semiring< NumericInterval > & srB )
	: super( true, true, false, false )
	, m_semiringNum( sr.clone() )
	, m_semiringBound( srB.clone() ) {
}

PolynomialFunction::PolynomialFunction( const VariableDomain & defaultDomain,
		const Semiring< Numeric > & sr,
		const Semiring< NumericInterval > & srB )
	: super( defaultDomain, true, true, false, false )
	, m_semiringNum( sr.clone() )
	, m_semiringBound( srB.clone() ) {
}

PolynomialFunction::PolynomialFunction( const string & filename,
		const Semiring< Numeric > & sr,
		const Semiring< NumericInterval > & srB )
	: super( true, true, false, false )
	, m_semiringNum( sr.clone() )
	, m_semiringBound( srB.clone() ) {

	load( filename );
}

PolynomialFunction::~PolynomialFunction() {
	delete m_semiringNum;
	delete m_semiringBound;
}



// load a polynomial from the specified file
bool PolynomialFunction::load( const string &file_in ) {
	// open the specified file
	fs::fstream ifs( file_in, std::ios_base::in );
	if ( !ifs.is_open() ) {
		std::cerr << "PolynomialFunction::Load: Failed to open file for reading " <<
				file_in << std::endl;
		return false;
	}

	VariableID unique_var_id = 0;
	string line;

	try {
		std::getline( ifs, line );
		string name, domain;

		while ( !ifs.eof() && !ifs.bad() ) {
			// ignore empty or comment lines
			if ( line.empty() || line[0] == '#' ) {
				std::getline( ifs, line );
				continue;
			}

			size_t equals_pos = line.find( '=' );
			if ( equals_pos != string::npos ) {
				readVariable( line, equals_pos, unique_var_id );
			} else {
				readFactor( line, unique_var_id );
			}

			std::getline( ifs, line );
		}

	} catch ( std::exception & e ) {
		std::cerr << "PolynomialFunction::load(" << file_in << ") caught exception: "
				<< std::endl << e.what() << std::endl;
	}

	// close file
	if ( ifs.is_open() ) ifs.close();

	init();

	return ( variables.size() > 0 && factors.size() > 0 );
}


// read in a variable from the specified line of a file
bool PolynomialFunction::readVariable( string const & line, size_t equals_pos,
		VariableID & unique_var_id ) {
	string name = BOOSTNS::trim_copy( line.substr( 0, equals_pos ) );
	string domain = BOOSTNS::trim_copy( line.substr( equals_pos + 1 ) );

	string lcname = BOOSTNS::algorithm::to_lower_copy( name );
	try {
		if ( lcname == "default" ) {
			defaultDomain = VariableDomain( domain );
		} else {
			addVariable( name, domain, unique_var_id );
		}
	} catch ( BOOSTNS::bad_lexical_cast & e ) {
		std::cout << "ERROR occurred while trying to add variable with name " <<
				lcname << " and domain string " << domain << ": " << e.what() << std::endl;
		throw e;
	}
	return true;
}


//// read in a factor from the specified line of a file
//bool PolynomialFunction::readFactor( string const & line, VariableID & unique_var_id ) {
//	// create a new factor and put it in the factor list
//	Factor * f = new SimpleSumFactor( factors.size() );
//	factors.push_back( f );
//
//	string name;
//	Variable * v( NULL );
//
//	// tokenize to get the factor's variables
//	BOOSTNS::tokenizer<> tok( line );
//	for ( BOOSTNS::tokenizer<>::iterator it = tok.begin();
//			it != tok.end(); ++it ) {
//
//		name = *it;
//		v = NULL;
//
//		// add this variable to this polynomial
//		v = addVariable( name, "", unique_var_id );
//		assert( (VariableCount) variables.size() == unique_var_id );
//
//		// cross reference the variable and the factor
//		f->addVariable( v );
//	}
//
//	return true;
//}

// read in a factor from the specified line of a file
bool PolynomialFunction::readFactor( string const & line,
									 VariableID & unique_var_id ) {
	// create a new factor and put it in the factor list
	NonlinearProductFactor * f = new NonlinearProductFactor( factors.size() );
	factors.push_back( f );

	string name;
	Numeric exponent = 1.0;
	Variable * v( NULL );

	std::vector< string > sv1, sv2;
	BOOSTNS::split( sv1, line, BOOSTNS::is_any_of( "," ) );
	assert( !sv1.empty() );

	// Iterate through the split up parts of the string -- these should be of
	// the form [ x0^e ], where x0 is the variable name and e is the exponent.
	// These can also just be a single constant with no variables, which is
	// then simply interpreted as [ k ], where k is the constant.
	for ( const string & s : sv1 ) {
		BOOSTNS::split( sv2, s, BOOSTNS::is_any_of( "^" ) );
		assert( sv2.size() == 1 || sv2.size() == 2 );
		if ( sv2.empty() || sv2.size() > 2 ) {
			std::cout << "ERROR: Invalid term when reading PolynomialFunction at '"
				<< s << "' in line: '" << line << "'" << std::endl;
			throw "Invalid split when reading PolynomialFunction";
		}

		name = BOOSTNS::algorithm::to_lower_copy( BOOSTNS::trim_copy( sv2[0] ) );
		exponent = 1.0;

		Variable * v = NULL;

		if ( sv2.size() == 1 ) {
			try {
				Numeric constant = BOOSTNS::lexical_cast< Numeric >( name );
				// if this doesn't throw, then this is simply a constant
				f->setCoeff( constant );
			} catch ( BOOSTNS::bad_lexical_cast & ) {
				// if this threw, then it's not a number
				v = addVariable( name, "", unique_var_id );
			}
		} else {
			assert( sv2.size() == 2 );
			v = addVariable( name, "", unique_var_id );
			try {
				exponent = BOOSTNS::lexical_cast< Numeric >(
						BOOSTNS::trim_copy( sv2[1] ) );
			} catch ( BOOSTNS::bad_lexical_cast & e ) {
				std::cout << "ERROR: Invalid exponent when reading PolynomialFunction at '"
						<< s << "' in line: '" << line << "'" << std::endl;
				throw e;
			}
		}

		if ( v != NULL ) f->addVariable( v, exponent, 0, false );
	}

	return true;
}


// write this polynomial to a file
bool PolynomialFunction::save( const string & file_out, State * state ) {

	// TODO: deal with PolynomialFunction::save when state != NULL (write out
	// the state at the end)
	assert( state == NULL );

	fs::fstream ofs( file_out, std::ios_base::out );
	if ( !ofs.is_open() ) {
		std::cerr << "PolynomialFunction::Save: Failed to open file for writing (" <<
				file_out << ")" << std::endl;
		return false;
	}
	bool result = true;

	string line;

	try {
		// write the default variable domain as "default: <domain>\n"
		ofs << "default: " << defaultDomain << std::endl;

		for ( VariablePtrVec::value_type & v : variables ) {
			// write as "<varname>: <domain>\n"
			ofs << v->getName() << ": " << v->getDomain() << std::endl;
		}

		ofs << std::endl;

		for ( Factor * f : factors ) {
			const VariablePtrVec & vars = f->getVariables();
			// for each variable in this factor, write the <varname>
			for ( VariablePtrVec::const_iterator it( vars.begin() );
					it != vars.end(); ++it ) {
				ofs << ( it != vars.begin() ? "," : "" ) << (*it)->getName();
			}
			ofs << std::endl;
		}

	} catch ( std::exception & e ) {
		std::cerr << "PolynomialFunction::save(" << file_out << ") caught exception: "
				<< std::endl << e.what() << std::endl;
		result = false;
	}

	if ( ofs.is_open() ) ofs.close();
	return result;
}

} // namespace rdis
