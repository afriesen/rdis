/*
 * OptimizableFunction.cpp
 *
 *  Created on: Jan 11, 2013
 *      Author: afriesen
 */

#include "common.h"
#include "OptimizableFunction.h"

#include BOOSTPATH/make_shared.hpp>
#include BOOSTPATH/lexical_cast.hpp>
#include BOOSTPATH/algorithm/string.hpp>

#include <string>

namespace fs = BOOSTNS::filesystem;


namespace rdis {

OptimizableFunction::OptimizableFunction( bool hasDerivatives, bool hasBounds,
		bool resetFactorsOnChange, bool hasBlockedVars )
	: numVars( 0 )
	, defaultDomain( "0:0" )
	, m_hasDerivatives( hasDerivatives )
	, m_hasBounds( hasBounds )
	, m_resetFactorsOnChange( resetFactorsOnChange )
	, m_hasBlockedVars( hasBlockedVars )
{}

OptimizableFunction::OptimizableFunction( const VariableDomain & defaultDomain,
		bool hasDerivatives, bool hasBounds, bool resetFactorsOnChange,
		bool hasBlockedVars )
	: numVars( 0 )
	, defaultDomain( defaultDomain )
	, m_hasDerivatives( hasDerivatives )
	, m_hasBounds( hasBounds )
	, m_resetFactorsOnChange( resetFactorsOnChange )
	, m_hasBlockedVars( hasBlockedVars )
{}

OptimizableFunction::~OptimizableFunction() {
	for ( Variable * v : variables ) {
		delete v;
	}

	for ( Factor * f : factors ) {
		delete f;
	}

	variables.clear();
	factors.clear();
}


Variable * OptimizableFunction::addVariable( string const & name,
		string const & domain, VariableID & varID ) {

	// if this variable was already added then just return it
	for ( Variable *v : variables ) {
		if ( v->getName() == name ) return v;
	}

	Variable *v = domain.empty() ?
				  new Variable( varID, name, defaultDomain ) :
				  new Variable( varID, name, domain );

	assert((VariableID) variables.size() == varID );
	variables.push_back( v );

	++varID;
	++numVars;

	return v;
}


void OptimizableFunction::init() {
	// initialize each factor
	for ( Factor * f : factors ) {
		f->init();
	}

	assert( numVars == (VariableCount) variables.size() );
}


Numeric OptimizableFunction::eval() const {
	Numeric ferr = 0.0;
	return evalFactors( factors, ferr, true );
}


Numeric OptimizableFunction::evalFactors( const FactorPtrVec & fctrs,
		Numeric & ferr, bool useCached ) const {

	Numeric feval = semiring().ProductIdentity(), newfeval = 0.0;
	ferr = semiring().ProductIdentity();

//	std::cout << "current state: " << std::endl;
//	for ( const Variable * v : variables ) {
//		std::cout << ( v->isAssigned() ?
//				BOOSTNS::lexical_cast< string >( v->eval() ) : "NA" ) << ", ";
//	}
//	std::cout << std::endl;

	for ( FactorPtrVec::const_iterator it( fctrs.begin() ), end( fctrs.end() );
			it != end; ++it ) {

		const Factor & f( **it );
		if ( !f.isAssigned() && !f.areAllVarsAssigned() ) continue;

		// eval the factor and keep track of the product and the total error

		newfeval = ( useCached ? f.eval() : f.evalNoCache() );
//		std::cout << "f" << f.getID() << " eval: " << newfeval << std::endl;
		feval = semiring().Product( feval, newfeval );

		ferr = semiring().Product( ferr, f.getSimplificationError() );

//		std::cout << f.getID() << " (" << f.eval() << "), ";

//		// eval the derivatives if we're using them
//		if ( useDerivatives ) {
//			tempPG.clear();
//			f.computeGradient( tempPG );
//			productGradient( vd.grad, tempPG );
//		}

		++instrumentation::g_numProducts;
	}

	return feval;
}


Numeric OptimizableFunction::eval( const State & x ) {

	State xprev( variables.size() );
	BOOSTNS::dynamic_bitset<> wasassigned( variables.size(), false );

	// get the current values, if they exist
	for ( Variable * v : variables ) {
		if ( v->isAssigned() ) {
			xprev[v->getID()] = v->eval();
			wasassigned.set( v->getID() );
		}
	}

	// set the new values
	for ( Variable * v : variables ) {
		v->assign( x[v->getID()] );
		onVarAssigned( v->getID(), x[v->getID()] );
	}

//	// eval, making sure not to cache anything
//	Numeric fv = semiring().ProductIdentity();
//	for ( Factor * f : factors ) {
//		fv = semiring().Product( fv, f->evalNoCache( ) );
//	}

	Numeric ferr = 0.0;
	const Numeric fv = evalFactors( factors, ferr, false );

	// set the vars back
	for ( Variable * v : variables ) {
		if ( wasassigned.test( v->getID() ) ) {
			v->assign( xprev[v->getID()] );
			onVarAssigned( v->getID(), xprev[v->getID()] );
		} else {
			v->unassign();
			onVarUnassigned( v->getID() );
		}
	}

	return fv;
}


NumericInterval OptimizableFunction::computeBounds() const {
	return computeBounds( factors, -1 );
}


NumericInterval OptimizableFunction::computeBounds( const FactorPtrVec & fctrs,
		VariableID assignedVID ) const {

	if ( !hasBounds() ) return NumericInterval( NumericMIN, NumericMAX );

	assert( hasBounds() );
	NumericInterval bounds( semiringB().ProductIdentity() );

	for ( const Factor * f : fctrs ) {

//		if ( assignedVID >= 0 &&
//				( f->isAssigned() && f->getAssignedKey() != assignedVID ) ) {
//			continue;
//		}
//		if ( assignedVID >= 0 && !f->containsVar( assignedVID ) ) continue;

		if ( assignedVID < 0
				|| ( f->isAssigned() && f->getAssignedKey() == assignedVID )
				|| !f->isAssigned()
//				|| ( !f->isAssigned() && f->containsVar( assignedVID ) )
				) {

			NumericInterval fb = f->computeBounds();
			bounds = semiringB().Product( bounds, fb );

//			++instrumentation::g_counter2;
		}
	}

	return bounds;
}


void OptimizableFunction::computeGradient( NumericVec & gradient,
		bool checkGrad ) const {
	gradient.clear();
	gradient.resize( getNumVars(), semiring().ProductIdentity() );

	PartialGradient pg = pgtemp1;

	computeGradient( factors, pg, checkGrad );

	for ( const auto & dfdv : pg ) {
		gradient[dfdv.first] = dfdv.second;
	}
}


void OptimizableFunction::computeGradient( const FactorPtrVec & facs,
		PartialGradient & gradient, bool checkGrad ) const {

	if ( semiring().Product( 2, 6 ) == 8 ) {
		computeGradientOfSum( facs, gradient );
	} else {
		assert( semiring().Product( 2, 6 ) == 12 );
		computeGradientOfProduct( facs, gradient );
	}

	if ( checkGrad ) checkGradient( facs, gradient );
}


void OptimizableFunction::computeGradientOfSum( const FactorPtrVec & facs,
		PartialGradient & gradient ) const {
	gradient.clear();
	gradient.reserve( getNumVars() );

	assert( semiring().Product( 2, 6 ) == 8 );

	PartialGradient pg;
	for ( const Factor * f : facs ) {
		pg.clear();
		f->computeGradient( pg );

		productGradient( gradient, pg, semiring() );
	}
}


void OptimizableFunction::computeGradientOfProduct( const FactorPtrVec & facs,
		PartialGradient & gradient ) const {

	gradient.clear();
	gradient.reserve( getNumVars() );

	assert( semiring().Product( 2, 6 ) == 12 );

	Numeric ferr = 0.0;
	PartialGradient g;
	bool iszero = false;
	Numeric funcEval = evalFactors( facs, ferr, true );

//	std::cout << "funceval: " << funcEval << std::endl;

	for ( const Factor * f : facs ) {
		g.clear();
		f->computeGradient( g );

		const Numeric facEval = f->eval();

		iszero = ( facEval == 0 ); //approxeq( facEval, 0.0, 1e-16 );

		if ( !iszero ) funcEval /= facEval;

//		std::cout << "fac " << f->getID() << ": " << funcEval << " ( / " <<
//				facEval << ")" << std::endl;

		for ( const auto & vdf : g ) {
			auto it = gradient.find( vdf.first );
			if ( it == gradient.end() ) {
				gradient[vdf.first] = funcEval * vdf.second;
			} else {
				it->second += funcEval * vdf.second;
			}
		}

		if ( !iszero ) funcEval *= facEval;
	}
}


void OptimizableFunction::print( std::ostream & os ) const {
	os << "# vars: " << getNumVars() << ", # terms: "
			<< factors.size() << std::endl;
	os << "domain: " << defaultDomain;
	if ( m_hasBounds ) os << ", bounds: " << computeBounds();
	else os << " (bounds not supported on this function type)";
	os << std::endl;

	os << /*std::endl <<*/ "terms: " << std::endl;
	for ( const Factor * fp : factors ) {
		os << "\t" << *fp << "\t\t-- (bounds " << fp->computeBounds() <<
				")" << std::endl;
	}
 }


std::ostream & operator<<( std::ostream & out, OptimizableFunction const & p ) {
	p.print( out );
	return out;
}


bool operator==( const OptimizableFunction & lhs,
		const OptimizableFunction & rhs ) {
	if ( lhs.numVars != rhs.numVars ||
			lhs.factors.size() != rhs.factors.size() ) {
		return false;
	}

	for ( size_t fi = 0; fi < lhs.factors.size(); ++fi ) {
		if ( lhs.factors[fi] != rhs.factors[fi] ) return false;
		if ( *lhs.factors[fi] != *rhs.factors[fi] ) return false;
		if ( lhs.factors[fi]->getVariables() !=
				rhs.factors[fi]->getVariables() ) {
			return false;
		}
	}

	return true;
}


bool operator!=( const OptimizableFunction & lhs,
		const OptimizableFunction & rhs ) {
	return !( lhs == rhs );
}


bool OptimizableFunction::readAndSplit( fs::fstream & ifs, SplitVector & sv ) {

	using namespace BOOSTNS;
	string line;

	sv.clear();
	while ( !ifs.eof() && !ifs.bad() ) {
		std::getline( ifs, line );
		if ( line.empty() || line[0] == '#' ) continue;

//		split( sv, line, is_space(), token_compress_on );
		split( sv, line, is_any_of( " \t,;" ), token_compress_on );

		return true;
	}

	return false;
}


void OptimizableFunction::checkGradient( const FactorPtrVec & facs,
		const PartialGradient & g, const Numeric h ) const {

//	const Numeric h = 1e-8;
	Numeric origeval = 0.0, ferr = 0.0;
	PartialGradient fd;
	fd.reserve( g.size() );

	Variable * v = NULL;

//	for ( Variable * v : variables ) {
	for ( const auto & pgv : g ) {
		v = variables[pgv.first];

		origeval = v->eval();

		v->assign( v->eval() - h / 2.0 );
		onVarAssigned( v->getID(), v->eval() );
		Numeric fl = evalFactors( facs, ferr, false );

		v->assign( v->eval() + h );
		onVarAssigned( v->getID(), v->eval() );
		Numeric fh = evalFactors( facs, ferr, false );

		// restore v
		v->assign( origeval );
		onVarAssigned( v->getID(), v->eval() );

		fd[v->getID()] = ( fh - fl ) / h;
//		std::cout << "fd for " << v->getID() << ": " << fh << " - " << fl <<
//				" / " << h << " = " << fd[v->getID()] << std::endl;
	}

	assert( fd.size() == g.size() );

	bool diff = false;

//	for ( Variable * v : variables ) {
	for ( const auto & pgv : g ) {
		v = variables[pgv.first];

		assert( g.find( v->getID() ) != g.end() );
		const Numeric fdv = fd.at( v->getID() );

		Numeric fdg = fdv - g.at( v->getID() );
		fdg = ( fdv == 0 ? fdg : fdg / fdv );

		if ( !approxeq( fdg, 0.0, 1e-3 ) ) {
			diff = true;

			std::cout << "func grad -- var " << v->getID() << ": df/dv = " <<
					g.at( v->getID() ) << ", fd = " << fd.at( v->getID() ) <<
					" -- diff = " << fd.at( v->getID() ) - g.at( v->getID() ) <<
					", fdg = " << fdg << std::endl;
		}
	}

//	assert( !diff );
}


} // namespace rdis
