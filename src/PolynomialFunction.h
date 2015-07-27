/*
 * PolynomialFunction.h
 *
 *  Created on: Apr 10, 2013
 *      Author: afriesen
 */

#ifndef POLYNOMIALFUNCTION_H_
#define POLYNOMIALFUNCTION_H_

#include "common.h"
#include "OptimizableFunction.h"

namespace rdis {

// this function is simply a sum of NonlinearProductFactor terms, which, when
// read from a file, evaluates to
//  	f(x) = \sum_{g_j \in terms} [ c_j * \prod_{x_i \in vars(g_j)} x_i^e_i ]
// where c_j is the (constant) factor coefficient and e_i are the exponents of
// the individual variables
class PolynomialFunction
		: public rdis::OptimizableFunction {

	typedef OptimizableFunction super;

protected:
	const Semiring< Numeric > * m_semiringNum;
	const Semiring< NumericInterval > * m_semiringBound;
public:
	virtual const Semiring< Numeric > & semiring() const { return *m_semiringNum; }
	virtual const Semiring< NumericInterval > & semiringB() const { return *m_semiringBound; }

public:
	PolynomialFunction(
			const Semiring< Numeric > & sr = MinSumSemiring< Numeric >(),
			const Semiring< NumericInterval > & srB =
					MinSumSemiring< NumericInterval >() );

	explicit PolynomialFunction( const VariableDomain & defaultDomain,
			const Semiring< Numeric > & sr = MinSumSemiring< Numeric >(),
			const Semiring< NumericInterval > & srB =
					MinSumSemiring< NumericInterval >() );

	explicit PolynomialFunction( const string & filename,
			const Semiring< Numeric > & sr = MinSumSemiring< Numeric >(),
			const Semiring< NumericInterval > & srB =
					MinSumSemiring< NumericInterval >() );

	virtual ~PolynomialFunction();

public:
	// load/save this polynomial from/to a file
	virtual bool load( const string &file );
	virtual bool save( const string &file, State * state = NULL  );

private:
	bool readVariable( string const & line, size_t equals_pos, VariableID & unqvid );
	bool readFactor( string const & line, VariableID & unique_var_id );

};

} // namespace rdis

#endif // POLYNOMIALFUNCTION_H_
