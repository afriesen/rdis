/*
 * CGDSubspaceOptimizer.h
 *
 *  Created on: Mar 25, 2014
 *      Author: afriesen
 */

#ifndef RDIS_CGDSUBSPACEOPTIMIZER_H_
#define RDIS_CGDSUBSPACEOPTIMIZER_H_

#include "common.h"
#include "SubspaceOptimizer.h"
#include "OptimizableFunction.h"

#include "external/include/minimize_nrc.h"


namespace rdis {

// conjugate gradient descent subspace optimizer
class CGDSubspaceOptimizer
		: public SubspaceOptimizer {

	typedef SubspaceOptimizer super;

public:
	CGDSubspaceOptimizer( OptimizableFunction & f_ );
	virtual ~CGDSubspaceOptimizer();

public:
	virtual Numeric optimize( const VariablePtrVec & vars,
			const FactorPtrVec & factors, NumericVec & xinit,
			Numeric & deltaFval, const bool printdbg );

protected:
	// temporary vector and gradient storage to reduce memory overhead
	NumericVec nvtmp;
	IntervalVec ivtmp;
	PartialGradient pgtmp;
	PartialGradient pgtmp2;


protected:
	// structure to use NRC minimization methods on a sub-function of the main
	// optimizable function (i.e., a subset of the variables and a subset of the
	// factors)
	struct SubfunctionFD {
		SubfunctionFD( OptimizableFunction & func_, bool doAscent,
				const VariablePtrVec & vars_, const FactorPtrVec & facs_,
				IntervalVec & ivtemp,
				PartialGradient & pgtemp1_, PartialGradient & pgtemp2_ );

		// evaluate this subfunction at the specified location (must be in same
		// order as <vars> vector)
		nrc::Doub operator()( nrc::VecDoub_I & x );

		// compute the derivative (gradient) at the specified location and put
		// the result in <deriv>
		void df( nrc::VecDoub_I & x, nrc::VecDoub_O & deriv );


		bool quickAssignVals( const NumericVec & xval );


		// function to update on assign
		OptimizableFunction & func;

		// coefficient used to switch between ascent and descent
		const Numeric coeff;

		// set of vars we're optimizing over (all the vectors must stay
		// consistent with this ordering)
		const VariablePtrVec & vars;

		// set of factors that define the function we're evaluating
		const FactorPtrVec & facs;

		// the bounds of the vars we're optimizing over
		IntervalVec & bounds;

		PartialGradient & pgtemp1;
		PartialGradient & pgtemp2;
	};
};

} // namespace rdis

#endif // RDIS_CGDSUBSPACEOPTIMIZER_H_
