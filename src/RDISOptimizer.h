/*
 * OptimizerDG.h
 *
 *  Created on: Jul 25, 2013
 *      Author: afriesen
 */

#ifndef RDIS_RDISOPTIMIZER_H_
#define RDIS_RDISOPTIMIZER_H_

#include "common.h"
#include "Optimizer.h"
#include "Component.h"
#include "SubspaceOptimizer.h"

#include BOOSTPATH/chrono.hpp>
#include BOOSTPATH/tuple/tuple.hpp>
#include BOOSTPATH/graph/adjacency_list.hpp>
#include BOOSTPATH/program_options.hpp>

#include BOOSTPATH/container/vector.hpp>
#include BOOSTPATH/function.hpp>
#include BOOSTPATH/random.hpp>

#include <list>
#include <queue>

namespace rdis {

typedef std::vector< VariableCount > VariableCountVec;
typedef BOOSTNS::uint32_t uint32;


class RDISOptimizer
		: public Optimizer {

typedef Optimizer super;

	enum VarSelHeur {
		All = 0,
		Random,
		HMetis,
		PaToH
	};

	enum ValSelType {
		Init = 0,
		AlternatingMin,
		Restart
	};

	struct RDISOptionValues {
		RDISOptionValues()
			: timeAsSeed( false )
			, useStaticDecomp( false )
			, desiredNumVars( 0 )
			, heur( PaToH )
			, minNumRestarts( 1 )
			, nRRperLvl( 2 )
			, nRRatTop( nRRperLvl )
			, maxNumAssignsUntilRR( 10 )
			, alwaysUseAllFactors( false )
			, noAssignLimitAtTop( true )
		{}

		// true to use the current time as a random seed
		bool timeAsSeed;

		// true to use a static decomposition (decompose once, re-use from then on)
		bool useStaticDecomp;

		// variable selection options
		uint32 desiredNumVars;
		VarSelHeur heur;

		// value selection options
		uint32 minNumRestarts;
		uint32 nRRperLvl;
		uint32 nRRatTop;
		uint32 maxNumAssignsUntilRR;
		bool alwaysUseAllFactors;
		bool noAssignLimitAtTop;
	};

public:
	RDISOptimizer( OptimizableFunction & of, SubspaceOptimizer & ssopt,
			BOOSTNS::program_options::variables_map options_ =
					BOOSTNS::program_options::variables_map(),
			size_t cacheSize = 100000 );
	virtual ~RDISOptimizer();

	static void addOptions(
			BOOSTNS::program_options::options_description & desc );

public:
	Numeric optimize( const State * xinitial = NULL,
			Numeric finit = NumericMAX,
			Clock::time_point timeout = Clock::time_point::max(),
			bool printInfo = true,
			Clock::time_point actualStartTime = Clock::time_point::max() );

protected:
	void doOptimization( Component & c );

	// pick a variable to assign from this component
	virtual void getVarsToAssign( Component & ccomp, VarDataPtrVec & vds );
	void ensureFactorWillBeAssigned( Component & ccomp, VarDataPtrVec & vars );
	void addBlockRemainders( Component & ccomp, VarDataPtrVec & vars,
							 const size_t initpos );
	virtual bool getVarsToAssignStatic( Component & ccomp, VarDataPtrVec & vds );
	virtual void getVarsToAssignMOMS( Component & ccomp, VarDataPtrVec & vds,
									  bool weighted );
	virtual void getVarsToAssignHMetis( Component & ccomp, VarDataPtrVec & vds );
	virtual void getVarsToAssignPaToH( Component & ccomp, VarDataPtrVec & vds );

	bool checkEmpty( Component & ccomp );
	bool checkUnassignedBound( Component & ccomp );
	bool checkAssignedBound( Component & ccomp );

	bool getFromCache( Component & ccomp );

	virtual bool getValueFromDomain( Component & ccomp,
			const VarDataPtrVec & vds, NumericVec & value );

	ValSelType getSSInitialVal( Component & ccomp,
			const VarDataPtrVec & vars, NumericVec & value,
			const bool doAlternatingMin, const bool forceRR,
			const bool wasAssigned );

	void quickAssignSSInitialVal( Component & ccomp, const VarDataPtrVec & vars,
			const NumericVec & value, const bool doAlternatingMin,
			const SubdomainCSP & sdprev, const bool useAllFactors );

	void quickUnassignSSInitialVal( Component & ccomp, const VarDataPtrVec & vars,
			const NumericVec & value, const bool doAlternatingMin,
			const SubdomainCSP & sdprev, const bool useAllFactors );

	Numeric sampleRandomState( const VarData & vd );

	void quickAssignVals( const VarDataPtrVec & vars, const NumericVec & xval );
	void quickUnassignVals( const VarDataPtrVec & vars );

	void assign( Component & ccomp, const VarDataPtrVec & vds,
			const NumericVec & xval );
	void unassign( Component & ccomp, const VarDataPtrVec & vds );
	void unassign( VariableID vid ); // TODO: remove this function? (brought from old Optimizer.cpp)

	// simplify/unsimplify factors and return true if any were (un)simplified
	bool simplifyFactors( Component & ccomp, const VarDataPtrVec & vds );
	bool unsimplifyFactors( Component & ccomp );

	// simplify/unsimplify factors at the top level and return the eval of the
	// simplified factors
	Numeric simplifyFactors( Numeric & ferr );
	void unsimplifyFactors();

	void productComponentVals( Component & parent, Component & child );

	virtual void updateDomain( Component & ccomp );
	bool computeCurrentFunctionEval( Component & ccompIn,
			const SubdomainCSP & newsd, Numeric & curFuncEval );

	void addToCache( Component & ccomp );
	void removeChildrenFromCache( Component & ccomp );

	// set the initial vals to assign to each var
	virtual void setInitialState( const State & xinit );
	void setInitialValFromChildren( const Subdomain & sd );

protected:
	// call to initialize the memory used for optimization
	virtual void init( size_t cacheSize );

	virtual void processOptions();

	// call before performing optimization
	virtual void prepareToOptimize();

	virtual bool checkTimedOut( Component & ccomp );

protected:
	void writeGraph() const;
	void writeFactorGraph( const Component & ccomp ) const;

	void printVars( const VarDataPtrVec & vds ) const;

public:
	virtual void printStats( std::ostream & os = std::cout,
			bool printInfo = true ) const;
	friend std::ostream & operator<<( std::ostream & out,
			RDISOptimizer const & p );

public:
	// takes the epsilon for the full function, not per factor
	Numeric getEpsilon() { return epsilon; }
	void setEpsilon( Numeric eps ) {
		epsilon = eps;
	}

public:
	virtual void setApproxFactors( bool approxf ) { approxFactors = approxf; }
	virtual void setUseCache( bool useCache_ ) { useCache = useCache_; }
	virtual void setUseBounds( bool useBounds_ ) { useBounds = useBounds_; }
	virtual void setUseLinearCache( bool useLinCache_ ) { useLinearCache = useLinCache_; }

protected:
	// a pool of subdomains (so we're not constantly creating/deleting these)
	SubdomainPool sdpool;

	// pool of reusable connected components
	ComponentPool ccpool;

	// random number generator
	BOOSTNS::random::mt19937 rng;

	// the subspace optimizer to use for choosing values of subsets of variables
	SubspaceOptimizer & ssopt;

	// the options specified on the command line / options file / wherever that
	// control the optimization
	const BOOSTNS::program_options::variables_map options;

	// struct to store the values of the options
	RDISOptionValues optionVals;

	// the variables and all the relevant data + structures about them that we
	// need to maintain and manipulate
	OptimizerVariableData vardata;

	// connectivity graph of the variables
	ConnectivityGraph connGraph;

	// the cache we use to store already-evaluated components
	ComponentCache cache;

	// the top ccomp used for this optimization -- stored to allow for static
	// decomposition across multiple external restarts
	ComponentIP topccomp;

	// maximum error to tolerate (per factor) when splitting variable domains
	Numeric epsilon;

	// true to indicate the connectivity graph needs to be initialized
	bool cgNeedsInit;

	// current number of simplified factors
	VariableCount nSimplifiedFactors;

	// the best / most recent full function eval found yet
	Numeric bestFullFunctionEval;
	Numeric prevFullFunctionEval;

	bool approxFactors; 	// if false, don't approximate factors
	bool useCache;			// if false, don't use the cache
	bool useLinearCache;	// if true, use a linear-sized cache
	bool useBounds;			// if false, don't use bounds (i.e., don't do BnB)

	VariableCount numassigned;
	VariableCountVec numassign, numcachehits, numboundhits;

	// records the order in which variables were first assigned - useful for DBG
	std::string assignmentOrder;
	bool printedAO;

	// temporary vector storage to reduce memory overhead
	VariableIDVec vidstmp;
	VariablePtrVec varstmp;
	FactorPtrVec factmp;
};

} // namespace rdis

#endif // RDIS_RDISOPTIMIZER_H_
