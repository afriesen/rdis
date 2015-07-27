//
// Created by Abe Friesen on 6/15/15.
//

#ifndef RDIS_COMPONENT_H
#define RDIS_COMPONENT_H

#include "IntrusivePtrPool.h"

#include "State.h"
#include "OptimizableFunction.h"
#include "ComponentCache.h"
#include "ConnectivityGraph.h"
#include "OptimizerVariableData.h"


namespace rdis {

class Component;

typedef BOOSTNS::intrusive_ptr< Component > ComponentIP;
typedef BOOSTNS::intrusive_ptr< const Component > ComponentCIP;

typedef IntrusivePtrPool< Component, OptimizableFunction &,
        const FactorPtrVec &, OptimizerVariableData *, ConnectivityGraph *,
        bool > ComponentPool;

// a comparator to order components for optimization
class ComponentComparator {
public:
    ComponentComparator( OptimizerVariableData * ovd_ )
            : ovd( ovd_ )
    {}
    bool operator()( const ComponentIP & lhs, const ComponentIP & rhs ) const;

private:
    OptimizerVariableData * ovd;
};

typedef BOOSTNS::container::flat_multiset< ComponentIP, ComponentComparator >
            Subcomponents;


// class representing a connected component of variables given some specific
// set of assignments and resulting graph state
class Component
        : public IntrusivePtrPoolObj< Component, ComponentPool > {

public:
    Component( OptimizableFunction & func_, const FactorPtrVec & allFctrs,
            OptimizerVariableData * vd, ConnectivityGraph * cg,
            bool staticDecomp );
    ~Component();

    // fully initialize this component
    void init( Numeric epsilon_, size_t level_, Component * parent_,
            ConnectivityGraph::IntraCCompIterator ccit );
    void init( Numeric epsilon_, size_t level_,
            const VariablePtrVec & vars_ );

    // re-initialize this component (parent must be same and nothing must
    // have changed in this component)
    void reinit( Numeric epsilon_, Component * parent_ );

    // reset this component so it can be re-used with diff vars, etc.
    void reset( bool fullreset );

public:
    // called by IntrusivePtrPool -- equivalent to ::reset(true)
    virtual void clear();

public:
    // called when a variable has been chosen to be assigned from this
    // component, and is then assigned
    void onVarsAssigned( const VarDataPtrVec & vds,
            const NumericVec & xval );

    // decompose this component into its subcomponents (based on the
    // assignment and resulting simplification of the factors and
    // connectivity graph)
    void decompose();

    // called after a child has been fully (or partially) evaluated
    // -- updates fmin estimates
    void onChildEvaluated( Component & child );

    // called when the variable chosen from this component is unassigned
    void onVarUnassigned();

    // compute and set the fmin for this component based on the fmin of its
    // parent and siblings
    void computeFMin();

    // this component's evaluation was retrieved from the cache so set the
    // retrieved value in this component -- return true if this was
    // successful, false if this should be computed anyway
    bool setCachedValue( const SubdomainCSP & sdopt );

    // helper function to perform quick unassignment for all children ccomps
    void quickUnassignChildren();

    // perform a quick (un)assignment of the variables in this component to
    // their optimal values (i.e., don't alter things like the conn. graph
    // or anything else -- just the vars themselves, so factor computations
    // are more accurate)
    void quickAssignOptVals( bool print = false );
    void quickUnassignOptVals();

public:
    void quickAssignFromSD( const Subdomain & sd, bool print = false );
    void quickUnassignFromSD( const Subdomain & sd, bool print = false );

public:
    Numeric getEpsilon() const { return epsilon; }
    size_t getLevel() const { return level; }
    const VariableIDVec & getVars() const { return vars; }
    const FactorPtrVec & getFactors() const { return factors; }

    // return the parent / children (subcomponents) of this component
    Component * getParent() const { return parent; }
    const Subcomponents & getChildren() const { return children; }

    bool containsVar( VariableID vid ) const {
        if ( vid < 0 || vid >= hasVar.size() ) return false;
        return hasVar.test( vid );
    }

    VarDataPtrVec & getVarsToAssign() { return varsToAssign; }
    NumericVec & getValToAssign() { return valToAssign; }

    // get the VarDataVec of the assigned variables from this component
    const VarDataPtrVec & getAssignedVDs() const { return assignedVDs; }
    const VariableIDVec & getAssignedVids() const { return assignedvids; }

    VarData * getRepresentativeVD() const {
        return ( assignedVDs.empty() ? NULL : assignedVDs.front() );
    }

    bool haveVarsBeenAssigned() const { return varsHaveBeenAssigned; }

    SortedFactorIDs & getSimplifiedFactors() { return simplifiedFactors; }
    const SortedFactorIDs & getSimplifiedFactors() const { return simplifiedFactors; }

    bool hasBeenEvald() const { return hasBeenEvaluated; }
    void setHasBeenEvald( bool hbe ) { hasBeenEvaluated = hbe; }

    bool wasLastEvalOpt() const { return lastEvalWasOpt; }
    void setLastEvalWasOpt( bool lewo ) { lastEvalWasOpt = lewo; }

    Numeric getPartialEval() const { return partialEval; }
    void setPartialEval( Numeric newpe ) { partialEval = newpe; }

    const SubdomainCSP & getPrevSD() const { return prevsd; }
    void setPrevSD( const SubdomainCSP & sd ) { prevsd = sd; }

    const SubdomainCSP getCurSD() const { return cursd; }
    void setCurSD( const SubdomainSP & sd ) { cursd = sd; }
    void updateCurSD( Numeric curfx, Numeric curferr, bool isbound,
            SubdomainCSP child = SubdomainCSP() ) {
        assert( cursd );
        cursd->fx = curfx;
        cursd->ferr = curferr;
        cursd->isBound = isbound;
        if ( child ) cursd->addChild( child );
    }

    const SubdomainCSP & getOptSD() const { return optsd; }
    void clearOptSD() { optsd.reset(); }
    void setOptSD( const SubdomainCSP & sd ) { optsd = sd; }
    void setOptSD( SubdomainSP sd, Numeric fullEval, bool isbound ) {
        if ( optsd ) sd->setAsCopy( *optsd );
        sd->fx = fullEval;
        sd->isBound = isbound;
        optsd = sd;
    }

    bool getFullEvalWasCached() const { return fullEvalWasCached; }

    bool dontCache() const { return dontCacheFullEval; }

    Numeric getFMin() const { assert( isFminSet ); return fmin; }
    void setFMin( Numeric newfmin ) { fmin = newfmin; isFminSet = true; }

    // get the cache key for this component
    CacheKey & getCacheKey() {
        if ( !isCacheKeySet ) createCacheKey();
        return ckey;
    }

    size_t getNumRandomRestarts() const { return numRandomRestarts; }
    void incrementNumRandomRestarts() {
        ++numRandomRestarts;
        numVAsinceLastRR = 0;
    }

    size_t getNumVAsinceLastRR() const { return numVAsinceLastRR; }

    bool hasNonmonotoneEvals() const { return noticedNonmonotoneEvals; }
    void setNonmonotoneEvals( bool nme ) { noticedNonmonotoneEvals = nme; }

    // get the upper+lower bounds over the evaluation of this component (over all
    // factors given that all variables in this component are unassigned)
    NumericInterval getUnassignedBounds();

    // get the lower+upper bounds over the evaluation of this component (over all
    // factors, given that the chosen variable in this comp. is assigned)
    NumericInterval getAssignedBounds();

    State & getCurFuncEvalState() { return curFuncEvalState; }

    bool getCurFEStateIsSet() const { return curFEStateIsSet; }
    void setCurFEStateIsSet( bool cfeis ) { curFEStateIsSet = cfeis; }

    bool areVarsQuickAssigned() const { return varsAreQuickAssigned; }

public:
    std::string toString() const;
    friend std::ostream & operator<<( std::ostream & os, const Component & c );

protected:
    // determine and create the subcomponents that are the children of this
    // component now that assignedVD has been assigned
    void createChildren();

    // create the cache key for this component
    void createCacheKey();

    // compute the current (i.e., based on current var+factor assignments)
    // upper+lower bound from the factors making up this component
    NumericInterval computeBounds();

protected:
    // see RDISOptimizer.h for descriptions of these
    OptimizableFunction & func;
    const FactorPtrVec & allFactors;
    OptimizerVariableData * vardata;
    ConnectivityGraph * connGraph;
    bool useStaticDecomp;
    Numeric epsilon;

    // the level (in the optimization) that this component is at; components
    // form a tree, so this corresponds to the depth in the tree (starts at 0).
    size_t level;

    // the list of variables in this component
    VariableIDVec vars;

    // the set of factors in this component
    FactorPtrVec factors;

    // the parent component (of which this is a subcomponent) and the
    // children (subcomponents) of this component
    Component * parent;
    Subcomponents children;

    // flags indicating whether or not a variable is present in this ccomp
    BOOSTNS::dynamic_bitset<> hasVar;

    // the variables and their values that will be assigned (mainly exists
    // here for memory efficiency)
    VarDataPtrVec varsToAssign;
    NumericVec valToAssign;

    // the variable(s) that were chosen to be assigned from this component
    // (thus resulting in the decomposition into the children subcomponents)
    VarDataPtrVec assignedVDs;
    VariableIDVec assignedvids;

    // set to true to indicate that a variable has been assigned from within
    // this component
    bool varsHaveBeenAssigned;

    // current list of factors that have were simplified when the assignedVD
    // was assigned its current value
    SortedFactorIDs simplifiedFactors;

    // set to true once this component has been evaluated (aka, either it
    // was optimized, it hit a bound, or it was retrieved from the cache)
    // and now we're done processing it
    bool hasBeenEvaluated;

    // set to true after evaluation of one value assignment if this assignment
    // became the new optimal value
    bool lastEvalWasOpt;

    // the evaluation of only those factors that became fully assigned when this
    // components chosen variables were assigned
    Numeric partialEval;

    // the most recent, the current (can be partial), and the optimal subdomain
    // found so far (they may be equal)
    SubdomainCSP prevsd, optsd;
    SubdomainSP cursd;

    // true to indicate that the full eval was retrieved from the cache
    bool fullEvalWasCached;

    // set to true to not put the full evaluation in the cache (necessary if
    // we're recomputing something that was in the cache in the first place)
    bool dontCacheFullEval;

    // flags to indicate whether cache key and unassigned and assigned
    // bounds have been computed
    bool isFminSet, isUABoundSet, isABoundSet, isCacheKeySet;

    // the minimum f(x) value in this subdomain that has either been observed
    // or was specified by the user / optimizer (which is why it exists
    // separately from optsd)
    Numeric fmin;

    // used for determining children's fmins, equal to:
    //  cfmin = fmin / ( partialEval * \Prod_{c \in chldrn} fminLB(c) )
    Numeric childFmin;

    // lower+upper bounds on this component's evaluation, given that all vars are
    // unassigned / some vars are assigned (but this component hasn't been split
    // into its subcomponents yet)
    NumericInterval unassignedBounds, assignedBounds;

    // key used for cache lookups for this component
    CacheKey ckey;

    // state used to compute the current function eval (cache it so we don't
    // have to recompute all the unchanging parts each time)
    State curFuncEvalState;
    bool curFEStateIsSet;

    // the number of times that vals were assigned based on a random restart
    // (for gradient descent version)
    size_t numRandomRestarts;

    // the number of var assigns since last random restart
    size_t numVAsinceLastRR;

    // true if we've noticed that function evals are non-monotonic (and thus
    // might not be converging) for this component
    bool noticedNonmonotoneEvals;

    // true if the vars in this component were quick-assigned and are still
    // in that state
    bool varsAreQuickAssigned;

    // store this for memory-reuse (efficiency) purposes -- clear() before use
    CGVertexIDSet cclabels;

}; // class Component

} // namespace rdis

#endif //RDIS_OPTIMIZERCOMPONENT_H
