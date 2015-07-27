//
// Created by Abe Friesen on 6/15/15.
//

#include "Component.h"

namespace rdis {

Component::Component( OptimizableFunction & func_,
            const FactorPtrVec & allFctrs, OptimizerVariableData * vd,
            ConnectivityGraph * cg, bool staticDecomp )
        : func( func_ )
        , allFactors( allFctrs )
        , vardata( vd )
        , connGraph( cg )
        , useStaticDecomp( staticDecomp )
        , epsilon( -1.0 )
        , level( -1 )
        , parent( NULL )
        , children( ComponentComparator( vd ) )
        , varsHaveBeenAssigned( false )
        , hasBeenEvaluated( false )
        , lastEvalWasOpt( false )
        , partialEval( NumericMAX )
        , fullEvalWasCached( false )
        , dontCacheFullEval( false )
        , isFminSet( false )
        , isCacheKeySet( false )
        , isUABoundSet( false )
        , isABoundSet( false )
        , fmin( NumericMAX )
        , childFmin( NumericMAX )
        , unassignedBounds( NumericMIN, NumericMAX )
        , assignedBounds( NumericMIN, NumericMAX )
        , curFEStateIsSet( false )
        , numRandomRestarts( 0 )
        , numVAsinceLastRR( 0 )
        , noticedNonmonotoneEvals( false )
        , varsAreQuickAssigned( false )
{
    hasVar.resize( vardata->numVars(), false );
}


Component::~Component() {
    children.clear();
}


void Component::init( Numeric epsilon_, size_t level_,
        Component * parent_, ConnectivityGraph::IntraCCompIterator ccit ) {

    epsilon = epsilon_;
    level = level_;

    reset( true );

    vars.reserve( 5 );
    factors.reserve( 5 );

    // create vars and factors vectors from ccit
    for ( ; !ccit.end(); ++ccit ) {
        ConnectivityGraph::VertexID vxid = *ccit;
        if ( connGraph->isVariableID( vxid ) ) {
            VariableID vid = connGraph->vxidToVID( vxid );
            assert( !vardata->isAssigned( vid ) );
            vars.push_back( vid );
            hasVar.set( vid );
        } else {
            factors.push_back( allFactors[connGraph->vxidToFID( vxid )] );
        }
    }

    parent = parent_;
    if ( parent->hasNonmonotoneEvals() ) setNonmonotoneEvals( true );

    assert( ckey.factors.empty() );
    std::sort( vars.begin(), vars.end() );
    std::sort( factors.begin(), factors.end(), FactorPtrComparator() );
}


void Component::init( Numeric epsilon_, size_t level_,
        const VariablePtrVec & varsIn ) {

    epsilon = epsilon_;
    level = level_;
    parent = NULL;

    assert( level == 0 ); // only call this at the top level

    reset( true );

    vars.reserve( varsIn.size() );
    for ( const Variable * v : varsIn ) {
        if ( !v->isAssigned() ) {
            vars.push_back( v->getID() );
            hasVar.set( v->getID() );
        }
    }

    assert( ckey.factors.empty() );
    std::sort( vars.begin(), vars.end() );

    // create the set of factors for this component
    for ( Factor * f : allFactors ) {
        for ( VariableID vid : vars ) {
            if ( vardata->var( vid )->isAssigned() ) continue;
            if ( f->containsVar( vid ) ) {
                factors.push_back( f );
                break;
            }
        }
    }

    std::sort( factors.begin(), factors.end(), FactorPtrComparator() );
}


void Component::reinit( Numeric epsilon_, Component * parent_ ) {

    epsilon = epsilon_;

    // make sure this has been cleaned if being reused

    assert( parent == parent_ );

    reset( false );
    parent = parent_;

    if ( parent != NULL && parent->hasNonmonotoneEvals() ) {
        setNonmonotoneEvals( true );
    }

    // if this is being reused and everything is the same, then don't recompute
    assert( !vars.empty() && !factors.empty() );
}


void Component::reset( bool fullreset ) {

    assert( !varsAreQuickAssigned );

    if ( fullreset ) {
        vars.clear();
        factors.clear();

        parent = NULL;
        children.clear();

        hasVar.reset();

        assignedVDs.clear();
        assignedvids.clear();

        isCacheKeySet = false;

        ckey.clear();
    }

    varsHaveBeenAssigned = false;
    hasBeenEvaluated = false;
    lastEvalWasOpt = false;

    simplifiedFactors.clear();

    partialEval = 0.0;
    prevsd.reset();
    cursd.reset();
    clearOptSD();
    fullEvalWasCached = false;

    // TODO: put fullEvalWasCached into subdomain

    dontCacheFullEval = false;

    isFminSet = false;
    isUABoundSet = false;
    isABoundSet = false;

    fmin = NumericMAX;
    childFmin = NumericMAX;

    numRandomRestarts = 0;
    numVAsinceLastRR = 0;

    noticedNonmonotoneEvals = false;

    unassignedBounds.assign( NumericMIN, NumericMAX );
    assignedBounds.assign( NumericMIN, NumericMAX );

    curFuncEvalState.clear();
    curFEStateIsSet = false;
}


void Component::clear() {
    reset( true );
}


void Component::onVarsAssigned( const VarDataPtrVec & vds,
        const NumericVec & /*xval*/ ) {

    assert( assignedVDs.empty() || assignedVDs == vds );
    assert( isFminSet );

    if ( assignedVDs.empty() ) {
        assignedVDs = vds;
        assert( assignedvids.empty() );
        assignedvids.clear();
        for ( const VarData * vd : vds ) {
            assignedvids.push_back( vd->v->getID() );
        }
    }

    varsHaveBeenAssigned = true;
    fullEvalWasCached = false;
    isABoundSet = false;

    ++numVAsinceLastRR;

    assert( !areVarsQuickAssigned() );
    for ( ComponentIP & c : children ) {
        assert( !c->areVarsQuickAssigned() );
    }
}


void Component::decompose() {

    createChildren();

    // compute childFmin = fmin(par) - ( partialEval(par) + \Sum_{x \in ch(par)}(fminLB(x)) )
    // which is used to efficiently set the fmins of the children later

    assert( isFminSet ); // we assume fmin has been set already
    childFmin = varsHaveBeenAssigned ? fmin - partialEval : fmin;

    for ( const ComponentIP & c : children ) {
        const Numeric cfLB = c->getUnassignedBounds().lower();
        childFmin -= cfLB;

//        if ( std::isinf( childFmin ) ) {
//            childFmin = childFmin < 0 ? NumericMIN : NumericMAX;
//        }
        assert( !std::isinf( childFmin ) );
    }
}


void Component::onChildEvaluated( Component & child ) {

    assert( child.hasBeenEvald() );
    assert( child.isUABoundSet );

    const Numeric tol = 1e-6; // use epsilon instead?

    const NumericInterval chuab = child.unassignedBounds;

#ifdef DEBUG
    const bool lbwrong = !approxleq( chuab.lower(), child.getOptSD()->fx, tol );
    const bool ubwrong = !approxgeq( chuab.upper(), child.getOptSD()->fx, tol );
    const bool wasCached = child.fullEvalWasCached;

    if ( ( lbwrong || ubwrong ) && !wasCached ) {

        SubdomainCSP sdopt;
        if ( child.haveVarsBeenAssigned() ) {
            sdopt = child.getOptSD();
        }

        std::cout << "child " << child << ": child.fullEval: "
                << child.getOptSD()->fx << ", diffLB: " <<
                ( child.getOptSD()->fx - chuab.lower() ) << ", diffUB: " <<
                ( chuab.upper() - child.getOptSD()->fx ) << " (UAB: " <<
                child.unassignedBounds << ") -- is bound? " <<
                child.getOptSD()->isBound << ": ";
        if ( sdopt ) std::cout << *sdopt;
        else std::cout << "N/A";
        std::cout << std::endl;

        for ( Factor * f : child.getFactors() ) {
            std::cout << f->getID() << " assigned? " << f->isAssigned() <<
                    ", bounds " << f->computeBounds( VariableIDVec( 1, -1 ) ) <<
                    ", allVarsAssigned? " << f->areAllVarsAssigned() <<
                    ", eval? " << ( f->areAllVarsAssigned() ||
                                    f->isAssigned() ? f->eval() : -111111 )
                    << ", NUA " << f->numUnassignedVars() << std::endl;
        }

        std::cout << "parent: " << *this << std::endl;

        child.quickAssignOptVals();
        for ( Factor * f : child.getFactors() ) {
            std::cout << f->getID() << " assigned? " << f->isAssigned() <<
                    ", bounds " << f->computeBounds( VariableIDVec( 1, -1 ) ) <<
                    ", allVarsAssigned? " << f->areAllVarsAssigned() <<
                    ", eval? " << ( f->areAllVarsAssigned() || f->isAssigned() ? f->eval() : -111111 )
                    << ", NUA " << f->numUnassignedVars() << std::endl;
        }
        child.quickUnassignOptVals();
    }

    assert( wasCached || ( !lbwrong && !ubwrong ) );
    assert( approxleq( chuab.lower(), child.getOptSD()->fx, tol ) || wasCached );
    assert( approxgeq( chuab.upper(), child.getOptSD()->fx, tol ) || wasCached );
    // these might assert sometimes...?
#endif // DEBUG

    // replace the LB in the child fmin estimator with the actual computed eval
//	if ( !child.fullEvalIsBound ) {
    childFmin += chuab.lower();
    childFmin -= child.getOptSD()->fx;
//	}

    if ( isABoundSet ) {
        assert( child.isUABoundSet );

        Numeric assignedLB = 0.0;
        Numeric assignedUB = 0.0;

        // TODO: this could be done more efficiently
        for ( const ComponentIP & c : children ) {
//			std::cout << *c << " ";
            if ( c->hasBeenEvald() ) {
//			if ( c->getFullEvalWasCached() ) {
                assignedLB += child.getOptSD()->fx;
                assignedUB += child.getOptSD()->fx;
//				std::cout << c->fullEval << " (" << assignedLB << "), ";
            } else {
                assert( c->isUABoundSet );
                assignedLB += c->unassignedBounds.lower();
                assignedUB += c->unassignedBounds.upper();
//				std::cout << c->originalUnassignedBounds.lower() << " ("
//                      << assignedLB << "), ";
            }
        }

        assert( assignedLB <= assignedUB );
        assignedBounds.assign( assignedLB, assignedUB );
    }
}


std::string Component::toString() const {
    std::string s = "{";
    for ( auto it = vars.begin(); it != vars.end(); ++it ) {
        s += ( it == vars.begin() ? "" : "," );
        s += BOOSTNS::lexical_cast< string >( *it );
    }

    s += " | ";
    s += BOOSTNS::lexical_cast< string >( factors.size() );
//	for ( Factor * f : factors ) {
//		s += ( f == factors.front() ? "" : "," );
//		s += BOOSTNS::lexical_cast< string >( f->getID() );
//	}

    s += "}";
    return s;

}


void Component::onVarUnassigned() {
    if ( hasBeenEvaluated ) return;

    assert( simplifiedFactors.empty() );

    hasBeenEvaluated = true;

//	std::cout << "component " << *this << " full eval " << fullEval <<
//			" (bound: " << ( fullEvalIsBound ? "yes" : "no" ) << ")" << std::endl;

//	if ( cachedEval ) {
//		std::cout << "\tcached eval: " << *cachedEval << std::endl;
//		std::cout << "\tfull eval: " << *optimum << std::endl;
//		std::cout << "\tcache key: " << ckey << ", " << ckey.lookup << std::endl;
//	}

    assert( !areVarsQuickAssigned() );
    for ( ComponentIP & c : children ) {
        assert( !c->areVarsQuickAssigned() );
    }
}


void Component::computeFMin() {
    const NumericInterval uab = getUnassignedBounds();
    fmin = parent->childFmin + uab.lower();
    assert( !std::isnan( fmin ) );
    isFminSet = true;
//	std::cout << "fmin for component " << *this << " is " << fmin << std::endl;
}


bool Component::setCachedValue( const SubdomainCSP & sdopt ) {

    assert( sdopt );

    // pretend this var was chosen to be assigned
    assignedVDs.clear();
    for ( VariableID vid : sdopt->vids ) {
        VarData & vd( vardata->get( vid ) );
        assignedVDs.push_back( &vd );
    }
//    clearOptSD();
    varsHaveBeenAssigned = true;

    // set the computed value of this component as the cache value
    assert( !getOptSD() );
    setOptSD( sdopt );
    fullEvalWasCached = true;
    isABoundSet = false;

    // mark this as having been evaluated
    hasBeenEvaluated = true;

    return true;
}


void Component::quickUnassignChildren() {
    for ( ComponentIP & c : children ) {
        if ( c->areVarsQuickAssigned() ) c->quickUnassignOptVals();
    }
}


void Component::quickAssignOptVals( bool print ) {
    assert( haveVarsBeenAssigned() && hasBeenEvald() /*&& !getFullEvalIsBound()*/ );
    const SubdomainCSP & fopt = getOptSD();

    assert( fopt );
    if ( !fopt ) return;
    if ( print ) std::cout << "quick assign of SD: " << *fopt << std::endl;
    assert( !fopt->isBound );

    quickAssignFromSD( *fopt, print );
    varsAreQuickAssigned = true;
}


void Component::quickUnassignOptVals() {
    assert( varsAreQuickAssigned );
    quickUnassignFromSD( *getOptSD() );
    varsAreQuickAssigned = false;
}


void Component::quickAssignFromSD( const Subdomain & sd, bool print ) {
    if ( print ) {
        std::cout << "quick assigning vars " << sd.vids << " vals " << sd.x
                << std::endl;
    }

    for ( size_t i = 0; i < sd.vids.size(); ++i ) {
        vardata->var( sd.vids[i] )->assign( sd.x[i] );
        func.onVarAssigned( sd.vids[i], sd.x[i] );
    }

    for ( const SubdomainCSP & c : sd.children ) {
        assert( c );
        quickAssignFromSD( *c, print );
    }
}

void Component::quickUnassignFromSD( const Subdomain & sd, bool print ) {
    if ( print ) std::cout << "quick unassigning vars " << sd.vids << std::endl;

    for ( VariableID vid : sd.vids ) {
        vardata->var( vid )->unassign();
        func.onVarUnassigned( vid );
    }

    for ( const SubdomainCSP & c : sd.children ) {
        assert( c );
        quickUnassignFromSD( *c );
    }
}


NumericInterval Component::getUnassignedBounds() {
    if ( !isUABoundSet ) {
        assert( assignedVDs.empty() || !getRepresentativeVD()->v->isAssigned() );
        unassignedBounds = computeBounds();
        isUABoundSet = true;
    }

    return unassignedBounds;
}


// get the lower bound over the evaluation of this component (over all
// factors, given that the chosen variable in this comp. is assigned)
NumericInterval Component::getAssignedBounds() {
    if ( !isABoundSet ) {
        assert( !assignedVDs.empty() && getRepresentativeVD()->v->isAssigned() );
        assignedBounds = computeBounds();
        isABoundSet = true;
    }

    return assignedBounds;
}


void Component::createChildren() {
//	assert( assignedVD != NULL && assignedVD->v->isAssigned() );

    if ( !useStaticDecomp ) {
        assert( children.empty() || (*children.begin())->hasPool() );
        children.clear();
    } else {
        if ( !children.empty() ) {
            for ( const ComponentIP & cc : children ) {
                cc->reinit( epsilon, this );
            }
            return;
        }
    }

    cclabels.clear();
    cclabels.reserve( vars.size() / 2 + 1 );

    // create the child subcomponents of this component

    for ( VariableID vid : vars ) {
//		if ( varsHaveBeenAssigned && vid == assignedVD->v->getID() ) continue;
        if ( varsHaveBeenAssigned && vardata->isAssigned( vid ) ) continue;
        cclabels.emplace( connGraph->getComponentLabel( vid ) );
    }

    for ( ConnectivityGraph::VertexID ccvxid : cclabels ) {
        ComponentIP child;
        if ( pool != NULL ) {
//			std::cout << "ccomp pool size: " << pool->size() << std::endl;
            child = pool->get( func, allFactors, vardata, connGraph,
                    useStaticDecomp );
        } else {
            assert( false );
            child.reset( new Component( func, allFactors, vardata, connGraph,
                    useStaticDecomp ) );
        }
        child->init( epsilon, level+1, this,
                     connGraph->getIntraCCompIterator( ccvxid ) );
        children.insert( child );
    }
}


void Component::createCacheKey() {
    assert( ckey.factors.empty() );
    assert( ckey.unassignedVars.empty() );
    assert( ckey.assignedVars.empty() );

    // add these vars to the cache key
    for ( auto & vid : vars ) {
        assert( vardata->isUnassigned( vid ) );
        ckey.unassignedVars.emplace( vid );
    }

    // set the ckey factors as the component's factors
    ckey.factors.clear();
    for ( Factor * f : factors ) {
        ckey.factors.emplace( f->getID() );
    }

    // add all assigned variables in factors to the cache key
    for ( Factor * f : factors ) {
        if ( f->numAssignedVars() > 0 ) {
            size_t nadded( 0 );

            for ( auto & v : f->getVariables() ) {
                if ( !v->isAssigned() ||
                     !f->wasVarAssignedBeforeSimpl( v->getID() ) ) {
                    continue;
                }

                const VarData & vd = vardata->get( v->getID() );
                ckey.assignedVars.emplace( vd.assignmentIndex, v->getID(), v  );

                if ( ++nadded >= f->numAssignedVars() ) break;
            }
        }
    }

    isCacheKeySet = true;
}


NumericInterval Component::computeBounds() {
    NumericInterval bounds = func.computeBounds( factors,
            haveVarsBeenAssigned() ? getRepresentativeVD()->v->getID() : -1 );

    assert( !std::isnan( bounds.lower() ) && !std::isnan( bounds.upper() ) );

    return bounds;
}



bool ComponentComparator::operator()( const ComponentIP & lhs,
        const ComponentIP & rhs ) const {
    // process smaller components first as this will give us better bounds
    // for the larger components
    return ( lhs->getVars().size() < rhs->getVars().size() );
}


std::ostream & operator<<( std::ostream & os, const Component & c ) {
    os << c.toString();
    return os;
}


} // namespace rdis
