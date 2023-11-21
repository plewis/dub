#pragma once

namespace proj {

    class Particle {
            
        public:
        
            struct CoalProposal {
                // This object is populated during a call to proposeCoalescence.
                
                // Stack of increments to species tree and all gene trees
                stack<double>               _dt;
                
                // Stack of ancstral species added in the course of the proposal
                stack<tuple<unsigned, unsigned, SMCGlobal::species_t> > _speciations;
                
                // The index of the gene in which the proposed coalescence occurred
                int                         _coal_gene;
                
                // The species in which the the proposed coalescence occurred
                SMCGlobal::species_t        _coal_spp;

                // The ancestral node added to create the proposed coalescence event
                Node *                      _coal_anc;
                
                void clear() {
                    _dt = {};
                    _speciations = {};
                    _coal_gene = -1;
                    _coal_spp = 0;
                    _coal_anc = nullptr;
                }
            };

            Particle();
            Particle(const Particle & other);
            ~Particle();
            
            void clear();
            void setData(Data::SharedPtr data);
            
            void resetSpeciesForest();
            void resetGeneForests(bool compute_partials);
            
            void incrementSpeciations() {_nspeciations++;}
            unsigned getSpeciations() const {return _nspeciations;}
            void clearSpeciations() {_nspeciations = 0;}
            
            void threadComputePartials(unsigned first, unsigned last);
            void computeAllPartials();
            void stowAllPartials();
            
            double calcTotalCoalRate(vector<SMCGlobal::species_tuple_t> & species_tuples, double speciation_increment);
            void advanceAllLineagesBy(double dt);
            double proposeCoalescence(unsigned seed, unsigned step, unsigned pindex, CoalProposal & proposal, bool compute_partial, bool make_permanent);
            void finalizeProposal(CoalProposal & proposal);
            void reverseProposal(CoalProposal & proposal);
            void advance(unsigned step, unsigned pindex, CoalProposal & proposal, bool compute_partial, bool make_permanent);
            
            unsigned getCount() const {return _count;}
            void setCount(unsigned cnt) {_count = cnt;}
            
#if defined(USING_MULTITHREADING)
            unsigned getXtra() const {return _xtra;}
            void setXtra(unsigned x) {_xtra = x;}
            
            unsigned getBeginIndex() const {return _begin_index;}
            void setBeginIndex(unsigned x) {_begin_index = x;}
#endif

#if defined(DEBUGGING_SANITY_CHECK)
            void debugStoreForests();
            void debugCheckForests();
#endif

            double getLogWeight() const {return _log_weight;}
            
            Lot::SharedPtr getLot() const {return _lot;}
            void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            
            enum last_event_t {
                LAST_EVENT_UNDEFINED = 0,
                LAST_EVENT_SPECIATION = 1,
                LAST_EVENT_COALESCENCE = 2
            };
            
            void setLastEvent(last_event_t last) {_last_event = last;}
            bool lastEventUndefined() const {return _last_event == LAST_EVENT_UNDEFINED;}
            bool lastEventSpeciation() const {return _last_event == LAST_EVENT_SPECIATION;}
            bool lastEventCoalescence() const {return _last_event == LAST_EVENT_COALESCENCE;}
            
            void refreshHeightsInternalsPreorders();
            
            SpeciesForest       & getSpeciesForest();
            const SpeciesForest & getSpeciesForest() const;

            vector<GeneForest>       & getGeneForests();
            const vector<GeneForest> & getGeneForests() const;
            
            GeneForest       & getGeneForest(unsigned gene);
            const GeneForest & getGeneForest(unsigned gene) const;
            
            unsigned debugCountNumCoalEvents() const;
            void debugCheckPartials() const;
            void debugShowAllGeneForests() const;
            void debugCheckAllPrevSpeciesStacksEmpty() const;
            
            void copyParticleFrom(const Particle & other);
            void operator=(const Particle & other);
                                
        protected:
                
            PartialStore::partial_t pullPartial(unsigned gene);
            void stowPartial(unsigned gene, Node * nd);

            Data::SharedPtr        _data;
            vector<GeneForest>     _gene_forests;
            SpeciesForest          _species_forest;
            unsigned               _nspeciations;
            last_event_t           _last_event;
            unsigned               _count;
#if defined(USING_MULTITHREADING)
            unsigned               _xtra;
            unsigned               _begin_index;
#endif

#if defined(DEBUGGING_SANITY_CHECK)
            string                  _debug_sfbefore;
            vector<string>          _debug_gfbefore;
#endif

            mutable double         _log_weight;
            mutable Lot::SharedPtr _lot;
    };
        
    inline Particle::Particle() {
        _lot.reset(new Lot());
        clear();
    }

    inline Particle::Particle(const Particle & other) {
        _lot.reset(new Lot());
        copyParticleFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _log_weight = 0.0;
        _gene_forests.clear();
        _species_forest.clear();
        _nspeciations = 0;
        _last_event = LAST_EVENT_UNDEFINED;
    }
    
    inline void Particle::setData(Data::SharedPtr data) {
        _data = data;
    }
            
    inline void Particle::resetSpeciesForest() {
        _species_forest.createTrivialForest();
    }

    inline void Particle::resetGeneForests(bool compute_partials) {
        assert(SMCGlobal::_ntaxa > 0);
        assert(SMCGlobal::_ngenes > 0);
        _gene_forests.clear();
        _gene_forests.resize(SMCGlobal::_ngenes);
        unsigned g = 0;
        for (auto & gf : _gene_forests) {
            gf.setData(_data);
            gf.setRelRate(SMCGlobal::_relrate_for_gene[g]);
            gf.setGeneIndex(g++);
            gf.createTrivialForest(compute_partials);
        }
    }

    inline void Particle::threadComputePartials(unsigned first, unsigned last) {
        for (unsigned k = first; k < last; k++) {
            _gene_forests[k].computeAllPartials();
        }
    }
    
    inline void Particle::computeAllPartials() {
//#if defined(USING_MULTITHREADING)
//        vector<thread> threads;
//        for (unsigned i = 0; i < SMCGlobal::_nthreads; i++) {
//            threads.push_back(thread(&Particle::threadComputePartials,
//                this,
//                SMCGlobal::_thread_first_gene[i],
//                SMCGlobal::_thread_last_gene[i])
//            );
//        }
//
//        // The join function causes this loop to pause until the ith thread finishes
//        for (unsigned i = 0; i < threads.size(); i++) {
//            threads[i].join();
//        }
//#else
        for (auto & gf : _gene_forests) {
            gf.computeAllPartials();
        }
//#endif
    }
    
    inline void Particle::stowAllPartials() {
        for (auto & gf : _gene_forests) {
            gf.stowAllPartials();
        }
    }
    
    inline PartialStore::partial_t Particle::pullPartial(unsigned gene) {
        assert(gene < _gene_forests.size());

        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
#if defined(USING_MULTITHREADING)
        {
            lock_guard<mutex> guard(SMCGlobal::_mutex);
            ptr = ps.getPartial(gene);
        }
#else
        ptr = ps.getPartial(gene);
#endif
        return ptr;
    }

    inline void Particle::stowPartial(unsigned gene, Node * nd) {
        assert(gene < _gene_forests.size());
        assert(nd);
        assert(nd->_partial);
#if defined(USING_MULTITHREADING)
        {
            lock_guard<mutex> guard(SMCGlobal::_mutex);
            ps.putPartial(gene, nd->_partial);

            // Decrement shared pointer reference count
            nd->_partial.reset();
        }
#else
        ps.putPartial(gene, nd->_partial);
        
        // Decrement shared pointer reference count
        nd->_partial.reset();
#endif
    }

    inline void Particle::finalizeProposal(CoalProposal & proposal) {
        //TODO: Particle::finalizeProposal
        // Clear _prev_species_stack for all nodes in species tree and all gene trees
        unsigned num_speciations = (unsigned)proposal._speciations.size();
        if (num_speciations > 0) {
            // Empty the _prev_species_stack for all nodes in the species tree
            for (auto & nd : _species_forest._nodes) {
                nd.emptyPrevSpeciesStack();
            }
            
            // Empty the _prev_species_stack for all nodes in all gene trees
            for (auto & gf : _gene_forests) {
                for (auto & nd : gf._nodes) {
                    nd.emptyPrevSpeciesStack();
                }
            }
        }
        
        debugCheckAllPrevSpeciesStacksEmpty();
    }
    
    inline void Particle::reverseProposal(CoalProposal & proposal) {
        //TODO: Particle::reverseProposal
        unsigned num_speciations = (unsigned)proposal._speciations.size();
        unsigned num_increments = (unsigned)proposal._dt.size();
        assert(num_increments = num_speciations + 1);
        
        // Identify the nodes involved in the coalescent event
        Node * anc = proposal._coal_anc;
        assert(anc);
        Node * lchild = anc->_left_child;
        assert(lchild);
        Node * rchild = lchild->_right_sib;
        assert(rchild);
        
        // Return partial to be reused next time
        stowPartial(proposal._coal_gene, anc);

        // Reverse the coalescence join
        GeneForest & gf = _gene_forests[proposal._coal_gene];
        gf.unjoinLineagePair(anc, lchild, rchild);
        gf.stowNode(anc);
        anc = nullptr;
        
        // Reverse the increment associated with the coalescent event
        double dt = proposal._dt.top();
        proposal._dt.pop();
        advanceAllLineagesBy(-dt);
        
        // Now reverse all speciation events, if there were any
        while (!proposal._speciations.empty()) {
            auto spptuple = proposal._speciations.top();
            proposal._speciations.pop();
            
            unsigned left_pos        = get<0>(spptuple);
            unsigned right_pos       = get<1>(spptuple);
            SMCGlobal::species_t spp = get<2>(spptuple);
            
            // Identify the nodes involved in the speciation event
            anc = _species_forest.findSpecies(spp);
            assert(anc);
            lchild = anc->_left_child;
            assert(lchild);
            rchild = lchild->_right_sib;
            assert(rchild);
            
            assert(fabs(anc->getEdgeLength()) < 0.00001);
            
            // Reverse the speciation join
            _species_forest.unjoinLineagePair(anc, lchild, rchild);

            // Update lineage vector
            _species_forest.addTwoRemoveOneAt(_species_forest._lineages, left_pos, lchild, right_pos, rchild, anc);
            
            // Return anc to the pool of unused nodes
            _species_forest.stowNode(anc);
            anc = nullptr;
            
            // Reverse the increment associated with the speciation event
            dt = proposal._dt.top();
            proposal._dt.pop();
            advanceAllLineagesBy(-dt);
        }
        
        //temporary!
        //if (SMCGlobal::_debugging) {
        //    output(format("species forest after reverting proposal:\n  %s\n") % _species_forest.makeNewick(9, true, false),2);
        //}

        if (num_speciations > 0) {
            // Revert all nodes in all gene trees for which _prev_species_stack is non-empty
            for (auto & gf : _gene_forests) {
                for (auto & nd : gf._nodes) {
                    nd.revertSpecies();
                }
            }
        }
        
        debugCheckAllPrevSpeciesStacksEmpty();
    }
        
    inline void Particle::debugCheckAllPrevSpeciesStacksEmpty() const {
#if defined(DEBUGGING_SANITY_CHECK)
        unsigned node_index = 0;
        
        // Ensure that all species tree nodes have an empty _prev_species_stack
        for (auto & nd : _species_forest._nodes) {
            if (nd.canRevertSpecies()) {
                cerr << str(format("Node %d in species forest has a non-empty _prev_species_stack") % node_index) << endl;
            }
            assert(!nd.canRevertSpecies());
            node_index++;
        }

        // Ensure that all gene tree nodes have an empty _prev_species_stack
        for (auto & gene_forest : _gene_forests) {
            node_index = 0;
            for (auto & nd : gene_forest._nodes) {
                if (nd.canRevertSpecies()) {
                    cerr << str(format("Node %d in gene forest %d has nonempty species stack") % node_index % gene_forest.getGeneIndex()) << endl;
                }
                assert(!nd.canRevertSpecies());
                node_index++;
            }
        }
#endif
    }
    
    inline double Particle::proposeCoalescence(unsigned seed, unsigned step, unsigned pindex, CoalProposal & proposal, bool compute_partial, bool make_permanent) {
        //TODO: Particle::proposeCoalescence
        setSeed(seed);
        
        advance(step, pindex, proposal, compute_partial, make_permanent);
        while (lastEventSpeciation()) {
            advance(step, pindex, proposal, compute_partial, make_permanent);
        }
        
        double log_weight = getLogWeight();
        return log_weight;
    }
    
    inline double Particle::calcTotalCoalRate(vector<SMCGlobal::species_tuple_t> & species_tuples, double speciation_increment) {
        double total_rate = 0.0;
        for (auto & gene_forest : _gene_forests) {
            total_rate += gene_forest.calcTotalRate(species_tuples, speciation_increment);
        }
        return total_rate;
    }
    
    inline void Particle::advanceAllLineagesBy(double dt) {
        _species_forest.advanceAllLineagesBy(dt);
        for (auto & gene_forest : _gene_forests) {
            gene_forest.advanceAllLineagesBy(dt);
        }
    }

    inline void Particle::advance(unsigned step, unsigned pindex, CoalProposal & proposal, bool compute_partial, bool make_permanent) {
        //TODO: Particle::advance
        // Create species_tuples vector. Each 3-tuple entry stores:
        //  1. number of lineages
        //  2. gene index
        //  3. species within gene (ignored if species tree)
        vector<SMCGlobal::species_tuple_t> species_tuples;
        
        // Draw a speciation increment Delta. Note: speciation_increment will
        // equal "infinity" if species tree is complete.
        auto incr_rate = _species_forest.drawIncrement(_lot);
        double Delta = incr_rate.first;
        double speciation_rate = incr_rate.second;
        
        // Visit each species within each locus, storing the number of lineages in
        // species_tuples and computing the total coalescence rate (total_rate).
        double total_rate = calcTotalCoalRate(species_tuples, Delta);
        
        // Draw coalescence increment delta ~ Exponential(total_rate)
        double delta = (total_rate > 0.0 ? -log(1.0 - _lot->uniform())/total_rate : SMCGlobal::_infinity);
        
        // If delta > Delta, then a speciation event happened; otherwise a coalescent event happened.
        bool is_speciation = (delta > Delta);
                
        //temporary!
        //if (SMCGlobal::_debugging && is_speciation) {
        //    output(format("species forest before speciation:\n  %s\n") % _species_forest.makeNewick(9, true, false),2);
        //}

        // Advance all forests by increment dt
        double dt = is_speciation ? Delta : delta;
        advanceAllLineagesBy(dt);
        
        // Record the increment
        proposal._dt.push(dt);

        if (is_speciation) {
            _last_event = Particle::LAST_EVENT_SPECIATION;
            
            // Create speciation event (species forest saves nodes involved to
            // a stack to allow reversion)
            SMCGlobal::species_t left_spp;
            SMCGlobal::species_t right_spp;
            SMCGlobal::species_t anc_spp;

            // Pull next available node
            Node * anc = _species_forest.pullNode();
            
            // The speciationEvent function joins two nodes chosen randomly and stores
            // the species involved in variables left_spp, right_spp, and anc_spp
            unsigned left_pos, right_pos;
            _species_forest.speciationEvent(_lot, anc, left_pos, right_pos, left_spp, right_spp, anc_spp);
            
            //temporary!
            //if (SMCGlobal::_debugging) {
            //    output(format("species forest after speciation:\n  %s\n") % _species_forest.makeNewick(9, true, false),2);
            //}
            
            // Record the species added to the species tree
            // left_pos and right_pos are needed to ensure that things get
            // put back exactly the way they were when speciations are reverted
            proposal._speciations.push(make_tuple(left_pos, right_pos, anc_spp));
                        
            // Advise all gene trees of the change in the species tree
            // Nodes that are reassigned save their previous state to
            // their _prev_species_stack to allow reversion
            for (auto & gene_forest : _gene_forests) {
                gene_forest.mergeSpecies(left_spp, right_spp, anc_spp);
            }
        }
        else {
            _last_event = Particle::LAST_EVENT_COALESCENCE;
            
            assert(SMCGlobal::_prior_prior);    // prior-post not yet ready
            if (SMCGlobal::_prior_prior) {
                // Choose which locus and species in which to have a coalescence

                // Create vector of coalescence rates by extracting first element
                // of each species tuple (number of lineages n) and constructing rate
                // as n(n-1)/theta. The theta is omitted because it is identical for
                // all rates.
#if 1
                // This is correct
                vector<double> probs(species_tuples.size());
                transform(species_tuples.begin(), species_tuples.end(), probs.begin(), [total_rate](SMCGlobal::species_tuple_t & spp_tuple){
                    double n = (double)get<0>(spp_tuple);
                    return n*(n - 1.0)/(SMCGlobal::_theta*total_rate);
                });
                
#if defined(DEBUGGING_SANITY_CHECK)
                double check = accumulate(probs.begin(), probs.end(), 0.0);
                assert(fabs(check - 1.0) < 0.0001);
#endif
                
#else
                // This is clearly incorrect
                vector<double> coal_rates(species_tuples.size());
                transform(species_tuples.begin(), species_tuples.end(), coal_rates.begin(), [](SMCGlobal::species_tuple_t & spp_tuple){
                    double n = (double)get<0>(spp_tuple);
                    return n;
                });
                                
                vector<double> probs(coal_rates.size());
                SMCGlobal::normalizeRates(coal_rates, probs);
#endif
                
                // Choose gene and species
                unsigned             which = SMCGlobal::multinomialDraw(_lot, probs);
                assert(which < probs.size());
                unsigned                 g = get<1>(species_tuples[which]);
                SMCGlobal::species_t   spp = get<2>(species_tuples[which]);

                GeneForest & gene_forest = _gene_forests[g];
                
                // Pull next available node
                Node * anc = gene_forest.pullNode();
                if (compute_partial)
                    anc->_partial = pullPartial(g);
                    
                // The coalescenceEvent function joins two nodes chosen randomly
                // from species spp in gene g. Returns calculated log weight if
                // compute_partial is true; otherwise returns 0.0 for log weight.
                _log_weight = gene_forest.coalescentEvent(_lot, spp, anc, compute_partial, make_permanent);
                
                // Save information about the coalescent event proposed
                proposal._coal_gene = g;
                proposal._coal_spp  = spp;
                proposal._coal_anc  = anc;
                
                if (compute_partial) {
                    // Adjust weight for fact that proposal differs from prior
                    if (_species_forest.getNumLineages() > 1) {
                        _log_weight -= speciation_rate*(delta - Delta);
                        _log_weight -= log(speciation_rate);
                    }
                    assert(!isnan(_log_weight));
                    assert(!isinf(_log_weight));
                }
            }
            else {
                // prior post code will go here when written
            }
            
            // Only speciation events should have log_weight 0.0, unless we're
            // simulating, in which case calculate_partial will be false
            assert(_log_weight != 0.0 || !compute_partial);
        }
    
        if (is_speciation) {
            incrementSpeciations();
        }
    }

#if defined(DEBUGGING_SANITY_CHECK)
    inline void Particle::debugStoreForests() {
        _debug_sfbefore = _species_forest.makeNewick(6,true,false);
        _debug_gfbefore.clear();
        for (auto gf : _gene_forests) {
            _debug_gfbefore.push_back(gf.makeNewick(6,true,false));
        }
    }
    
    inline void Particle::debugCheckForests() {
        string sfafter = _species_forest.makeNewick(6,true,false);
        assert(_debug_sfbefore == sfafter);
        unsigned g = 0;
        for (auto gf : _gene_forests) {
            string gfafter = gf.makeNewick(6,true,false);
            assert(gfafter == _debug_gfbefore[g++]);
        }
    }
#endif

    inline void Particle::refreshHeightsInternalsPreorders() {
        _species_forest.heightsInternalsPreorders();
        for (auto & gf : _gene_forests) {
            gf.heightsInternalsPreorders();
        }
    }
    
    inline SpeciesForest & Particle::getSpeciesForest() {
        return _species_forest;
    }
    
    inline const SpeciesForest & Particle::getSpeciesForest() const {
        return _species_forest;
    }
    
    inline GeneForest & Particle::getGeneForest(unsigned gene) {
        assert(gene < SMCGlobal::_ngenes);
        return _gene_forests[gene];
    }

    inline vector<GeneForest> & Particle::getGeneForests() {
        return _gene_forests;
    }
    
    inline const vector<GeneForest> & Particle::getGeneForests() const {
        return _gene_forests;
    }
    
    inline const GeneForest & Particle::getGeneForest(unsigned gene) const {
        assert(gene < SMCGlobal::_ngenes);
        return _gene_forests[gene];
    }

    inline unsigned Particle::debugCountNumCoalEvents() const {
        // Returns number of coalescent events over all gene trees
        unsigned total = 0;
        for (auto & gf : _gene_forests) {
            unsigned n = gf.getNumLineages();
            total += SMCGlobal::_ntaxa - n;
        }
        return total;
    }
            
    inline void Particle::debugCheckPartials() const {
        for (auto & gf : _gene_forests) {
            gf.debugCheckPartials();
        }
    }
    
    inline void Particle::debugShowAllGeneForests() const {
        for (auto & gf : _gene_forests) {
            cerr << gf.makeNewick(0, true, false) << endl;
        }
    }
    
    inline void Particle::copyParticleFrom(const Particle & other) {
        clear();

        // Performs a deep copy of other to this particle
        _data = other._data;
        
        // Copy gene forests
        assert(SMCGlobal::_ngenes == other._gene_forests.size());
        _gene_forests.resize(SMCGlobal::_ngenes);
        for (unsigned i = 0; i < SMCGlobal::_ngenes; i++)
            _gene_forests[i] = other._gene_forests[i];

        // Copy species forest
        _species_forest = other._species_forest;
        
        _nspeciations = other._nspeciations;
        _last_event = other._last_event;
        _count = other._count;
#if defined(USING_MULTITHREADING)
        _xtra = other._xtra;
        _begin_index = other._begin_index;
#endif
        
        // No need to copy _log_weight
        // No need to copy _lot
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
}
