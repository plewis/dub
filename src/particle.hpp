#pragma once

namespace proj {

    class Particle {
            
        public:
        
            Particle();
            Particle(const Particle & other);
            ~Particle();
            
            void clear();
            void setData(Data::SharedPtr data);
            
            void resetSpeciesForest();
            void resetGeneForests(bool compute_partials);
            
            //void incrementSpeciations() {_nspeciations++;}
            //unsigned getSpeciations() const {return _nspeciations;}
            //void clearSpeciations() {_nspeciations = 0;}
            
            void threadComputePartials(unsigned first, unsigned last);
            void computeAllPartials();
            void stowAllPartials();
            
            double calcLogLikelihood();
            
            double calcTotalCoalRate(double speciation_increment);
            void clearMarkAllForests();
            void revertToMarkAllForests();
            void advanceAllLineagesBy(double dt);
            double proposeCoalescence(unsigned seed, unsigned step, unsigned pindex, bool compute_partial, bool make_permanent);
            void finalizeProposalAllForests();
            double priorPrior(unsigned step, unsigned pindex, double total_rate, bool compute_partial);
            double priorPost(unsigned step, unsigned pindex, bool make_permanent);
            void advance(unsigned step, unsigned pindex, bool compute_partial, bool make_permanent);
            
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
            
            void setLastProposedGene(unsigned g);
            unsigned getLastProposedGene();
            
            void setLastProposedSpecies(SMCGlobal::species_t s);
            SMCGlobal::species_t getLastProposedSpecies();
            
            void setLastProposedFirstIndex(unsigned f);
            unsigned getLastProposedFirstIndex();
            
            void setLastProposedSecondIndex(unsigned s);
            unsigned getLastProposedSecondIndex();
                                
        protected:
                
            PartialStore::partial_t pullPartial(unsigned gene);
            void stowPartial(unsigned gene, Node * nd);

            Data::SharedPtr        _data;
            vector<GeneForest>     _gene_forests;
            SpeciesForest          _species_forest;
            //unsigned               _nspeciations;
            last_event_t           _last_event;
            unsigned               _count;
            
#if defined(USING_MULTITHREADING)
            unsigned               _xtra;
            unsigned               _begin_index;
#endif

            vector<Node::species_tuple_t> _species_tuples;
            
            unsigned                _last_proposed_gene;
            SMCGlobal::species_t    _last_proposed_spp;
            unsigned                _last_proposed_first_index;
            unsigned                _last_proposed_second_index;

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
        //_nspeciations = 0;
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
    
    inline double Particle::calcLogLikelihood() {
        double log_likelihood = 0.0;
        for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
            log_likelihood += _gene_forests[g].calcLogLikelihood();
        }
        return log_likelihood;
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
    
    inline void Particle::setLastProposedGene(unsigned g) {
        _last_proposed_gene = g;
    }

    inline unsigned Particle::getLastProposedGene() {
        return _last_proposed_gene;
    }
    
    inline void Particle::setLastProposedSpecies(SMCGlobal::species_t s) {
        _last_proposed_spp = s;
    }

    inline SMCGlobal::species_t Particle::getLastProposedSpecies() {
        return _last_proposed_spp;
    }
    
    inline void Particle::setLastProposedFirstIndex(unsigned f) {
        _last_proposed_first_index = f;
    }

    inline unsigned Particle::getLastProposedFirstIndex() {
        return _last_proposed_first_index;
    }
    
    inline void Particle::setLastProposedSecondIndex(unsigned s) { _last_proposed_second_index = s;
    }
    
    inline unsigned Particle::getLastProposedSecondIndex() {return _last_proposed_second_index;}

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

    inline void Particle::finalizeProposalAllForests() {
        //TODO: Particle::finalizeProposal
        _species_forest.finalizeProposal();
        for (auto & gf : _gene_forests) {
            gf.finalizeProposal();
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
    
    inline double Particle::proposeCoalescence(unsigned seed, unsigned step, unsigned pindex, bool compute_partial, bool make_permanent) {
        //TODO: Particle::proposeCoalescence
        setSeed(seed);
        
        clearMarkAllForests();
        
        advance(step, pindex, compute_partial, make_permanent);
        while (lastEventSpeciation()) {
            advance(step, pindex, compute_partial, make_permanent);
        }
        
        if (make_permanent) {
            finalizeProposalAllForests();
        }
        else
            revertToMarkAllForests();

        double log_weight = getLogWeight();
        return log_weight;
    }
    
    inline double Particle::calcTotalCoalRate(double speciation_increment) {
        double total_rate = 0.0;
        for (auto & gene_forest : _gene_forests) {
            total_rate += gene_forest.calcTotalRate(_species_tuples, speciation_increment);
        }
        return total_rate;
    }
    
    inline void Particle::clearMarkAllForests() {
        _species_forest.clearMark();
        for (auto & gf : _gene_forests) {
            gf.clearMark();
        }
    }
        
    inline void Particle::revertToMarkAllForests() {
        _species_forest.revertToMark();
        for (auto & gf : _gene_forests) {
            gf.revertToMark();
        }
    }
    
    inline void Particle::advanceAllLineagesBy(double dt) {
        _species_forest.advanceAllLineagesBy(dt);
        for (auto & gene_forest : _gene_forests) {
            gene_forest.advanceAllLineagesBy(dt);
        }
    }
    
    inline double Particle::priorPrior(unsigned step, unsigned pindex, double total_rate, bool compute_partial) {
        double log_weight = 0.0;
        
        // Prior-prior chooses a locus and species in which to have a coalescence
        // from the prior. The probability of coalescence in one particular
        // gene-species pair i is
        //   p(i) = r_i/total_rate
        // where
        //   r_i = n_i*(n_i-1)/theta
        //   total_rate = sum_i r_i
        //   n_i = number of lineages in pair i.
        vector<double> probs(_species_tuples.size());
        transform(_species_tuples.begin(), _species_tuples.end(), probs.begin(), [total_rate](Node::species_tuple_t & spp_tuple){
            double n = (double)get<0>(spp_tuple);
            return n*(n - 1.0)/(SMCGlobal::_theta*total_rate);
        });
                
#if defined(DEBUGGING_SANITY_CHECK)
        double check = accumulate(probs.begin(), probs.end(), 0.0);
        assert(fabs(check - 1.0) < 0.0001);
#endif
                                
        // Choose which gene-species pair in which the coalescence will happen
        unsigned which = SMCGlobal::multinomialDraw(_lot, probs);
        assert(which < probs.size());
        
        // Get a reference to the relevant gene forest
        unsigned g = get<1>(_species_tuples[which]);
        GeneForest & gf = _gene_forests[g];
        
        // Get the species involved
        SMCGlobal::species_t spp = get<2>(_species_tuples[which]);
            
        // Pull next available node
        Node * anc_node = gf.pullNode();
        if (compute_partial)
            anc_node->_partial = pullPartial(g);
                
        // Create variables in which to store indexes (relative
        // to spp) of nodes joined
        unsigned first_index = 0;
        unsigned second_index = 0;
                
        // coalescentEventPriorPrior joins two nodes chosen randomly
        // from species spp in gene forest g. Returns calculated log weight
        // if compute_partial is true; otherwise returns 0.0 for log weight.
        //_log_weight = gf.coalescentEventPriorPrior(_lot, spp, anc, first_index, second_index, compute_partial);

        // Get vector of nodes in the specified species spp
        unsigned i = 0;
        unsigned j = 0;
        Node * first_node = nullptr;
        Node * second_node = nullptr;
        try {
            auto & node_vect = gf._lineages_within_species.at(spp);
            unsigned n = (unsigned)node_vect.size();
            assert(n > 1);
            
            // Choose a random pair of lineages to join
            pair<unsigned,unsigned> chosen_pair = _lot->nchoose2(n);
            i = chosen_pair.first;
            j = chosen_pair.second;
            
            // Return the chosen pair in the supplied reference variables
            first_index = i;
            second_index = j;
            
            // Join the chosen pair of lineages
            first_node  = node_vect[i];
            second_node = node_vect[j];
            gf.joinLineagePair(anc_node, first_node, second_node);
        }
        catch (const out_of_range &) {
            throw XProj(gf.lineagesWithinSpeciesKeyError(spp));
        }

        anc_node->setSpecies(spp);
        
        assert(first_node->getSpecies() == spp);
        assert(second_node->getSpecies() == spp);
        
        if (compute_partial) {
            // Compute partial likelihood array of ancestral node
            assert(_data);
            assert(anc_node->_left_child);
            assert(anc_node->_left_child->_right_sib);
            assert(anc_node->_partial);

            log_weight = gf.calcPartialArray(anc_node);
            assert(!isnan(log_weight));
            assert(!isinf(log_weight));
        }
        
        gf._mark_anc_nodes.push(anc_node);
        gf._mark_left_right_pos.push(make_pair(first_index, second_index));
        
        _last_proposed_gene = g;
        _last_proposed_spp = spp;
        _last_proposed_first_index = first_index;
        _last_proposed_second_index = second_index;
        
        return log_weight;
    }

    inline double Particle::priorPost(unsigned step, unsigned pindex, bool make_permanent) {
        // Prior-post tries every possible join within each locus-species pair,
        // then chooses one of these possibilities using a multinomial draw from
        // the normalized weights. The log weight of the proposal is the log of
        // the normalizing constant for the multinomial draw.
        double log_weight = 0.0;
        
        if (make_permanent) {
            unsigned                       g = _last_proposed_gene;
            SMCGlobal::species_t         spp = _last_proposed_spp;
            unsigned             first_index = _last_proposed_first_index;
            unsigned            second_index = _last_proposed_second_index;
            
            // Get reference to correct gene forest
            GeneForest & gf = _gene_forests[g];
            
            // Pull next available node
            Node * anc_node = gf.pullNode();
            anc_node->_partial = pullPartial(g);
            
            // Join the chosen pair of lineages
            try {
                auto & node_vect = gf._lineages_within_species.at(spp);
                Node * first_node  = node_vect[first_index];
                Node * second_node = node_vect[second_index];
                gf.joinLineagePair(anc_node, first_node, second_node);
            }
            catch (const out_of_range &) {
                throw XProj(gf.lineagesWithinSpeciesKeyError(spp));
            }
            
            anc_node->setSpecies(spp);

            // Recalculate partial for anc_node
            double logw = gf.calcPartialArray(anc_node);
            assert(!isnan(logw));
            assert(!isinf(logw));

            gf._mark_anc_nodes.push(anc_node);
            gf._mark_left_right_pos.push(make_pair(first_index, second_index));
        }
        else {
            // Create a vector to hold the log weight from every possible join.
            vector<double> log_weights;
            
            // Create a map relating an index into log_weights (key) to a tuple (value)
            // containing:
            //  <0> the index into _species_tuples
            //  <1> the first index of the join, and
            //  <2> the second index of the join.
            map<unsigned, tuple<unsigned, unsigned, unsigned> > index_map;
            
            // The 4th element of _species_tuples provides pointers to the root nodes
            // of each lineage in that locus-species combination
            unsigned log_weight_index = 0;
            unsigned species_tuples_index = 0;
            Node * first_node = nullptr;
            Node * second_node = nullptr;
            for (auto t : _species_tuples) {
                // Get number of lineages
                unsigned n = get<0>(t);
                
                // Get locus
                unsigned g = get<1>(t);
                
                // Get species
                SMCGlobal::species_t spp = get<2>(t);
                
                // Get vector of nodes in the specified species spp
                Node::ptr_vect_t & node_vect = get<3>(t);
                assert(n == node_vect.size());
                assert(n > 1);
                
                // Get reference to correct gene forest
                GeneForest & gf = _gene_forests[g];
                
                // Pull next available node
                Node * anc_node = gf.pullNode();
                anc_node->_partial = pullPartial(g);

                // Visit each possible pair of nodes to join and compute
                // log weight of that join
                for (unsigned i = 0; i < n - 1; i++) {
                    for (unsigned j = i + 1; j < n; j++) {
                        index_map[log_weight_index] = make_tuple(species_tuples_index, i, j);
                        
                        // Join the chosen pair of lineages
                        first_node  = node_vect[i];
                        second_node = node_vect[j];

                        // Coalescent events should not cross species boundaries
                        assert(first_node->getSpecies() == spp);
                        assert(second_node->getSpecies() == spp);
                        
                        gf.joinLineagePair(anc_node, first_node, second_node);

                        // Compute partial likelihood array of ancestral node
                        assert(anc_node->_left_child);
                        assert(anc_node->_left_child->_right_sib);
                        assert(anc_node->_partial);

                        double logw = gf.calcPartialArray(anc_node);
                        assert(!isnan(logw));
                        assert(!isinf(logw));
                        log_weights.push_back(logw);
                        assert(log_weights.size() == log_weight_index + 1);
                        
                        // Unjoin the chosen pair of lineages
                        gf.unjoinLineagePair(anc_node, first_node, second_node);
                        
                        log_weight_index++;
                    }   // j loop
                } // i loop
                
                // Return partial
                gf.stowPartial(anc_node);

                // Return anc node
                gf.stowNode(anc_node);
                anc_node = nullptr;
                
                species_tuples_index++;
            } // species_tuples loop
            
            // Compute log of sum of weights
            double sum_log_weights = SMCGlobal::calcLogSum(log_weights);
                              
            // Create multinomial probability distribution by normalizing
            // the log weights
            vector<double> probs(log_weights.size(), 0.0);
            transform(log_weights.begin(), log_weights.end(), probs.begin(), [sum_log_weights](double logw){return exp(logw - sum_log_weights);});
            
            // Choose one pair
            unsigned which = SMCGlobal::multinomialDraw(_lot, probs);
            auto chosen = index_map[which];
            species_tuples_index  = get<0>(chosen);
            unsigned first_index  = get<1>(chosen);
            unsigned second_index = get<2>(chosen);

            unsigned             g         = get<1>(_species_tuples[species_tuples_index]);
            SMCGlobal::species_t spp       = get<2>(_species_tuples[species_tuples_index]);
            Node::ptr_vect_t &   node_vect = get<3>(_species_tuples[species_tuples_index]);

            log_weight = sum_log_weights;

            // Get reference to correct gene forest
            GeneForest & gf = _gene_forests[g];
            
            // Pull next available node
            Node * anc_node = gf.pullNode();
            anc_node->_partial = pullPartial(g);
            
            // Join the chosen pair of lineages
            first_node  = node_vect[first_index];
            second_node = node_vect[second_index];
            gf.joinLineagePair(anc_node, first_node, second_node);
            
            anc_node->setSpecies(spp);

            //TODO: necessary?
            double logw = gf.calcPartialArray(anc_node);
            assert(!isnan(logw));
            assert(!isinf(logw));

            gf._mark_anc_nodes.push(anc_node);
            gf._mark_left_right_pos.push(make_pair(first_index, second_index));

            _last_proposed_gene = g;
            _last_proposed_spp = spp;
            _last_proposed_first_index = first_index;
            _last_proposed_second_index = second_index;
        }
        return log_weight;
    }

    inline void Particle::advance(unsigned step, unsigned pindex, bool compute_partial, bool make_permanent) {
        //TODO: Particle::advance
        
        // Clear species_tuples vector. Each 4-tuple entry stores:
        //  1. number of lineages
        //  2. gene index
        //  3. species within gene
        //  4. vector<Node *> holding lineage roots for gene/species combination
        _species_tuples.clear();
        
        // Draw a speciation increment Delta. Note: speciation_increment will
        // equal "infinity" if species tree is complete.
        auto incr_rate = _species_forest.drawIncrement(_lot);
        double Delta = incr_rate.first;
        double speciation_rate = incr_rate.second;
        
        // Visit each species within each locus, storing the number of lineages in
        // _species_tuples and computing the total coalescence rate (total_rate).
        double total_rate = calcTotalCoalRate(Delta);
        
        // Draw coalescence increment delta ~ Exponential(total_rate)
        double delta = (total_rate > 0.0 ? -log(1.0 - _lot->uniform())/total_rate : SMCGlobal::_infinity);
        
        // If delta > Delta, then a speciation event happened; otherwise a coalescent event happened.
        bool is_speciation = (delta > Delta);
                
        // Advance all forests by increment dt
        double dt = is_speciation ? Delta : delta;

        // Forests save dt to allow reversion
        advanceAllLineagesBy(dt);
        
        if (is_speciation) {
            _last_event = Particle::LAST_EVENT_SPECIATION;
            
            // Create speciation event
            SMCGlobal::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(_lot, left_spp, right_spp, anc_spp);
                                    
            // Advise all gene trees of the change in the species tree
            // Nodes that are reassigned save their previous state
            // to allow reversion
            for (auto & gene_forest : _gene_forests) {
                gene_forest.mergeSpecies(left_spp, right_spp, anc_spp);
            }
        }
        else {
            _last_event = Particle::LAST_EVENT_COALESCENCE;
            
            if (SMCGlobal::_prior_post) {
                _log_weight = priorPost(step, pindex, make_permanent);
            }
            else {
                _log_weight = priorPrior(step, pindex, total_rate, compute_partial);
            }
            
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
        
        //_nspeciations = other._nspeciations;
        _last_event = other._last_event;
        _count = other._count;
#if defined(USING_MULTITHREADING)
        _xtra = other._xtra;
        _begin_index = other._begin_index;
#endif
        // no need to copy these as they are just temporary workspaces
        _last_proposed_gene = 0;
        _last_proposed_spp = 0;
        _last_proposed_first_index   = 0;
        _last_proposed_second_index  = 0;
        
        // No need to copy _log_weight
        // No need to copy _lot
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
}
