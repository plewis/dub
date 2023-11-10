#pragma once

extern proj::Lot rng;
extern proj::PartialStore ps;

namespace proj {

    class Particle;
    
    class GeneForest : public Forest {
    
        friend class Particle;
        
        public:
        
            GeneForest();
            ~GeneForest();
                                
            void setData(Data::SharedPtr d);
            void setRelRate(double r);
            pair<bool,double> advanceGeneForest(unsigned step,
                                                unsigned particle,
                                                bool simulating);
            void simulateData(Data::SharedPtr data,
                                unsigned starting_site,
                                unsigned nsites);
            void setGeneIndex(unsigned i);
            void simulateGeneTree(unsigned gene);
            
            double getLastLogLikelihood() const {return _prev_log_likelihood;}
            double calcLogLikelihood() const;
            static void computeLeafPartials(unsigned gene, Data::SharedPtr data);
            static void releaseLeafPartials(unsigned gene);
            
            void setPriorPost(bool use_prior_post);
            
            typedef tuple<double, unsigned, SMCGlobal::species_t, unsigned, unsigned> coal_tuple_t;

            //void possibleCoalescentEvents(vector<coal_tuple_t> & coals);
            //void coalescePair(SMCGlobal::species_t spp, unsigned i, unsigned j);
            double calcTotalRate(vector<SMCGlobal::species_tuple_t> & species_tuples, double speciation_increment);
            double coalescentEvent(Lot::SharedPtr lot, SMCGlobal::species_t spp, bool compute_partial);
            void mergeSpecies(SMCGlobal::species_t left_species, SMCGlobal::species_t right_species, SMCGlobal::species_t anc_species);
            
            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void createTrivialForest(bool compute_partials = true);
            bool isSpeciesForest() const {return false;}
            
            static void clearLeafPartials();
            
            void operator=(const GeneForest & other);
            
        protected:
                  
            double calcTransitionProbability(unsigned from, unsigned to, double edge_length);
            void debugComputeLeafPartials(unsigned gene, int number, PartialStore::partial_t partial);
            double calcPartialArray(Node * new_nd);
            void computeAllPartials();
            void buildLineagesWithinSpeciesMap();
            double computeCoalRatesForSpecies(vector<SMCGlobal::species_t> & species, vector<double> & rates);
            
            static PartialStore::leaf_partials_t _leaf_partials;
            
            // NOTE: any variables added must be copied in operator=
            
            // key is species index, value is vector of Node pointers
            map<SMCGlobal::species_t, Node::ptr_vect_t > _lineages_within_species;
            
            Data::SharedPtr _data;
            double _relrate;
            unsigned _gene_index;
            bool _prior_post;
    };
    
    inline GeneForest::GeneForest() {
        _relrate = 1.0;
        clear();
    }

    inline GeneForest::~GeneForest() {
        clear();
    }
    
    inline void GeneForest::clear() {
        //_data.reset();
        //_relrate = 1.0;
        
        // Return partials to PartialStore
        for (auto & nd : _nodes) {
            if (nd._partial) {
                nd._partial.reset();
            }
        }
        
        Forest::clear();
        _lineages_within_species.clear();
    }

    inline void GeneForest::clearLeafPartials() {
        for (unsigned g = 0; g < _leaf_partials.size(); ++g) {
            _leaf_partials[g].clear();
        }
        _leaf_partials.clear();
    }
    
    inline void GeneForest::setData(Data::SharedPtr d) {
        _data = d;
    }
    
    inline void GeneForest::setRelRate(double r) {
        _relrate = r;
    }
    
    inline void GeneForest::setGeneIndex(unsigned i) {
        assert(i < SMCGlobal::_ngenes);
        _gene_index = i;
    }
    
    inline double GeneForest::calcTotalRate(vector<SMCGlobal::species_tuple_t> & species_tuples, double speciation_increment) {
        double total_rate = 0.0;
        
        // Build _lineages_within_species, a map that provides a vector
        // of Node pointers for each species
        buildLineagesWithinSpeciesMap();

        for (auto & kvpair : _lineages_within_species) {
            // Get number of lineages in this species
            double n = (double)kvpair.second.size();
            if (n > 1) {
                total_rate += n*(n-1)/SMCGlobal::_theta;
                species_tuples.push_back(make_tuple(n, _gene_index, kvpair.first));
            }
        }
        
        return total_rate;
    }

#if 0
//    inline void GeneForest::possibleCoalescentEvents(vector<coal_tuple_t> & coals) {
//        for (auto & kv : _lineages_within_species) {
//            // Get species
//            SMCGlobal::species_t spp = kv.first;
//
//            // Get vector of nodes in species spp
//            auto & node_vect = _lineages_within_species.at(spp);
//
//            // Continue if this species has at least 2 nodes
//            unsigned n = (unsigned)node_vect.size();
//            if (n > 1) {
//                auto partial = ps.getPartial(_gene_index);
//                for (unsigned i = 0; i < n - 1; ++i) {
//                    for (unsigned j = i + 1; j < n; ++j) {
//                        // Join the i,j pair of lineages
//                        Node * first_node  = node_vect[i];
//                        Node * second_node = node_vect[j];
//                        Node * anc_node    = joinLineagePair(first_node, second_node);
//                        anc_node->setSpecies(spp);
//
//                        // Compute partial likelihood array of ancestral node
//                        assert(_data);
//                        assert(anc_node->_left_child);
//                        assert(anc_node->_left_child->_right_sib);
//                        assert(anc_node->_partial == nullptr);
//                        anc_node->_partial = partial;
//                        double log_weight = calcPartialArray(anc_node);
//                        assert(!isnan(log_weight));
//                        assert(!isinf(log_weight));
//
//                        // vector<coal_tuple_t> & coals
//                        // typedef tuple<
//                        //      double,
//                        //      GeneForest &,
//                        //      SMCGlobal::species_t,
//                        //      unsigned,
//                        //      unsigned
//                        // > coal_tuple_t;
//                        coals.push_back(make_tuple(log_weight, _gene_index, spp, i, j));
//
//                        // Unjoin lineage pair
//                        unjoinLineagePair(anc_node, first_node, second_node);
//                        anc_node->_partial.reset();
//                    } // j loop
//                } // i loop
//                ps.putPartial(_gene_index, partial);
//                partial.reset();
//            }   // more than 1 lineage in species
//        } // loop across species
//    }
#endif

#if 0
//    inline void GeneForest::coalescePair(SMCGlobal::species_t spp, unsigned i, unsigned j) {
//        auto & node_vect = _lineages_within_species.at(spp);
//        Node * first_node  = node_vect[i];
//        Node * second_node = node_vect[j];
//        Node * anc_node    = joinLineagePair(first_node, second_node);
//        anc_node->setSpecies(spp);
//
//        // Compute partial likelihood array of ancestral node
//        assert(_data);
//        assert(anc_node->_left_child);
//        assert(anc_node->_left_child->_right_sib);
//        assert(anc_node->_partial == nullptr);
//        anc_node->_partial = ps.getPartial(_gene_index);
//        calcPartialArray(anc_node);
//        
//        // Update lineage vector
//        removeTwoAddOne(_lineages_within_species.at(spp), first_node, second_node, anc_node);
//        removeTwoAddOne(_lineages, first_node, second_node, anc_node);
//    }
#endif

    inline double GeneForest::coalescentEvent(Lot::SharedPtr lot, SMCGlobal::species_t spp, bool compute_partial) {
        double log_weight = 0.0;
        
        // Get vector of nodes in the specified species spp
        auto & node_vect = _lineages_within_species.at(spp);
        unsigned n = (unsigned)node_vect.size();
        assert(n > 1);
        
        if (compute_partial) {
#if defined(PRIOR_POST)
//            // Prior-post
//            //unsigned npairs = n*(n-1)/2;
//            vector< pair<unsigned, unsigned> > ijpair;
//            vector<double> log_weights;
//            auto partial = ps.getPartial(_gene_index);
//            for (unsigned i = 0; i < n - 1; ++i) {
//                for (unsigned j = i + 1; j < n; ++j) {
//                    // Join the i,j pair of lineages
//                    Node * first_node  = node_vect[i];
//                    Node * second_node = node_vect[j];
//                    Node * anc_node    = joinLineagePair(first_node, second_node);
//                    anc_node->setSpecies(spp);
//
//                    // Compute partial likelihood array of ancestral node
//                    assert(_data);
//                    assert(anc_node->_left_child);
//                    assert(anc_node->_left_child->_right_sib);
//                    assert(anc_node->_partial == nullptr);
//                    anc_node->_partial = partial;
//                    log_weight = calcPartialArray(anc_node);
//                    assert(!isnan(log_weight));
//                    assert(!isinf(log_weight));
//
//                    log_weights.push_back(log_weight);
//                    ijpair.push_back(make_pair(i,j));
//
//                    // Unjoin lineage pair
//                    unjoinLineagePair(anc_node, first_node, second_node);
//                    anc_node->_partial.reset();
//                }
//            }
//
//            // Compute sum of weights on log scale
//            double log_sum_weights = SMCGlobal::calcLogSum(log_weights);
//
//            // Normalize log weights to create a discrete probability distribution
//            vector<double> probs(log_weights.size());
//            transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
//
//            assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
//
//            // Choose one pair to join
//            unsigned which_pair = SMCGlobal::multinomialDraw(probs);
//            pair<unsigned, unsigned> chosen_pair = ijpair[which_pair];
//
//            Node * first_node  = node_vect[chosen_pair.first];
//            Node * second_node = node_vect[chosen_pair.second];
//            Node * anc_node    = joinLineagePair(first_node, second_node);
//            anc_node->setSpecies(spp);
//
//            // Compute partial likelihood array of ancestral node
//            assert(_data);
//            assert(anc_node->_left_child);
//            assert(anc_node->_left_child->_right_sib);
//            assert(anc_node->_partial == nullptr);
//            anc_node->_partial = partial;
//            log_weight = calcPartialArray(anc_node);
//            assert(!isnan(log_weight));
//            assert(!isinf(log_weight));
//
//            // Update lineage vector
//            removeTwoAddOne(_lineages_within_species.at(spp), first_node, second_node, anc_node);
//            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
//
//            // The log_weight equals log_sum_weights for prior-post
//            log_weight = log_sum_weights;
#else // not PRIOR_POST
            // Choose a random pair of lineages to join
            pair<unsigned,unsigned> chosen_pair = lot->nchoose2(n);
            unsigned i = chosen_pair.first;
            unsigned j = chosen_pair.second;
            
            // Join the chosen pair of lineages
            Node * first_node  = node_vect[i];
            Node * second_node = node_vect[j];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            anc_node->setSpecies(spp);
            
            // Compute partial likelihood array of ancestral node
            assert(_data);
            assert(anc_node->_left_child);
            assert(anc_node->_left_child->_right_sib);
            assert(anc_node->_partial == nullptr);

            // Grab one partial from partial storage
#if defined(USING_MULTITHREADING)
            {
                lock_guard<mutex> guard(SMCGlobal::_mutex);
                anc_node->_partial = ps.getPartial(_gene_index);
            }
#else
            anc_node->_partial = ps.getPartial(_gene_index);
#endif

            log_weight = calcPartialArray(anc_node);
            assert(!isnan(log_weight));
            assert(!isinf(log_weight));

            // Update lineage vector
            removeTwoAddOne(_lineages_within_species.at(spp), first_node, second_node, anc_node);
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
#endif
        }
        else {
            // Simulating gene tree
            // Choose a random pair of lineages to join
            pair<unsigned,unsigned> chosen_pair = rng.nchoose2(n);
            unsigned i = chosen_pair.first;
            unsigned j = chosen_pair.second;
            
            // Join the chosen pair of lineages
            Node * first_node  = node_vect[i];
            Node * second_node = node_vect[j];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            anc_node->setSpecies(spp);
            
            // Update lineage vector
            removeTwoAddOne(_lineages_within_species.at(spp), first_node, second_node, anc_node);
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
        }

// 0 below is because this section will need to be reworked if LOW_MEM is defined
#if 0 && defined(LOW_MEM)
//        // Stow partials from the two children of anc_node (if they are internals)
//        // as they will never be needed again. Don't stow leaf partials because
//        // those have copies in _leaf_partials
//        if (first_node->_left_child) {
//            // Internal node partials should point to just one object
//            assert(first_node->_partial.use_count() == 1);
//            ps.putPartial(_gene_index, first_node->_partial);
//        }
//        first_node->_partial.reset();
//        if (first_node->_left_child) {
//            // Internal node partials should point to just one object
//            assert(second_node->_partial.use_count() == 1);
//            ps.putPartial(_gene_index, second_node->_partial);
//        }
//        second_node->_partial.reset();
#endif
        
        return log_weight;
    }
    
    inline void GeneForest::mergeSpecies(SMCGlobal::species_t left_species, SMCGlobal::species_t right_species, SMCGlobal::species_t anc_species) {
        // Every lineage previously assigned to left_species or right_species
        // should be reassigned to anc_species
                    
        // Create a functor that assigns anc_species to the supplied nd if it is currently
        // in either left_species or right_species
        auto reassign = [&left_species, &right_species, &anc_species, this](Node * nd) {
            SMCGlobal::species_t & ndspp = nd->getSpecies();
            if (ndspp == left_species || ndspp == right_species) {
                nd->setSpecies(anc_species);
            }
        };
        
        // Apply functor reassign to each node in _lineages
        for_each(_lineages.begin(), _lineages.end(), reassign);
    }

    inline pair<bool,double> GeneForest::advanceGeneForest(unsigned step, unsigned particle, bool simulating) {
        return make_pair(true, 0.1);
#if 0
//        // Returns true if a coalescent event was realized, false if chosen coalescence was
//        // deep or if all species have only one lineage.
//
//        // If coalescence_occurred, log_sum_weights is the log weight returned
//        double log_sum_weights = 0.0;
//
//        // Assumes that _epochs is a vector of pointers to SpeciationEpoch objects, sorted
//        // from smallest to largest height.
//        assert(checkEpochs(_epochs, Epoch::init_epoch | Epoch::coalescent_epoch | Epoch::speciation_epoch));
//
//        // Create vector of lineage counts within species to use in
//        // creating CoalescentEpoch entries
//        Epoch::lineage_counts_t species_lineage_counts;
//
//        // Build _lineages_within_species, a map that provides a vector
//        // of Node pointers for each species
//        buildLineagesWithinSpeciesMap();
//
//        debugShowEpochs(_epochs);
//
//        // Get iterator (sit) to next speciation epoch
//        double h = _forest_height;
//        assert(h != Forest::_infinity);
//        auto sit = find_if(_epochs.begin(), _epochs.end(), [h](const Epoch & e){
//            return (e.isSpeciationEpoch() && e._valid && e._height > h);
//        });
//        bool speciation_found = (bool)(sit != _epochs.end());
//        double time_to_next_speciation = speciation_found ? (sit->_height - _forest_height) : Forest::_infinity;
//
//        // Should be at least one species with lineages left to coalesce; if not,
//        // it means we have not calculated the number of steps needed correctly
//        unsigned sz = (unsigned)_lineages_within_species.size();
//        assert(sz > 0);
//
//        // Compute a vector of coalescence rates for each species
//        // The sum of these is the total coalescence rate used to determine
//        // the time to the next coalescence event.
//        vector<Node::species_t> species(sz);
//        vector<double> rates(sz);
//        double total_rate = computeCoalRatesForSpecies(species, rates);
//
//        // A negative value for coalecence time t means all species are
//        // represented by just one lineage and thus coalescence is not
//        // possible within this epoch.
//        double t = -1.0;
//
//        // If total_rate is not zero, then not all species have been reduced
//        // to a single lineage and coalescence is possible (but may be deep)
//        if (total_rate > 0.0) {
//            t = -log(1 - rng.uniform())/total_rate;
//        }
//
//        // Four possibilities:
//        //  case    t  speciation_found   description
//        //  ----  ---  ----------------   ---------------------------------------------------------
//        //    1    -1       true          No coalescence (each species just 1 lineage)
//        //    2    >0       true          a. Coalescence if t < time_to_next_speciation
//        //                                b. No coalescence if t > time_to_next_speciation
//        //    3    -1       false         Should never happen (multiple species in final ancestor)
//        //    4    >0       false         Coalescence in the ancestral population
//        assert(t > -1.0 || speciation_found);   // check for case 3
//
//        // coalescence_occurred is true if case 2a or 4, false if case 1 or 2b
//        bool coalescence_occurred = t > 0.0 && t < time_to_next_speciation;
//        if (coalescence_occurred) {
//            // Extend all lineages by t, updating _forest_height accordingly
//            advanceAllLineagesBy(t);
//
//            // Normalize the rates to create a vector of probabilities
//            vector<double> probs(rates.size(), 0.0);
//            transform(rates.begin(), rates.end(), probs.begin(), [total_rate](double r){return r/total_rate;});
//
//            // Choose in which species the next coalescence event occurs
//            unsigned k = Forest::multinomialDraw(probs);
//            Node::species_t s = species[k];
//
//            // Get number of lineages in the chosen species
//            unsigned n = (unsigned)_lineages_within_species[s].size();
//
//            Node * first_node  = nullptr;
//            Node * second_node = nullptr;
//            Node * anc_node    = nullptr;
//            pair<unsigned,unsigned> chosen_pair;
//
//            if (!SMCGlobal::_prior_post) {
//                // Perform prior-prior (randomly join two lineages in species s)
//                // Choose a random pair of lineages to join
//                chosen_pair = rng.nchoose2(n);
//                unsigned i = chosen_pair.first;
//                unsigned j = chosen_pair.second;
//
//                // Join the chosen pair of lineages
//                auto & node_vect = _lineages_within_species.at(s);
//                first_node  = node_vect[i];
//                second_node = node_vect[j];
//                anc_node    = joinLineagePair(first_node, second_node);
//                anc_node->setSpecies(s);
//
//                // Compute partial likelihood array of ancestral node
//                assert(anc_node->_left_child);
//                assert(anc_node->_left_child->_right_sib);
//                assert(anc_node->_partial == nullptr);
//                anc_node->_partial = ps.getPartial(_gene_index);
//                if (!simulating) {
//                    assert(_data);
//                    calcPartialArray(anc_node);
//                }
//
//                double log_likelihood = 0.0;
//                double log_weight = 0.0;
//                if (!simulating) {
//                    // Compute log Felsenstein likelihood
//                    log_likelihood = calcLogLikelihood();
//                    log_weight = log_likelihood - _prev_log_likelihood;
//                }
//
//                // Function returns log_sum_weights; name is holdover from when
//                // prior-post was only thing done
//                log_sum_weights = log_weight;
//
//                // Update previous log-likelihood and log-coalescent-likelihood values
//                _prev_log_likelihood = log_likelihood;
//            }
//            else {
//                // Perform prior-post (multinomial draw from all n-choose-2 possible joins in species s)
//                unsigned npairs = n*(n-1)/2;
//                vector<double> log_weights;
//                vector<pair<unsigned,unsigned> > pairs;
//                vector<PartialStore::partial_t> partials;
//                vector<double> log_likelihoods;
//
//                for (unsigned i = 0; i < n - 1; i++) {
//                    for (unsigned j = i + 1; j < n; j++) {
//                        // Get pointers to nodes i and j within species s
//                        auto & node_vect = _lineages_within_species.at(s);
//                        first_node  = node_vect[i];
//                        second_node = node_vect[j];
//                        anc_node    = joinLineagePair(first_node, second_node);
//                        anc_node->setSpecies(s);
//
//                        // Temporarily remove first_node and second_node from _lineages, adding anc_node
//                        removeTwoAddOne(_lineages, first_node, second_node, anc_node);
//
//                        // Update map relating species (keys) to counts of lineages (values)
//                        updateLineageCountsVector(species_lineage_counts);
//
//                        // Create a coalescence epoch for this join
//                        Epoch cepoch(Epoch::coalescent_epoch, _forest_height);
//                        cepoch._gene = _gene_index;
//                        cepoch._lineage_counts = species_lineage_counts;
//                        cepoch._species = anc_node->getSpecies();
//                        auto eit = speciation_found ? insertEpochBefore(_epochs, cepoch, sit) : pushBackEpoch(_epochs, cepoch);
//
//                        // Compute partial likelihood array
//                        assert(anc_node->_partial == nullptr);
//                        anc_node->_partial = ps.getPartial(_gene_index);
//                        assert(anc_node->_left_child);
//                        assert(anc_node->_left_child->_right_sib);
//                        if (!simulating) {
//                            assert(_data);
//                            calcPartialArray(anc_node);
//                        }
//
//                        // Compute log_weight corresponding to this proposed join
//                        double log_likelihood = 0.0;
//                        double log_coalescent_likelihood = 0.0;
//                        double log_weight = 0.0;
//                        if (!simulating) {
//                            // Compute log Felsenstein likelihood
//                            log_likelihood = calcLogLikelihood();
//                            log_weight = log_likelihood - _prev_log_likelihood;
//
//                            // Compute log coalescent likelihood
//                            log_coalescent_likelihood = calcLogCoalescentLikelihood(_epochs, _gene_index);
//                        }
//
//                        // Store the pair in pairs vector
//                        pairs.push_back(make_pair(i,j));
//
//                        // Store the ancestor's partial in the partials vector
//                        partials.push_back(anc_node->_partial);
//                        anc_node->_partial.reset();
//
//                        // Store the log_weight in the log_weights vector
//                        log_weights.push_back(log_weight);
//
//                        // Store the log_likelihood in the log_likelihoods vector
//                        log_likelihoods.push_back(log_likelihood);
//
//                        // Unjoin lineage pair
//                        unjoinLineagePair(anc_node, first_node, second_node);
//
//                        // Add first_node and second_node back into _lineages, removing anc_node
//                        addTwoRemoveOne(_lineages, first_node, second_node, anc_node);
//
//                        // Remove the coalescent epoch added for this pair
//                        removeEpochAt(_epochs, eit);
//                    }
//                }
//
//                // Compute sum of weights on log scale (this is *the* log_weight for prior-post)
//                log_sum_weights = calcLogSum(log_weights);
//
//                // Normalize log weights to create a discrete probability distribution
//                probs.resize(log_weights.size());
//                transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
//
//                // Choose one pair to join
//                unsigned which_pair = Forest::multinomialDraw(probs);
//                chosen_pair = pairs[which_pair];
//
//                // Update previous log-likelihood
//                _prev_log_likelihood = log_likelihoods[which_pair];
//
//                auto & node_vect = _lineages_within_species.at(s);
//                first_node  = node_vect[chosen_pair.first];
//                second_node = node_vect[chosen_pair.second];
//                anc_node    = joinLineagePair(first_node, second_node);
//
//                // Ancestral node is in the same species as the nodes that were joined
//                anc_node->setSpecies(s);
//
//                // Set ancestral node partial likelihood array to the one already calculated
//                if (!simulating) {
//                    assert(_data);
//                    assert(anc_node->_partial == nullptr);
//                    assert(anc_node->_left_child);
//                    assert(anc_node->_left_child->_right_sib);
//                    anc_node->_partial = partials[which_pair];
//                    partials[which_pair].reset();
//                }
//
//                // Clear the temporary partials vector
//                partials.clear();
//            }
//
//            // Clear partials for first_node and second_node as they will no longer be needed
//            first_node->_partial.reset();
//            second_node->_partial.reset();
//
//            // Update lineage vector
//            removeTwoAddOne(_lineages_within_species.at(s), first_node, second_node, anc_node);
//            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
//
//            // Clear partial for final lineage if we are down to just one
//            if (_lineages.size() == 1) {
//                anc_node->_partial.reset();
//            }
//
//            // Update species_lineage_counts to reflect updated _lineages_within_species
//            updateLineageCountsVector(species_lineage_counts);
//
//            // Create entry in _epochs for this coalescent event
//            Epoch cepoch(Epoch::coalescent_epoch, _forest_height);
//            cepoch._gene = _gene_index;
//            cepoch._lineage_counts = species_lineage_counts;
//            cepoch._species = anc_node->getSpecies();
//            if (speciation_found)
//                insertEpochBefore(_epochs, cepoch, sit);
//            else
//                pushBackEpoch(_epochs, cepoch);
//        }
//        else {
//            // No gene tree coalescence before the end of the species tree epoch
//            // Extend all lineages by time_to_next_speciation
//            assert(speciation_found);
//            assert(time_to_next_speciation != Forest::_infinity);
//            advanceAllLineagesBy(time_to_next_speciation);
//
//            // Ensure that this speciation epoch will not be used again until it is reset
//            // by a call to Epoch::resetAllEpochs. This is necessary because roundoff error
//            // can sometimes result in an already-used speciation epoch being reused unless
//            // it is flagged as having already been consumed.
//            sit->_valid = false;
//
//            // Every lineage previously assigned to left_species or right_species
//            // should be reassigned to anc_species
//
//            // Begin by getting the left, right, and anc species from the iterator (sit)
//            // sitting on the next speciation epoch
//            Node::species_t left_species  = sit->_left_species;
//            Node::species_t right_species = sit->_right_species;
//            Node::species_t anc_species   = sit->_anc_species;
//
//            // Create a functor that assigns anc_species to the supplied nd if it is currently
//            // in either left_species or right_species
//            auto reassign = [&left_species, &right_species, &anc_species, this](Node * nd) {
//                Node::species_t & ndspp = nd->getSpecies();
//                if (ndspp == left_species || ndspp == right_species) {
//                    nd->setSpecies(anc_species);
//                }
//            };
//
//            // Apply functor reassign to each node in _lineages
//            for_each(_lineages.begin(), _lineages.end(), reassign);
//        }
//        return make_pair(coalescence_occurred, log_sum_weights);
#endif
    }
    
    inline void GeneForest::simulateData(Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        
#if !defined(USE_JUKE_CANTOR_MODEL)
#   error Must define USE_JUKE_CANTOR_MODEL as that is the only model currently implemented
#endif

        // Create vector of states for each node in the tree
        unsigned nnodes = (unsigned)_nodes.size();
        vector< vector<unsigned> > sequences(nnodes);
        for (unsigned i = 0; i < nnodes; i++) {
            sequences[i].resize(nsites, 4);
        }
        
        // Walk through tree in preorder sequence, simulating all sites as we go
        //    DNA   state      state
        //         (binary)  (decimal)
        //    A      0001        1
        //    C      0010        2
        //    G      0100        4
        //    T      1000        8
        //    ?      1111       15
        //    R      0101        5
        //    Y      1010       10
        
        // Simulate starting sequence at the root node
        Node * nd = *(_lineages.begin());
        unsigned ndnum = nd->_number;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = rng.randint(0,3);
        }
        
        nd = findNextPreorder(nd);
        while (nd) {
            ndnum = nd->_number;
            assert(ndnum < nnodes);

            // Get reference to parent sequence
            assert(nd->_parent);
            unsigned parnum = nd->_parent->_number;
            assert(parnum < nnodes);
            
            // Evolve nd's sequence given parent's sequence and edge length
            for (unsigned i = 0; i < nsites; i++) {
                unsigned from_state = sequences[parnum][i];
                double cum_prob = 0.0;
                double u = rng.uniform();
                for (unsigned to_state = 0; to_state < 4; to_state++) {
                    cum_prob += calcTransitionProbability(from_state, to_state, nd->_edge_length);
                    if (u < cum_prob) {
                        sequences[ndnum][i] = to_state;
                        break;
                    }
                }
                assert(sequences[ndnum][i] < 4);
            }
            
            // Move to next node in preorder sequence
            nd = findNextPreorder(nd);
        }

        assert(data);
        Data::data_matrix_t & dm = data->getDataMatrixNonConst();
        
        // Copy sequences to data object
        for (unsigned t = 0; t < SMCGlobal::_ntaxa; t++) {
            // Allocate row t of _data's _data_matrix data member
            dm[t].resize(starting_site + nsites);
            
            // Get reference to nd's sequence
            unsigned ndnum = _nodes[t]._number;
            
            // Translate to state codes and copy
            for (unsigned i = 0; i < nsites; i++) {
                dm[t][starting_site + i] = (Data::state_t)1 << sequences[ndnum][i];
            }
        }
    }
    
    // This is an override of the abstract base class function
    inline void GeneForest::createTrivialForest(bool compute_partials) {
        assert(SMCGlobal::_ntaxa > 0);
        assert(SMCGlobal::_ntaxa == SMCGlobal::_taxon_names.size());
        clear();
        _nodes.resize(2*SMCGlobal::_ntaxa - 1);
        for (unsigned i = 0; i < SMCGlobal::_ntaxa; i++) {
            string taxon_name = SMCGlobal::_taxon_names[i];
            _nodes[i]._number = i;
            _nodes[i]._name = taxon_name;
            _nodes[i].setEdgeLength(0.0);
            _nodes[i]._height = 0.0;
            try {
                Node::setSpeciesBit(_nodes[i]._species, SMCGlobal::_taxon_to_species.at(taxon_name), /*init_to_zero_first*/true);
            } catch(out_of_range) {
                throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % taxon_name));
            }
            if (compute_partials) {
                assert(_data);
                assert(_leaf_partials[_gene_index].size() > i);
                _nodes[i]._partial = _leaf_partials[_gene_index][i];
            }
            _lineages.push_back(&_nodes[i]);
        }
        _forest_height = 0.0;
        _next_node_index = SMCGlobal::_ntaxa;
        _next_node_number = SMCGlobal::_ntaxa;
        _prev_log_likelihood = 0.0;
    }
    
    inline void GeneForest::debugComputeLeafPartials(unsigned gene, int number, PartialStore::partial_t partial) {
        assert(_data);

        unsigned npatterns = _data->getNumPatternsInSubset(gene);
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(gene);
        unsigned first_pattern = be.first;

        // Set each partial according to the observed data for this leaf
        auto data_matrix=_data->getDataMatrix();
        for (unsigned p = 0; p < npatterns; p++) {
            unsigned pp = first_pattern + p;
            for (unsigned s = 0; s < SMCGlobal::_nstates; s++) {
                Data::state_t state = (Data::state_t)1 << s;
                Data::state_t d = data_matrix[number][pp];
                double result = state & d;
                partial->_v[p*SMCGlobal::_nstates + s] = (result == 0.0 ? 0.0 : 1.0);
            }
            // Ensure that the _nstates partials add up to at least 1
            assert(accumulate(partial->_v.begin() + p*SMCGlobal::_nstates, partial->_v.begin() + p*SMCGlobal::_nstates + SMCGlobal::_nstates, 0.0) >= 1.0);
        }
    }
    
    inline double GeneForest::calcTransitionProbability(unsigned from, unsigned to, double edge_length) {
        assert(_relrate > 0.0);
        double transition_prob = 0.0;
#if defined(USE_JUKE_CANTOR_MODEL)
        if (from == to) {
            transition_prob = 0.25 + 0.75*exp(-4.0*_relrate*edge_length/3.0);
        }
        
        else {
            transition_prob = 0.25 - 0.25*exp(-4.0*_relrate*edge_length/3.0);
        }
#else   // HKY with hard-coded empirical frequencies specific to *beast tutorial example!
        double pi[] = {0.25, 0.25, 0.25, 0.25};
#       error need to create _index if this code is used
        assert(_index >= 0 && _index <= 3);
        if (_index == 1) {
            pi[0] = 0.40783;
            pi[1] = 0.20913;
            pi[2] = 0.19049;
            pi[3] = 0.19255;
        }
        else if (_index == 2) {
            pi[0] = 0.24551;
            pi[1] = 0.21386;
            pi[2] = 0.23698;
            pi[3] = 0.30364;
        }
        else if (_index == 3) {
            pi[0] = 0.20185;
            pi[1] = 0.21339;
            pi[2] = 0.21208;
            pi[3] = 0.37268;
        }
        double Pi[] = {pi[0] + pi[2], pi[1] + pi[3], pi[0] + pi[2], pi[1] + pi[3]};
        bool is_transition = (from == 0 && to == 2) || (from == 1 && to == 3) || (from == 2 && to == 0) || (from == 3 && to == 1);
        bool is_same = (from == 0 && to == 0) || (from == 1 && to == 1) | (from == 2 && to == 2) | (from == 3 && to == 3);
        bool is_transversion = !(is_same || is_transition);

        // JC expected number of substitutions per site
        // v = betat*(AC + AG + AG + CA + CG + CT + GA + GC + GT + TA + TC + TG)
        //   = betat*12*(1/16)
        //   = (3/4)*betat
        // betat = (4/3)*v

        // HKY expected number of substitutions per site
        //  v = betat*(AC + AT + CA + CG + GC + GT + TA + TG) + kappa*betat*(AG + CT + GA + TC)
        //    = 2*betat*(AC + AT + CG + GT + kappa(AG + CT))
        //    = 2*betat*((A + G)*(C + T) + kappa(AG + CT))
        //  betat = v/[2*( (A + G)(C + T) + kappa*(AG + CT) )]
        double kappa = 4.0349882; // (614.0*4.6444 + 601.0*4.0624 + 819.0*3.558)/(614.0 + 601.0 + 819.0);
        //double kappa = 1.0;
        double betat = 0.5*_relrate*edge_length/((pi[0] + pi[2])*(pi[1] + pi[3]) + kappa*(pi[0]*pi[2] + pi[1]*pi[3]));
        
        if (is_transition) {
            double pi_j = pi[to];
            double Pi_j = Pi[to];
            transition_prob = pi_j*(1.0 + (1.0 - Pi_j)*exp(-betat)/Pi_j - exp(-betat*(kappa*Pi_j + 1.0 - Pi_j))/Pi_j);
        }
        else if (is_transversion) {
            double pi_j = pi[to];
            transition_prob = pi_j*(1.0 - exp(-betat));
        }
        else {
            double pi_j = pi[to];
            double Pi_j = Pi[to];
            transition_prob = pi_j*(1.0 + (1.0 - Pi_j)*exp(-betat)/Pi_j) + (Pi_j - pi_j)*exp(-betat*(kappa*Pi_j + 1.0 - Pi_j))/Pi_j;
        }
#endif
        return transition_prob;
    }

    inline double GeneForest::calcPartialArray(Node * new_nd) {
        // Computes the partial array for new_nd and returns the difference in
        // log likelihood due to the addition of new_nd
        //char base[] = {'A','C','G','T'};
        
        // Get pattern counts
        auto counts = _data->getPatternCounts();

        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_gene_index);
        unsigned first_pattern = be.first;
                
        // Ensure tree is dichotomous
        assert(new_nd->_left_child);
        assert(new_nd->_left_child->_right_sib);
        assert(!new_nd->_left_child->_right_sib->_right_sib);
        
        auto & parent_partial_array = new_nd->_partial->_v;
        unsigned npatterns = _data->getNumPatternsInSubset(_gene_index);
        for (Node * child = new_nd->_left_child; child; child = child->_right_sib) {
            assert(child->_partial);
            auto & child_partial_array = child->_partial->_v;

            for (unsigned p = 0; p < npatterns; p++) {
                //unsigned pp = first_pattern + p;

                for (unsigned s = 0; s < SMCGlobal::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < SMCGlobal::_nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbability(s, s_child, child->_edge_length);
                        double child_partial = child_partial_array[p*SMCGlobal::_nstates + s_child];
                                                
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    
                    if (child == new_nd->_left_child)
                        parent_partial_array[p*SMCGlobal::_nstates + s] = sum_over_child_states;
                    else {
                        parent_partial_array[p*SMCGlobal::_nstates + s] *= sum_over_child_states;
                    }
                }   // parent state loop
            }   // pattern loop
        }   // child loop

        // Compute the ratio of after to before likelihoods
        //TODO: make more efficient
        double prev_loglike = 0.0;
        double curr_loglike = 0.0;
        auto & newnd_partial_array = new_nd->_partial->_v;
        auto & lchild_partial_array = new_nd->_left_child->_partial->_v;
        auto & rchild_partial_array = new_nd->_left_child->_right_sib->_partial->_v;
        for (unsigned p = 0; p < npatterns; p++) {
            unsigned pp = first_pattern + p;
            //unsigned count = counts[pp];
            double left_sitelike = 0.0;
            double right_sitelike = 0.0;
            double newnd_sitelike = 0.0;
            for (unsigned s = 0; s < SMCGlobal::_nstates; s++) {
                left_sitelike += 0.25*lchild_partial_array[p*SMCGlobal::_nstates + s];
                right_sitelike += 0.25*rchild_partial_array[p*SMCGlobal::_nstates + s];
                newnd_sitelike += 0.25*newnd_partial_array[p*SMCGlobal::_nstates + s];
            }
            prev_loglike += log(left_sitelike)*counts[pp];
            prev_loglike += log(right_sitelike)*counts[pp];
            curr_loglike += log(newnd_sitelike)*counts[pp];
        }
        
        return curr_loglike - prev_loglike;
    }
    
    inline double GeneForest::calcLogLikelihood() const {
        // Computes the log of the Felsenstein likelihood. Note that this function just
        // carries out the final summation. It assumes that all nodes (leaf nodes included)
        // have partial likelihood arrays already computed.
        if (!_data)
            return 0.0;
            
        // Compute log likelihood of every lineage
        double total_log_likelihood = 0.0;
        
        // Get the number of patterns
        unsigned npatterns = _data->getNumPatternsInSubset(_gene_index);
        
        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_gene_index);
        unsigned first_pattern = be.first;
        
        // Get the name of the gene (data subset)
        string gene_name = _data->getSubsetName(_gene_index);

        // Get pattern counts
        auto counts = _data->getPatternCounts();
        
        // Sum log likelihood across all lineages. If forest is a tree, there will
        // be just one lineage, which represents the root of the gene tree.
        unsigned tmp = 0;
        for (auto nd : _lineages) {
        
            double log_like = 0.0;
            for (unsigned p = 0; p < npatterns; p++) {
                unsigned pp = first_pattern + p;
                double site_like = 0.0;
                for (unsigned s = 0; s < SMCGlobal::_nstates; s++) {
                    double child_partial = nd->_partial->_v[p*SMCGlobal::_nstates+s];
                    site_like += 0.25*child_partial;
                }
                
                log_like += log(site_like)*counts[pp];
            }
            total_log_likelihood += log_like;
            tmp++;
        }

        return total_log_likelihood;
    }

    inline void GeneForest::buildLineagesWithinSpeciesMap() {
        // Assumes every node in _lineages is correctly assigned to a species
        
        _lineages_within_species.clear();
        for (auto nd : _lineages) {
            // Add nd to the vector of nodes belonging to this species
            SMCGlobal::species_t spp = nd->getSpecies();
            _lineages_within_species[spp].push_back(nd);
        }
    }
    
    inline void GeneForest::simulateGeneTree(unsigned gene) {
        createTrivialForest(false);
        setGeneIndex(gene);
        
        unsigned nsteps = SMCGlobal::_ntaxa - 1;
        for (unsigned step = 0; step < nsteps; ++step) {
            //output(format("  step %d of %d") % step % nsteps, 2);
            bool coalescent_event = false;
            while (!coalescent_event) {
                auto result = advanceGeneForest(step, 0, true);
                coalescent_event = result.first;
            }
        }
    }
    
    inline double GeneForest::computeCoalRatesForSpecies(vector<SMCGlobal::species_t> & species, vector<double> & rates) {
        unsigned i = 0;
        for (auto & x : _lineages_within_species) {
            // Key is the species s (a set containing one or more integer values)
            SMCGlobal::species_t s = x.first;
            assert(species.size() > i);
            species[i] = s;
            
            // Value is vector of node pointers belonging to species s
            unsigned n = (unsigned)x.second.size();
            assert(n > 0);
            
            // Rate of coalescence in species s is n choose 2 divided by theta/2
            // or, equivalently, n*(n-1)/theta
            double coal_rate = float(n)*(n-1)/SMCGlobal::_theta;
            rates[i] = coal_rate;
            
            ++i;
        }

        // The total coalescence rate is the sum of individual species-specific rates
        double total_rate = accumulate(rates.begin(), rates.end(), 0.0);
        return total_rate;
    }
    
    inline void GeneForest::setPriorPost(bool use_prior_post) {
        _prior_post = use_prior_post;
    }
    
    inline void GeneForest::operator=(const GeneForest & other) {
        Forest::operator=(other);
        _data = other._data;
        _gene_index = other._gene_index;
        _prior_post = other._prior_post;

        // Node _partial data members not copied in base class because partials are
        // only relevant for gene forests (and gene index must be known)
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            // Find out whether this and other have partials
            bool this_partial_exists = (_nodes[i]._partial != nullptr);
            bool other_partial_exists = (other._nodes[i]._partial != nullptr);
            
            if (this_partial_exists && other_partial_exists) {
                // Sanity check: make sure _partials are both of the correct length
                assert(other._nodes[i]._partial->_v.size() == ps.getNElements(_gene_index));
                
#if defined(LOW_MEM)
                // Copy array, not just the shared pointer
                copy(other._nodes[i]._partial->_v.begin(), other._nodes[i]._partial->_v.end(), _nodes[i]._partial->_v.begin());
#else
                // Just copy the shared pointer
                _nodes[i]._partial.reset();
                _nodes[i]._partial = other._nodes[i]._partial;
#endif
            }
            else if (this_partial_exists && !other_partial_exists) {
                // Sanity check: make sure this _partial is of the correct length
                assert(_nodes[i]._partial->_v.size() == ps.getNElements(_gene_index));
                
                // OK to set partial to null
                _nodes[i]._partial.reset();
            }
            else if (other_partial_exists && !this_partial_exists) {
                // Sanity check: make sure other _partial is of the correct length
                assert(other._nodes[i]._partial->_v.size() == ps.getNElements(_gene_index));
                
                // OK to copy
#if defined(LOW_MEM)
                // Copy array, not just the shared pointer
                _nodes[i]._partial = ps.getPartial(_gene_index);
                copy(other._nodes[i]._partial->_v.begin(), other._nodes[i]._partial->_v.end(), _nodes[i]._partial->_v.begin());
#else
                // Just copy the shared pointer
                _nodes[i]._partial = other._nodes[i]._partial;
#endif
            }
        }
        
        // No need to copy _lineages_within_species because it is
        // only used in GeneForest::advanceGeneForest and is rebuilt
        // every time it is used
    }
    
    inline void GeneForest::releaseLeafPartials(unsigned gene) {
        assert(_leaf_partials.size() == SMCGlobal::_ngenes);
        _leaf_partials[gene].clear();
    }
    
    inline void GeneForest::computeLeafPartials(unsigned gene, Data::SharedPtr data) {
        assert(data);
        assert(_leaf_partials.size() == 0 || _leaf_partials.size() == SMCGlobal::_ngenes);
        assert(SMCGlobal::_ngenes > 0);
        assert(SMCGlobal::_ntaxa > 0);
        assert(SMCGlobal::_nstates == 4);
        
        // Allocate a vector of leaf partials for each gene
        _leaf_partials.resize(SMCGlobal::_ngenes);
                
        // Get reference to raw data matrix, which is a vector of vector<state_t>
        auto data_matrix = data->getDataMatrix();

        // Create vector of leaf partials
        
        // Get number of patterns and first pattern index for gene g
        unsigned npatterns = data->getNumPatternsInSubset(gene);
        Data::begin_end_pair_t be = data->getSubsetBeginEnd(gene);
        unsigned first_pattern = be.first;
        
        _leaf_partials[gene].resize(SMCGlobal::_ntaxa);

        for (unsigned t = 0; t < SMCGlobal::_ntaxa; ++t) {
            PartialStore::partial_t partial_ptr = ps.getPartial(gene);
            vector<double> & leaf_partial = partial_ptr->_v;
            
            // Set each partial according to the observed data for leaf t
            for (unsigned p = 0; p < npatterns; p++) {
                unsigned pp = first_pattern + p;
                for (unsigned s = 0; s < SMCGlobal::_nstates; s++) {
                    Data::state_t state = (Data::state_t)1 << s;
                    Data::state_t d = data_matrix[t][pp];
                    double result = state & d;
                    leaf_partial[p*SMCGlobal::_nstates + s] = (result == 0.0 ? 0.0 : 1.0);
                }
                
                // Ensure that the _nstates partials add up to at least 1
                assert(accumulate(partial_ptr->_v.begin() + p*SMCGlobal::_nstates, partial_ptr->_v.begin() + p*SMCGlobal::_nstates + SMCGlobal::_nstates, 0.0) >= 1.0);
            }

            _leaf_partials[gene][t] = partial_ptr;
        }
    }

    inline void GeneForest::computeAllPartials() {
        // Assumes _leaf_partials have been computed but that every node in the tree
        // has _partial equal to nullptr.
        assert(_data);
                
        if (_preorders.size() == 0) {
            refreshAllPreorders();
        }
        
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                assert(nd->_partial == nullptr);
                if (nd->_left_child) {
                    nd->_partial = ps.getPartial(_gene_index);
                    calcPartialArray(nd);
                }
                else {
                    nd->_partial = _leaf_partials[_gene_index][nd->_number];
                }
            }
        }
    }
}
