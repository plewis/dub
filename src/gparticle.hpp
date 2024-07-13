#pragma once

extern void doof(string msg, proj::G::verbosity_t verb);
extern void doof(format & fmt, proj::G::verbosity_t level);
extern proj::PartialStore         ps;
extern proj::Partition::SharedPtr partition;
extern proj::Data::SharedPtr      data;

namespace proj {

    class Bundle;

    // Gene tree particle
    class GParticle : public Particle {
    
        friend class Bundle;
        
        public:
            GParticle() : _locus_index(-1), _log_weight(G::_negative_infinity), _species_tree(nullptr) {}
            ~GParticle() {}
            
            // Getters and setters
            double getLogWeight() const {return _log_weight;}
            void setLogWeight(double logw) {_log_weight = logw;}
            
            unsigned getLocusIndex() const {return _locus_index;}
            void setLocusIndex(unsigned i) {_locus_index = i;}
                        
            SParticle * getSpeciesTree() {return _species_tree;}
            const SParticle * getSpeciesTree() const {return _species_tree;}
            void setSpeciesTree(SParticle * s) {_species_tree = s;}
            
            double calcTransitionProbability(unsigned from, unsigned to, double edge_length, double relrate);
            double calcPartialArray(Node * new_nd);
            Node * coalesce(Lot::SharedPtr lot);
            void simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites);
            
#if defined(STOW_UNUSED_PARTIALS)
            void deleteStowedPartials();
            void enumerateUnusedPartials(vector<map<PartialStore::partial_t, int> > & unused) const;
#endif

            // Overrides of base class abstract virtual functions
            pair<double, double> drawIncrement(Lot::SharedPtr lot);
            Node * joinRandomPair(Lot::SharedPtr lot);
            double calcLogLikelihood() const;
            void createTrivialForest();
            string info() const {
                //return str(format("  G-%d-%d") % getLocusIndex() % getIndex());
                return str(format("%.5f") % getLogWeight());
            }
            
            static void clearLeafPartials();
            static void computeLeafPartials();
            
            void operator=(const GParticle & other);
            
        protected:
            int          _locus_index;
            double       _log_weight;
            SParticle *  _species_tree;
            
            static PartialStore::leaf_partials_t _leaf_partials;
    };

    inline double GParticle::calcLogLikelihood() const {
        // Computes the log of the Felsenstein likelihood.
        // Note that this function just carries out the final summation.
        // It assumes that all nodes (leaf nodes included) have
        // partial likelihood arrays already computed.
        if (!::data)
            return 0.0;
            
        // Compute log likelihood of every lineage
        double total_log_likelihood = 0.0;
        
        // Get the number of patterns
        unsigned npatterns = ::data->getNumPatternsInSubset(_locus_index);
        
        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = ::data->getSubsetBeginEnd(_locus_index);
        unsigned first_pattern = be.first;
        
        // Get the name of the gene (data subset)
        string gene_name = ::data->getSubsetName(_locus_index);

        // Get pattern counts
        auto counts = ::data->getPatternCounts();
        
        // Sum log likelihood across all lineages
        unsigned tmp = 0;
        for (auto nd : _lineages) {
            double log_like = 0.0;
            for (unsigned p = 0; p < npatterns; p++) {
                unsigned pp = first_pattern + p;
                double site_like = 0.0;
                for (unsigned s = 0; s < G::_nstates; s++) {
                    double child_partial = nd->_partial->_v[p*G::_nstates+s];
                    site_like += 0.25*child_partial;
                }
                
                log_like += log(site_like)*counts[pp];
            }

            total_log_likelihood += log_like;
            tmp++;
        }

        return total_log_likelihood;
    }
    
    inline void GParticle::createTrivialForest() {
        assert(::data);
        _nodes.clear();
        _nodes.resize(2*G::_ntaxa - 1);
        _lineages.clear();
        _lineages.resize(G::_ntaxa);
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            _lineages[i] = &_nodes[i];
            string nm = G::_taxon_names[i];
            _lineages[i]->setName(nm);
            _lineages[i]->setNumber(i);
            if (!::data->empty()) {
                assert(GParticle::_leaf_partials[_locus_index].size() > i);
                _nodes[i].setPartial(GParticle::_leaf_partials[_locus_index][i]);
            }
            G::setSpeciesBit(_lineages[i]->getSpecies(), G::_taxon_to_species[nm], /*init_to_zero_first*/true);
        }
        _next_node_number = G::_ntaxa;
        _returned_node_numbers.clear();
        _nleaves = G::_ntaxa;
        _is_species_tree = false;
        refreshSplits();
    }
    
    inline pair<double, double> GParticle::drawIncrement(Lot::SharedPtr lot) {
        // Choose an increment at random from the Exp(rate) prior
        // where rate is the total rate over all species currently
        // represented among the lineages
        
        // Create map in which keys are species and values are counts of lineages
        map<G::species_t, unsigned > species_counts;
        for (auto nd : _lineages) {
            species_counts[nd->getSpecies()]++;
        }

        // Calculate total rate based on theta and numbers of lineages in each species
        output("\n", G::VDEBUG);
        double total_rate = 0.0;
        for (auto & p : species_counts) {
            output(format("~~> species %d contains %d lineages\n") % p.first % p.second, G::VDEBUG);

            // Get theta value for current species (p.first)
            const SParticle::theta_map_t & theta_map = _species_tree->getThetaMapConst();
            double theta = theta_map.at(p.first);

            // Get number of lineages in this species
            unsigned       n = p.second;
            
            // Calculate coalescent rate within this species
            double rate = 1.0*n*(n-1)/theta;
            
            total_rate += rate;
        }
        
        // Draw increment from Exponential with rate total_rate
        double incr = G::_infinity;
        if (total_rate > 0.0) {
            // At least one species has more than one lineage, so
            // coalescence is possible in at least one species
            incr = lot->gamma(1.0, 1.0/total_rate);
        }
        output(format("~~> total rate = %g, increment = %g\n") % total_rate % incr, G::VDEBUG);

        return make_pair(incr, total_rate);
    }

    inline Node * GParticle::joinRandomPair(Lot::SharedPtr lot) {
        // Create map in which keys are species and values are vectors of Node *
        //TODO: inefficient, nearly same thing done in drawIncrement
        map<G::species_t, vector<Node *> > species_nodevector;
        for (auto nd : _lineages) {
            species_nodevector[nd->getSpecies()].push_back(nd);
        }

        // Calcuate rates vector
        //TODO: inefficient, should reuse calculation in drawIncrement
        vector<double> rates;
        double total_rate = 0.0;
        for (auto & p : species_nodevector) {
            // Get theta value for current species (p.first)
            const SParticle::theta_map_t & theta_map = _species_tree->getThetaMapConst();
            double theta = theta_map.at(p.first);
            
            // Get number of lineages in this species
            unsigned n = (unsigned)p.second.size();
            
            // Calculate coalescent rate within this species
            double rate = 1.0*n*(n-1)/theta;
            
            total_rate += rate;
            rates.push_back(rate);
        }
        
        // Normalize rates vector
        transform(rates.begin(), rates.end(), rates.begin(), [total_rate](double r){return r/total_rate;});
        
        // Determine species in which coalescence occurred
        unsigned which = G::multinomialDraw(lot, rates);
        auto it = species_nodevector.begin();
        advance(it, which);
        assert(it != species_nodevector.end());
        vector<Node *> & choices = it->second;
        
        G::species_t chosen_species = it->first;
        output(format("joined two lineages in species %d\n") % chosen_species, G::VDEBUG);

        // Choose two lineages within species spp to join
        pair<unsigned, unsigned> p = lot->nchoose2((unsigned)choices.size());
        Node * lchild = choices[p.first];
        Node * rchild = choices[p.second];

        // Grab next node to use as ancestor
        unsigned node_number = _next_node_number;
        size_t n = _returned_node_numbers.size();
        if (n > 0) {
            node_number = _returned_node_numbers[n-1];
            _returned_node_numbers.pop_back();
        }
        else {
            _next_node_number++;
        }
        assert(_nodes.size() > node_number);
        Node * anc = &_nodes[node_number];
        anc->setNumber(node_number);
        
        makeAnc(anc, lchild, rchild);

        if (!::data->empty()) {
            // Calculate partial for anc and store _log_weight
            PartialStore::partial_t part = ps.pullPartial(_locus_index);
            anc->setPartial(part);
            _log_weight = calcPartialArray(anc);
        }
        
        return anc;
    }

    inline Node * GParticle::coalesce(Lot::SharedPtr lot) {
        assert(_species_tree);
        vector<G::join_info_t> spec_rec;
        _species_tree->recordJoinInfo(spec_rec);
        
        Node * anc = nullptr;
        bool done = false;
        while (!done) {
            // Determine the height of the next species tree join
            G::join_info_t ceiling = _species_tree->heightOfNextNodeAbove(_height, spec_rec);
            
            // height (+inf if species tree complete)
            double upper_bound = get<0>(ceiling);

            auto incr_rate = drawIncrement(lot);
            
            if (_height + incr_rate.first < upper_bound) {
                // Increment does not take us to the next speciation event, so
                // add a coalescence event and we're done
                output(format("DONE! height + incr = %g < %g = upper_bound\n") % (_height + incr_rate.first) % upper_bound, G::VDEBUG);
                extendAllLineagesBy(incr_rate.first);
                anc = joinRandomPair(lot);
                done = true;
            }
            else {
                // Increment takes us beyond the next speciation event, so
                // advance all lineages to upper_bound and try again
                output(format("height + incr = %g > %g = upper_bound\n") % (_height + incr_rate.first) % upper_bound, G::VDEBUG);
                extendAllLineagesBy(upper_bound - _height);
                
                G::species_t spp = get<2>(ceiling) | get<3>(ceiling);
                if (spp == G::_species_zero) {
                    // Species tree needs to be extended before continuing
                    // Note: joinThenIncrement returns ancestral node, but
                    // we are ignoring that return value here.
                    _species_tree->joinThenIncrement(lot);
                    
                    // Refresh speciation info
                    _species_tree->recordJoinInfo(spec_rec);
                    
                    // Penultimate entry in spec_rec should now have height _height
                    // and has info about which species were joined
                    unsigned n = (unsigned)spec_rec.size();
                    assert(n > 1);
                    G::join_info_t latest_join = spec_rec[n - 2];
                    assert(fabs(get<0>(latest_join) - _height) < G::_small_enough);
                    G::species_t lspp = get<2>(latest_join);
                    G::species_t rspp = get<3>(latest_join);
                    
                    // Every lineage in species lspp or rspp should now be
                    // in species lspp | rspp
                    for (auto nd : _lineages) {
                        G::species_t spp = nd->getSpecies();
                        if (spp == lspp || spp == rspp) {
                            nd->setSpecies(lspp | rspp);
                        }
                    }
                }
                else {
                    // Species tree does not need to be extended, and ceiling
                    // contains info about species that should be merged
                    G::species_t lspp = get<2>(ceiling);
                    G::species_t rspp = get<3>(ceiling);
                    
                    // Every lineage in species lspp or rspp should now be
                    // in species lspp | rspp
                    for (auto nd : _lineages) {
                        G::species_t spp = nd->getSpecies();
                        if (spp == lspp || spp == rspp) {
                            nd->setSpecies(lspp | rspp);
                        }
                    }
                }   // extension unnecessary
            }   // _height + incr_rate.first >= upper_bound
        } // while (!done)
        
        return anc;
    }

    inline void GParticle::clearLeafPartials() {
        _leaf_partials.clear();
    }
    
    inline void GParticle::computeLeafPartials() {
        assert(::data);
        assert(_leaf_partials.size() == 0);
        assert(G::_nloci > 0);
        assert(G::_ntaxa > 0);
        assert(G::_nstates == 4);
        
        for (unsigned locus = 0; locus < G::_nloci; locus++) {
                
            // Allocate a vector of leaf partials for each locus
            _leaf_partials.resize(G::_nloci);
                    
            // Get reference to raw data matrix, which is a vector of vector<state_t>
            auto data_matrix = ::data->getDataMatrix();

            // Create vector of leaf partials
            
            // Get number of patterns and first pattern index for locus
            unsigned npatterns = ::data->getNumPatternsInSubset(locus);
            Data::begin_end_pair_t be = ::data->getSubsetBeginEnd(locus);
            unsigned first_pattern = be.first;
            
            _leaf_partials[locus].resize(G::_ntaxa);

            for (unsigned t = 0; t < G::_ntaxa; ++t) {
                PartialStore::partial_t partial_ptr = ps.pullPartial(locus);
                
                vector<double> & leaf_partial = partial_ptr->_v;
                
                // Set each partial according to the observed data for leaf t
                for (unsigned p = 0; p < npatterns; p++) {
                    unsigned pp = first_pattern + p;
                    for (unsigned s = 0; s < G::_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[t][pp];
                        double result = state & d;
                        leaf_partial[p*G::_nstates + s] = (result == 0.0 ? 0.0 : 1.0);
                    }
                    
                    // Ensure that the _nstates partials add up to at least 1
                    assert(accumulate(partial_ptr->_v.begin() + p*G::_nstates, partial_ptr->_v.begin() + p*G::_nstates + G::_nstates, 0.0) >= 1.0);
                }

                _leaf_partials[locus][t] = partial_ptr;
            }   // loop over taxa
        }  // loop over loci
    }
    
    inline double GParticle::calcTransitionProbability(unsigned from, unsigned to, double edge_length, double relrate) {
        assert(relrate > 0.0);
        double transition_prob = 0.0;
#if defined(USE_JUKE_CANTOR_MODEL)
        if (from == to) {
            transition_prob = 0.25 + 0.75*exp(-4.0*relrate*edge_length/3.0);
        }
        
        else {
            transition_prob = 0.25 - 0.25*exp(-4.0*relrate*edge_length/3.0);
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

    inline double GParticle::calcPartialArray(Node * new_nd) {
        // Computes the partial array for new_nd and returns the difference in
        // log likelihood due to the addition of new_nd
        
        // Get pattern counts
        auto counts = ::data->getPatternCounts();

        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = ::data->getSubsetBeginEnd(_locus_index);
        unsigned first_pattern = be.first;
                
        // Ensure tree is dichotomous
        assert(new_nd->_left_child);
        assert(new_nd->_left_child->_right_sib);
        assert(!new_nd->_left_child->_right_sib->_right_sib);
        assert(new_nd->_left_child->_partial);
        assert(new_nd->_left_child->_right_sib->_partial);
        
        auto & parent_partial_array = new_nd->_partial->_v;
        unsigned npatterns = ::data->getNumPatternsInSubset(_locus_index);
        for (Node * child = new_nd->_left_child; child; child = child->_right_sib) {
            assert(child->_partial);
            
#if defined(STOW_UNUSED_PARTIALS)
            if (child->_partial->_in_storage) {
                string s = G::memoryAddressAsString(child->_partial.get());
                output(format("attempting to access stowed partial: %s (locus %d)\n") % s % _locus_index, G::VSTANDARD);
                assert(false);
            }
#endif
            
            auto & child_partial_array = child->_partial->_v;

            for (unsigned p = 0; p < npatterns; p++) {
                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbability(s, s_child, child->_edge_length, /*relrate*/1.0);
                        double child_partial = child_partial_array[p*G::_nstates + s_child];
                                                
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    
                    if (child == new_nd->_left_child)
                        parent_partial_array[p*G::_nstates + s] = sum_over_child_states;
                    else {
                        parent_partial_array[p*G::_nstates + s] *= sum_over_child_states;
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
            for (unsigned s = 0; s < G::_nstates; s++) {
                left_sitelike += 0.25*lchild_partial_array[p*G::_nstates + s];
                right_sitelike += 0.25*rchild_partial_array[p*G::_nstates + s];
                newnd_sitelike += 0.25*newnd_partial_array[p*G::_nstates + s];
            }
            prev_loglike += log(left_sitelike)*counts[pp];
            prev_loglike += log(right_sitelike)*counts[pp];
            curr_loglike += log(newnd_sitelike)*counts[pp];
        }
        
        return curr_loglike - prev_loglike;
    }

    inline void GParticle::operator=(const GParticle & other) {
        Particle::operator=(other);
        _locus_index = other._locus_index;
        _log_weight = other._log_weight;
        // _species_tree does not need to be copied

        // Node _partial data members not copied in base class because partials are
        // only relevant for gene trees (and locus index must be known)
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            // Find out whether this and other have partials
            bool this_partial_exists = (_nodes[i]._partial != nullptr);
            bool other_partial_exists = (other._nodes[i]._partial != nullptr);
            
            if (this_partial_exists && other_partial_exists) {
                // Sanity check: make sure _partials are both of the correct length
                assert(other._nodes[i]._partial->_v.size() == ps.getNElements(_locus_index));
                
                // Just copy the shared pointer
                _nodes[i]._partial = other._nodes[i]._partial;
            }
            else if (this_partial_exists && !other_partial_exists) {
                // Sanity check: make sure this _partial is of the correct length
                assert(_nodes[i]._partial->_v.size() == ps.getNElements(_locus_index));
                
                // OK to set partial to null
                _nodes[i]._partial.reset();
            }
            else if (other_partial_exists && !this_partial_exists) {
                // Sanity check: make sure other _partial is of the correct length
                assert(other._nodes[i]._partial->_v.size() == ps.getNElements(_locus_index));
                
                // Just copy the shared pointer
                _nodes[i]._partial = other._nodes[i]._partial;
            }
        }
    }

#if defined(STOW_UNUSED_PARTIALS)
    inline void GParticle::deleteStowedPartials() {
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        for (auto & pre : preorders) {
            for (auto nd : pre) {
                if (nd->_partial && nd->_partial->_in_storage) {
                    nd->_partial.reset();
                }
            }
        }
    }
#endif
    
#if defined(STOW_UNUSED_PARTIALS)
    inline void GParticle::enumerateUnusedPartials(vector<map<PartialStore::partial_t, int> > & unused) const {
        // Set count for any root or leaf node partials to -1.
        // Increment count for partials above the root node of each lineage
        // as long as count is not -1.
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        for (auto & pre : preorders) {
            for (auto nd : pre) {
                if (nd->_partial) {
                    bool is_root = !nd->_parent;
                    bool is_leaf = !nd->_left_child;
                        
                    if (is_root || is_leaf) {
                        unused[_locus_index][nd->_partial] = -1;
                    }
                    else {
                        // If not a root or a leaf, must be an internal
                        assert(nd->_left_child);
                        bool taboo = unused[_locus_index][nd->_partial] < 0;
                        if (!taboo) {
                            unused[_locus_index][nd->_partial]++;
                        }
                    }
                }
            }
        }
    }
#endif

    inline void GParticle::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        
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
            sequences[ndnum][i] = lot->randint(0,3);
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
                double u = lot->uniform();
                for (unsigned to_state = 0; to_state < 4; to_state++) {
                    cum_prob += calcTransitionProbability(from_state, to_state, nd->_edge_length, /*relrate*/1.0);
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
        for (unsigned t = 0; t < G::_ntaxa; t++) {
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
}

