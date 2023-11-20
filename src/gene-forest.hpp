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
            
            unsigned getGeneIndex() const;
            void setGeneIndex(unsigned i);
            
            void simulateGeneTree(unsigned gene);
            
            double getLastLogLikelihood() const {return _prev_log_likelihood;}
            double calcLogLikelihood() const;
            static void computeLeafPartials(unsigned gene, Data::SharedPtr data);
            static void releaseLeafPartials(unsigned gene);
            
            void setPriorPost(bool use_prior_post);
            
            typedef tuple<double, unsigned, SMCGlobal::species_t, unsigned, unsigned> coal_tuple_t;

            string lineagesWithinSpeciesKeyError(SMCGlobal::species_t spp);

            double calcTotalRate(vector<SMCGlobal::species_tuple_t> & species_tuples, double speciation_increment);
            double coalescentEvent(Lot::SharedPtr lot, SMCGlobal::species_t spp, Node * anc, bool compute_partial, bool make_permanent);
            void mergeSpecies(SMCGlobal::species_t left_species, SMCGlobal::species_t right_species, SMCGlobal::species_t anc_species);
            
            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void createTrivialForest(bool compute_partials = true);
            bool isSpeciesForest() const {return false;}
            
            void debugCheckPartials(bool verbose = false) const;
            static void clearLeafPartials();
            
            void operator=(const GeneForest & other);
            
        protected:
                  
            double calcTransitionProbability(unsigned from, unsigned to, double edge_length);
            void debugComputeLeafPartials(unsigned gene, int number, PartialStore::partial_t partial);
            double calcPartialArray(Node * new_nd);
            void computeAllPartials();
            PartialStore::partial_t pullPartial();
            void stowPartial(Node * nd);
            void stowAllPartials();
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
        // Reset partial pointers for all nodes to decrement
        // their use counts. Note that any given partial may
        // still being used in some other particle, so we
        // do not want to stow the partial back to PartialStore
        // as that would invalidate the partial everywhere.
        for (auto & nd : _nodes) {
            nd._partial.reset();
        }
        
        Forest::clear();
        _lineages_within_species.clear();
    }

    inline void GeneForest::debugCheckPartials(bool verbose) const {
        double tmp = 0.0;
        if (_preorders.size() == 0) {
            throw XProj(format("GeneForest::debugCheckPartials: gene %d has no lineages!") % _gene_index);
        }
        for (auto & preorder : _preorders) {
            for (auto nd : preorder) {
                // Check whether _partial exists
                if (!nd->_partial) {
                    throw XProj(format("GeneForest::debugCheckPartials: node %d from gene %d has no partial") % nd->_number % _gene_index);
                }
                else if (verbose) {
                    cerr << str(format("nd->_number = %d") % nd->_number) << endl;
                    cerr << str(format("  nd->_name = %s") % nd->_name) << endl;
                    cerr << str(format("  is leaf? %s") % (nd->_left_child ? "no" : "yes")) << endl;
                    cerr << str(format("  nd->_partial.get()) = %s") % SMCGlobal::memoryAddressAsString(nd->_partial.get())) << endl;
                    cerr << str(format("  nd->_partial->_v[0] = %.5f") % nd->_partial->_v[0]) << endl;
                    cerr << str(format("  nd->_partial->_v[1] = %.5f") % nd->_partial->_v[1]) << endl;
                    cerr << str(format("  nd->_partial->_v[2] = %.5f") % nd->_partial->_v[2]) << endl;
                    cerr << str(format("  nd->_partial->_v[3] = %.5f") % nd->_partial->_v[3]) << endl;
                }
                
                // Check whether _partial has zeros for every pattern
                vector<double> & v = nd->_partial->_v;
                tmp = accumulate(v.begin(), v.end(), 0.0);
                if (tmp == 0.0) {
                    throw XProj(format("GeneForest::debugCheckPartials: node %d from gene %d has empty partial") % nd->_number % _gene_index);
                }
                
                // Check whether _partial has at least one pattern with all-zero partials
                unsigned nstates = SMCGlobal::_nstates;
                unsigned npatterns = (unsigned)nd->_partial->_v.size()/nstates;
                for (unsigned pat = 0; pat < npatterns; pat++) {
                    double sum_partials = 0.0;
                    for (unsigned s = 0; s < nstates; s++) {
                        double vpat = v[pat*nstates + s];
                        assert(!isnan(vpat));
                        assert(!isinf(vpat));
                        sum_partials += vpat;
                    }
                    if (sum_partials  == 0.0) {
                        throw XProj(format("GeneForest::debugCheckPartials: node %d from gene %d has partial with pattern %d all zeros") % nd->_number % _gene_index % pat);
                    }
                }
            }
        }
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
    
    inline unsigned GeneForest::getGeneIndex() const {
        return _gene_index;
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
    
    inline string GeneForest::lineagesWithinSpeciesKeyError(SMCGlobal::species_t spp) {
        string msg = str(format("GeneForest::coalescentEvent species %d not found in _lineages_within_species map for gene %d\n") % spp % _gene_index);
        if (_lineages_within_species.size() == 0) {
            msg += "There are actually no entries in _lineages_within_species!\n";
        }
        else {
            msg += str(format("Here are the %d species that are keys in the map:\n") % _lineages_within_species.size());
            for (auto kv : _lineages_within_species) {
                msg += str(format("  species %d contains %d lineages\n") % kv.first % kv.second.size());
            }
        }
        return msg;
    }

    inline double GeneForest::coalescentEvent(Lot::SharedPtr lot, SMCGlobal::species_t spp, Node * anc_node, bool compute_partial, bool make_permanent) {
        double log_weight = 0.0;
        
        // Get vector of nodes in the specified species spp
        unsigned i = 0;
        unsigned j = 0;
        Node * first_node = nullptr;
        Node * second_node = nullptr;
        try {
            auto & node_vect = _lineages_within_species.at(spp);
            unsigned n = (unsigned)node_vect.size();
            assert(n > 1);
            
            // Choose a random pair of lineages to join
            pair<unsigned,unsigned> chosen_pair = lot->nchoose2(n);
            i = chosen_pair.first;
            j = chosen_pair.second;
            
            // Join the chosen pair of lineages
            first_node  = node_vect[i];
            second_node = node_vect[j];
            joinLineagePair(anc_node, first_node, second_node);
        }
        catch (const out_of_range & oor) {
            throw XProj(lineagesWithinSpeciesKeyError(spp));
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

            log_weight = calcPartialArray(anc_node);
            assert(!isnan(log_weight));
            assert(!isinf(log_weight));
        }
        if (make_permanent) {
            // Empty anc_node's _prev_species_stack since this will be a permanent change
            anc_node->emptyPrevSpeciesStack();
        
            // Update lineage vector since this will be a permanent change
            try {
                removeTwoAddOne(_lineages_within_species.at(spp), first_node, second_node, anc_node);
            } catch(const out_of_range & oor) {
                throw XProj(lineagesWithinSpeciesKeyError(spp));
            }
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            refreshAllPreorders();
            
#if defined(DEBUGGING_SANITY_CHECK)
            debugCheckPartials();
#endif
        }
        
        return log_weight;
    }
    
    inline void GeneForest::mergeSpecies(SMCGlobal::species_t left_species, SMCGlobal::species_t right_species, SMCGlobal::species_t anc_species) {
        // Every lineage previously assigned to left_species or right_species
        // should be reassigned to anc_species

        // Create a functor that assigns anc_species to the supplied nd if it is currently
        // in either left_species or right_species
        auto reassign = [left_species, right_species, anc_species, this](Node * nd) {
            SMCGlobal::species_t ndspp = nd->getSpecies();
            if (ndspp == left_species || ndspp == right_species) {
                nd->setSpecies(anc_species);
            }
        };
        
        // Apply functor reassign to each node in _lineages
        for_each(_lineages.begin(), _lineages.end(), reassign);
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
            } catch(const out_of_range & oor) {
                throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % taxon_name));
            }
            if (compute_partials) {
                assert(_data);
                assert(_leaf_partials[_gene_index].size() > i);
                _nodes[i]._partial = _leaf_partials[_gene_index][i];
            }
            _lineages.push_back(&_nodes[i]);
        }
        refreshAllPreorders();
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
        assert(new_nd->_left_child->_partial);
        assert(new_nd->_left_child->_right_sib->_partial);
        
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
    
    inline void GeneForest::releaseLeafPartials(unsigned gene) {
        assert(_leaf_partials.size() == SMCGlobal::_ngenes);
        _leaf_partials[gene].clear();
    }
    
    inline PartialStore::partial_t GeneForest::pullPartial() {
        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
#if defined(USING_MULTITHREADING)
        {
            lock_guard<mutex> guard(SMCGlobal::_mutex);
            ptr = ps.getPartial(_gene_index);
        }
#else
        ptr = ps.getPartial(_gene_index);
#endif
        return ptr;
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

    inline void GeneForest::stowPartial(Node * nd) {
        assert(nd);
        assert(nd->_partial);
#if defined(USING_MULTITHREADING)
        {
            lock_guard<mutex> guard(SMCGlobal::_mutex);
            ps.putPartial(_gene_index, nd->_partial);
        }
#else
        ps.putPartial(_gene_index, nd->_partial);
#endif
        // Decrement shared pointer reference count
        nd->_partial.reset();
    }

    inline void GeneForest::stowAllPartials() {
        for (auto & nd : _nodes) {
            // Stow partials belonging to internal nodes
            // Partials for leaf nodes should be reset but not stowed
            // because they are copies of _leaf_partials and were not
            // obtained via PartialStore::getPartial.
            if (nd._left_child && nd._partial) {
                stowPartial(&nd);
            }
            nd._partial.reset();
        }
    }
    
    inline void GeneForest::computeAllPartials() {
        // Assumes _leaf_partials have been computed but that every node in the tree
        // has _partial equal to nullptr.
        assert(_data);
                
        if (_preorders.size() == 0) {
            refreshAllPreorders();
        }
#if defined(DEBUGGING_SANITY_CHECK)
        else {
            debugCheckAllPreorders();
        }
#endif
        
        for (auto & preorder : _preorders) {
            // Visit nodes in the subtree rooted at preorder in post-order sequence
            for (auto nd : boost::adaptors::reverse(preorder)) {
                assert(nd->_partial == nullptr);
                if (nd->_left_child) {
                    assert(nd->_left_child->_partial != nullptr);
                    assert(nd->_left_child->_right_sib->_partial != nullptr);
                    nd->_partial = pullPartial();
                    calcPartialArray(nd);
                }
                else {
                    nd->_partial = _leaf_partials[_gene_index][nd->_number];
                }
            }
        }
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
                
                // Just copy the shared pointer
                _nodes[i]._partial.reset();
                _nodes[i]._partial = other._nodes[i]._partial;
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
                
                // Just copy the shared pointer
                _nodes[i]._partial = other._nodes[i]._partial;
            }
        }
        
        // No need to copy _lineages_within_species because it is
        // only used in GeneForest::advanceGeneForest and is rebuilt
        // every time it is used
    }
        
}
