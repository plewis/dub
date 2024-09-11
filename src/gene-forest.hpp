#pragma once

extern proj::Lot::SharedPtr rng;
extern proj::PartialStore ps;

namespace proj {
    
    class GeneForest : public Forest {
    
        friend class Particle;
        
        public:
        
            GeneForest();
            ~GeneForest();
            
            unsigned getGeneLength() const;
                                
            void setParticle(Particle * p);
            void setData(Data::SharedPtr d);
            void setRelRate(double r);
            pair<bool,double> advanceGeneForest(unsigned step,
                                                unsigned particle,
                                                bool simulating);
            void simulateData(Lot::SharedPtr lot, Data::SharedPtr data,
                                unsigned starting_site,
                                unsigned nsites);
            
            unsigned getGeneIndex() const;
            void setGeneIndex(unsigned i);
            
            void simulateGeneTree(unsigned gene);
            
            void buildCoalInfoVect();
            
#if defined(UPGMA_WEIGHTS)
            void debugShowDistanceMatrix(const vector<double> & d) const;
#   if defined(UPGMA_CONSTRAINED)
            void constructUPGMA(const Forest & species_forest);
#   else
            void constructUPGMA();
#   endif
            void destroyUPGMA();
#endif

            void resetPrevLogLikelihood() {_prev_log_likelihood = _log_likelihood;}
            double getPrevLogLikelihood() const {return _prev_log_likelihood;}
            void computeAllPartials();
#if defined(UPGMA_CONSTRAINED)
            double calcLogLikelihood(const Forest & species_forest);
#else
            double calcLogLikelihood();
#endif
            static void computeLeafPartials(unsigned gene, Data::SharedPtr data);
            static void releaseLeafPartials(unsigned gene);
            
            void setPriorPost(bool use_prior_post);
            
            typedef tuple<double, unsigned, G::species_t, unsigned, unsigned> coal_tuple_t;

            string lineagesWithinSpeciesKeyError(G::species_t spp);

            double calcTotalRate(vector<Node::species_tuple_t> & species_tuples);

            void mergeSpecies(double height, G::species_t left_species, G::species_t right_species, G::species_t anc_species);
            
            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void createTrivialForest(bool compute_partials = true);
            bool isSpeciesForest() const {return false;}
            void setSpeciesFromNodeName(Node * nd);
            
            void debugCheckPartials(bool verbose = false) const;
            static void clearLeafPartials();
            
            void addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient);
            
            void saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap = false) const;
            void recordHeights(vector<double> & height_vect) const;
            
            //static pair<double,double> calcTreeDistances(GeneForest & ref, GeneForest & test);

            double calcPartialArray(Node * new_nd);
            
            void operator=(const GeneForest & other);
            
            typedef shared_ptr<GeneForest> SharedPtr;
            
        protected:
                  
            double calcTransitionProbability(unsigned from, unsigned to, double edge_length);
            void debugComputeLeafPartials(unsigned gene, int number, PartialStore::partial_t partial);
            PartialStore::partial_t pullPartial();
            void stowPartial(Node * nd);
            void stowAllPartials();
            void buildLineagesWithinSpeciesMap();
            //double computeCoalRatesForSpecies(vector<G::species_t> & species, vector<double> & rates);
            
            static PartialStore::leaf_partials_t _leaf_partials;
            
            // NOTE: any variables added must be copied in operator=
            
            // key is species index, value is vector of Node pointers
            map<G::species_t, Node::ptr_vect_t > _lineages_within_species;

            Particle * _particle;
            Data::SharedPtr _data;
            double _relrate;
            unsigned _gene_index;
            bool _prior_post;

#if defined(UPGMA_WEIGHTS)
            stack<Node *> _upgma_additions;
            map<Node *, double> _upgma_starting_edgelen;
#endif
    };
    
    inline GeneForest::GeneForest() {
        _relrate = 1.0;
        clear();
    }

    inline GeneForest::~GeneForest() {
        clear();
    }
    
    inline void GeneForest::clear() {
        _particle = nullptr;
        
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

    inline void GeneForest::buildCoalInfoVect() {
        // Assumes heights of all nodes are accurate
        //
        _coalinfo.clear();
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != G::_infinity);
                    assert(nd->_left_child->_right_sib);
                    assert(nd->_left_child->_right_sib->_right_sib == nullptr);
                    nd->_species = (nd->_left_child->_species | nd->_left_child->_right_sib->_species);
                    addCoalInfoElem(nd, _coalinfo);
                }
                else {
                    // nd is a leaf node
                    unsigned spp_index = G::_taxon_to_species[nd->_name];
                    nd->_species = (G::species_t)1 << spp_index;
                }
            }
        }
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
                    cerr << str(format("  nd->_partial.get()) = %s") % G::memoryAddressAsString(nd->_partial.get())) << endl;
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
                unsigned nstates = G::_nstates;
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
    
    inline void GeneForest::setParticle(Particle * p) {
        _particle = p;
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
        assert(i < G::_nloci);
        _gene_index = i;
    }

    inline unsigned GeneForest::getGeneLength() const {
        assert(_data);
        assert(_data->_partition);
        return _data->_partition->numSitesInSubset(_gene_index);
    }

    //inline double GeneForest::calcTotalRate(vector<Node::species_tuple_t> & species_tuples, double speciation_increment) {
    inline double GeneForest::calcTotalRate(vector<Node::species_tuple_t> & species_tuples) {
        double total_rate = 0.0;
        
        // Build _lineages_within_species, a map that provides a vector
        // of Node pointers for each species
        buildLineagesWithinSpeciesMap();

        for (auto & kvpair : _lineages_within_species) {
            // kvpair.first is the species
            // kvpair.second is a vector of lineages within that species
            // Get number of lineages in this species
            unsigned n = (unsigned)kvpair.second.size();
            if (n > 1) {
                total_rate += 1.0*n*(n-1)/G::_theta;
                Node::species_tuple_t x = make_tuple(n, _gene_index, kvpair.first, kvpair.second);
                species_tuples.push_back(x);
            }
        }
        
        return total_rate;
    }
    
    inline string GeneForest::lineagesWithinSpeciesKeyError(G::species_t spp) {
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
        
    inline void GeneForest::mergeSpecies(double height, G::species_t left_species, G::species_t right_species, G::species_t anc_species) {
        // Every node previously assigned to left_species
        // or right_species should be reassigned to anc_species
        // above the specified height. Assumes nd->_height is correct.

        // Create a functor that assigns anc_species to the
        // supplied nd if it is currently in either left_species
        // or right_species
        auto reassign = [height, left_species, right_species, anc_species, this](Node * nd) {
            double h = nd->_height;
            double l = nd->_edge_length;
            G::species_t ndspp = nd->getSpecies();
            if (h + l > height && (ndspp == left_species || ndspp == right_species)) {
                nd->setSpecies(anc_species);
            }
        };
        
        // Apply functor reassign to each node in _lineages
        //for_each(_lineages.begin(), _lineages.end(), reassign);
        for (auto preorder : _preorders) {
            for_each(preorder.begin(), preorder.end(), reassign);
        }
    }
    
    inline void GeneForest::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        
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
    
    // This is an override of an abstract base class function
    inline void GeneForest::setSpeciesFromNodeName(Node * nd) {
        if (G::_taxon_to_species.count(nd->_name) == 0)
            throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % nd->_name));
        else {
            Node::setSpeciesBit(nd->_species, G::_taxon_to_species.at(nd->_name), /*init_to_zero_first*/true);
        }
    }
            
    // This is an override of the abstract base class function
    inline void GeneForest::createTrivialForest(bool compute_partials) {
        assert(G::_ntaxa > 0);
        assert(G::_ntaxa == G::_taxon_names.size());
        clear();
        unsigned nnodes = 2*G::_ntaxa - 1;
        _nodes.resize(nnodes);
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            string taxon_name = G::_taxon_names[i];
            _nodes[i]._number = (int)i;
            _nodes[i]._my_index = (int)i;
            _nodes[i]._name = taxon_name;
            _nodes[i].setEdgeLength(0.0);
            _nodes[i]._height = 0.0;
            if (G::_taxon_to_species.count(taxon_name) == 0)
                throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % taxon_name));
            else {
                Node::setSpeciesBit(_nodes[i]._species, G::_taxon_to_species.at(taxon_name), /*init_to_zero_first*/true);
            }
            if (compute_partials) {
                assert(_data);
                assert(_leaf_partials[_gene_index].size() > i);
                _nodes[i]._partial = _leaf_partials[_gene_index][i];
            }
            _lineages.push_back(&_nodes[i]);
        }
        
        // Add all remaining nodes to _unused_nodes vector
        _unused_nodes.clear();
        for (unsigned i = G::_ntaxa; i < nnodes; i++) {
            _nodes[i]._my_index = (int)i;
            _nodes[i]._number = -1;
            _unused_nodes.push_back(i);
        }
        
        refreshAllPreorders();
        _forest_height = 0.0;
        //_next_node_index = G::_ntaxa;
        //_next_node_number = G::_ntaxa;
        _log_likelihood = 0.0;
        _prev_log_likelihood = 0.0;
        _log_prior = 0.0;
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
            for (unsigned s = 0; s < G::_nstates; s++) {
                Data::state_t state = (Data::state_t)1 << s;
                Data::state_t d = data_matrix[number][pp];
                double result = state & d;
                partial->_v[p*G::_nstates + s] = (result == 0.0 ? 0.0 : 1.0);
            }
            // Ensure that the _nstates partials add up to at least 1
            assert(accumulate(partial->_v.begin() + p*G::_nstates, partial->_v.begin() + p*G::_nstates + G::_nstates, 0.0) >= 1.0);
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

                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbability(s, s_child, child->_edge_length);
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
    
#if defined(UPGMA_WEIGHTS)
struct negLogLikeDist {
    negLogLikeDist(unsigned npatterns, unsigned first, const Data::pattern_counts_t & counts, const vector<double> & same, const vector<double> & diff, double v0)
        : _npatterns(npatterns), _first(first), _counts(counts), _same(same), _diff(diff), _v0(v0) {}
    
    double operator()(double const & v) {
        double edgelen = v + _v0;
        double tprob_same = 0.25 + 0.75*exp(-4.0*edgelen/3.0);
        double tprob_diff = 0.25 - 0.25*exp(-4.0*edgelen/3.0);

        double log_like = 0.0;
        for (unsigned p = 0; p < _npatterns; p++) {
            double site_like = 0.25 * (tprob_same * _same[p] + tprob_diff * _diff[p]);
            log_like += log(site_like) * _counts[_first + p];
        }
        
        return -log_like;
    }
    
    private:
        unsigned _npatterns;
        unsigned _first;
        const Data::pattern_counts_t & _counts;
        const vector<double> & _same;
        const vector<double> & _diff;
        double _v0;
};
#endif

#if defined(UPGMA_WEIGHTS)
    inline void GeneForest::debugShowDistanceMatrix(const vector<double> & d) const {
        // d is a 1-dimensional vector that stores the lower triangle of a square matrix
        // (not including diagonals) in row order
        //
        // For example, for a 4x4 matrix (- means non-applicable):
        //
        //       0  1  2  3
        //     +-----------
        //  0  | -  -  -  -
        //  1  | 0  -  -  -
        //  2  | 1  2  -  -
        //  3  | 3  4  5  -
        //
        // For this example, d = {0, 1, 2, 3, 4 ,5}
        //
        // See this explanation for how to index d:
        //   https://math.stackexchange.com/questions/646117/how-to-find-a-function-mapping-matrix-indices
        //
        // In short, d[k] is the (i,j)th element, where k = i(i-1)/2 + j
        //       i   j   k = i*(i-1)/2 + j
        //       1   0   0 = 1*0/2 + 0
        //       2   0   1 = 2*1/2 + 0
        //       2   1   2 = 2*1/2 + 1
        //       3   0   3 = 3*2/2 + 0
        //       3   1   4 = 3*2/2 + 1
        //       3   2   5 = 3*2/2 + 2
        //
        // Number of elements in d is n(n-1)/2
        // Solving for n, and letting x = d.size(),
        //  x = n(n-1)/2
        //  2x = n^2 - n
        //  0 = a n^2 + b n + c, where a = 1, b = -1, c = -2x
        //  n = (-b += sqrt(b^2 - 4ac))/(2a)
        //    = (1 + sqrt(1 + 8x))/2
        double x = (double)d.size();
        double dbln = (1.0 + sqrt(1.0 + 8.0*x))/2.0;
        unsigned n = (unsigned)dbln;
        
        output(format("\nDistance matrix (%d x %d):\n") % n % n, 0);

        // Column headers
        output(format("%12d") % " ", 0);
        for (unsigned j = 0; j < n; j++) {
            output(format("%12d") % j, 0);
        }
        output("\n", 0);
        
        unsigned k = 0;
        for (unsigned i = 0; i < n; i++) {
            output(format("%12d") % i, 0);
            for (unsigned j = 0; j < n; j++) {
                if (j < i) {
                    double v = d[k++];
                    if (v == G::_infinity)
                        output("         inf", 0);
                    else
                        output(format("%12.5f") % v, 0);
                }
                else {
                    output("         inf", 0);
                }
            }
            output("\n", 0);
        }
        output("\n", 0);
    }
#endif
 
#if defined(UPGMA_WEIGHTS)
#if defined(UPGMA_CONSTRAINED)
    inline void GeneForest::constructUPGMA(const Forest & species_forest) {
#else
    inline void GeneForest::constructUPGMA() {
#endif
        // Get the name of the gene (data subset)
        string gene_name = _data->getSubsetName(_gene_index);

        // debugging output
        // output(format("\nGene forest for locus \"%s\" before UPGMA:\n%s\n") % gene_name % makeNewick(9, /*use_names*/true, /*coalunits*/false), 0);
        // output(format("  Height before UPGMA = %g\n") % _forest_height, 0);
                
        // Get the number of patterns
        unsigned npatterns = _data->getNumPatternsInSubset(_gene_index);

        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_gene_index);
        unsigned first_pattern = be.first;
        
        // Get pattern counts
        auto counts = _data->getPatternCounts();
        
        // Create vectors to store products of same-state and different-state partials
#if !defined(USE_JUKE_CANTOR_MODEL)
        throw XProj("GeneForest::constructUPGMA function assumes JC69 but USE_JUKE_CANTOR_MODEL was not #defined");
#endif
        vector<double> same_state(npatterns, 0.0);
        vector<double> diff_state(npatterns, 0.0);
        
        // Create a map relating position in dij vector to row,col in distance matrix
        map<unsigned, pair<unsigned, unsigned>> dij_row_col;
        
        // Create distance matrix dij and workspace dij2 used to build next dij
        // Both dij and dij2 are 1-dimensional vectors that store only the
        // lower diagonal of the distance matrix (excluding diagonal elements)
        unsigned n = (unsigned)_lineages.size();
        vector<double> dij(n*(n-1)/2, G::_infinity);
        vector<double> dij2;
        
        // Calculate distances between all pairs of lineages
        for (unsigned i = 1; i < n; i++) {
            for (unsigned j = 0; j < i; j++) {
                Node * lnode = _lineages[i];
                Node * rnode = _lineages[j];
                
                // Fill same_state and diff_state vectors
                same_state.assign(npatterns, 0.0);
                diff_state.assign(npatterns, 0.0);
                for (unsigned p = 0; p < npatterns; p++) {
                    for (unsigned lstate = 0; lstate < G::_nstates; lstate++) {
                        double lpartial = lnode->_partial->_v[p*G::_nstates + lstate];
                        for (unsigned rstate = 0; rstate < G::_nstates; rstate++) {
                            double rpartial = rnode->_partial->_v[p*G::_nstates + rstate];
                            if (lstate == rstate)
                                same_state[p] += lpartial*rpartial;
                            else
                                diff_state[p] += lpartial*rpartial;
                        }
                    }
                }
                
                double min_dist = 0.0;
                double max_dist = min_dist + 5.0; //TODO: replace arbitrary value 5.0
                
#if defined(UPGMA_CONSTRAINED)
                // Determine minimum distance based on species tree
                G::species_t lspp = lnode->getSpecies();
                G::species_t rspp = rnode->getSpecies();
                double min_height = 0.0;
                if (lspp != rspp) {
                    min_height = species_forest.mrcaHeight(lspp, rspp);
                }
                min_dist = 2.0*min_height;
                max_dist = min_dist + 5.0; //TODO: replace arbitrary value 5.0
#endif

                // Optimize edge length using black-box maximizer
                double v0 = lnode->getEdgeLength() + rnode->getEdgeLength();
                negLogLikeDist f(npatterns, first_pattern, counts, same_state, diff_state, v0);
                auto r = boost::math::tools::brent_find_minima(f, min_dist, max_dist, std::numeric_limits<double>::digits);
                double maximized_log_likelihood = -r.second;
                unsigned k = i*(i-1)/2 + j;
                dij[k] = r.first;
                dij_row_col[k] = make_pair(i,j);
                
                // output(format("d[%d] = %.5f (i = %d, j = %d, logL = %.5f)\n") % k % dij[k] % i % j % maximized_log_likelihood, 0);
            }
        }

        // debugShowDistanceMatrix(dij);

        // Create a map relating nodes in _lineages to rows of dij
        // Also save starting edge lengths so they can be restored in destroyUPGMA()
        map<Node *, unsigned> row;
        _upgma_starting_edgelen.clear();
        for (unsigned i = 0; i < n; i++) {
            Node * nd = _lineages[i];
            _upgma_starting_edgelen[nd] = nd->_edge_length;
            row[nd] = i;
        }

        // Build UPGMA tree on top of existing forest
        assert(_upgma_additions.empty());
        unsigned nsteps = n - 1;
        while (nsteps > 0) {
            // Find smallest entry in d
            auto it = min_element(dij.begin(), dij.end());
            unsigned offset = (unsigned)distance(dij.begin(), it);
            auto p = dij_row_col.at(offset);
            unsigned i = p.first;
            unsigned j = p.second;
            
            // Update all leading edge lengths
            double v = *it;
            for (auto nd : _lineages) {
                nd->_edge_length += 0.5*v;
            }
            
            //debugShowLineages();

            // Join lineages i and j
            Node * anc = pullNode();
            Node * lnode = _lineages[i];
            Node * rnode = _lineages[j];
            anc->_left_child = lnode;
            anc->_right_sib = nullptr;
            anc->_parent = nullptr;
            lnode->_right_sib = rnode;
            lnode->_parent = anc;
            rnode->_parent = anc;
            
            // Nodes added to _upgma_additions will be removed in destroyUPGMA()
            _upgma_additions.push(anc);
            
            // Remove lnode and rnode from _lineages and add anc at the end
            removeTwoAddOne(_lineages, lnode, rnode, anc);
            row[anc] = i;
                        
            //debugShowLineages();
            // output(format("\nJoining lineages %d and %d\n") % i % j, 0);

            anc->_partial = pullPartial();
            calcPartialArray(anc);
            
            // Update distance matrix
            for (unsigned k = 0; k < n; k++) {
                if (k != i && k != j) {
                    unsigned ik = (i > k) ? (i*(i-1)/2 + k) : (k*(k-1)/2 + i);
                    unsigned jk = (j > k) ? (j*(j-1)/2 + k) : (k*(k-1)/2 + j);
                    double a = dij[ik];
                    double b = dij[jk];
                    dij[ik] = 0.5*(a + b);
                    dij[jk] = G::_infinity;
                }
            }
            
            // Sanity check
            for (auto nd : _lineages) {
                assert(!nd->_right_sib);
                assert(!nd->_parent);
            }
            
            // Build new distance matrix
            unsigned n2 = (unsigned)_lineages.size();
            assert(n2 == n - 1);
            unsigned dim2 = n2*(n2-1)/2;
            dij2.resize(dim2);
            dij2.assign(dim2, G::_infinity);
            
            // Calculate distances between all pairs of lineages
            dij_row_col.clear();
            for (unsigned i2 = 1; i2 < n2; i2++) {
                for (unsigned j2 = 0; j2 < i2; j2++) {
                    Node * lnode2 = _lineages[i2];
                    Node * rnode2 = _lineages[j2];
                    unsigned i = row[lnode2];
                    unsigned j = row[rnode2];
                    unsigned k2 = i2*(i2-1)/2 + j2;
                    unsigned k = i*(i-1)/2 + j;
                    if (j > i) {
                        k = j*(j-1)/2 + i;
                    }
                    dij2[k2] = dij[k];
                    dij_row_col[k2] = make_pair(i2,j2);
                }
            }
            
            // debugShowDistanceMatrix(dij2);
            
            // Set up for next iteration
            dij = dij2;
            n = n2;
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                row[nd] = i;
            }
            
            --nsteps;
        }
        
        // debugging output
        // output(format("\nGene forest for locus \"%s\" after UPGMA:\n%s\n") % gene_name % makeNewick(9, /*use_names*/true, /*coalunits*/false), 0);
        // output(format("  Height after UPGMA = %g\n") % _forest_height, 0);
    }
#endif
    
#if defined(UPGMA_WEIGHTS)
    inline void GeneForest::destroyUPGMA() {
        while (!_upgma_additions.empty()) {
            Node * anc = _upgma_additions.top();
            Node * lnode = anc->_left_child;
            assert(lnode);
            Node * rnode = lnode->_right_sib;
            assert(rnode);
            assert(!rnode->_right_sib);
            lnode->_right_sib = nullptr;
            rnode->_right_sib = nullptr;
            lnode->_parent = nullptr;
            rnode->_parent = nullptr;
            addTwoRemoveOne(_lineages, lnode, rnode, anc);
            stowNode(anc);
            _upgma_additions.pop();
        }
        
        // Restore starting edge lengths
        for (auto nd : _lineages) {
            nd->_edge_length = _upgma_starting_edgelen.at(nd);
        }
        
        // output("\nIn GeneForest::destroyUPGMA:\n", 0);
        // output(format("  Height before refreshAllHeightsAndPreorders = %g\n") % _forest_height, 0);
        // refreshAllHeightsAndPreorders();
        // output(format("  newick = %s\n") % makeNewick(9, /*use_names*/true, /*coalunits*/false), 0);
        // output(format("  Height after refreshAllHeightsAndPreorders = %g\n") % _forest_height, 0);
        // output("\n", 0);
    }
#endif
    
#if defined(UPGMA_CONSTRAINED)
    inline double GeneForest::calcLogLikelihood(const Forest & species_forest) {
#else
    inline double GeneForest::calcLogLikelihood() {
#endif
        // Computes the log of the Felsenstein likelihood. Note that this function just
        // carries out the final summation. It assumes that all nodes (leaf nodes included)
        // have partial likelihood arrays already computed.
        if (!_data)
            return 0.0;
            
        
            
        // Compute log likelihood of every lineage
        double total_log_likelihood = 0.0;
        
#if defined(UPGMA_WEIGHTS)
        // Build remainder of the tree using UPGMA if forest has non-zero height
        // If height is zero, then we do not want to use UPGMA completion
        bool trivial_forest = (getHeight() == 0.0);
#   if defined(UPGMA_CONSTRAINED)
        if (!trivial_forest)
            constructUPGMA(species_forest);
#   else
        if (!trivial_forest)
            constructUPGMA();
#   endif
#endif
        
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
                for (unsigned s = 0; s < G::_nstates; s++) {
                    double child_partial = nd->_partial->_v[p*G::_nstates+s];
                    site_like += 0.25*child_partial;
                }
                
                log_like += log(site_like)*counts[pp];
            }

            total_log_likelihood += log_like;
            tmp++;
        }

        // output(format("GeneForest::calcLogLikelihood for locus \"%s\"\n") % gene_name, 0);
        // output(format("  newick = \"%s\"\n") % makeNewick(9, true, false), 0);
        // output(format("  total_log_likelihood = %.9f\n") % total_log_likelihood, 0);

#if defined(UPGMA_WEIGHTS)
        if (!trivial_forest)
            destroyUPGMA();
#endif
        
        _log_likelihood = total_log_likelihood;

        return total_log_likelihood;
    }
    
    inline void GeneForest::buildLineagesWithinSpeciesMap() {
        // Assumes every node in _lineages is correctly assigned to a species
        _lineages_within_species.clear();
        for (auto nd : _lineages) {
            // Add nd to the vector of nodes belonging to this species
            G::species_t spp = nd->getSpecies();
            _lineages_within_species[spp].push_back(nd);
        }
    }
    
    inline void GeneForest::simulateGeneTree(unsigned gene) {
        createTrivialForest(false);
        setGeneIndex(gene);
        
        unsigned nsteps = G::_ntaxa - 1;
        for (unsigned step = 0; step < nsteps; ++step) {
            //output(format("  step %d of %d") % step % nsteps, 2);
            bool coalescent_event = false;
            while (!coalescent_event) {
                auto result = advanceGeneForest(step, 0, true);
                coalescent_event = result.first;
            }
        }
    }
    
    inline void GeneForest::setPriorPost(bool use_prior_post) {
        _prior_post = use_prior_post;
    }
    
    inline void GeneForest::releaseLeafPartials(unsigned gene) {
        assert(_leaf_partials.size() == G::_nloci);
        _leaf_partials[gene].clear();
    }
    
    inline PartialStore::partial_t GeneForest::pullPartial() {
        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
        ptr = ps.getPartial(_gene_index);
        return ptr;
    }

    inline void GeneForest::computeLeafPartials(unsigned gene, Data::SharedPtr data) {
        assert(data);
        assert(_leaf_partials.size() == 0 || _leaf_partials.size() == G::_nloci);
        assert(G::_nloci > 0);
        assert(G::_ntaxa > 0);
        assert(G::_nstates == 4);
        
        // Allocate a vector of leaf partials for each gene
        _leaf_partials.resize(G::_nloci);
                
        // Get reference to raw data matrix, which is a vector of vector<state_t>
        auto data_matrix = data->getDataMatrix();

        // Create vector of leaf partials
        
        // Get number of patterns and first pattern index for gene g
        unsigned npatterns = data->getNumPatternsInSubset(gene);
        Data::begin_end_pair_t be = data->getSubsetBeginEnd(gene);
        unsigned first_pattern = be.first;
        
        _leaf_partials[gene].resize(G::_ntaxa);

        for (unsigned t = 0; t < G::_ntaxa; ++t) {
            PartialStore::partial_t partial_ptr = ps.getPartial(gene);
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

            _leaf_partials[gene][t] = partial_ptr;
        }
    }

    inline void GeneForest::stowPartial(Node * nd) {
        assert(nd);
        assert(nd->_partial);
        ps.putPartial(_gene_index, nd->_partial);

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
        
    inline void GeneForest::addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient) {
        // Assumes nd is an internal node
        assert(nd->_left_child);
        
        recipient.push_back(
            make_tuple(
                nd->_height,
                _gene_index + 1,
                vector<G::species_t>({
                    nd->_left_child->_species,
                    nd->_left_child->_right_sib->_species,
                })
            )
        );
    }

    inline void GeneForest::saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap) const {
        // Appends to coalinfo_vect; clear before calling if desired
        // GeneForest version ignores cap argument.
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this
        
        // coalinfo_t is a tuple with these elements:
        // - height of node
        // - 1-offset gene index (0 means speciation)
        // - vector of child species
        
        // Should only be called for complete gene trees
        assert(_lineages.size() == 1);

        // Copy tuples stored in _coalinfo to end of coalinfo_vect
        coalinfo_vect.insert(coalinfo_vect.end(), _coalinfo.begin(), _coalinfo.end());
    }
    
    void GeneForest::recordHeights(vector<double> & height_vect) const {
        // Appends to height_vect; clear before calling if desired
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this

        // Should only be called for complete gene trees
        assert(_lineages.size() == 1);
                        
        for (auto nd : boost::adaptors::reverse(_preorders[0])) {
            if (nd->_left_child) {
                // internal
                height_vect.push_back(nd->_height);
            }
        }
        
    }
    
}
