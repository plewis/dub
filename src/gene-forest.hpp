#pragma once

extern proj::Lot::SharedPtr rng;
extern proj::PartialStore ps;

namespace proj {

#if defined(LAZY_COPYING)
    class GeneForestExtension;
#endif
    
    class GeneForest : public Forest {
    
        friend class Particle;
        
#if defined(LAZY_COPYING)
        friend class GeneForestExtension;
#endif
        
        public:
        
            typedef shared_ptr<GeneForest> SharedPtr;
            typedef shared_ptr<const GeneForest> ConstSharedPtr;
            
            GeneForest();
            ~GeneForest();
            
#if defined(LAZY_COPYING)
            void addIncrAndJoin(double incr, const Split & lsplit, const Split & rsplit, GeneForestExtension & gfx);
            void copyLineageSpecies(vector<G::species_t> & species_of_lineages) const;
            void debugCheckBleedingEdge(string msg, double anc_height) const;
#endif
                        
            unsigned getGeneLength() const;
                                            
            Data::SharedPtr getData();
            void setData(Data::SharedPtr d);
            void setRelRate(double r);
            pair<bool,double> advanceGeneForest(unsigned step,
                                                unsigned particle,
                                                bool simulating);
            void simulateData(Lot::SharedPtr lot, Data::SharedPtr data,
                                unsigned starting_site,
                                unsigned nsites);
            
            double getRelRate() const;
            
            int getGeneIndex() const;
            void setGeneIndex(unsigned i);
                        
            void setLogLikelihood(double lnl);
            double getLogLikelihood() const;
            
            void simulateGeneTree(unsigned gene);
            
            void buildCoalInfoVect();
            
            double calcLogLikelihood();

            void computeAllPartials();
            static void computeLeafPartials(unsigned gene, Data::SharedPtr data);
            static void releaseLeafPartials(unsigned gene);
            
            void setPriorPost(bool use_prior_post);
            
            typedef tuple<double, unsigned, G::species_t, unsigned, unsigned> coal_tuple_t;

            string lineagesWithinSpeciesKeyError(G::species_t spp);

#if defined(LAZY_COPYING)
            void mergeSpecies(G::species_t left_species, G::species_t right_species, G::species_t anc_species);
#else
            double calcTotalRate(vector<Node::species_tuple_t> & species_tuples) const;
            void mergeSpecies(double height, G::species_t left_species, G::species_t right_species, G::species_t anc_species);
#endif
            
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

            
#if defined(LAZY_COPYING)
            double calcPartialArray(Node * new_nd, const Node * lchild, const Node * rchild) const;
#else
            double calcPartialArray(Node * new_nd);
#endif
            
            void operator=(const GeneForest & other);
            
        protected:
                  
            double calcTransitionProbability(unsigned from, unsigned to, double edge_length) const;
            double calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) const;
            void debugComputeLeafPartials(unsigned gene, int number, PartialStore::partial_t partial);
            PartialStore::partial_t pullPartial();
            void stowPartial(Node * nd);
            void stowAllPartials();
            void buildLineagesWithinSpeciesMap() const;
            //double computeCoalRatesForSpecies(vector<G::species_t> & species, vector<double> & rates);
            
            static PartialStore::leaf_partials_t _leaf_partials;
            
            // NOTE: any variables added must be copied in operator=
            
            // key is species index, value is vector of Node pointers
            mutable map<G::species_t, Node::ptr_vect_t > _lineages_within_species;

            Data::SharedPtr _data;
            int _gene_index;    // -1 (not yet set), 0, 1, ..., G::_nloci - 1
            double _relrate;
            bool _prior_post;
            
            double _log_likelihood;
    };
    
    inline GeneForest::GeneForest() :
            _relrate(1.0),
            _gene_index(-1)
    {
        clear();
    }

    inline GeneForest::~GeneForest() {
        clear();
    }
    
    inline void GeneForest::clear() {
        _log_likelihood = 0.0;
        stowAllPartials();
        Forest::clear();
        _lineages_within_species.clear();
    }
    
    inline void GeneForest::setLogLikelihood(double lnl) {
        _log_likelihood = lnl;
    }

    inline double GeneForest::getLogLikelihood() const {
        return _log_likelihood;
    }

    inline void GeneForest::buildCoalInfoVect() {
        // Assumes heights of all nodes are accurate
        
#if defined(LAZY_COPYING)
        // Assumes this is not a gene forest extension
        assert(getBoundaryExtension().first < 0);
#endif
        
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
                    unsigned spp_index = G::_taxon_to_species.at(nd->_name);
                    nd->_species = (G::species_t)1 << spp_index;
                }
            }
        }
    }

    inline void GeneForest::debugCheckPartials(bool verbose) const {
        assert(_gene_index >= 0);
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
    
    inline Data::SharedPtr GeneForest::getData() {
        return _data;
    }
    
    inline void GeneForest::setData(Data::SharedPtr d) {
        _data = d;
    }
    
    inline void GeneForest::setRelRate(double r) {
        _relrate = r;
    }
    
    inline double GeneForest::getRelRate() const {
        return _relrate;
    }
    
    inline int GeneForest::getGeneIndex() const {
        return _gene_index;
    }
    
    inline void GeneForest::setGeneIndex(unsigned i) {
        assert(i < G::_nloci);
        _gene_index = i;
    }

    inline unsigned GeneForest::getGeneLength() const {
        assert(_data);
        assert(_data->_partition);
        assert(_gene_index >= 0);
        return _data->_partition->numSitesInSubset(_gene_index);
    }

#if defined(LAZY_COPYING)
#else
    //inline double GeneForest::calcTotalRate(vector<Node::species_tuple_t> & species_tuples, double speciation_increment) {
    inline double GeneForest::calcTotalRate(vector<Node::species_tuple_t> & species_tuples) const {
        assert(_gene_index >= 0);
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

                // species_tuple_t elements:
                //  0. number of lineages (unsigned)
                //  1. gene index (unsigned)
                //  2. species within gene (G::species_t)
                //  3. lineage roots (Node::ptr_vect_t)

                Node::species_tuple_t x = make_tuple(n, _gene_index, kvpair.first, kvpair.second);
                species_tuples.push_back(x);
            }
        }
        
        return total_rate;
    }
#endif
    
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
        
#if defined(LAZY_COPYING)
    inline void GeneForest::mergeSpecies(G::species_t left_species, G::species_t right_species, G::species_t anc_species) {
        // Every node previously assigned to left_species
        // or right_species should be reassigned to anc_species
        
        // Create a functor that assigns anc_species to the
        // supplied nd if it is currently in either left_species
        // or right_species
        auto reassign = [left_species, right_species, anc_species, this](Node * nd) {
            G::species_t ndspp = nd->getSpecies();
            if (ndspp == left_species || ndspp == right_species) {
                nd->setSpecies(anc_species);
            }
        };
        
        // Apply functor reassign to each node in _lineages
        for_each(_lineages.begin(), _lineages.end(), reassign);
    }
#else
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
            if (h + l + G::_small_enough > height && (ndspp == left_species || ndspp == right_species)) {
                nd->setSpecies(anc_species);
            }
        };
        
        // Apply functor reassign to each node in _lineages
        for_each(_lineages.begin(), _lineages.end(), reassign);
        
        //for (auto preorder : _preorders) {
        //    for_each(preorder.begin(), preorder.end(), reassign);
        //}
    }
#endif
    
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
        
        // Draw equilibrium base frequencies from Dirichlet
        // having parameter G::_comphet
        vector<double> basefreq = {0.25, 0.25, 0.25, 0.25};
        if (G::_comphet != G::_infinity) {
            // Draw 4 Gamma(G::_comphet, 1) variates
            double A = lot->gamma(G::_comphet, 1.0);
            double C = lot->gamma(G::_comphet, 1.0);
            double G = lot->gamma(G::_comphet, 1.0);
            double T = lot->gamma(G::_comphet, 1.0);
            double total = A + C + G + T;
            basefreq[0] = A/total;
            basefreq[1] = C/total;
            basefreq[2] = G/total;
            basefreq[3] = T/total;
        }
        
        // Simulate starting sequence at the root node
        Node * nd = *(_lineages.begin());
        unsigned ndnum = nd->_number;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = G::multinomialDraw(lot, basefreq);
        }
        
        nd = findNextPreorder(nd);
        while (nd) {
            ndnum = nd->_number;
            assert(ndnum < nnodes);

            // Get reference to parent sequence
            assert(nd->_parent);
            unsigned parnum = nd->_parent->_number;
            assert(parnum < nnodes);
            
            // Choose relative rate for this node's edge
            // Lognormal m=mean v=variance
            //   sigma^2 = log(v/m^2 + 1)
            //   mu = log m - sigma^2/2
            // Mean m=1 because these are relative rates, so
            //   sigma^2 = log(v + 1)
            //   mu = -sigma^2/2
            double edge_relrate = 1.0;
            if (G::_edge_rate_variance > 0.0) {
                double sigma2 = log(1.0 + G::_edge_rate_variance);
                double sigma = sqrt(sigma2);
                double mu = -0.5*sigma2;
                double normal_variate = sigma*lot->normal() + mu;
                edge_relrate = exp(normal_variate);
            }
            
            // Evolve nd's sequence given parent's sequence and edge length
            for (unsigned i = 0; i < nsites; i++) {
                // Choose relative rate for this site
                double site_relrate = 1.0;
                if (G::_asrv_shape != G::_infinity)
                    site_relrate = lot->gamma(G::_asrv_shape, 1.0/G::_asrv_shape);
                unsigned from_state = sequences[parnum][i];
                double cum_prob = 0.0;
                double u = lot->uniform();
                for (unsigned to_state = 0; to_state < 4; to_state++) {
                    cum_prob += calcSimTransitionProbability(from_state, to_state, basefreq, site_relrate*edge_relrate*nd->_edge_length);
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
        assert(_gene_index >= 0);
        assert(G::_ntaxa > 0);
        assert(G::_ntaxa == G::_taxon_names.size());
        clear();
        _nleaves = G::_ntaxa;
        unsigned nnodes = 2*_nleaves - 1;
        _nodes.resize(nnodes);
        for (unsigned i = 0; i < _nleaves; i++) {
            string taxon_name = G::_taxon_names[i];
            _nodes[i]._number = (int)i;
            _nodes[i]._name = taxon_name;
            _nodes[i].setEdgeLength(0.0);
            _nodes[i]._height = 0.0;
            _nodes[i]._split.resize(_nleaves);
            _nodes[i]._split.setBitAt(i);
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
        for (unsigned i = _nleaves; i < nnodes; i++) {
            _nodes[i]._number = (int)i;
            _unused_nodes.push_back(i);
        }
        
        refreshAllPreorders();
        _forest_height = 0.0;
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
    
    inline double GeneForest::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) const {
        assert(pi.size() == 4);
        assert(fabs(accumulate(pi.begin(), pi.end(), 0.0) - 1.0) < G::_small_enough);
        assert(_relrate > 0.0);
        double transition_prob = 0.0;
        
        // F81 transition probabilities
        double Pi[] = {pi[0] + pi[2], pi[1] + pi[3], pi[0] + pi[2], pi[1] + pi[3]};
        bool is_transition = (from == 0 && to == 2) || (from == 1 && to == 3) || (from == 2 && to == 0) || (from == 3 && to == 1);
        bool is_same = (from == 0 && to == 0) || (from == 1 && to == 1) | (from == 2 && to == 2) | (from == 3 && to == 3);
        bool is_transversion = !(is_same || is_transition);

        // HKY expected number of substitutions per site
        //  v = betat*(AC + AT + CA + CG + GC + GT + TA + TG) + kappa*betat*(AG + CT + GA + TC)
        //    = 2*betat*(AC + AT + CG + GT + kappa(AG + CT))
        //    = 2*betat*((A + G)*(C + T) + kappa(AG + CT))
        //  betat = v/[2*( (A + G)(C + T) + kappa*(AG + CT) )]
        double kappa = 1.0;
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
        return transition_prob;
    }
    
    inline double GeneForest::calcTransitionProbability(unsigned from, unsigned to, double edge_length) const {
        assert(_relrate > 0.0);
        double transition_prob = 0.0;
#if defined(USE_JUKE_CANTOR_MODEL)
        if (from == to) {
            transition_prob = 0.25 + 0.75*exp(-4.0*_relrate*edge_length/3.0);
        }
        
        else {
            transition_prob = 0.25 - 0.25*exp(-4.0*_relrate*edge_length/3.0);
        }
#else
        // HKY with hard-coded empirical frequencies specific to *beast tutorial example!
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

#if defined(LAZY_COPYING)
    inline double GeneForest::calcPartialArray(
            Node * new_nd,
            const Node * lchild,
            const Node * rchild) const {
        assert(_gene_index >= 0);

        // Computes the partial array for new_nd and returns the difference in
        // log likelihood due to the addition of new_nd
        //char base[] = {'A','C','G','T'};
        
        // Get pattern counts
        auto counts = _data->getPatternCounts();

        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_gene_index);
        unsigned first_pattern = be.first;
                
        auto & parent_partial_array = new_nd->_partial->_v;
        unsigned npatterns = _data->getNumPatternsInSubset(_gene_index);
#if 1
        // Determine if there is an edge length extension (this would be the
        // case if new_nd comes from a gene forest extension)
        double lchild_stem_height = lchild->_height + lchild->_edge_length;
        double rchild_stem_height = rchild->_height + rchild->_edge_length;
        assert(fabs(lchild_stem_height - rchild_stem_height) < G::_small_enough);
        
        // Calculate the edge length extension
        double edgelen_extension = new_nd->_height - lchild_stem_height;
        
        // Edge length extension may be slightly negative due to roundoff
        assert(edgelen_extension >= -G::_small_enough);
        if (edgelen_extension < 0.0)
            edgelen_extension = 0.0;
            
        for (const Node * child : {lchild, rchild})  {
            assert(child->_partial);
            auto & child_partial_array = child->_partial->_v;
                
            double pr_same = calcTransitionProbability(0, 0, child->_edge_length + edgelen_extension);
            double pr_diff = calcTransitionProbability(0, 1, child->_edge_length + edgelen_extension);
            for (unsigned p = 0; p < npatterns; p++) {
                //unsigned pp = first_pattern + p;

                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? pr_same : pr_diff);
                        double child_partial = child_partial_array[p*G::_nstates + s_child];
                                                
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    
                    if (child == lchild)
                        parent_partial_array[p*G::_nstates + s] = sum_over_child_states;
                    else {
                        parent_partial_array[p*G::_nstates + s] *= sum_over_child_states;
                    }
                }   // parent state loop
            }   // pattern loop
        }
#else
        for (Node * child = new_nd->_left_child; child; child = child->_right_sib) {
            assert(child->_partial);
            auto & child_partial_array = child->_partial->_v;

            // If this gene forest is an extension, check to see if edge length
            // needs to be extended to account for both the edge length in the
            // parent forest as well as the delta accumulated in the extension
            double edgelen_extension = 0.0;
            bool straddler = (new_nd->_height > _starting_height) && (child->_height < _starting_height);
            if (_is_extension && straddler) {
                edgelen_extension = _proposed_delta;
            }
                
            double pr_same = calcTransitionProbability(0, 0, child->_edge_length + edgelen_extension);
            double pr_diff = calcTransitionProbability(0, 1, child->_edge_length + edgelen_extension);
            for (unsigned p = 0; p < npatterns; p++) {
                //unsigned pp = first_pattern + p;

                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? pr_same : pr_diff);
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
#endif

        // Compute the ratio of after to before likelihoods
        //TODO: make more efficient
        double prev_loglike = 0.0;
        double curr_loglike = 0.0;
        auto & newnd_partial_array = new_nd->_partial->_v;
        auto & lchild_partial_array = lchild->_partial->_v;
        auto & rchild_partial_array = rchild->_partial->_v;
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
        
        G::_npartials_calculated++;
        return curr_loglike - prev_loglike;
    }
#else
    inline double GeneForest::calcPartialArray(Node * new_nd) {
        assert(_gene_index >= 0);

        // Computes the partial array for new_nd and returns the difference in
        // log likelihood due to the addition of new_nd
        //char base[] = {'A','C','G','T'};
        
        // Get pattern counts
        auto counts = _data->getPatternCounts();

        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_gene_index);
        unsigned first_pattern = be.first;
                
        auto & parent_partial_array = new_nd->_partial->_v;
        unsigned npatterns = _data->getNumPatternsInSubset(_gene_index);
        for (Node * child = new_nd->_left_child; child; child = child->_right_sib) {
            assert(child->_partial);
            auto & child_partial_array = child->_partial->_v;
                
            double pr_same = calcTransitionProbability(0, 0, child->_edge_length);
            double pr_diff = calcTransitionProbability(0, 1, child->_edge_length);
            for (unsigned p = 0; p < npatterns; p++) {
                //unsigned pp = first_pattern + p;

                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? pr_same : pr_diff);
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
        
        G::_npartials_calculated++;
        return curr_loglike - prev_loglike;
    }
#endif
    
    inline double GeneForest::calcLogLikelihood()
    {
        assert(_gene_index >= 0);

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
                for (unsigned s = 0; s < G::_nstates; s++) {
                    double child_partial = nd->_partial->_v[p*G::_nstates+s];
                    site_like += 0.25*child_partial;
                }
                
                log_like += log(site_like)*counts[pp];
            }

            total_log_likelihood += log_like;
            tmp++;
        }

        // output(format("GeneForest::calcLogLikelihood for locus \"%s\"\n") % gene_name);
        // output(format("  newick = \"%s\"\n") % makeNewick(9, true, false));
        // output(format("  total_log_likelihood = %.9f\n") % total_log_likelihood);

        _log_likelihood = total_log_likelihood;

        return total_log_likelihood;
    }
    
    inline void GeneForest::buildLineagesWithinSpeciesMap() const {
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
            //output(format("  step %d of %d") % step % nsteps);
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
#if defined(USING_MULTITHREADING)
        lock_guard<mutex> guard(G::_mutex);
#endif
        assert(_gene_index >= 0);
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
#if defined(USING_MULTITHREADING)
        lock_guard<mutex> guard(G::_mutex);
#endif
        assert(_gene_index >= 0);

        if (nd && nd->_left_child && nd->_partial) {
            // Nothing to do if nd or nd->_partial is null
            // Only internal partials are every reset/stowed;
            // leaf partials are calculated at the beginning
            // and never changed.
            if (nd->_partial.use_count() == 1) {
                // Partial is not being used by any other node, so it is
                // safe to stow it in PartialStore for reuse later
                ps.putPartial(_gene_index, nd->_partial);
            }
        
            nd->_partial.reset();
        }
    }
    
    inline void GeneForest::stowAllPartials() {
        for (auto & nd : _nodes) {
            stowPartial(&nd);
            nd._partial.reset();
        }
    }
    
    inline void GeneForest::computeAllPartials() {
        // Assumes _leaf_partials have been computed but that every node in the tree
        // has _partial equal to nullptr.
        assert(_data);
        assert(_gene_index >= 0);
                
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
#if defined(LAZY_COPYING)
                    calcPartialArray(nd, nd->_left_child, nd->_left_child->_right_sib);
#else
                    calcPartialArray(nd);
#endif
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
        _relrate = other._relrate;
        _prior_post = other._prior_post;
        _log_likelihood = other._log_likelihood;
        
        assert(_gene_index >= 0);

        // Node _partial data members not copied in base class because partials are
        // only relevant for gene forests (and gene index must be known)
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            // Find out whether this and other have partials
            bool this_partial_exists = (_nodes[i]._partial != nullptr);
            bool other_partial_exists = (other._nodes[i]._partial != nullptr);
            
            if (this_partial_exists && other_partial_exists) {
#if defined(USING_MULTITHREADING)
                // not thread safe to access PartialStore function
#else
                // Sanity check: make sure _partials are both of the correct length
                assert(other._nodes[i]._partial->_v.size() == ps.getNElements(_gene_index));
#endif
                
                // Just copy the shared pointer
                _nodes[i]._partial.reset();
                _nodes[i]._partial = other._nodes[i]._partial;
            }
            else if (this_partial_exists && !other_partial_exists) {
#if defined(USING_MULTITHREADING)
                // not thread safe to access PartialStore function
#else
                // Sanity check: make sure this _partial is of the correct length
                assert(_nodes[i]._partial->_v.size() == ps.getNElements(_gene_index));
#endif
                
                // OK to set partial to null
                _nodes[i]._partial.reset();
            }
            else if (other_partial_exists && !this_partial_exists) {
#if defined(USING_MULTITHREADING)
                // not thread safe to access PartialStore function
#else
                // Sanity check: make sure other _partial is of the correct length
                assert(other._nodes[i]._partial->_v.size() == ps.getNElements(_gene_index));
#endif
                
                // Just copy the shared pointer
                _nodes[i]._partial = other._nodes[i]._partial;
            }
        }
        
        // No need to copy _lineages_within_species because it is
        // only used in GeneForest::advanceGeneForest and is rebuilt
        // every time it is used
    }
        
    inline void GeneForest::addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient)
     {
        assert(_gene_index >= 0);

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
    
#if defined(LAZY_COPYING)
    void GeneForest::copyLineageSpecies(vector<G::species_t> & species_of_lineages) const {
        species_of_lineages.resize(_lineages.size());
        unsigned i = 0;
        for (auto nd : _lineages) {
            species_of_lineages[i++] = nd->_species;
        }
    }
#endif

#if defined(LAZY_COPYING)
    inline void GeneForest::debugCheckBleedingEdge(string msg, double anc_height) const {
        output(format("\n|~~~> %s...\n") % msg, G::LogCateg::DEBUGGING);
        
        // Find maximum height of all nodes in _lineages vector
        double maxh = 0.0;
        for (auto nd : _lineages) {
            if (nd->_height > maxh)
                maxh = nd->_height;
        }

        // Find maximum height + edgelen of all nodes in _lineages vector
        double maxhplus = 0.0;
        for (auto nd : _lineages) {
            if (nd->_height + nd->_edge_length > maxhplus)
                maxhplus = nd->_height + nd->_edge_length;
        }
        
        output(format("\n%12.9f = maxh\n%12.9f = maxhplus\n%12.9f = forest height\n%12.9f = anc height\n\n") % maxh % maxhplus % _forest_height % anc_height, G::LogCateg::DEBUGGING);

        double diff = fabs(maxhplus - _forest_height);
        //if (diff > G::_small_enough) {
        //    cerr << endl;
        //}
        assert(diff <= G::_small_enough);
        
        double diff2 = fabs(_forest_height - anc_height);
        //if (diff2 > G::_small_enough) {
        //    cerr << endl;
        //}
        assert(diff2 <= G::_small_enough);
    }
#endif

}
