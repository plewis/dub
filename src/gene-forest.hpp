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
            void simulateData(Data::SharedPtr data, unsigned starting_site, unsigned nsites);
            void setGeneIndex(unsigned i);
            void simulateGeneTree(unsigned gene, epoch_list_t & epochs);
            pair<bool,double> advanceGeneForest(unsigned step, unsigned particle, bool simulating = false);
            
            double calcLogLikelihood() const;
            
            static void computeLeafPartials(Data::SharedPtr data);

            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void createTrivialForest(bool compute_partials = true);
            bool isSpeciesForest() const {return false;}
            
            void operator=(const GeneForest & other);
            
        protected:
                  
            Epoch & createInitEpoch();
            epoch_list_t & digest();
            double calcTransitionProbability(unsigned from, unsigned to, double edge_length);
            void debugComputeLeafPartials(unsigned gene, int number, PartialStore::partial_t partial);
            void calcPartialArray(Node * new_nd);
            void computeAllPartials();
            void initSpeciesLineageCountsVector(Epoch::lineage_counts_t & species_lineage_counts) const;
            void buildLineagesWithinSpeciesMap();
            double computeCoalRatesForSpecies(vector<Node::species_t> & species, vector<double> & rates);
            void updateLineageCountsVector(Epoch::lineage_counts_t & slc);
            
            static PartialStore::leaf_partials_t _leaf_partials;
            
            // NOTE: any variables added must be copied in operator=
            
            // key is species index, value is vector of Node pointers
            map<Node::species_t, Node::ptr_vect_t > _lineages_within_species;
            
            Data::SharedPtr _data;
            unsigned _gene_index;
    };
    
    inline GeneForest::GeneForest() {
    }

    inline GeneForest::~GeneForest() {
        clear();
    }
    
    inline void GeneForest::clear() {
        // Return partials to PartialStore
        for (auto & nd : _nodes) {
            if (nd._partial) {
                ps.stowPartial(nd._partial, _gene_index);
                nd._partial.reset();
            }
        }
        
        Forest::clear();
        _lineages_within_species.clear();
    }

    inline void GeneForest::setData(Data::SharedPtr d) {
        _data = d;
    }
    
    inline void GeneForest::setGeneIndex(unsigned i) {
        assert(i < Forest::_ngenes);
        _gene_index = i;
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
        for (unsigned t = 0; t < _ntaxa; t++) {
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
        assert(Forest::_ntaxa > 0);
        assert(Forest::_ntaxa == Forest::_taxon_names.size());
        clear();
        _nodes.resize(2*Forest::_ntaxa - 1);
        for (unsigned i = 0; i < Forest::_ntaxa; i++) {
            string taxon_name = Forest::_taxon_names[i];
            _nodes[i]._number = i;
            _nodes[i]._name = taxon_name;
            _nodes[i].setEdgeLength(0.0);
            _nodes[i]._height = 0.0;
            _nodes[i]._species = {Forest::_taxon_to_species[taxon_name]};
            if (compute_partials) {
                assert(_data);
                assert(_leaf_partials[_gene_index].size() > i);
                _nodes[i]._partial = ps.getPartial(_gene_index);
                (*_nodes[i]._partial) = (*_leaf_partials[_gene_index][i]);
            }
            _lineages.push_back(&_nodes[i]);
        }
        _forest_height = 0.0;
        _next_node_index = Forest::_ntaxa;
        _next_node_number = Forest::_ntaxa;
        
        _prev_log_likelihood = 0.0;
        
        _prev_log_coalescent_likelihood = 0.0;
    }
    
    inline Epoch & GeneForest::createInitEpoch() {
        // Create an init entry in _epochs for this gene that just stores the starting
        // species_lineage_counts vector
        refreshAllPreorders();

        Epoch iepoch(Epoch::init_epoch, -1.0);
        iepoch._gene = _gene_index;
        initSpeciesLineageCountsVector(iepoch._lineage_counts);
        auto it = pushFrontEpoch(_epochs, iepoch);
        
        return *it;
    }
    
    inline epoch_list_t & GeneForest::digest() {
        // Assumes gene forest is a complete tree
        assert(_lineages.size() == 1);
        refreshAllPreorders();
        
        // Begin with a clean slate
        _epochs.clear();
        
        // Create an init entry in _epochs for this gene that just stores the starting
        // species_lineage_counts vector

        // Initialize vector of lineage counts for each species
        // For example, if there are 3 species, then species_lineage_counts
        // will have 3 + 2 = 5 elements and only the first 3 will have non-zero
        // counts upon initialization. E.g. [2,3,2,0,0]. This says that species
        // 0 has 2 lineages, species 1 has 3 lineages, and species 2 has 2 lineages.
        Epoch & init_epoch = createInitEpoch();
        Epoch::lineage_counts_t species_lineage_counts = init_epoch._lineage_counts;

        // Record coalescent events in the gene forest using a post-order
        // traversal, calculating heights of nodes along the way
        Node::ptr_vect_t & preorder = _preorders[0];
        for (Node * nd : boost::adaptors::reverse(preorder)) {
            Node * lchild = nd->getLeftChild();
            if (lchild) {
                // nd is an internal
                Node * rchild = lchild->getRightSib();
                
                // Assume tree has no polytomies
                assert(!rchild->getRightSib());

                // Gather info needed to compute height of nd
                double height_left  = lchild->getHeight();
                double brlen_left   = lchild->getEdgeLength();
                double height_right = rchild->getHeight();
                double brlen_right  = rchild->getEdgeLength();

                // Height calculated through left and right children
                // should be identical if this is an ultrametric time tree
                double hleft        = height_left + brlen_left;
                double hright       = height_right + brlen_right;
                assert(fabs(hleft - hright) < Forest::_small_enough);
                
                // Set height of this internal node so that later nodes can use it
                double h = 0.5*(hleft + hright);
                nd->setHeight(h);

                // Internal node's species is the union of its childrens' species
                nd->setSpeciesToUnion(lchild->getSpecies(), rchild->getSpecies());
                
                //unsigned ndnum = nd->getNumber();
                
                // Create entry in _epochs for this coalescent event
                Epoch cepoch(Epoch::coalescent_epoch, h);
                cepoch._gene = _gene_index;
                cepoch._coalescence_node = nd;
                cepoch._lineage_counts = species_lineage_counts;
                cepoch._species = nd->getSpecies();
                pushBackEpoch(_epochs, cepoch);
            }
            else {
                // nd is a leaf
                nd->setHeight(0.0);
            }
        }

        // Sort epochs by height
        _epochs.sort(epochLess);
        
        return _epochs;
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
            for (unsigned s = 0; s < Forest::_nstates; s++) {
                Data::state_t state = (Data::state_t)1 << s;
                Data::state_t d = data_matrix[number][pp];
                double result = state & d;
                (*partial)[p*Forest::_nstates + s] = (result == 0.0 ? 0.0 : 1.0);
            }
            // Ensure that the _nstates partials add up to at least 1
            assert(accumulate(partial->begin()+p*Forest::_nstates, partial->begin()+p*Forest::_nstates + Forest::_nstates, 0.0) >= 1.0);
        }
    }
    
    inline double GeneForest::calcTransitionProbability(unsigned from, unsigned to, double edge_length){
        double transition_prob = 0.0;
#if defined(USE_JUKE_CANTOR_MODEL)
        if (from == to) {
            transition_prob = 0.25 + 0.75*exp(-4.0*edge_length/3.0);
        }
        
        else {
            transition_prob = 0.25 - 0.25*exp(-4.0*edge_length/3.0);
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
        double betat = 0.5*edge_length/((pi[0] + pi[2])*(pi[1] + pi[3]) + kappa*(pi[0]*pi[2] + pi[1]*pi[3]));
        
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

    inline void GeneForest::calcPartialArray(Node * new_nd) {
        auto & parent_partial_array = *(new_nd->_partial);
        unsigned npatterns = _data->getNumPatternsInSubset(_gene_index);
        for (Node * child = new_nd->_left_child; child; child = child->_right_sib) {
            auto & child_partial_array = *(child->_partial);

            for (unsigned p = 0; p < npatterns; p++) {

                for (unsigned s = 0; s < Forest::_nstates; s++) {

                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < Forest::_nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbability(s, s_child, child->_edge_length);
                        double child_partial = child_partial_array[p*Forest::_nstates + s_child];
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    
                    if (child == new_nd->_left_child)
                        parent_partial_array[p*Forest::_nstates + s] = sum_over_child_states;
                    else
                        parent_partial_array[p*Forest::_nstates + s] *= sum_over_child_states;
                }   // parent state loop
            }   // pattern loop
        }   // child loop
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
        //unsigned last_pattern = be.second;
        
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
                for (unsigned s = 0; s < Forest::_nstates; s++) {
                    double child_partial = (*nd->_partial)[p*Forest::_nstates+s];
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
            _lineages_within_species[nd->getSpecies()].push_back(nd);
        }
    }
    
    inline void GeneForest::initSpeciesLineageCountsVector(Epoch::lineage_counts_t & species_lineage_counts) const {
        // Assumes _preorders is up-to-date.
        // Assumes forest is trivial.
        // Assumes every leaf node is assigned to the correct species
        // Resizes supplied species_lineage_counts vector to 2*nspecies - 1 elements
        // and zeros it. Each leaf node visited increases the count of that leaf's species by one.
        //unsigned sz = (unsigned)(2*Forest::_nspecies - 1);
        //species_lineage_counts.resize(sz);
        //species_lineage_counts.assign(sz, 0);
        
        // Start out each species with count of zero
        species_lineage_counts.clear();
        for (unsigned i = 0; i < Forest::_nspecies; i++) {
            Node::species_t s = {i};
            species_lineage_counts[s] = 0;
        }
        
        for (auto preorder : _preorders) {
            for (auto nd : preorder) {
                if (!nd->_left_child) {
                    if (nd->_species.size() > 1) {
                        // Leaf nodes are supposed to belong to single species
                        // so if _species_set includes more than one species,
                        // then this must be left over from when the gene forest
                        // was built conditional on a previous species tree
                        unsigned spp = Forest::_taxon_to_species[nd->_name];
                        nd->setSpeciesFromUnsigned(spp);
                    }
                    species_lineage_counts[nd->getSpecies()]++;
                }
            }
        }
    }
    
    inline void GeneForest::simulateGeneTree(unsigned gene, epoch_list_t & epochs) {
        createTrivialForest(false);
        //setData(_data);
        setGeneIndex(gene);
        copyEpochsFrom(epochs);
        createInitEpoch();
        resetAllEpochs(_epochs);
        debugShowEpochs(_epochs);
        
        unsigned nsteps = Forest::_ntaxa - 1;
        for (unsigned step = 0; step < nsteps; ++step) {
            debugShow(format("  step %d of %d") % step % nsteps);
            bool coalescent_event = false;
            while (!coalescent_event) {
                auto result = advanceGeneForest(step, 0, true);
                coalescent_event = result.first;
            }
        }
    }
    
    //inline void GeneForest::updateLineageCountsVector(Epoch::lineage_counts_t & slc, unsigned slc_size) {
    inline void GeneForest::updateLineageCountsVector(Epoch::lineage_counts_t & slc) {
        //slc.assign(slc_size, 0);
        slc.clear();
        for (auto & x : _lineages_within_species) {
            //assert(x.first.size() == 1);
            //unsigned s = *(x.first.begin());
            unsigned n = (unsigned)x.second.size();
            //slc[s] = n;
            slc[x.first] = n;
        }
    }
    
    inline double GeneForest::computeCoalRatesForSpecies(vector<Node::species_t> & species, vector<double> & rates) {
        unsigned i = 0;
        for (auto & x : _lineages_within_species) {
            // Key is the species s (a set containing one or more integer values)
            Node::species_t s = x.first;
            species[i] = s;
            
            // Value is vector of node pointers belonging to species s
            unsigned n = (unsigned)x.second.size();
            assert(n > 0);
            
            // Rate of coalescence in species s is n choose 2 divided by theta/2
            // or, equivalently, n*(n-1)/theta
            double coal_rate = float(n)*(n-1)/Forest::_theta;
            rates[i] = coal_rate;
            
            ++i;
        }

        // The total coalescence rate is the sum of individual species-specific rates
        double total_rate = accumulate(rates.begin(), rates.end(), 0.0);
        return total_rate;
    }
    
    inline pair<bool,double> GeneForest::advanceGeneForest(unsigned step, unsigned particle, bool simulating) {
        // Returns true if a coalescent event was realized, false if chosen coalescence was
        // deep or if all species have only one lineaged.
        
        // If coalescence_occurred, log_sum_weights is the log weight returned
        double log_sum_weights = 0.0;
    
        // Assumes that _epochs is a vector of pointers to SpeciationEpoch objects, sorted
        // from smallest to largest height.
        assert(checkEpochs(_epochs, Epoch::init_epoch | Epoch::coalescent_epoch | Epoch::speciation_epoch));
        
        // Create vector of lineage counts within species to use in
        // creating CoalescentEpoch entries
        Epoch::lineage_counts_t species_lineage_counts;
        
        // Build _lineages_within_species, a map that provides a vector
        // of Node pointers for each species
        buildLineagesWithinSpeciesMap();
        
        debugShowEpochs(_epochs);
        
        // Get iterator (sit) to next speciation epoch
        double h = _forest_height;
        assert(h != Forest::_infinity);
        auto sit = find_if(_epochs.begin(), _epochs.end(), [h](const Epoch & e){
            return (e.isSpeciationEpoch() && e._valid && e._height > h);
        });
        bool speciation_found = (bool)(sit != _epochs.end());
        double time_to_next_speciation = speciation_found ? (sit->_height - _forest_height) : Forest::_infinity;
        
        // Should be at least one species with lineages left to coalesce; if not,
        // it means we have not calculated the number of steps needed correctly
        unsigned sz = (unsigned)_lineages_within_species.size();
        assert(sz > 0);
        
        // Compute a vector of coalescence rates for each species
        // The sum of these is the total coalescence rate used to determine
        // the time to the next coalescence event.
        vector<Node::species_t> species(sz);
        vector<double> rates(sz);
        double total_rate = computeCoalRatesForSpecies(species, rates);
        
        // A negative value for coalecence time t means all species are
        // represented by just one lineage and thus coalescence is not
        // possible within this epoch.
        double t = -1.0;
        
        // If total_rate is not zero, then not all species have been reduced
        // to a single lineage and coalescence is possible (but may be deep)
        if (total_rate > 0.0) {
            t = -log(1 - rng.uniform())/total_rate;
        }
        
        // Four possibilities:
        //  case    t  speciation_found   description
        //  ----  ---  ----------------   ---------------------------------------------------------
        //    1    -1       true          No coalescence (each species just 1 lineage)
        //    2    >0       true          a. Coalescence if t < time_to_next_speciation
        //                                b. No coalescence if t > time_to_next_speciation
        //    3    -1       false         Should never happen (multiple species in final ancestor)
        //    4    >0       false         Coalescence in the ancestral population
        assert(t > -1.0 || speciation_found);   // check for case 3

        // coalescence_occurred is true if case 2a or 4, false if case 1 or 2b
        bool coalescence_occurred = t > 0.0 && t < time_to_next_speciation;
        if (coalescence_occurred) {
            // Extend all lineages by t
            advanceAllLineagesBy(t);
            
            // Normalize the rates to create a vector of probabilities
            vector<double> probs(rates.size(), 0.0);
            transform(rates.begin(), rates.end(), probs.begin(), [total_rate](double r){return r/total_rate;});
            
            // Choose in which species the next coalescence event occurs
            unsigned k = Forest::multinomialDraw(probs);
            Node::species_t s = species[k];
            
            // Get number of lineages in the chosen species
            unsigned n = (unsigned)_lineages_within_species[s].size();
            
            // Perform prior-post (multinomial draw from all n-choose-2 possible joins in species s)
            unsigned npairs = n*(n-1)/2;
            vector<double> log_weights;
            vector<pair<unsigned,unsigned> > pairs;
            vector<PartialStore::partial_t> partials(npairs);
            
            vector<double> log_likelihoods;
            vector<double> log_coalescent_likelihoods;
            unsigned the_pair = 0;
            for (unsigned i = 0; i < n - 1; i++) {
                for (unsigned j = i + 1; j < n; j++) {
                    Node * first_node  = _lineages_within_species[s][i];
                    Node * second_node = _lineages_within_species[s][j];
                    Node * anc_node    = joinLineagePair(first_node, second_node);
                    removeTwoAddOne(_lineages, first_node, second_node, anc_node);
                    anc_node->setSpecies(s);
                    
                    updateLineageCountsVector(species_lineage_counts);
                    Epoch cepoch(Epoch::coalescent_epoch, _forest_height);
                    cepoch._gene = _gene_index;
                    cepoch._coalescence_node = anc_node;
                    cepoch._lineage_counts = species_lineage_counts;
                    cepoch._species = anc_node->getSpecies();
                    auto eit = speciation_found ? insertEpochBefore(_epochs, cepoch, sit) : pushBackEpoch(_epochs, cepoch);
                    
                    // Compute partial likelihood array
                    assert(anc_node->_partial == nullptr);
                    anc_node->_partial = ps.getPartial(_gene_index);
                                                            
                    assert(anc_node->_left_child);
                    assert(anc_node->_left_child->_right_sib);
                    if (!simulating) {
                        assert(_data);
                        calcPartialArray(anc_node);
                    }
                    
                    double log_likelihood = 0.0;
                    double log_coalescent_likelihood = 0.0;
                    double log_weight = 0.0;
                    if (!simulating) {
                        // Compute log Felsenstein likelihood
                        log_likelihood = calcLogLikelihood();
                        log_weight = log_likelihood - _prev_log_likelihood;
                
                        // Compute log coalescent likelihood
                        log_coalescent_likelihood = calcLogCoalescentLikelihood(_epochs, _gene_index);
                        log_weight += log_coalescent_likelihood - _prev_log_coalescent_likelihood;
                    }
                    
                    // Store all relevant information about this pair
                    pairs.push_back(make_pair(i,j));
                    
                    partials[the_pair++] = anc_node->_partial;
                    anc_node->_partial.reset();

                    log_weights.push_back(log_weight);
                    log_likelihoods.push_back(log_likelihood);
                    log_coalescent_likelihoods.push_back(log_coalescent_likelihood);
                    
                    // Unjoin lineage pair
                    unjoinLineagePair(anc_node, first_node, second_node);
                    addTwoRemoveOne(_lineages, first_node, second_node, anc_node);
                                        
                    removeEpochAt(_epochs, eit);
                }
            }

            // Compute sum of weights on log scale
            log_sum_weights = calcLogSum(log_weights);
            
            // Normalize log weights to create a discrete probability distribution
            probs.resize(log_weights.size());
            transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
            
            // Choose one pair to join
            unsigned which_pair = Forest::multinomialDraw(probs);
            auto chosen_pair = pairs[which_pair];
            
            // Update previous log-likelihood and log-coalescent-likelihood values
            _prev_log_likelihood            = log_likelihoods[which_pair];
            _prev_log_coalescent_likelihood = log_coalescent_likelihoods[which_pair];
            
            Node * first_node  = _lineages_within_species[s][chosen_pair.first];
            Node * second_node = _lineages_within_species[s][chosen_pair.second];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            
            // Ancestral node is in the same species as the nodes that were joined
            anc_node->setSpecies(s);
            
            // Set ancestral node partial likelihood array to the one already calculated
            if (!simulating) {
                assert(_data);
                assert(anc_node->_partial == nullptr);
                //anc_node->_partial = ps.getPartial(_gene_index); //bug: just transfer existing, no need to allocate anew
                assert(anc_node->_left_child);
                assert(anc_node->_left_child->_right_sib);
                anc_node->_partial = partials[which_pair];
                partials[which_pair].reset();
                            
                // Store unused partials
                for (unsigned p = 0; p < partials.size(); p++) {
                    if (p != which_pair) {
                        ps.stowPartial(partials[p], _gene_index);
                        partials[p].reset();
                    }
                }
            }
            
            // Update lineage vector
            removeTwoAddOne(_lineages_within_species[s], first_node, second_node, anc_node);
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            
            // Create entry in _epochs for this coalescent event
            updateLineageCountsVector(species_lineage_counts);
            Epoch cepoch(Epoch::coalescent_epoch, _forest_height);
            cepoch._gene = _gene_index;
            cepoch._coalescence_node = anc_node;
            cepoch._lineage_counts = species_lineage_counts;
            cepoch._species = anc_node->getSpecies();
            if (speciation_found)
                insertEpochBefore(_epochs, cepoch, sit);
            else
                pushBackEpoch(_epochs, cepoch);

            debugShow(format("    coalesced %d and %d in species %d at height %.5f\n    %s") % chosen_pair.first % chosen_pair.second % *(s.begin()) % _forest_height % makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
        }
        else {
            // No gene tree coalescence before the end of the species tree epoch
            // Extend all lineages by time_to_next_speciation
            assert(speciation_found);
            assert(time_to_next_speciation != Forest::_infinity);
            advanceAllLineagesBy(time_to_next_speciation);
            
            // Ensure that this speciation epoch will not be used again until it is reset
            sit->_valid = false;
            
            // Every lineage previously assigned to left_species or right_species
            // should be reassigned to anc_species
            Node::species_t left_species  = sit->_left_species;
            Node::species_t right_species = sit->_right_species;
            Node::species_t anc_species   = sit->_anc_species;
            auto reassign = [&left_species, &right_species, &anc_species](Node * nd) {
                Node::species_t & ndsppset = nd->getSpecies();
                if (ndsppset == left_species || ndsppset == right_species)
                    nd->setSpecies(anc_species);
            };
            for_each(_lineages.begin(), _lineages.end(), reassign);
            
            debugShow(format("    hit speciation event at height %.5f\n    %s") % _forest_height % makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
        }

        return make_pair(coalescence_occurred, log_sum_weights);
    }

    inline void GeneForest::operator=(const GeneForest & other) {
        Forest::operator=(other);
        _data = other._data;
        _gene_index = other._gene_index;

        // Node _partial data members not copied in base class because
        // partials are only relevant for gene forests (and gene index
        // must be known
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            // Find out whether this and other have partials
            bool this_partial_exists = (_nodes[i]._partial != nullptr);
            bool other_partial_exists = (other._nodes[i]._partial != nullptr);
            
            if (this_partial_exists && other_partial_exists) {
                // Sanity check: make sure _partials are both of the correct length
                assert(other._nodes[i]._partial->size() == ps.getNElements(_gene_index));
                assert(_nodes[i]._partial->size() == ps.getNElements(_gene_index));
                
                // OK to copy
                *(_nodes[i]._partial)  = *(other._nodes[i]._partial);
            }
            else if (this_partial_exists && !other_partial_exists) {
                // Sanity check: make sure this _partial is of the correct length
                assert(_nodes[i]._partial->size() == ps.getNElements(_gene_index));
                
                // OK to stow existing partial
                ps.stowPartial(_nodes[i]._partial, _gene_index);
                
                // OK to set partial to null
                _nodes[i]._partial = nullptr;
            }
            else if (other_partial_exists && !this_partial_exists) {
                // Sanity check: make sure other _partial is of the correct length
                assert(other._nodes[i]._partial->size() == ps.getNElements(_gene_index));
                
                // Grab a partial of the correct length
                _nodes[i]._partial     = ps.getPartial(_gene_index);
                assert(_nodes[i]._partial->size() == ps.getNElements(_gene_index));
                
                // OK to copy
                *(_nodes[i]._partial)  = *(other._nodes[i]._partial);
            }
        }
        
        // No need to copy _lineages_within_species because it is
        // only used in GeneForest::advanceGeneForest and is rebuilt
        // every time it is used
    }
    
    inline void GeneForest::computeLeafPartials(Data::SharedPtr data) {
        assert(data);
        assert(_leaf_partials.size() == 0);
        assert(Forest::_ngenes > 0);
        assert(Forest::_ntaxa > 0);
        assert(Forest::_nstates == 4);
        
        // Allocate a vector of leaf partials for each gene
        _leaf_partials.resize(Forest::_ngenes);
                
        // Get reference to raw data matrix, which is a vector of vector<state_t>
        auto data_matrix = data->getDataMatrix();

        // Create vector of leaf partials
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
            // Get number of patterns and first pattern index for gene g
            unsigned npatterns = data->getNumPatternsInSubset(g);
            Data::begin_end_pair_t be = data->getSubsetBeginEnd(g);
            unsigned first_pattern = be.first;
            
            // Allocate a partial array for each taxon for gene g
            // _leaf_partials is vector< vector<partial_t> >
            // _leaf_partials[g] is vector<partial_t>
            // _leaf_partials[g][t] is partial_t (i.e. shared pointer to partial array)
            _leaf_partials[g].resize(Forest::_ntaxa);

            for (unsigned t = 0; t < Forest::_ntaxa; ++t) {
                PartialStore::partial_t partial_ptr = ps.getPartial(g);
                vector<double> & leaf_partial = *partial_ptr;
                
                // Set each partial according to the observed data for leaf t
                for (unsigned p = 0; p < npatterns; p++) {
                    unsigned pp = first_pattern + p;
                    for (unsigned s = 0; s < Forest::_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[t][pp];
                        double result = state & d;
                        leaf_partial[p*Forest::_nstates + s] = (result == 0.0 ? 0.0 : 1.0);
                    }
                    
                    // Ensure that the _nstates partials add up to at least 1
                    assert(accumulate(partial_ptr->begin()+p*Forest::_nstates, partial_ptr->begin()+p*Forest::_nstates + Forest::_nstates, 0.0) >= 1.0);
                }

                _leaf_partials[g][t] = partial_ptr;
            }
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
                nd->_partial = ps.getPartial(_gene_index);
                if (nd->_left_child) {
                    calcPartialArray(nd);
                }
                else {
                    *(nd->_partial) = *(_leaf_partials[_gene_index][nd->_number]);
                }
            }
        }
    }
}
