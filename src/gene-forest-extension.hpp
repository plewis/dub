#pragma once

#if defined(LAZY_COPYING)

namespace proj {

    class GeneForestExtension {
        public:
            
                            GeneForestExtension();
                                
            void            dock(const GeneForest::SharedPtr gf, PartialStore::partial_t partial);
            void            undock();

            G::uint_pair_t  chooseNodesToJoin(const vector<unsigned> & node_indices) const;
            G::species_t    chooseSpecies(double total_rate, double theta) const;
            unsigned        speciesCount(G::species_t spp) const;
            double          getProposedDelta() const;
            double          getHeight() const;
            double          getLogWeight() const;
            const Node *    getProposedAnc() const;
            const Node *    getProposedLChild() const;
            const Node *    getProposedRChild() const;
            
            const G::merge_vect_t & getMergers() const;
            
            double          calcTotalRate(double theta);
            void            mergeSpecies(G::species_t left_spp, G::species_t right_spp);
            void            advanceProposedDeltaBy(double dt);
            void            coalesce(double total_rate);
            void            debugCheckSpeciesVect() const;
            
            PartialStore::partial_t getPartial();

        private:
        
            // key: species ==> value: vector of indices into docked forest _lineages
            typedef unordered_map<G::species_t, vector<unsigned> > sppmap_t;

            GeneForest::ConstSharedPtr    _docked_gene_forest;
            double                        _log_weight;
            double                        _proposed_delta;
            Node                          _proposed_anc;
            const Node *                  _proposed_lchild;
            const Node *                  _proposed_rchild;
            vector<G::species_t>          _species_vect;
            G::merge_vect_t               _mergers;
            sppmap_t                      _sppmap;
    };
    
    inline GeneForestExtension::GeneForestExtension() {
        undock();
    }
    
    inline void GeneForestExtension::dock(const GeneForest::SharedPtr gf, PartialStore::partial_t partial) {
        // Check to make sure this extension was previously undocked
        assert(gf);
        assert(_docked_gene_forest == nullptr);
        assert(_species_vect.empty());

        // Attach GeneForest
        _docked_gene_forest = gf;
        
        // Create vector of species corresponding to lineages in _docked_gene_forest
        _docked_gene_forest->copyLineageSpecies(_species_vect);
        
        _proposed_delta = 0.0;
        _proposed_anc._height = gf->getForestHeight();
        _proposed_anc._edge_length = 0.0;
        _proposed_anc._partial = partial;
    }
    
    inline void GeneForestExtension::undock() {
        _docked_gene_forest.reset();
        _species_vect.clear();
        _mergers.clear();
        _sppmap.clear();
        _log_weight = 0.0;
        _proposed_delta = 0.0;
        _proposed_anc.clearPointers();
        _proposed_anc._number = -2;
        _proposed_anc._name = "fake";
        _proposed_anc._height = 0.0;
        _proposed_anc._partial.reset();
        _proposed_anc._edge_length = 0.0;
        _proposed_anc._flags = 0;
        _proposed_anc._species = 0;
        _proposed_anc._split.clear();
        _proposed_lchild = nullptr;
        _proposed_rchild = nullptr;
    }

    inline unsigned GeneForestExtension::speciesCount(G::species_t spp) const {
        return (unsigned)count(_species_vect.begin(), _species_vect.end(), spp);
    }

    inline const Node * GeneForestExtension::getProposedAnc() const {
        return &_proposed_anc;
    }
    
    inline const Node * GeneForestExtension::getProposedLChild() const {
        return _proposed_lchild;
    }
    
    inline const Node * GeneForestExtension::getProposedRChild() const {
        return _proposed_rchild;
    }
    
    inline double GeneForestExtension::getProposedDelta() const {
        return _proposed_delta;
    }
    
    inline double GeneForestExtension::getHeight() const {
        return _proposed_anc._height;
    }
    
    inline double GeneForestExtension::getLogWeight() const {
        return _log_weight;
    }
    
    inline const G::merge_vect_t & GeneForestExtension::getMergers() const {
        return _mergers;
    }
    
    inline void GeneForestExtension::advanceProposedDeltaBy(double dt) {
        assert(dt > 0.0);
        _proposed_delta += dt;
        _proposed_anc._height += dt;
    }
    
    inline void GeneForestExtension::mergeSpecies(G::species_t left_spp, G::species_t right_spp) {
        // Assumes height of _proposed_anc is the level at which the merger occurs
        
        // Add element to _mergers vector
        _mergers.push_back(make_tuple(_proposed_anc._height, left_spp, right_spp));
            
        // Merge species in _species_vect
        G::species_t anc_spp = (left_spp | right_spp);
        for (auto & s : _species_vect) {
            if (s == left_spp || s == right_spp) {
                s = anc_spp;
            }
        }
    }
    
    inline double GeneForestExtension::calcTotalRate(double theta) {
        // Enumerate lineages belonging to each current species
        _sppmap.clear();
        unsigned i = 0;
        for (auto s : _species_vect) {
            _sppmap[s].push_back(i++);
        }
        
        // Calculate total rate and return it
        double total_rate = 0.0;
        for (auto & sv : _sppmap) {
            unsigned n = (unsigned)sv.second.size();
            double r = 1.0*n*(n-1)/theta;
            total_rate += r;
        }
        return total_rate;
    }
    
    inline G::species_t GeneForestExtension::chooseSpecies(double total_rate, double theta) const {
        unsigned nspecies = (unsigned)_sppmap.size();
        assert(nspecies > 0);
        
        vector<double> probs(nspecies, 0.0);
        vector<G::species_t> species_vect;
        species_vect.reserve(nspecies);
        
        // If n_i = number of lineages in species i, the probability of coalescence is
        //   r_i / total_rate = (n_i*(n_i-1)/theta) / (sum_i r_i)
        unsigned i = 0;
        for (auto & sv : _sppmap) {
            species_vect.push_back(sv.first);
            double n = (double)sv.second.size();
            probs[i++] = n*(n - 1.0)/(theta*total_rate);
        }
        
        // Choose a species using a multinomial draw from probs
        unsigned which = G::multinomialDraw(::rng, probs);
        assert(which < probs.size());
        G::species_t spp = species_vect[which];
        return spp;
    }
    
    inline G::uint_pair_t GeneForestExtension::chooseNodesToJoin(const vector<unsigned>  & node_indices) const {
        unsigned n = (unsigned)node_indices.size();
        assert(n > 1);
        
        // Choose a random pair of lineages to join
        pair<unsigned,unsigned> chosen_pair = ::rng->nchoose2(n);
        unsigned i = node_indices[chosen_pair.first];
        unsigned j = node_indices[chosen_pair.second];
        
        // Sanity checks
        assert(i != j);
        assert(i < _docked_gene_forest->_lineages.size());
        assert(j < _docked_gene_forest->_lineages.size());

        return make_pair(i,j);
    }

    inline void GeneForestExtension::coalesce(double total_rate) {
        // Choose the species in which the coalescence will happen
        G::species_t spp = chooseSpecies(total_rate, G::_theta);
        _proposed_anc.setSpecies(spp);
                    
        // Choose the two nodes to join
        G::uint_pair_t chosen = chooseNodesToJoin(_sppmap[spp]);
        
        // Get pointers to the two nodes to join
        _proposed_lchild = _docked_gene_forest->_lineages[chosen.first];
        _proposed_rchild = _docked_gene_forest->_lineages[chosen.second];
        
        // Set _proposed_anc's split to union of the two child splits
        _proposed_anc._split.resize(G::_ntaxa);
        _proposed_anc._split += _proposed_lchild->_split;
        _proposed_anc._split += _proposed_rchild->_split;

        // Compute partial likelihood array of ancestral node
        _log_weight = _docked_gene_forest->calcPartialArray(&_proposed_anc, _proposed_lchild, _proposed_rchild);

        assert(!isnan(_log_weight));
        assert(!isinf(_log_weight));
    }
    
    inline PartialStore::partial_t GeneForestExtension::getPartial() {
        return _proposed_anc._partial;
    }

    inline void GeneForestExtension::debugCheckSpeciesVect() const {
        unsigned n = (unsigned)_species_vect.size();
        assert(_docked_gene_forest->_lineages.size() == n);
        for (unsigned i = 0; i < n; i++) {
            assert(_docked_gene_forest->_lineages[i]->_species == _species_vect[i]);
        }
    }

}
#endif

