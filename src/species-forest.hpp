#pragma once

class Particle;

namespace proj {

    class SpeciesForest : public Forest {
    
        friend class Particle;
        
        public:
        
            SpeciesForest();
             ~SpeciesForest();
             
            void speciationEvent(Lot::SharedPtr lot, SMCGlobal::species_t & left, SMCGlobal::species_t & right, SMCGlobal::species_t & anc);
            
            Node * findSpecies(SMCGlobal::species_t spp);
            
            static pair<double,double> calcTreeDistances(SpeciesForest & ref, SpeciesForest & test);

            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void clearMark();
            void revertToMark();
            void finalizeProposal();
            void createTrivialForest(bool compute_partials = false);
            bool isSpeciesForest() const {return true;}

            pair<double,double> drawIncrement(Lot::SharedPtr lot);
            double calcLogSpeciesTreeDensity(double lambda) const;

            void saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect) const;
            void recordHeights(vector<double> & height_vect) const;
            
            void operator=(const SpeciesForest & other);

        protected:
            
            double advanceSpeciesForest(unsigned particle, unsigned step);
            
            // NOTE: any variables added must be copied in operator=
    };
    
    inline SpeciesForest::SpeciesForest() {
    }

    inline SpeciesForest::~SpeciesForest() {
        clear();
    }

    inline void SpeciesForest::clear() {
        Forest::clear();
    }

    inline pair<double,double> SpeciesForest::calcTreeDistances(SpeciesForest & ref, SpeciesForest & test) {
        // Store splits from reference tree
        Split::treeid_t ref_splits;
        ref.storeSplits(ref_splits);
        
        // Store edge lengths from reference tree
        map<Split, double> ref_edgelen_map;
        ref.storeEdgelensBySplit(ref_edgelen_map);
                
        // Store splits from test tree
        Split::treeid_t test_splits;
        test.storeSplits(test_splits);
        
        // Store edge lengths from test tree
        map<Split, double> test_edgelen_map;
        test.storeEdgelensBySplit(test_edgelen_map);
                
        // Now calculate squares for leaf nodes, storing in KLleaves
        std::vector<double> KLleaves(SMCGlobal::_nspecies);
        Split s;
        s.resize(SMCGlobal::_nspecies);
        Split sroot;
        sroot.resize(SMCGlobal::_nspecies);
        for (unsigned i = 0; i < SMCGlobal::_nspecies; i++) {
            s.clear();
            s.setBitAt(i);
            sroot.setBitAt(i);
            assert(ref_edgelen_map.count(s) == 1);
            assert(test_edgelen_map.count(s) == 1);
            double ref_edge_length  = ref_edgelen_map[s];
            double test_edge_length = test_edgelen_map[s];
            double square = pow(test_edge_length - ref_edge_length, 2.0);
            KLleaves[i] = square;
        }
        
        // Store union of refsplits and testsplits in allsplits
        Split::treeid_t all_splits;
        set_union(
            ref_splits.begin(), ref_splits.end(),
            test_splits.begin(), test_splits.end(),
            std::inserter(all_splits, all_splits.begin()));
        
        // Traverse allsplits, storing squared branch length differences in KLinternals
        std::vector<double> KLinternals(all_splits.size());
        double RFdist = 0.0;
        unsigned i = 0;
        for (auto s : all_splits) {
            if (s == sroot)
                continue;
            bool s_in_ref  = ref_edgelen_map.count(s) == 1;
            bool s_in_test = test_edgelen_map.count(s) == 1;
            assert(s_in_ref || s_in_test);
            if (!s_in_ref) {
                double test_edge_length = test_edgelen_map[s];
                double test_square = pow(test_edge_length, 2.0);
                KLinternals[i++] = test_square;
                RFdist += 1.0;
            }
            else if (!s_in_test) {
                double ref_edge_length = ref_edgelen_map[s];
                double ref_square = pow(ref_edge_length, 2.0);
                KLinternals[i++] = ref_square;
                RFdist += 1.0;
            }
            else {
                double test_edge_length = test_edgelen_map[s];
                double ref_edge_length = ref_edgelen_map[s];
                double square = pow(test_edge_length - ref_edge_length, 2.0);
                KLinternals[i++] = square;
            }
        }
            
        // Calculate KL distance
        double KFSS = 0.0;
        for (auto square : KLinternals) {
            KFSS += square;
        }
        for (auto square : KLleaves) {
            KFSS += square;
        }
        assert(KFSS >= 0.0);
        double KFdist = sqrt(KFSS);
        return make_pair(KFdist, RFdist);
    }
    
    inline Node * SpeciesForest::findSpecies(SMCGlobal::species_t spp) {
        Node * the_node = nullptr;
        for (auto nd : _lineages) {
            if (nd->_species == spp) {
                the_node = nd;
                break;
            }
        }
        assert(the_node);
        return the_node;
    }

    inline pair<double,double> SpeciesForest::drawIncrement(Lot::SharedPtr lot) {
        unsigned n = (unsigned)_lineages.size();
        double incr = SMCGlobal::_infinity;
        double rate = 0.0;
        if (n > 1) {
            rate = SMCGlobal::_lambda*n;
            incr = -log(1.0 - lot->uniform())/rate;
        }
        return make_pair(incr, rate);
    }

    // This is an override of an abstract base class function
    inline void SpeciesForest::createTrivialForest(bool compute_partials) {
        assert(!compute_partials); // catch non-sensical true value
        clear();
        Forest::_nodes.resize(2*SMCGlobal::_nspecies - 1);
        for (unsigned i = 0; i < SMCGlobal::_nspecies; i++) {
            _nodes[i]._number = (int)i;
            _nodes[i]._name = SMCGlobal::_species_names[i];
            _nodes[i]._edge_length = 0.0;
            _nodes[i]._height = 0.0;
            Node::setSpeciesBit(_nodes[i]._species, i, /*init_to_zero_first*/true);
            _lineages.push_back(&_nodes[i]);
        }
        refreshAllPreorders();
        _forest_height = 0.0;
        _next_node_index = SMCGlobal::_nspecies;
        _next_node_number = SMCGlobal::_nspecies;
        _prev_log_likelihood = 0.0;
        //_prev_log_coalescent_likelihood = 0.0;
    }
    
    inline double SpeciesForest::calcLogSpeciesTreeDensity(double lambda) const {
        // Assume that this species forest is fully resolved
        assert(_preorders.size() == 1);
                
        // Build vector of internal node heights
        vector<double> internal_heights;
        for (const Node * nd : _preorders[0]) {
            if (nd->getLeftChild()) {
                internal_heights.push_back(nd->getHeight());
            }
        }
        
        // Number of internal nodes should be _nspecies - 1
        assert(internal_heights.size() == SMCGlobal::_nspecies - 1);
        
        // Sort heights
        sort(internal_heights.begin(), internal_heights.end());
        
        double log_prob_density = 0.0;
        unsigned n = SMCGlobal::_nspecies;
        double h0 = 0.0;
        for (auto it = internal_heights.begin(); it != internal_heights.end(); ++it) {
            double h = *it;
            double r = lambda*n;
            double logr = log(r);
            double t = h - h0;
            double log_exponential_density = logr - r*t;
            log_prob_density += log_exponential_density;
            h0 = h;
            n--;
        }
        
        return log_prob_density;
    }
            
    inline void SpeciesForest::clearMark() {
        // Clear the stacks so that future speciation events can be recorded
        _mark_increments = {};
        _mark_anc_nodes = {};
        _mark_left_right_pos = {};
    }
    
    inline void SpeciesForest::revertToMark() {
        // Should be one increment corresponding to the coalescent event
        assert(!_mark_increments.empty());
        double dt = _mark_increments.top();
        _mark_increments.pop();
        advanceAllLineagesBy(-dt);
        
        assert(_mark_left_right_pos.size() == _mark_anc_nodes.size());
        
        while (!_mark_left_right_pos.empty()) {
            pair<unsigned, unsigned> p = _mark_left_right_pos.top();
            _mark_left_right_pos.pop();
            unsigned left_pos  = p.first;
            unsigned right_pos = p.second;
            
            // Get pointer to ancestral node representing the speciation event
            Node * anc = _mark_anc_nodes.top();
            _mark_anc_nodes.pop();
            
            assert(anc);
            assert(fabs(anc->getEdgeLength()) < 0.00001);
            
            // Get pointer to ancestral node's left child
            Node * lchild = anc->getLeftChild();
            assert(lchild);

            // Get pointer to ancestral node's right child
            Node * rchild = lchild->getRightSib();
            rchild = lchild->_right_sib;
            assert(rchild);
            
            // Reverse the speciation join
            unjoinLineagePair(anc, lchild, rchild);
            
            // Update lineage vector
            addTwoRemoveOneAt(_lineages, left_pos, lchild, right_pos, rchild, anc);

            // Return anc to the pool of unused nodes
            stowNode(anc);
            anc = nullptr;
            
            // Remove increment associated with this speciation event
            double dt = _mark_increments.top();
            _mark_increments.pop();
            advanceAllLineagesBy(-dt);
        }
        assert(_mark_increments.empty());
    }
    
    inline void SpeciesForest::finalizeProposal() {
        // Empty the _prev_species_stack for all nodes in the species tree
        for (auto & nd : _nodes) {
            nd.emptyPrevSpeciesStack();
        }
        clearMark();
    }
    
    inline void SpeciesForest::speciationEvent(Lot::SharedPtr lot, SMCGlobal::species_t & left_spp, SMCGlobal::species_t & right_spp, SMCGlobal::species_t & anc_spp) {
        //BOOKMARK: speciationEvent
        unsigned nlineages = (unsigned)_lineages.size();
        
        // Choose two lineages to join
        assert(nlineages > 1);
        auto chosen_pair = lot->nchoose2(nlineages);
        unsigned left_pos = chosen_pair.first;
        unsigned right_pos = chosen_pair.second;
        Node * first_node  = _lineages[left_pos];
        Node * second_node = _lineages[right_pos];
        
        // Get species for the two lineages to join
        left_spp = first_node->getSpecies();
        right_spp = second_node->getSpecies();

        // Create ancestral node
        Node * anc_node = pullNode();
        joinLineagePair(anc_node, first_node, second_node);
        anc_node->setSpeciesToUnion(left_spp, right_spp);

        // Get species for the new ancestral node
        anc_spp = anc_node->getSpecies();
        
        // Update _lineages vector
        removeTwoAddOne(_lineages, first_node, second_node, anc_node);
        refreshAllPreorders();
        
        // Save info needed to revert
        _mark_anc_nodes.push(anc_node);
        _mark_left_right_pos.push(make_pair(left_pos, right_pos));
    }
        
    inline void SpeciesForest::operator=(const SpeciesForest & other) {
        Forest::operator=(other);
        
        // Should not be anything in the mark stacks if copying
        assert(_mark_increments.empty());
        assert(_mark_anc_nodes.empty());
        assert(_mark_left_right_pos.empty());

#if defined(DEBUGGING_SANITY_CHECK)
        // Sanity check: make sure that _partials do not exist for either this nor other
        for (unsigned i = 0; i < _nodes.size(); i++) {
            assert((_nodes[i]._partial == nullptr && other._nodes[i]._partial == nullptr));
        }
#endif
    }

    inline void SpeciesForest::saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect) const {
        // Appends to params and splits; clear these before calling if desired
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this

        // Should only be called for complete species trees
        assert(_lineages.size() == 1);
                        
        // Save all internal node heights
        vector< tuple<double, SMCGlobal::species_t, SMCGlobal::species_t> > height_spp_tuples;
        for (auto nd : boost::adaptors::reverse(_preorders[0])) {
            if (nd->_left_child) {
                // internal
                //string split_repr = nd->_split.createPatternRepresentation();
                height_spp_tuples.push_back(make_tuple(
                    nd->_height,
                    nd->_left_child->_species,
                    nd->_left_child->_right_sib->_species));
            }
        }
        
        // Sort heights from smallest to largest
        sort(height_spp_tuples.begin(), height_spp_tuples.end());
        
        // Compute increments between coalescent events and store in params
        unsigned n = SMCGlobal::_nspecies;
        assert(height_spp_tuples.size() == SMCGlobal::_nspecies - 1);
        for (auto h : height_spp_tuples) {
            coalinfo_vect.push_back(make_tuple(get<0>(h), 0, get<1>(h), get<2>(h)));
            n--;
        }
    }
    
    void SpeciesForest::recordHeights(vector<double> & height_vect) const {
        // Appends to height_vect; clear before calling if desired
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this

        // Should only be called for complete species trees
        assert(_lineages.size() == 1);
                        
        for (auto nd : boost::adaptors::reverse(_preorders[0])) {
            if (nd->_left_child) {
                // internal
                height_vect.push_back(nd->_height);
            }
        }
    }
}
