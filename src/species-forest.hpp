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
        _log_species_tree_prior = 0.0;
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
            SMCGlobal::species_t spp = anc->getSpecies();
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
        //TODO: SpeciesForest::speciationEvent
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

}
