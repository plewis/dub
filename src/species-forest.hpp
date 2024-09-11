#pragma once

class Particle;

namespace proj {

    class SpeciesForest : public Forest {
    
        friend class Particle;
        
        public:
        
            SpeciesForest();
             ~SpeciesForest();
             
            void speciationEvent(Lot::SharedPtr lot, G::species_t & left, G::species_t & right, G::species_t & anc);
            
            Node * findSpecies(G::species_t spp);
                                                
            void fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const;
            tuple<double,double,double> drawTruncatedIncrement(Lot::SharedPtr lot, double truncate_at);
            pair<double,double> drawIncrement(Lot::SharedPtr lot);
            double calcLogSpeciesTreeDensity(double lambda) const;

            // Overrides of base class functions
            void clear();
            void operator=(const SpeciesForest & other);
            
            // Overrides of abstract base class functions
            void finalizeProposal();
            void createTrivialForest(bool compute_partials = false);
            bool isSpeciesForest() const {return true;}
            void setSpeciesFromNodeName(Node * nd);
            void addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient);
            
            void buildCoalInfoVect();
            void saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap = false) const;
            void recordHeights(vector<double> & height_vect) const;
            
            tuple<double, G::species_t, G::species_t, G::species_t> findNextSpeciationEvent(double starting_height) const;
            void rebuildStartingFromHeight(double starting_height);
            
            typedef shared_ptr<SpeciesForest> SharedPtr;

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
    
    inline Node * SpeciesForest::findSpecies(G::species_t spp) {
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

    inline tuple<double,double,double> SpeciesForest::drawTruncatedIncrement(Lot::SharedPtr lot, double truncate_at) {
        // Draw from exponential distribution truncated at the specified value
        unsigned n = (unsigned)_lineages.size();
        double incr = G::_infinity;
        double rate = 0.0;
        double cum_prob_upper_bound = 1.0;
        if (n > 1) {
            rate = G::_lambda*n;
            double cum_prob_upper_bound = 1.0 - exp(-rate*truncate_at);
            incr = -log(1.0 - cum_prob_upper_bound*lot->uniform())/rate;
        }
        return make_tuple(incr, rate, cum_prob_upper_bound);
    }
    
    inline pair<double,double> SpeciesForest::drawIncrement(Lot::SharedPtr lot) {
        unsigned n = (unsigned)_lineages.size();
        double incr = G::_infinity;
        double rate = 0.0;
        if (n > 1) {
            rate = G::_lambda*n;
            incr = -log(1.0 - lot->uniform())/rate;
        }
        return make_pair(incr, rate);
    }

    // This is an override of an abstract base class function
    inline void SpeciesForest::setSpeciesFromNodeName(Node * nd) {
        auto it = find(G::_species_names.begin(), G::_species_names.end(), nd->_name);
        if (it == G::_species_names.end())
            throw XProj(str(format("Could not find an index for the species name \"%s\"") % nd->_name));
        else {
            unsigned i = (unsigned)distance(G::_species_names.begin(), it);
            Node::setSpeciesBit(nd->_species, i, /*init_to_zero_first*/true);
        }
    }
            
    // This is an override of an abstract base class function
    inline void SpeciesForest::createTrivialForest(bool compute_partials) {
        assert(!compute_partials); // catch non-sensical true value
        clear();
        unsigned nnodes = 2*G::_nspecies - 1;
        Forest::_nodes.resize(nnodes);
        for (unsigned i = 0; i < G::_nspecies; i++) {
            string species_name = G::_species_names[i];
            _nodes[i]._number = (int)i;
            _nodes[i]._my_index = (int)i;
            _nodes[i]._name = species_name;
            _nodes[i]._edge_length = 0.0;
            _nodes[i]._height = 0.0;
            Node::setSpeciesBit(_nodes[i]._species, i, /*init_to_zero_first*/true);
            _lineages.push_back(&_nodes[i]);
        }

        // Add all remaining nodes to _unused_nodes vector
        _unused_nodes.clear();
        for (unsigned i = G::_nspecies; i < nnodes; i++) {
            _nodes[i]._my_index = (int)i;
            _nodes[i]._number = -1;
            _unused_nodes.push_back(i);
        }
        
        refreshAllPreorders();
        _forest_height = 0.0;
        //_next_node_index = G::_nspecies;
        //_next_node_number = G::_nspecies;
        _log_likelihood = 0.0;
        _prev_log_likelihood = 0.0;
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
        assert(internal_heights.size() == G::_nspecies - 1);
        
        // Sort heights
        sort(internal_heights.begin(), internal_heights.end());
        
        double log_prob_density = 0.0;
        unsigned n = G::_nspecies;
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
    
    inline void SpeciesForest::fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const {
        // No fixing up to do if there are no species tree joins
        if (sppinfo_vect.empty())
            return;
                    
        debugCheckCoalInfoSorted(coalinfo_vect);
        debugCheckCoalInfoSorted(sppinfo_vect);
        
        // Example:
        //          -- before --           -- after ---
        //  height   left  right  species    left   right
        // -------  -----  -----  -------  -----  -----
        // 0.01167      E      E              E       E
        // 0.01389      B      B              B       B
        // 0.02047      B      B              B       B
        // 0.02079      D      D              D       D
        // 0.02150      D      D              D       D
        // 0.03438      C      C              C       C
        // 0.08638  ------------       CD  ------------ <-- currently considering
        // 0.11725  ------------      ACD  ------------
        // 0.13886      C      A             CD       A
        // 0.13942      D      A             CD       A
        // 0.17349     AC     AD             AC      AD
        // 0.39425    ACD      E            ACD       E
        // 0.41254      B   ACDE              B    ACDE
        //
        //          -- before --           -- after ---
        //  height   left  right  species    left   right
        // -------  -----  -----  -------  -----  -----
        // 0.01167      E      E              E       E
        // 0.01389      B      B              B       B
        // 0.02047      B      B              B       B
        // 0.02079      D      D              D       D
        // 0.02150      D      D              D       D
        // 0.03438      C      C              C       C
        // 0.08638  ------------       CD  ------------
        // 0.11725  ------------      ACD  ------------ <-- currently considering
        // 0.13886     CD      A            ACD     ACD
        // 0.13942     CD      A            ACD     ACD
        // 0.17349     AC     AD            ACD     ACD
        // 0.39425    ACD      E            ACD       E
        // 0.41254      B   ACDE              B    ACDE

        unsigned jstart = 0;
        for (auto & si : sppinfo_vect) {
            // Get handles for current sppinfo_vect element
            double                 h0 = get<0>(si);
            vector<G::species_t> & v0 = get<2>(si);
            
            // Can't assume v0 has size 2 because this species tree join may be the fake
            // join that combines all remaining species into an ancestral species
            // for incomplete species forests
            G::species_t           s0 = accumulate(v0.begin(), v0.end(), (G::species_t)0, [](const G::species_t & next, const G::species_t & cum){return (cum | next);});
            
            // Advance through coalinfo_vect to the first event that might be
            // affected by the current species tree join
            unsigned j = jstart;
            assert(j < coalinfo_vect.size());
            double h = get<0>(coalinfo_vect[j]);
            while (h < h0) {
                j++;
                assert(j < coalinfo_vect.size());
                h = get<0>(coalinfo_vect[j]);
            }
            
            // Next time can start at this point in the coalinfo_vect
            jstart = j;
            
            // Perform replacements in all subsequent gene tree coalinfo elements
            while (j < coalinfo_vect.size()) {
                h = get<0>(coalinfo_vect[j]);
                vector<G::species_t> & v = get<2>(coalinfo_vect[j]);
                assert(v.size() == 2);
                if (subsumed(v[0], s0)) {
                    v[0] = s0;
                }
                if (subsumed(v[1], s0)) {
                    v[1] = s0;
                }
                j++;
            }
        }
    }
        
    inline void SpeciesForest::speciationEvent(Lot::SharedPtr lot, G::species_t & left_spp, G::species_t & right_spp, G::species_t & anc_spp) {
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
    }
        
    inline void SpeciesForest::operator=(const SpeciesForest & other) {
        Forest::operator=(other);
    }
    
    inline void SpeciesForest::addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient) {
        // Assumes nd is an internal node
        assert(nd->_left_child);
        
        recipient.push_back(
            make_tuple(
                nd->_height,
                0,
                vector<G::species_t>({
                    nd->_left_child->_species,
                    nd->_left_child->_right_sib->_species,
                })
            )
        );
    }

    inline void SpeciesForest::saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap) const {
        // Appends to coalinfo_vect; clear coalinfo_vect before calling if desired
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this

        // Dump _coalinfo into coalinfo_vect
        if (!_coalinfo.empty()) {
            coalinfo_vect.insert(coalinfo_vect.begin(), _coalinfo.begin(), _coalinfo.end());
        }
        
        if (cap) {
            // Create an entry pooling the remaining species
            if (_lineages.size() > 1) {
                vector<G::species_t> sppvect;
                for (auto nd : _lineages) {
                    sppvect.push_back(nd->_species);
                }
                coalinfo_vect.push_back(make_tuple(
                    _forest_height,
                    0,
                    sppvect)
                );
            }
        }
    }
    
    inline void SpeciesForest::recordHeights(vector<double> & height_vect) const {
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
    
    inline void SpeciesForest::buildCoalInfoVect() {
        //TODO: GeneForest has same function: move to base class Forest?
        // Should only be called for complete species trees
        assert(_lineages.size() == 1);

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
    
    inline tuple<double, G::species_t, G::species_t, G::species_t> SpeciesForest::findNextSpeciationEvent(double starting_height) const {
        // Assumes that species forest is a tree with just one subtree
        assert(_preorders.size() == 1);
        assert(_lineages.size() == 1);
        
        double speciation_height = G::_infinity;
        G::species_t left_spp    = (G::species_t)0;
        G::species_t right_spp   = (G::species_t)0;
        G::species_t anc_spp     = (G::species_t)0;
        
        Node * root = _lineages[0];
        double root_height = root->getHeight();
        assert(fabs(root_height - getHeight()) < G::_small_enough);

        if (root_height > starting_height) {
                
            // Find node with smallest height that is nevertheless higher than starting_height
            pair<double, Node *> best = make_pair(G::_infinity, nullptr);
            for (auto nd : _preorders[0]) {
                if (nd->getLeftChild()) {  // Only internal nodes specify speciation events
                    double h = nd->getHeight();
                    if (h > starting_height && h < best.first) {
                        best = make_pair(h, nd);
                    }
                }
            }
            assert(best.first < G::_infinity);
            speciation_height = best.first;
            left_spp  = best.second->getLeftChild()->getSpecies();
            right_spp = best.second->getLeftChild()->getRightSib()->getSpecies();
            anc_spp   = best.second->getSpecies();
        }

        return make_tuple(speciation_height, left_spp, right_spp, anc_spp);
    }
    
    void SpeciesForest::rebuildStartingFromHeight(double starting_height) {
        refreshAllHeightsAndPreorders();
        
        if (_lineages.size() == 1) {
            // Species forest is a complete tree
            
            // Start at root (highest node)
            Node * nd = _lineages[0];
            if (nd->getHeight() <= starting_height) {
                // Nothing to be done
                return;
            }
        
            bool done = false;
            while (!done) {
                vector<Node *> nodes_to_remove;
                vector<Node *> nodes_to_add;
                bool found = false;
                for (auto nd : _lineages) {
                    // Remove nd if its height is greater than starting_height
                    if (nd->getHeight() > starting_height) {
                        found = true;
                        
                        // Identify left and right child
                        Node * lchild = nd->getLeftChild();
                        Node * rchild = lchild->getRightSib();
                        assert(lchild);
                        assert(rchild);
                        assert(!rchild->getRightSib());
                        
                        // Store nodes to add or remove (cannot do it now because we
                        // are currently iterating through _lineages
                        nodes_to_remove.push_back(nd);
                        nodes_to_add.push_back(lchild);
                        nodes_to_add.push_back(rchild);
                        
                        // Detach nd from its children and its children
                        // from the rest of the tree (the children are
                        // now independent subtrees in the forest)
                        lchild->setParent(nullptr);
                        rchild->setParent(nullptr);
                        lchild->setRightSib(nullptr);
                        rchild->setRightSib(nullptr);
                        nd->setLeftChild(nullptr);
                        assert(!nd->getRightSib());
                    }
                }
                
                for (auto nd : nodes_to_remove) {
                    // Get iterator to node to be deleted and remove from _lineages
                    auto it = find(_lineages.begin(), _lineages.end(), nd);
                    assert(it != _lineages.end());
                    _lineages.erase(it);
                    stowNode(nd);
                }
                
                for (auto nd : nodes_to_add) {
                    _lineages.push_back(nd);
                }
                
                if (!found)
                    done = true;
            }
        }
        else {
            // If species forest is not a complete tree, then it should
            // still be a trivial forest that needs to be constructed
            assert(_lineages.size() == G::_nspecies);
        }
        
        // Trim edge lengths of nodes remaining in _lineages so that
        // species tree height is exactly starting_height
        for (auto nd : _lineages) {
            double node_height = nd->getHeight();
            double node_edgelen = nd->getEdgeLength();
            double overlap = node_height + node_edgelen - starting_height;
            assert(overlap >= 0.0);
            nd->setEdgeLength(node_edgelen - overlap);
        }
        
        // Now build up complete species tree again staring with current height
        while (_lineages.size() > 1) {

            // Draw a speciation increment Delta. Note: speciation_increment will
            // equal "infinity" if species tree is complete.
            auto incr_rate = drawIncrement(::rng);
            double Delta = incr_rate.first;
            
            // advance all lineages by an amount Delta
            advanceAllLineagesBy(Delta);
            
            // Create speciation event
            G::species_t left_spp, right_spp, anc_spp;
            speciationEvent(::rng, left_spp, right_spp, anc_spp);
        }
        
        refreshAllHeightsAndPreorders();
        
        // debugging
        //output(format("  Rebuilt species tree starting from %g: %s\n") % starting_height % makeNewick(9, true, false), 0);
    }
    
}
