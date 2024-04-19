#pragma once

class Particle;

namespace proj {

    class SpeciesForest : public Forest {
    
        friend class Particle;
        
        public:
        
            SpeciesForest();
             ~SpeciesForest();
             
            void speciationEvent(Lot::SharedPtr lot, G::species_t & left, G::species_t & right, G::species_t & anc, bool mark);
            
            Node * findSpecies(G::species_t spp);
            
#if defined(EST_THETA)
            void drawLineageSpecificThetas();
            void drawThetaMean(double exponential_prior_rate);
            void setThetaMean(double thetamean) {_theta_mean = thetamean;}
            double getThetaMean() const {return _theta_mean;}
            double thetaForSpecies(G::species_t s) const;
#endif
            
            void setJointEstimation(bool is_joint_estimation) {_joint_estimation = is_joint_estimation;}
                        
            void fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect);
            pair<double,double> drawTruncatedIncrement(Lot::SharedPtr lot, double truncate_at);
            pair<double,double> drawIncrement(Lot::SharedPtr lot);
            double calcLogSpeciesTreeDensity(double lambda) const;

            // Overrides of base class functions
            void clear();
            void operator=(const SpeciesForest & other);
            
            // Overrides of abstract base class functions
            void revertToMark();
            void finalizeProposal();
            void createTrivialForest(bool compute_partials = false);
            bool isSpeciesForest() const {return true;}
            void setSpeciesFromNodeName(Node * nd);
            void recordCoalInfoAndClearMark();
            void addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient);
            
            void saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap = false) const;
            void recordHeights(vector<double> & height_vect) const;
            
            typedef shared_ptr<SpeciesForest> SharedPtr;

        protected:
            
            double advanceSpeciesForest(unsigned particle, unsigned step);
            
            // NOTE: any variables added must be copied in operator=
            bool _joint_estimation;
#if defined(EST_THETA)
            map<G::species_t, double> _theta_map;
            double _theta_mean;
#endif
    };
    
    inline SpeciesForest::SpeciesForest() {
    }

    inline SpeciesForest::~SpeciesForest() {
        clear();
    }

    inline void SpeciesForest::clear() {
        Forest::clear();
        _joint_estimation = true;
    }

#if defined(EST_THETA)
    inline void SpeciesForest::drawLineageSpecificThetas() {
        // Draw a theta value for each leaf species
        _theta_map.clear();
        for (auto nd : _lineages) {
            G::species_t s = nd->_species;
            if (G::_theta_mean_fixed > 0.0 && G::_theta_mean_frozen)
                _theta_map[s] = G::_theta_mean_fixed;
            else {
                double b = (G::_invgamma_shape - 1.0)*_theta_mean;
                _theta_map[s] = G::inverseGammaVariate(G::_invgamma_shape, b);
            }
        }
    }
    
    inline void SpeciesForest::drawThetaMean(double exponential_prior_rate) {
        // Should only be called for trivial forests
        assert(_lineages.size() == G::_nspecies);
        
        // Draw _theta_mean from Exponential prior
        _theta_mean = -log(1.0 - rng.uniform())/exponential_prior_rate;
    }
#endif

//    inline pair<double,double> SpeciesForest::calcTreeDistances(SpeciesForest & ref, SpeciesForest & test) {
//        // Store splits from reference tree
//        Split::treeid_t ref_splits;
//        ref.storeSplits(ref_splits);
//        
//        // Store edge lengths from reference tree
//        map<Split, double> ref_edgelen_map;
//        ref.storeEdgelensBySplit(ref_edgelen_map);
//                
//        // Store splits from test tree
//        Split::treeid_t test_splits;
//        test.storeSplits(test_splits);
//        
//        // Store edge lengths from test tree
//        map<Split, double> test_edgelen_map;
//        test.storeEdgelensBySplit(test_edgelen_map);
//                
//        // Now calculate squares for leaf nodes, storing in KLleaves
//        std::vector<double> KLleaves(G::_nspecies);
//        Split s;
//        s.resize(G::_nspecies);
//        Split sroot;
//        sroot.resize(G::_nspecies);
//        for (unsigned i = 0; i < G::_nspecies; i++) {
//            s.clear();
//            s.setBitAt(i);
//            sroot.setBitAt(i);
//            assert(ref_edgelen_map.count(s) == 1);
//            assert(test_edgelen_map.count(s) == 1);
//            double ref_edge_length  = ref_edgelen_map[s];
//            double test_edge_length = test_edgelen_map[s];
//            double square = pow(test_edge_length - ref_edge_length, 2.0);
//            KLleaves[i] = square;
//        }
//        
//        // Store union of refsplits and testsplits in allsplits
//        Split::treeid_t all_splits;
//        set_union(
//            ref_splits.begin(), ref_splits.end(),
//            test_splits.begin(), test_splits.end(),
//            std::inserter(all_splits, all_splits.begin()));
//        
//        // Traverse allsplits, storing squared branch length differences in KLinternals
//        std::vector<double> KLinternals(all_splits.size());
//        double RFdist = 0.0;
//        unsigned i = 0;
//        for (auto s : all_splits) {
//            if (s == sroot)
//                continue;
//            bool s_in_ref  = ref_edgelen_map.count(s) == 1;
//            bool s_in_test = test_edgelen_map.count(s) == 1;
//            assert(s_in_ref || s_in_test);
//            if (!s_in_ref) {
//                double test_edge_length = test_edgelen_map[s];
//                double test_square = pow(test_edge_length, 2.0);
//                KLinternals[i++] = test_square;
//                RFdist += 1.0;
//            }
//            else if (!s_in_test) {
//                double ref_edge_length = ref_edgelen_map[s];
//                double ref_square = pow(ref_edge_length, 2.0);
//                KLinternals[i++] = ref_square;
//                RFdist += 1.0;
//            }
//            else {
//                double test_edge_length = test_edgelen_map[s];
//                double ref_edge_length = ref_edgelen_map[s];
//                double square = pow(test_edge_length - ref_edge_length, 2.0);
//                KLinternals[i++] = square;
//            }
//        }
//            
//        // Calculate KL distance
//        double KFSS = 0.0;
//        for (auto square : KLinternals) {
//            KFSS += square;
//        }
//        for (auto square : KLleaves) {
//            KFSS += square;
//        }
//        assert(KFSS >= 0.0);
//        double KFdist = sqrt(KFSS);
//        return make_pair(KFdist, RFdist);
//    }
    
#if defined(EST_THETA)
    inline double SpeciesForest::thetaForSpecies(G::species_t s) const {
        if (_theta_map.count(s) == 0) {
            throw XProj(format("Could not find theta for species %d") % s);
        }
        return _theta_map.at(s);
    }
#endif

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

    inline pair<double,double> SpeciesForest::drawTruncatedIncrement(Lot::SharedPtr lot, double truncate_at) {
        // Draw from exponential distribution truncted at the specified value
        double upper_bound = truncate_at - _forest_height;
        unsigned n = (unsigned)_lineages.size();
        double incr = G::_infinity;
        double rate = 0.0;
        if (n > 1) {
            rate = G::_lambda*n;
            double cum_prob_upper_bound = 1.0 - exp(-rate*upper_bound);
            incr = -log(1.0 - cum_prob_upper_bound*lot->uniform())/rate;
        }
        return make_pair(incr, rate);
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
        Forest::_nodes.resize(2*G::_nspecies - 1);
        for (unsigned i = 0; i < G::_nspecies; i++) {
            string species_name = G::_species_names[i];
            _nodes[i]._number = (int)i;
            _nodes[i]._name = species_name;
            _nodes[i]._edge_length = 0.0;
            _nodes[i]._height = 0.0;
            Node::setSpeciesBit(_nodes[i]._species, i, /*init_to_zero_first*/true);
            _lineages.push_back(&_nodes[i]);
        }
        refreshAllPreorders();
        _forest_height = 0.0;
        _next_node_index = G::_nspecies;
        _next_node_number = G::_nspecies;
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
    
    inline void SpeciesForest::fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) {
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
    
    inline void SpeciesForest::recordCoalInfoAndClearMark() {
        double species_forest_height = _forest_height;
        unsigned nanc   = (unsigned)_mark_anc_nodes.size();
        unsigned nincr  = (unsigned)_mark_increments.size();
        unsigned nlrpos = (unsigned)_mark_left_right_pos.size();
        assert(nincr > 0);
        assert(nlrpos == nanc);
        assert(nincr == nlrpos || nincr == nlrpos + 1);
        
        if (nincr == nlrpos + 1) {
            // Remove increment corresponding to the coalescent event
            double dt = _mark_increments.top();
            species_forest_height -= dt;
            _mark_increments.pop();
            nincr = (unsigned)_mark_increments.size();
        }
                
        while (!_mark_left_right_pos.empty()) {
            //pair<unsigned, unsigned> p = _mark_left_right_pos.top();
            _mark_left_right_pos.pop();
            //unsigned left_pos  = p.first;
            //unsigned right_pos = p.second;
            
            // Get pointer to ancestral node representing the speciation event
            Node * anc = _mark_anc_nodes.top();
            assert(anc);
            //anc->_height = species_forest_height;
            _mark_anc_nodes.pop();
            
            // Get pointer to ancestral node's left child
            Node * lchild = anc->getLeftChild();
            assert(lchild);

            // Get pointer to ancestral node's right child
            Node * rchild = lchild->getRightSib();
            rchild = lchild->_right_sib;
            assert(rchild);
            
            // Update _coalinfo vector
            //_coalinfo.push_back(make_tuple(anc->_height, 0, vector<G::species_t>({lchild->_species, rchild->_species})));

            // Remove increment associated with this speciation event
            double dt = _mark_increments.top();
            species_forest_height -= dt;
            _mark_increments.pop();
        }
        
        assert(_mark_increments.empty());
        assert(_mark_anc_nodes.empty());
        assert(_mark_left_right_pos.empty());
    }
    
    inline void SpeciesForest::revertToMark() {
        unsigned nanc  = (unsigned)_mark_anc_nodes.size();
        unsigned nincr = (unsigned)_mark_increments.size();
        unsigned nlrpos = (unsigned)_mark_left_right_pos.size();
        if (_lineages.size() > 1) {
            assert(nincr > 0);
            assert(nincr == nlrpos || nincr == nlrpos + 1);
        }
        assert(nlrpos == nanc);
        
        if (nincr == nlrpos + 1) {
            // There is one more increment than the number of joins.
            // This is because either:
            // 1) _joint_estimation == true and the extra increment is from the
            //    coalescence event that ends the step; or
            // 2) _joint_estimation == false and the extra increment is from
            //    the first step in which there is only an increment and no join.
            // Either way, to revert we must first remove this extra (or only)
            // increment.
            double dt = _mark_increments.top();
            _mark_increments.pop();
            nincr = (unsigned)_mark_increments.size();
            advanceAllLineagesBy(-dt, /*mark*/false);
        }
        
        while (!_mark_left_right_pos.empty()) {
            if (_lineages.size() > 1 && !_joint_estimation) {
                // If species forests are being built conditional on
                // full gene trees, then each speciation event (except
                // the final one) is followed by an increment, so
                // reversion must remove the increment first.
                double dt = _mark_increments.top();
                _mark_increments.pop();
                advanceAllLineagesBy(-dt, /*mark*/false);
            }
            
            pair<unsigned, unsigned> p = _mark_left_right_pos.top();
            _mark_left_right_pos.pop();
            unsigned left_pos  = p.first;
            unsigned right_pos = p.second;
            
            // Get pointer to ancestral node representing the speciation event
            Node * anc = _mark_anc_nodes.top();
            _mark_anc_nodes.pop();
            assert(anc);
            assert(fabs(anc->getEdgeLength()) < 0.00001);
            
#if defined(EST_THETA)
            if (_joint_estimation) {
                // If not joint estimation, then we are integrating out
                // species-specific theta values, but if is joint estimation,
                // then each proposed ancestral species has its own theta
                // Erase theta corresponding to ancestral species
                G::species_t s = anc->getSpecies();
                
                // //temporary!
                // cerr << "erasing " << s << " from _theta_map" << endl;
                // cerr << "_theta_map.count(" << s << ") = " << // _theta_map.count(s) << endl;
                // cerr << "_theta_map before:" << endl;
                // for (const auto & z : _theta_map) {
                //     cerr << "  " << z.first << ": " << z.second << endl;
                // }
                
                unsigned n = (unsigned)_theta_map.erase(s);
                assert(n > 0);
                
                // //temporary!
                // cerr << n << " element(s) erased" << endl;
                // cerr << "_theta_map after:" << endl;
                // for (const auto & z : _theta_map) {
                //     cerr << "  " << z.first << ": " << z.second << endl;
                // }
            }
#endif

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

            if (_joint_estimation) {
                // If species forests are being built jointly with gene forests, then
                // each speciation event is preceded by an increment, so
                // reversion must remove the increment last
                double dt = _mark_increments.top();
                _mark_increments.pop();
                advanceAllLineagesBy(-dt, /*mark*/false);
            }
        }
        assert(_mark_increments.empty());
        assert(_mark_anc_nodes.empty());
        assert(_mark_left_right_pos.empty());
        _mark_coalinfo.clear();
    }
    
    inline void SpeciesForest::finalizeProposal() {
        // Empty the _prev_species_stack for all nodes in the species tree
        for (auto & nd : _nodes) {
            nd.emptyPrevSpeciesStack();
        }

        //recordCoalInfoAndClearMark();
        debugCheckMarkEmpty();
    }
    
    inline void SpeciesForest::speciationEvent(Lot::SharedPtr lot, G::species_t & left_spp, G::species_t & right_spp, G::species_t & anc_spp, bool mark) {
        //MARK: speciationEvent
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
        
#if defined(EST_THETA)
        assert(_theta_map.count(anc_spp) == 0);
        _theta_map[anc_spp] = G::inverseGammaVariate(G::_invgamma_shape, _theta_mean);
#endif
        
        // Update _lineages vector
        removeTwoAddOne(_lineages, first_node, second_node, anc_node);
        refreshAllPreorders();
        
        if (mark) {
            addCoalInfoElem(anc_node, _mark_coalinfo);
            _mark_anc_nodes.push(anc_node);
            _mark_left_right_pos.push(make_pair(left_pos, right_pos));
        }
        else {
            addCoalInfoElem(anc_node, _coalinfo);
        }
    }
        
    inline void SpeciesForest::operator=(const SpeciesForest & other) {
        Forest::operator=(other);
        _joint_estimation = other._joint_estimation;
        
#if defined(EST_THETA)
        _theta_mean = other._theta_mean;
        _theta_map = other._theta_map;
#endif

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
        
        // Dump _mark_coalinfo into coalinfo_vect
        if (!_mark_coalinfo.empty()) {
            coalinfo_vect.insert(coalinfo_vect.begin(), _mark_coalinfo.begin(), _mark_coalinfo.end());
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
