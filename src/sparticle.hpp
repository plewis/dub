#pragma once

extern void output(string msg, proj::G::verbosity_t verb);
extern void output(format & fmt, proj::G::verbosity_t level);
// extern proj::PartialStore         ps;
// extern proj::StopWatch            stopwatch;
extern proj::Lot::SharedPtr       rng;
// extern proj::Partition::SharedPtr partition;
// extern proj::Data::SharedPtr      data;

namespace proj {

    class Bundle;

    // Species tree particle
    class SParticle : public Particle {
        friend class Bundle;

        public:
            SParticle() {_is_species_tree = true; createTrivialForest();}
            ~SParticle() {}
            
            G::join_info_t heightOfNextNodeAbove(double h, const vector<G::join_info_t> & spec_info) const;
            
            void sanityCheck(unsigned step, unsigned bundle) const;
            void joinThenIncrement();
            void trimToHeight(double h0);
            
            // Overrides of base class virtual functions
            virtual void recordJoinInfo(vector<G::join_info_t> & join_info) const;

            // Overrides of base class abstract virtual functions
            pair<double, double> drawIncrement();
            void joinRandomPair();
            double calcLogLikelihood() const;
            void createTrivialForest();
            string info() const {
                return str(format("  S-%d (h=%g)") % _index % getHeight());
            }
    };
            
    inline double SParticle::calcLogLikelihood() const {
        // Returns log coalescent likelihood
        return 0.0;
    }
    
    inline void SParticle::createTrivialForest() {
        _nodes.clear();
        _nodes.resize(2*G::_nspecies - 1);
        _lineages.clear();
        _lineages.resize(G::_nspecies);
        for (unsigned i = 0; i < G::_nspecies; i++) {
            _lineages[i] = &_nodes[i];
            _lineages[i]->setName(G::_species_names[i]);
            _lineages[i]->setNumber(i);
            G::setSpeciesBit(_lineages[i]->getSpecies(), i, /*init_to_zero_first*/true);
        }
        _next_node_number = G::_nspecies;
        _returned_node_numbers.clear();
        _nleaves = G::_nspecies;
        refreshSplits();
    }
    
    inline void SParticle::joinRandomPair() {
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
        anc->setEdgeLength(0.0);
        
        // Choose two existing lineages at random to join
        pair<unsigned, unsigned> lr = rng->nchoose2((unsigned)_lineages.size());
        Node * lchild = _lineages[lr.first];
        Node * rchild = _lineages[lr.second];
        
        // Connect chosen lineages to ancestor
        makeAnc(anc, lchild, rchild);
        
        output(format("joining species %d and %d to form %d\n") % lchild->_species % rchild->_species % anc->_species, G::VTEMP);
    }

    inline pair<double, double> SParticle::drawIncrement() {
        // Choose an increment at random from the Exp(n*lambda) prior
        // where n is the current number of lineages
        unsigned n = (unsigned)_lineages.size();
        double rate = G::_lambda*n;
        double scale = 1.0/rate;
        double incr = rng->gamma(1.0, scale);
        return make_pair(incr, rate);
    }

    inline void SParticle::sanityCheck(unsigned step, unsigned bundle) const {
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        // Check whether all lineages are at the same height
        double href = -1.0;
        for (auto pre : preorders) {
            for (auto nd : boost::adaptors::reverse(pre)) {
                if (href < 0 && !nd->_left_child) {
                    // nd is first leaf encountered
                    href = 0.0;
                    while (nd->_parent) {
                        href += nd->_edge_length;
                        nd = nd->_parent;
                    }
                    
                    // Need to add root's edge length in case
                    // species tree not yet complete
                    href += nd->_edge_length;
                }
                else if (!nd->_left_child) {
                    // nd is a leaf, but not the first leaf encountered
                    double h = 0.0;
                    while (nd->_parent) {
                        h += nd->_edge_length;
                        nd = nd->_parent;
                    }

                    // Need to add root's edge length in case
                    // species tree not yet complete
                    h += nd->_edge_length;

                    if (fabs(href - h) > G::_epsilon) {
                        throw XProj(format("In step %d, bundle %d, root-to-leaf distances in species tree not all the same: %g vs %g") % step % bundle % href % h);
                    }
                }
            }   // postorder traversal
        }   // lineage loop
        
        // Check whether species tree height equals height just calculated
        if (fabs(_height - href) > G::_epsilon) {
            throw XProj(format("In step %d, bundle %d, height of species tree (%g) does not match calculated height (%g)") % step % bundle % _height % href);
        }
    }

    inline void SParticle::joinThenIncrement() {
        joinRandomPair();
        auto incr_rate = drawIncrement();
        extendAllLineagesBy(incr_rate.first);
    }

    inline void SParticle::trimToHeight(double h0) {
        // //temporary!
        // unsigned step = G::_step;
        // ofstream tmpf1("before-trimming.tre", ios::out | ios::app);
        // tmpf1 << "  tree t" << step << " = [&R]" << makeNewick(9, true) << ";\n";
        // tmpf1.close();

        //  Assume there are 3 lineages in the species tree
        //  After filtering gene trees from all loci, the
        //  height of the deepest coalescent event in any
        //  gene tree is shown as a dotted line
        //
        //   |   |   |   |   |  |   |   |
        //   |   |   |   |   |  |   +-+-+
        //   +-+-+   |   |   |  |     |
        //   ..|.....|...|...|..|.....|.. <-- highest gene tree
        //     |     |   +-+-+  |     |
        //     +--+--+     |    |     |
        //        |        |    +--+--+
        //        |        |       |
        //
        //  Trim back species tree to dotted line.
        //
        //   |   |   |   |   |  |   |   |
        //   |   |   |   |   |  |   +-+-+
        //   +-+-+   |   |   |  |     |
        //   ..|.....|...|...|..|.....|.. <-- highest gene tree
        //
        //  Add new increment to species tree (but no join).
        //
        //   |   |   |   |   |  |   |   |
        //   |   |   |   |   |  |   +-+-+
        //   +-+-+   |   |   |  |     |
        //   ..|.....|...|...|..|.....|.. <-- highest gene tree
        //     |     |   |   |  |     |
        //     |     |   |   |  |     |
        //
        assert(h0 > 0.0);
        assert(_lineages.size() > 0);
        double new_height = _height;
        bool done = false;
        while (!done) {
            // Determine highest node in _lineages
            Node * highest = _lineages[0];
            for (auto nd : _lineages) {
                if (nd->_height > highest->_height)
                    highest = nd;
            }
            
            done = true;
            if (highest->_height > h0) {
                // Detach children from highest
                assert(highest->_left_child);
                assert(highest->_left_child->_right_sib);
                assert(!highest->_right_sib);
                Node * lchild = highest->_left_child;
                Node * rchild = lchild->_right_sib;
                lchild->_parent = nullptr;
                rchild->_parent = nullptr;
                lchild->_right_sib = nullptr;
                highest->_left_child = nullptr;
                                
                // Remove highest from _lineages
                _returned_node_numbers.push_back(highest->_number);
                highest->_edge_length = 0.0;
                highest->_number = -1;  // flag node as being unused
                auto it = find(_lineages.begin(), _lineages.end(), highest);
                assert(it != _lineages.end());
                _lineages.erase(it);
                                
                // Add lchild and rchild to lineages
                _lineages.push_back(lchild);
                _lineages.push_back(rchild);
                
                done = false;
            }
        }
        
        // There are no longer any nodes in the species tree higher than h0
        // Choose a new species tree increment and adjust edge lengths of all
        // nodes in _lineages so that height of the species tree is h0 plus
        // the new increment.
        pair<double,double> incr_rate = drawIncrement();
        double incr = incr_rate.first;
        for (auto nd : _lineages) {
            //  a   b   c   d   e  f   g   h
            //  |   |   |   |   |  |   |   |
            //  |   |   |   |   |  |   +-j-+
            //  +-i-+   |   |   |  |     |
            //  ..|.....|...|...|..|.....|.. <-- highest gene tree (h0)
            //    |     |   |   |  |     |
            //    |     |          |     |
            //                     |     |
            //
            //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ <-- incr
            //
            double new_edge_length = h0 - nd->_height;
            new_edge_length += incr;
            nd->_edge_length = new_edge_length;
        }
        
        // Adjust species tree height
        _height = h0 + incr;
                
        // //temporary!
        // ofstream tmpf2("after-trimming.tre", ios::out | ios::app);
        // tmpf2 << "  tree t" << step << " = [&R]" << makeNewick(9, true) << ";\n";
        // tmpf2.close();
    }

    inline void SParticle::recordJoinInfo(vector<G::join_info_t> & join_info) const {
        Particle::recordJoinInfo(join_info);
        
        if (_lineages.size() == 1) {
            // If species tree is complete, add final element at infinity
            // with all zeros for species. This indicates that anything goes
            // with respect to coalescence because we are in the ancestral species
            join_info.push_back({G::_infinity, true, G::_species_zero, G::_species_zero});
        }
        else {
            // If species tree is not complete, add final element at _height
            // with all zeros for species representing the point at which
            // the next join will occur
            join_info.push_back({_height, true, G::_species_zero, G::_species_zero});
        }
    }
    
    inline G::join_info_t SParticle::heightOfNextNodeAbove(double h, const vector<G::join_info_t> & spec_info) const {
        // If there are no elements in spec_info, then species tree has not yet experienced a join,
        // in which case the height returned should be _height, which is the height at which the
        // next join will occur when it happens.
        G::join_info_t next_speciation = make_tuple(_height, true, G::_species_zero, G::_species_zero);
        for (auto & si : spec_info) {
            double ht = get<0>(si);
            if (ht > h) {
                next_speciation = si;
                break;
            }
        }
        return next_speciation;
    }

}

