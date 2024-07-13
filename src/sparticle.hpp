#pragma once

extern void output(string msg, proj::G::verbosity_t verb);
extern void output(format & fmt, proj::G::verbosity_t level);
extern proj::Lot::SharedPtr       rng;

namespace proj {

    class Bundle;

    // Species tree particle
    class SParticle : public Particle {
        friend class Bundle;

        public:
            typedef map<G::species_t, double> theta_map_t;
            
            SParticle() {_is_species_tree = true; createTrivialForest();}
            ~SParticle() {}
            
            // Getters and setters
            double getThetaMean() const {return _theta_mean;}
            void setThetaMean(double theta_mean) {_theta_mean = theta_mean; }
            
            theta_map_t & getThetaMap() {return _theta_map;}
            const theta_map_t & getThetaMapConst() const {return _theta_map;}
            
            bool getNewickTheta() const {return _newick_theta;}
            void setNewickTheta(bool save_thetas_in_newick) const {_newick_theta = save_thetas_in_newick;}
            
            G::join_info_t heightOfNextNodeAbove(double h, const vector<G::join_info_t> & spec_info) const;
            
            void        sanityCheck(unsigned step, unsigned bundle) const;
            Node *      joinThenIncrement(Lot::SharedPtr lot);
            void        trimToHeight(double h0, Lot::SharedPtr lot);
            
            double      drawThetaForSpecies(Lot::SharedPtr lot);
            void        initThetaMap(Lot::SharedPtr lot);

            // Overrides of base class virtual functions
            virtual void    recordJoinInfo(vector<G::join_info_t> & join_info) const;
            virtual string  makeNewick(unsigned precision, bool use_names) const;

            // Overrides of base class abstract virtual functions
            pair<double, double> drawIncrement(Lot::SharedPtr lot);
            Node * joinRandomPair(Lot::SharedPtr lot);
            double calcLogLikelihood() const;
            void createTrivialForest();
            string info() const {
                return str(format("  S-%d (h=%g)") % _index % getHeight());
            }
            
        private:
        
            double          _theta_mean;
            theta_map_t     _theta_map;
            
            mutable bool    _newick_theta;
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
    
    inline Node * SParticle::joinRandomPair(Lot::SharedPtr lot) {
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
        pair<unsigned, unsigned> lr = lot->nchoose2((unsigned)_lineages.size());
        Node * lchild = _lineages[lr.first];
        Node * rchild = _lineages[lr.second];
        
        // Connect chosen lineages to ancestor
        makeAnc(anc, lchild, rchild);
        
        output(format("joining species %d and %d to form %d\n") % lchild->_species % rchild->_species % anc->_species, G::VDEBUG);
        
        return anc;
    }

    inline pair<double, double> SParticle::drawIncrement(Lot::SharedPtr lot) {
        // Choose an increment at random from the Exp(n*lambda) prior
        // where n is the current number of lineages
        unsigned n = (unsigned)_lineages.size();
        double rate = G::_lambda*n;
        double scale = 1.0/rate;
        double incr = lot->gamma(1.0, scale);
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

    inline Node * SParticle::joinThenIncrement(Lot::SharedPtr lot) {
        Node * anc = joinRandomPair(lot);
        auto incr_rate = drawIncrement(lot);
        extendAllLineagesBy(incr_rate.first);
        return anc;
    }

    inline void SParticle::trimToHeight(double h0, Lot::SharedPtr lot) {
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
        pair<double,double> incr_rate = drawIncrement(lot);
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
    }

    inline void SParticle::recordJoinInfo(vector<G::join_info_t> & join_info) const {
        Particle::recordJoinInfo(join_info);
        
        if (_lineages.size() == 1) {
            // If species tree is complete, add final element at infinity
            // with all zeros for species. This indicates that anything goes
            // with respect to coalescence because we are in the ancestral species
            join_info.push_back({
                G::_infinity,      /*height*/
                true,              /*species tree*/
                G::_species_zero,  /*left child species*/
                G::_species_zero   /*right child species*/
            });
        }
        else {
            // If species tree is not complete, add final element at _height
            // with all zeros for species representing the point at which
            // the next join will occur
            join_info.push_back({
                _height,          /*height*/
                true,             /*species tree*/
                G::_species_zero, /*left child species*/
                G::_species_zero  /*right child species*/
            });
        }
    }
    
    inline G::join_info_t SParticle::heightOfNextNodeAbove(double h, const vector<G::join_info_t> & spec_info) const {
        // If there are no elements in spec_info, then species tree
        // has not yet experienced a join, in which case the height
        // returned should be _height, which is the height at which
        // the next join will occur when it happens. If the species
        // tree is complete, _height = G::_infinity.
        G::join_info_t next_speciation = make_tuple(
            _height,           /*height*/
            true,              /*is species tree*/
            G::_species_zero,  /*left child species*/
            G::_species_zero   /*right child species*/
        );
        for (auto & si : spec_info) {
            double ht = get<0>(si);
            if (ht > h) {
                next_speciation = si;
                break;
            }
        }
        
        return next_speciation;
    }

    inline string SParticle::makeNewick(unsigned precision, bool use_names) const {
        // Place basal polytomy (if there is one) at a height 10% greater than
        // _height
        double basal_polytomy_height = 0.0; //_height*(1.1);
        bool is_complete_tree = (bool)(_lineages.size() == 1);
        
        const format name_theta_edgelen( str(format("%%s[&theta=%%.%df]:%%.%df") % precision % precision) );
        const format name_edgelen( str(format("%%s:%%.%df") % precision) );
        
        const format name_theta( str(format("%%s[&theta=%%.%df]") % precision) );
        const format name_only("%%s");
        
        const format number_theta_edgelen( str(format("%%d[&theta=%%.%df]:%%.%df") % precision % precision) );
        const format number_edgelen( str(format("%%d:%%.%df") % precision) );
        
        const format number_theta( str(format("%%d[&theta=%%.%df]") % precision) );
        const format number_only("%%d");
        
        const format theta_edgelen( str(format(")[&theta=%%.%df]:%%.%df") % precision % precision) );
        const format edgelen_only( str(format("):%%.%df") % precision) );
        
        const format theta_only( str(format(")[&theta=%%.%df]") % precision) );
        const format right_paren(")");
        
        // Build preorders vector
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        vector<string> subtree_newicks;
        for (auto & preorder : preorders) {
            string subtree_newick;
            stack<Node *> node_stack;
            for (auto nd : preorder) {
                if (nd->_left_child) {
                    subtree_newick += "(";
                    node_stack.push(nd);
                }
                else {
                    double edge_length = nd->_edge_length;
                    if (!nd->_parent) {
                        // Subtree consists of just this one leaf node
                        edge_length += basal_polytomy_height;
                    }
                    if (use_names) {
                        if (precision > 0) {
                            if (_newick_theta)
                                subtree_newick += str(format(name_theta_edgelen) % nd->_name % _theta_map.at(nd->_species) % edge_length);
                            else
                                subtree_newick += str(format(name_edgelen) % nd->_name % edge_length);
                        }
                        else {
                            if (_newick_theta)
                                subtree_newick += str(format(name_theta) % nd->_name % _theta_map.at(nd->_species));
                            else
                                subtree_newick += str(format(name_only) % nd->_name );
                        }
                    } else {
                        if (precision > 0) {
                            if (_newick_theta)
                                subtree_newick += str(format(number_theta_edgelen) % (nd->_number + 1) % _theta_map.at(nd->_species) % edge_length);
                            else
                                subtree_newick += str(format(number_edgelen) % (nd->_number + 1) % edge_length);
                        }
                        else {
                            if (_newick_theta) {
                                subtree_newick += str(format(number_theta) % (nd->_number + 1) % _theta_map.at(nd->_species));
                            }
                            else {
                                subtree_newick += str(format(number_only) % (nd->_number + 1));
                            }
                        }
                    }
                    if (nd->_right_sib)
                        subtree_newick += ",";
                    else if (nd->_parent) {
                        Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                        double popped_edge_length = popped->_edge_length;
                        while (popped && !popped->_right_sib) {
                            node_stack.pop();
                            if (node_stack.empty()) {
                                if (is_complete_tree) {
                                    subtree_newick += ")";
                                }
                                else {
                                    // This is the root of one of several subtrees, so
                                    // it is important to preserve its edge length
                                    popped_edge_length += basal_polytomy_height;
                                    if (precision > 0) {
                                        if (_newick_theta)
                                            subtree_newick += str(format(theta_edgelen) % _theta_map.at(nd->_species) % popped_edge_length);
                                        else
                                            subtree_newick += str(format(edgelen_only) %  popped_edge_length);
                                    }
                                    else {
                                        if (_newick_theta)
                                            subtree_newick += str(format(theta_only) % _theta_map.at(nd->_species));
                                        else
                                            subtree_newick += str(format(right_paren));
                                    }
                                }
                                popped = 0;
                            }
                            else {
                                if (precision > 0) {
                                    if (_newick_theta)
                                        subtree_newick += str(format(theta_edgelen) % _theta_map.at(nd->_species) % popped_edge_length);
                                    else
                                        subtree_newick += str(format(edgelen_only) % popped_edge_length);
                                }
                                else {
                                    if (_newick_theta)
                                        subtree_newick += str(format(theta_only) % _theta_map.at(nd->_species));
                                    else
                                        subtree_newick += str(format(right_paren));
                                }
                                popped = node_stack.top();
                                popped_edge_length = popped->_edge_length;
                            }
                        }
                        if (popped && popped->_right_sib) {
                            node_stack.pop();
                            if (precision > 0) {
                                if (_newick_theta)
                                    subtree_newick += str(format(theta_edgelen) % _theta_map.at(nd->_species) % popped_edge_length);
                                else
                                    subtree_newick += str(format(edgelen_only) % popped_edge_length);
                                subtree_newick += ",";
                            }
                            else {
                                if (_newick_theta)
                                    subtree_newick += str(format(theta_only) % _theta_map.at(nd->_species));
                                else
                                    subtree_newick += str(format(right_paren));
                                subtree_newick += ",";
                            }
                        }
                    }
                }
            }
            subtree_newicks.push_back(subtree_newick);
        }
        
        string newick;
        if (is_complete_tree)
            newick = subtree_newicks[0];
        else {
            string separator = str(format(":%.5f,") % basal_polytomy_height);
            string insides = boost::join(subtree_newicks, ",");
            newick = str(format("(%s)") % insides);
        }

        return newick;
    }

    inline double SParticle::drawThetaForSpecies(Lot::SharedPtr lot) {
        double theta = 0.0;
        if (G::_fixed_theta) {
            theta = _theta_mean;
        }
        else {
            double b = (G::_invgamma_shape - 1.0)*_theta_mean;
            theta = G::inverseGammaVariate(G::_invgamma_shape, b, lot);
        }
        return theta;
    }
    
    inline void SParticle::initThetaMap(Lot::SharedPtr lot) {
        _theta_map.clear();
        for (unsigned i = 0; i < G::_nspecies; i++) {
            G::species_t s;
            G::setSpeciesBit(s, i, true);
            _theta_map[s] = drawThetaForSpecies(lot);
        }
    }
    
}

