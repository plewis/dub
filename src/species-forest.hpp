#pragma once

class Particle;

namespace proj {

    class SpeciesForest : public Forest {
    
        friend class Particle;
        
        public:
        
            SpeciesForest();
             ~SpeciesForest();
            
            void simulateSpeciesTree(epoch_list_t & epochs);

            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void    createTrivialForest(bool compute_partials = false);
            bool    isSpeciesForest() const {return true;}
            
            double calcLogSpeciesTreeDensity(double lambda) const;

            void operator=(const SpeciesForest & other);

        protected:
            
            epoch_list_t & digest();
            epoch_list_t::iterator findEpochGreaterThan(epoch_list_t & epochs, double h);
            epoch_list_t::iterator reconcileEpochsStartingAt(epoch_list_t & epochs, epoch_list_t::iterator start_iter, Node::species_t species1, Node::species_t species2, Node::species_t new_species);
            void insertSpeciationEpochAt(epoch_list_t & epochs, double h, Node::species_t species1, Node::species_t species2, Node::species_t new_species);
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

    // This is an override of an abstract base class function
    inline void SpeciesForest::createTrivialForest(bool compute_partials) {
        clear();
        Forest::_nodes.resize(2*Forest::_nspecies - 1);
        for (unsigned i = 0; i < Forest::_nspecies; i++) {
            _nodes[i]._number = (int)i;
            _nodes[i]._name = Forest::_species_names[i];
            _nodes[i]._edge_length = 0.0;
            _nodes[i]._height = 0.0;
            _nodes[i]._species = {i};
            _lineages.push_back(&_nodes[i]);
        }
        _forest_height = 0.0;
        _next_node_index = Forest::_nspecies;
        _next_node_number = Forest::_nspecies;
        _prev_log_likelihood = 0.0;
        _prev_log_coalescent_likelihood = 0.0;
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
        assert(internal_heights.size() == Forest::_nspecies - 1);
        
        // Sort heights
        sort(internal_heights.begin(), internal_heights.end());
        
        double log_prob_density = 0.0;
        unsigned n = Forest::_nspecies;
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
    
    inline epoch_list_t & SpeciesForest::digest() {
        // Populate _epochs based on species tree
        _epochs.clear();
        
        // Assumes _species_forest is a complete tree
        assert(_lineages.size() == 1);
        refreshAllPreorders();
        
        // Record speciation events in the species forest using a post-order
        // traversal, calculating heights of nodes along the way
        Node::ptr_vect_t & preorder = _preorders[0];
        for (Node * nd : boost::adaptors::reverse(preorder)) {
            Node * lchild = nd->getLeftChild();
            if (lchild) {
                // nd is an internal
                Node * rchild = lchild->getRightSib();
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
                
                // Create entry in _epochs for this speciation event
                //unsigned lnum = lchild->getNumber();
                //unsigned rnum = rchild->getNumber();
                //unsigned ndnum = nd->getNumber();
                nd->setSpeciesToUnion(lchild->getSpecies(), rchild->getSpecies());

                Epoch sepoch(Epoch::speciation_epoch, h);
                sepoch._left_species  = lchild->getSpecies();
                sepoch._right_species = rchild->getSpecies();
                sepoch._anc_species   = nd->getSpecies();
                //sepoch._left_species = lnum;
                //sepoch._right_species = rnum;
                //sepoch._anc_species = ndnum;
                pushBackEpoch(_epochs, sepoch);
                //_epochs.push_back(new SpeciationEpoch(h, lnum, rnum, ndnum));
            }
            else {
                // nd is a leaf
                nd->setHeight(0.0);
                nd->setSpeciesFromUnsigned((unsigned)nd->_number);
            }
        }

        // Sort epochs by height
        _epochs.sort(epochLess);
        
        return _epochs;
    }
    
    //inline void SpeciesForest::createFromTree(Tree::SharedPtr t) {
    //    // Assumes leaf numbers in t are indices into the GeneForest::_species_names vector
    //    createTrivialForest();
    //    for (Node * nd : boost::adaptors::reverse(t->_preorder)) {
    //        unsigned num = nd->getNumber();
    //        Node * lchild = nd->getLeftChild();
    //        unsigned lnum = lchild->getNumber();
    //        if (lchild) {
    //            // internal
    //            Node * rchild = lchild->getRightSib();
    //            unsigned rnum = rchild->getNumber();
    //            _nodes[num]._left_child = &_nodes[lnum];
    //            _nodes[lnum]._right_sib = &_nodes[rnum];
    //            _nodes[lnum]._parent = &_nodes[num];
    //            _nodes[rnum]._parent = &_nodes[num];
    //            _nodes[lnum]._edge_length = lchild->getEdgeLength();
    //            _nodes[rnum]._edge_length = rchild->getEdgeLength();
    //        }
    //        else {
    //            // leaf
    //            assert(nd->getName() == GeneForest::_species_names[num]);
    //        }
    //    }
    //}
    
    inline void SpeciesForest::simulateSpeciesTree(epoch_list_t & epochs) {
        createTrivialForest();
        unsigned nsteps = Forest::_nspecies - 1;
        for (unsigned i = 0; i < nsteps; ++i) {
            unsigned nlineages = (unsigned)_lineages.size();
            assert(nlineages > 1);
            
            // Choose waiting time until the next speciation event
            double u = rng.uniform();
            double r = Forest::_lambda*nlineages;
            double t = -log(1.0 - u)/r;
            //cout << str(format("  t = %g (u = %g, nlineages = %d, r = %g)\n") % t % u % nlineages % r);
                
            // Increment height of forest
            _forest_height += t;

            // Update edge lengths
            for (auto nd : _lineages)
                nd->_edge_length += t;

            // Choose two lineages to join
            auto chosen_pair = rng.nchoose2(nlineages);
            
            Node * first_node  = _lineages[chosen_pair.first];
            Node * second_node = _lineages[chosen_pair.second];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            anc_node->setSpeciesToUnion(first_node->getSpecies(), second_node->getSpecies());
            
            // Update lineage vector
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            
            Epoch epoch(Epoch::speciation_epoch, _forest_height);
            epoch._left_species  = first_node->getSpecies();
            epoch._right_species = second_node->getSpecies();
            epoch._anc_species   = anc_node->getSpecies();
            pushBackEpoch(epochs, epoch);
        }
        
        heightsInternalsPreorders();
    }
    
    inline epoch_list_t::iterator SpeciesForest::findEpochGreaterThan(epoch_list_t & epochs, double h) {
        // Find iterator pointing to first epoch with height > h
        auto it = epochs.begin();
        for (; it != epochs.end(); ++it) {
            if (it->_height < h)
                continue;
            else
                break;
        }
        return it;
    }
    
    inline epoch_list_t::iterator SpeciesForest::reconcileEpochsStartingAt(epoch_list_t & epochs, epoch_list_t::iterator start_iter, Node::species_t species1, Node::species_t species2, Node::species_t new_species) {
        auto it = start_iter;
        
        // Fix the subsequent coalescent epoch elements, replacing species1 and species2
        // everywhere with new_species until encountering either the end or the next speciation epoch
        for (; it != epochs.end(); ++it) {
            if (it->isSpeciationEpoch())
                continue;
                         
            // Let ss be an alias for the epoch's species
            Node::species_t & ss = it->_species;

            // Remove all elements of species1 from the set if found
            unsigned n1 = 0;
            for (auto it1 = species1.begin(); it1 != species1.end(); ++it1) {
                auto tmp = ss.find(*it1);
                if (tmp != ss.end()) {
                    ss.erase(tmp);
                    n1++;
                }
            }
                
            // Check to make sure that if one element of species1 was found then they all were found
            assert(n1 == 0 || n1 == species1.size());

            // Remove all elements of species2 from the set if found
            unsigned n2 = 0;
            for (auto it2 = species2.begin(); it2 != species2.end(); ++it2) {
                auto tmp = ss.find(*it2);
                if (tmp != ss.end()) {
                    ss.erase(tmp);
                    n2++;
                }
            }
                
            // Check to make sure that if one element of species2 was found then they all were found
            assert(n2 == 0 || n2 == species2.size());

            // If either species1 or species1 was found in species (and removed), need
            // to add new_species to the set
            if (n1 > 0 || n2 > 0) {
                ss.insert(new_species.begin(), new_species.end());
            }
            
            it->_lineage_counts[new_species] += it->_lineage_counts[species1];
            it->_lineage_counts[new_species] += it->_lineage_counts[species2];
            it->_lineage_counts[species1] = 0;
            it->_lineage_counts[species2] = 0;
        }
        
        return it;
    }
    
    inline void SpeciesForest::insertSpeciationEpochAt(epoch_list_t & epochs, double h, Node::species_t species1, Node::species_t species2, Node::species_t new_species) {
        // Find iterator pointing to first epoch with height > h
        auto it = findEpochGreaterThan(epochs, h);
        
        // Insert SpeciationEpoch just before it
        Epoch spp_epoch(Epoch::speciation_epoch, h);
        spp_epoch._left_species  = species1;
        spp_epoch._right_species = species2;
        spp_epoch._anc_species   = new_species;
        it = epochs.insert(it, spp_epoch);
        
        it = reconcileEpochsStartingAt(epochs, ++it, species1, species2, new_species);
        assert(it == epochs.end());
        debugShowEpochs(epochs);
    }
    
    inline double SpeciesForest::advanceSpeciesForest(unsigned particle, unsigned step) {
        // Step 0
        //    <------ gene tree 0 ---->      <------ gene tree 1 ---->
        //
        //    0   0   1   1   1   2   2      0   0   1   1   1   2   2  0   _
        //    |   |   |   |   |   |   |      |   |   |   |   |   |   |  1   |
        //    |   |   |   |   |   |   |      +-0-+   |   |   |   |   |  2   | first speciation event must occur
        //    +-0-+   |   |   |   |   |        |     +-1-+   |   |   |  3   | in this interval, regardless of
        //      |     |   |   |   |   |        |       |     |   |   |  4   |
        //      |     |   +-1-+   |   |        |       |     |   |   |  5   | the species joined
        //      |     |     |     +-2-+        |       |     +1,2+   |  6   -
        //      |     +--1--+       |          |       |       |     |  7
        //      |        |          |          |       +--1,2--+     |  8
        //      |        |          |          |           |         |  9
        //      |        |          |          |           |         | 10
        //      +---0,1--+          |          |           |         | 11
        //           |              |          +-- 0,1,2---+         | 12
        //           |              |                |               | 13
        //           |              |                |               | 14
        //           +-----0,1,2----+                |               | 15
        //                                           |               | 16
        //                                           +-----0,1,2-----+ 17
        //  _species_tree._epochs:
        //     type   height  gene    species
        //       I       -1      0              <-- these entries are used to
        //       I       -1      1              <-- initialize lineage counts
        //       C        2      1    0
        //       C        3      0    0
        //       C        3      1    1
        //       C        5      0    1
        //       C        6      0    2
        //       C        6      1    1,2       <-- this entry sets bound
        //       C        7      0    1
        //       C        8      1    1,2
        //       C       11      0    0,1
        //       C       12      1    0,1,2
        //       C       15      0    0,1,2
        //       C       17      1    0,1,2
        //
        //    Increment chosen at height 4.5, but no join performed first step.
        //    Coalescent likelihood returned after step 0 considers all coalescent
        //    events with heights > 4.5 to be in a single population.
        //
        // Step 1
        //    <------ gene tree 0 ---->      <------ gene tree 1 ---->
        //
        //    0   0   1   1   1   2   2      0   0   1   1   1   2   2  0
        //    |   |   |   |   |   |   |      |   |   |   |   |   |   |  1
        //    |   |   |   |   |   |   |      +-0-+   |   |   |   |   |  2
        //    +-0-+   |   |   |   |   |        |     +-1-+   |   |   |  3
        //      |     |   |   |   |   |        |       |     |   |   |  4
        //    --------------------------------------------------------- 4.5 -
        //      |     |   +-01+   |   |        |       |     |   |   |  5   | Second speciation event joining
        //      |     |     |     +-2-+        |       |     01,2+   |  6   - species 01 and 2 must occur in
        //      |     +--01-+       |          |       |       |     |  7     this interval
        //      |        |          |          |       +-01,2--+     |  8
        //      |        |          |          |           |         |  9
        //      |        |          |          |           |         | 10
        //      +----01--+          |          |           |         | 11
        //           |              |          +----01,2---+         | 12
        //           |              |                |               | 13
        //           |              |                |               | 14
        //           +------01,2-----+               |               | 15
        //                                           |               | 16
        //                                           +------01,2-----+ 17
        //
        //  _species_tree._epochs:
        //     type   height  gene    species
        //       I       -1      0    -        <-- these entries are used to
        //       I       -1      1    -        <-- initialize lineage counts
        //       C        2      1    0
        //       C        3      0    0
        //       C        3      1    1
        //       S       4.5     -    -        <-- species 0 and 1 joined to create species 01
        //       C        5      0    01
        //       C        6      0    2
        //       C        6      1    01,2     <-- this entry sets bound
        //       C        7      0    01
        //       C        8      1    01,2
        //       C       11      0    01,1
        //       C       12      1    01,2
        //       C       15      0    01,2
        //       C       17      1    01,2
        //
        //    Species 0 and 1 joined to create species 01.
        //    All epochs downstream replace species 0 and species 1 with 01.
        //    New increment chosen taking us to height 5.5.
        //    Coalescent likelihood returned after step 1 considers all coalescent
        //    events with heights > 5.5 to be in a single population.
        //
        // Step 2
        //    <------ gene tree 0 ---->      <------ gene tree 1 ---->
        //
        //    0   0   1   1   1   2   2      0   0   1   1   1   2   2  0
        //    |   |   |   |   |   |   |      |   |   |   |   |   |   |  1
        //    |   |   |   |   |   |   |      +-0-+   |   |   |   |   |  2
        //    +-0-+   |   |   |   |   |        |     +-1-+   |   |   |  3
        //      |     |   |   |   |   |        |       |     |   |   |  4
        //    --------------------------------------------------------- 4.5
        //      |     |   +-01+   |   |        |       |     |   |   |  5
        //    --------------------------------------------------------- 5.5
        //      |     |     |     +012+        |       |     +012+   |  6
        //      |     +-012-+       |          |       |       |     |  7
        //      |        |          |          |       +--012--+     |  8
        //      |        |          |          |           |         |  9
        //      |        |          |          |           |         | 10
        //      +---012--+          |          |           |         | 11
        //           |              |          +----012----+         | 12
        //           |              |                |               | 13
        //           |              |                |               | 14
        //           +------012-----+                |               | 15
        //                                           |               | 16
        //                                           +------012------+ 17
        //
        //  _species_tree._epochs:
        //     type   height  gene      species
        //       I       -1      0      -         <-- these entries are used to
        //       I       -1      1      -         <-- initialize lineage counts
        //       C        2      1      0
        //       C        3      0      0
        //       C        3      1      1
        //       S       4.5     -      -         <-- species 0 and 1 joined to create species 01
        //       C        5      0      01
        //       S       5.5     -      -         <-- species 01 and 2 joined to create species 012
        //       C        6      0      012
        //       C        6      1      012
        //       C        7      0      012
        //       C        8      1      012
        //       C       11      0      012
        //       C       12      1      012
        //       C       15      0      012
        //       C       17      1      012
        //
        //    Species 01 and 2 are joined at height 5.5 to create species 012.
        //    All epochs downstream replace species 01 and species 2 with 012.
        //    No more increments or joins are possible, so finish calculating the
        //    coalescent likelihood by considering the remaining coalescent epochs.
        
        // This is the return value
        double log_weight = 0.0;
        
        unsigned nlineages = (unsigned)_lineages.size();
        assert(nlineages > 1);
        
        if (step > 0) {
            pair<unsigned,unsigned> chosen_pair = rng.nchoose2(nlineages);
            Node * first_node  = _lineages[chosen_pair.first];
            Node * second_node = _lineages[chosen_pair.second];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            anc_node->setSpeciesToUnion(first_node->getSpecies(), second_node->getSpecies());
            
            // Insert a speciation epoch and revise subsequent coalescent epochs to replace the two
            // old species with the new species
            Node::species_t oldspp1 = first_node->getSpecies();
            Node::species_t oldspp2 = second_node->getSpecies();
            Node::species_t newspp  = anc_node->getSpecies();
            insertSpeciationEpochAt(_epochs, _forest_height, oldspp1, oldspp2, newspp);

            // Update lineage vector
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            --nlineages;
        }

        // Create a set of all species that currently (looking backward in time)
        // exist in species tree lineages
        set<Node::species_t> current_species;
        for (auto lit = _lineages.begin(); lit != _lineages.end(); lit++) {
            Node * nd = *lit;
            current_species.insert(nd->_species);
        }
        assert(current_species.size() == nlineages);
        
        // Find next coalescent event (over all gene trees) involving more than one species
        // That sets an upper bound on the time of the next speciation event, regardless of
        // which two species ultimately end up being joined.
        double maxT = calcMaxT(_epochs, _forest_height, current_species);
        double max_waiting_time = maxT - _forest_height;
            
        // Choose waiting time until the next speciation event (conditional on maxT)
        double u = rng.uniform();
        double r = Forest::_lambda*nlineages;
        double t = -log(1.0 - u*(1.0 - exp(-r*max_waiting_time)))/r;
            
        // Increment height of forest
        advanceAllLineagesBy(t);
        
        if (nlineages == 2) {
            // Down to final 2 lineages; go ahead and join them to complete the species tree
            Node * first_node  = _lineages[0];
            Node * second_node = _lineages[1];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            anc_node->setSpeciesToUnion(first_node->getSpecies(), second_node->getSpecies());
            
            // Insert a speciation epoch and revise subsequent coalescent epochs to replace the two
            // old species with the new species
            Node::species_t oldspp1 = first_node->getSpecies();
            Node::species_t oldspp2 = second_node->getSpecies();
            Node::species_t newspp  = anc_node->getSpecies();
            insertSpeciationEpochAt(_epochs, _forest_height, oldspp1, oldspp2, newspp);

            // Update lineage vector
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            --nlineages;
        }
                
        // Compute coalescent likelihood
        double log_coalescent_likelihood = 0.0;
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            log_coalescent_likelihood += calcLogCoalescentLikelihood(_epochs, g);
        }
        
        log_weight = log_coalescent_likelihood - _prev_log_coalescent_likelihood;
        _prev_log_coalescent_likelihood = log_coalescent_likelihood;

        debugShowEpochs(_epochs);
        debugShow(format("Species forest in particle %d after step %d (height = %.5f)\n    %s") % particle % step % _forest_height % makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
        
        return log_weight;
    }
    
    inline void SpeciesForest::operator=(const SpeciesForest & other) {
        Forest::operator=(other);
        
#if !defined(NDEBUG)
        // Sanity check: make sure that _partials do not exist for either this nor other
        for (unsigned i = 0; i < _nodes.size(); i++) {
            assert((_nodes[i]._partial == nullptr && other._nodes[i]._partial == nullptr));
        }
#endif
    }

}
