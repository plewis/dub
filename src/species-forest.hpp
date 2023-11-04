#pragma once

class Particle;

namespace proj {

    class SpeciesForest : public Forest {
    
        friend class Particle;
        
        public:
        
            SpeciesForest();
             ~SpeciesForest();
            
            void speciationEvent(SMCGlobal::species_t & left, SMCGlobal::species_t & right, SMCGlobal::species_t & anc);
            //void simulateSpeciesTree();

            // Overrides of base class functions
            void clear();
            
            // Overrides of abstract base class functions
            void    createTrivialForest(bool compute_partials = false);
            bool    isSpeciesForest() const {return true;}

#if defined(NEWWAY)
            double drawIncrement();
#else
            double calcRate(vector<SMCGlobal::rate_tuple_t> & rates);
#endif
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

#if defined(NEWWAY)
    inline double SpeciesForest::drawIncrement() {
        unsigned n = (unsigned)_lineages.size();
        double incr = SMCGlobal::_infinity;
        if (n > 1) {
            double r = SMCGlobal::_lambda*n;
            incr = -log(1.0 - rng.uniform())/r;
        }
        return incr;
    }
#else
    inline double SpeciesForest::calcRate(vector<SMCGlobal::rate_tuple_t> & rates) {
        unsigned n = (unsigned)_lineages.size();
        double r = 0.0;
        if (n > 1) {
            r = SMCGlobal::_lambda*n;
        
            // rates entry stores: (1) rate; (2) gene index (0 because species tree); and
            // (3) species (ignored for species tree, so just entering 0)
            rates.push_back(make_tuple(r, 0, (SMCGlobal::species_t)0));
        }
        return r;
    }
#endif

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
    
    inline void SpeciesForest::speciationEvent(SMCGlobal::species_t & left, SMCGlobal::species_t & right, SMCGlobal::species_t & anc) {
        unsigned nlineages = (unsigned)_lineages.size();
        
        // Choose two lineages to join
        assert(nlineages > 1);
        auto chosen_pair = rng.nchoose2(nlineages);
        Node * first_node  = _lineages[chosen_pair.first];
        Node * second_node = _lineages[chosen_pair.second];
        
        // Create ancestral node
        Node * anc_node = joinLineagePair(first_node, second_node);
        anc_node->setSpeciesToUnion(first_node->getSpecies(), second_node->getSpecies());
        
        // Update lineage vector
        removeTwoAddOne(_lineages, first_node, second_node, anc_node);
        
        // Return species joined in supplied reference variables
        left  = first_node->getSpecies();
        right = second_node->getSpecies();
        anc   = anc_node->getSpecies();
    }
    
//    inline void SpeciesForest::simulateSpeciesTree() {
//        createTrivialForest();
//        unsigned nsteps = SMCGlobal::_nspecies - 1;
//        for (unsigned i = 0; i < nsteps; ++i) {
//            unsigned nlineages = (unsigned)_lineages.size();
//            assert(nlineages > 1);
//
//            // Choose waiting time until the next speciation event
//            double u = rng.uniform();
//            double r = SMCGlobal::_lambda*nlineages;
//            double t = -log(1.0 - u)/r;
//
//            // Increment height of forest
//            _forest_height += t;
//
//            // Update edge lengths
//            for (auto nd : _lineages)
//                nd->_edge_length += t;
//
//            // Choose two lineages to join
//            auto chosen_pair = rng.nchoose2(nlineages);
//
//            Node * first_node  = _lineages[chosen_pair.first];
//            Node * second_node = _lineages[chosen_pair.second];
//            Node * anc_node    = joinLineagePair(first_node, second_node);
//            anc_node->setSpeciesToUnion(first_node->getSpecies(), second_node->getSpecies());
//
//            // Update lineage vector
//            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
//        }
//
//        heightsInternalsPreorders();
//    }
    
    inline double SpeciesForest::advanceSpeciesForest(unsigned particle, unsigned step) {
#if 0
        // Advances species forest by one join given the current gene trees.
        // Also adds to _log_species_tree_prior.
        //
        // Step 0
        //    <------ gene tree 0 ---->      <------ gene tree 1 ---->
        //
        //    0   0   1   1   1   2   2      0   0   1   1   1   2   2  0   _
        //    |   |   |   |   |   |   |      |   |   |   |   |   |   |  1   |
        //    |   |   |   |   |   |   |      +-0-+   |   |   |   |   |  2   | first speciation event must
        //    +-0-+   |   |   |   |   |        |     +-1-+   |   |   |  3   | occur in this interval,
        //      |     |   |   |   |   |        |       |     |   |   |  4   | regardless of the species
        //      |     |   +-1-+   |   |        |       |     |   |   |  5   | joined
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
        //    0   0   1   1   1   2   2      0   0   1   1   1   2   2  0   Begin by choosing two species
        //    |   |   |   |   |   |   |      |   |   |   |   |   |   |  1   to join: choose to join 0 and 1
        //    |   |   |   |   |   |   |      +-0-+   |   |   |   |   |  2   to create species 01. All 0 and 1
        //    +-0-+   |   |   |   |   |        |     +-1-+   |   |   |  3   species are converted to 01 below
        //      |     |   |   |   |   |        |       |     |   |   |  4   speciation event.
        //    --------------------------------------------------------- 4.5 -
        //      |     |   +-01+   |   |        |       |     |   |   |  5   | Second speciation event
        //      |     |     |     +-2-+        |       |     01,2+   |  6   - must occur in this interval
        //      |     +--01-+       |          |       |       |     |  7
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
            // Only join species after first step because a join in a species tree
            // only affects downstream (further into the past) coalescent events.
            // Joining first and then choosing increment allows filtering to take
            // account of the consequences of the join.
            pair<unsigned,unsigned> chosen_pair = rng.nchoose2(nlineages);
            Node * first_node  = _lineages[chosen_pair.first];
            Node * second_node = _lineages[chosen_pair.second];
            Node * anc_node    = joinLineagePair(first_node, second_node);
            anc_node->setSpeciesToUnion(first_node->getSpecies(), second_node->getSpecies());
            
            // Insert a speciation epoch and revise subsequent coalescent epochs to replace the two
            // old species with the new species
            //SMCGlobal::species_t oldspp1 = first_node->getSpecies();
            //SMCGlobal::species_t oldspp2 = second_node->getSpecies();
            //SMCGlobal::species_t newspp  = anc_node->getSpecies();

            // Update lineage vector
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            --nlineages;
        }

        // Create a set of all species that currently (looking backward in time) exist
        set<SMCGlobal::species_t> current_species;
        for (auto lit = _lineages.begin(); lit != _lineages.end(); lit++) {
            Node * nd = *lit;
            current_species.insert(nd->_species);
        }
        assert(current_species.size() == nlineages);
        
        // Find next coalescent event (over all gene trees) involving more than one species
        // That sets an upper bound on the time of the next speciation event, regardless of
        // which two species ultimately end up being joined.
        //double maxT = calcMaxT(_epochs, _forest_height, current_species);
        //double max_waiting_time = maxT - _forest_height;
            
        // Choose waiting time until the next speciation event (conditional on maxT)
        double u = rng.uniform();
        double r = SMCGlobal::_lambda*nlineages;
        //double t = -log(1.0 - u*(1.0 - exp(-r*max_waiting_time)))/r;
        double t = -log(1.0 - u)/r;
        
        // Update log of the species tree prior (prior is Exponential(r) but
        // conditioned on t < max_waiting_time)
        _log_species_tree_prior += log(r) - r*t;
        //_log_species_tree_prior -= log(1.0 - exp(-r*max_waiting_time));
            
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
            //SMCGlobal::species_t oldspp1 = first_node->getSpecies();
            //SMCGlobal::species_t oldspp2 = second_node->getSpecies();
            //SMCGlobal::species_t newspp  = anc_node->getSpecies();

            // Update lineage vector
            removeTwoAddOne(_lineages, first_node, second_node, anc_node);
            --nlineages;
        }
                
        // Compute coalescent likelihood
        //double log_coalescent_likelihood = 0.0;
        //for (unsigned g = 0; g < Forest::_ngenes; g++) {
        //    log_coalescent_likelihood += calcLogCoalescentLikelihood(_epochs, g);
        //}
        
        //log_weight = log_coalescent_likelihood - _prev_log_coalescent_likelihood;
        //_prev_log_coalescent_likelihood = log_coalescent_likelihood;

        //output(format("Species forest in particle %d after step %d (height = %.5f)\n    %s") % particle % step % _forest_height % makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false), 2);
        
        return log_weight;
#endif
        return 0.0;
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
