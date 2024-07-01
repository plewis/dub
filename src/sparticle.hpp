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
            SParticle() {createTrivialForest();}
            ~SParticle() {}
            
            void recordSpeciationInfo(vector<G::speciation_info_t> & spec_info) const;
            G::speciation_info_t heightOfNextNodeAbove(double h, const vector<G::speciation_info_t> & spec_info) const;
            
            void joinThenIncrement();
            
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
        _next_node = G::_nspecies;
        _nleaves = G::_nspecies;
        refreshSplits();
    }
    
    inline void SParticle::joinRandomPair() {
        // Grab next node to use as ancestor
        assert(_nodes.size() > _next_node);
        Node * anc = &_nodes[_next_node];
        anc->setNumber(_next_node++);
        
        // Choose two existing lineages at random to join
        pair<unsigned, unsigned> lr = rng->nchoose2(_lineages.size());
        Node * lchild = _lineages[lr.first];
        Node * rchild = _lineages[lr.second];
        
        // Connect chosen lineages to ancestor
        makeAnc(anc, lchild, rchild);
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

    inline void SParticle::joinThenIncrement() {
        joinRandomPair();
        auto incr_rate = drawIncrement();
        extendAllLineagesBy(incr_rate.first);
    }

    inline void SParticle::recordSpeciationInfo(vector<G::speciation_info_t> & spec_info) const {
        // Visit all internal nodes, for each recording the tuple:
        // height, node's species, left child's species, right child's species.
        spec_info.clear();

        // Build preorders vector
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);

        for (auto & preorder : preorders) {
            for (auto nd : preorder) {
                if (nd->getLeftChild()) {
                    // internal node
                    spec_info.push_back({
                        nd->getHeight(),
                        nd->getLeftChild()->getSpecies(),
                        nd->getLeftChild()->getRightSib()->getSpecies()
                    });
                }
            }
        }
        
        if (_lineages.size() == 1) {
            // If species tree is complete, add final element at infinity
            // with all zeros for species. This indicates that anything goes
            // with respect to coalescence because we are in the ancestral species
            spec_info.push_back({G::_infinity, G::_species_zero, G::_species_zero});
        }
        else {
            // If species tree is not complete, add final element at _height
            // with all zeros for species representing the point at which
            // the next join will occur
            spec_info.push_back({_height, G::_species_zero, G::_species_zero});
        }
        
        // Sort by increasing height above leaf level
        sort(spec_info.begin(), spec_info.end());
    }
    
    inline G::speciation_info_t SParticle::heightOfNextNodeAbove(double h, const vector<G::speciation_info_t> & spec_info) const {
        G::speciation_info_t next_speciation = make_tuple(G::_infinity, G::_species_zero, G::_species_zero);
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

