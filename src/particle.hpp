#pragma once

namespace proj {

    class Particle {
            
        public:
        
            Particle();
            Particle(const Particle & other);
            ~Particle();
            
            void clear();
            void setData(Data::SharedPtr data);
            
            void resetSpeciesForest();
            void resetGeneForests(bool compute_partials);
            
            void incrementSpeciations() {_nspeciations++;}
            unsigned getSpeciations() const {return _nspeciations;}
            void clearSpeciations() {_nspeciations = 0;}
            
            void advance(unsigned step, unsigned pindex, bool compute_partial);
            
            double getLogWeight() const {return _log_weight;}
            
            Lot::SharedPtr getLot() const {return _lot;}
            void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            
            enum last_event_t {
                LAST_EVENT_UNDEFINED = 0,
                LAST_EVENT_SPECIATION = 1,
                LAST_EVENT_COALESCENCE = 2
            };
            
            void setLastEvent(last_event_t last) {_last_event = last;}
            bool lastEventUndefined() const {return _last_event == LAST_EVENT_UNDEFINED;}
            bool lastEventSpeciation() const {return _last_event == LAST_EVENT_SPECIATION;}
            bool lastEventCoalescence() const {return _last_event == LAST_EVENT_COALESCENCE;}
            
            void refreshHeightsInternalsPreorders();
            
            SpeciesForest       & getSpeciesForest();
            const SpeciesForest & getSpeciesForest() const;

            vector<GeneForest>       & getGeneForests();
            const vector<GeneForest> & getGeneForests() const;
            
            GeneForest       & getGeneForest(unsigned gene);
            const GeneForest & getGeneForest(unsigned gene) const;
                                                
            void copyParticleFrom(const Particle & other);
            void operator=(const Particle & other);
                                
        protected:
                
            Data::SharedPtr _data;
            vector<GeneForest> _gene_forests;
            SpeciesForest _species_forest;
            unsigned _nspeciations;
            last_event_t _last_event;
            
            mutable double         _log_weight;
            mutable Lot::SharedPtr _lot;
    };
        
    inline Particle::Particle() {
        _lot.reset(new Lot());
        clear();
    }

    inline Particle::Particle(const Particle & other) {
        copyParticleFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _log_weight = 0.0;
        _gene_forests.clear();
        _species_forest.clear();
        _nspeciations = 0;
        _last_event = LAST_EVENT_UNDEFINED;
    }
    
    inline void Particle::setData(Data::SharedPtr data) {
        _data = data;
    }
            
    inline void Particle::resetSpeciesForest() {
        _species_forest.createTrivialForest();
    }

    inline void Particle::resetGeneForests(bool compute_partials) {
        assert(SMCGlobal::_ntaxa > 0);
        assert(SMCGlobal::_ngenes > 0);
        _gene_forests.clear();
        _gene_forests.resize(SMCGlobal::_ngenes);
        unsigned g = 0;
        for (auto & gf : _gene_forests) {
            gf.setData(_data);
            gf.setRelRate(SMCGlobal::_relrate_for_gene[g]);
            gf.setGeneIndex(g++);
            gf.createTrivialForest(compute_partials);
        }
    }

    inline void Particle::advance(unsigned step, unsigned pindex, bool compute_partial) {
        // Create species_tuples vector. Each 3-tuple entry stores:
        //  1. number of lineages
        //  2. gene index
        //  3. species within gene (ignored if species tree)
        vector<SMCGlobal::species_tuple_t> species_tuples;
        
        // Draw a speciation increment Delta. Note: speciation_increment will
        // equal "infinity" if species tree is complete.
        auto incr_rate = _species_forest.drawIncrement(_lot);
        double Delta = incr_rate.first;
        double speciation_rate = incr_rate.second;
        
        // Visit each species within each locus, storing the number of lineages in
        // species_tuples and computing the total coalescence rate (total_rate).
        // Next increment will be drawn from an Exponential distribution with rate
        // total_rate.
        double total_rate = 0.0;
        for (auto & gene_forest : _gene_forests) {
            total_rate += gene_forest.calcTotalRate(species_tuples, Delta);
        }

        // Draw a coalescence increment delta
        double delta = -log(1.0 - _lot->uniform())/total_rate;
        
        // If delta > Delta, then a speciation event happened; otherwise a coalescent event happened.
        bool is_speciation = (delta > Delta);
                
        // Advance all forests
        double dt = is_speciation ? Delta : delta;
        
        unsigned num_species_lineages = 1;
        if (speciation_rate > 0.0) {
            num_species_lineages = _species_forest.advanceAllLineagesBy(dt);
            if (fabs(SMCGlobal::_lambda*num_species_lineages - speciation_rate) > 0.0001) {
                throw XProj(format("%d | %.5f | %.5f | %.5f | %.5f > 0.0001") % num_species_lineages % SMCGlobal::_lambda % (SMCGlobal::_lambda*num_species_lineages) % speciation_rate % fabs(SMCGlobal::_lambda*num_species_lineages - speciation_rate));
            }
        }
        for (auto & gene_forest : _gene_forests) {
            gene_forest.advanceAllLineagesBy(dt);
        }

        if (is_speciation) {
            _last_event = Particle::LAST_EVENT_SPECIATION;
            
            // Create speciation event
            SMCGlobal::species_t left;
            SMCGlobal::species_t right;
            SMCGlobal::species_t anc;
            _species_forest.speciationEvent(_lot, left, right, anc);
                        
            // Advise all gene trees of the change in the species tree
            for (auto & gene_forest : _gene_forests) {
                gene_forest.mergeSpecies(left, right, anc);
            }
        }
        else {
            _last_event = Particle::LAST_EVENT_COALESCENCE;
            
            assert(SMCGlobal::_prior_prior);    // prior-post not yet ready
            if (SMCGlobal::_prior_prior) {
                // Choose which locus and species in which to have a coalescence

                // Create vector of lineage counts
                vector<double> lineage_counts(species_tuples.size());
                transform(species_tuples.begin(), species_tuples.end(), lineage_counts.begin(), [](SMCGlobal::species_tuple_t & spp_tuple){
                    return get<0>(spp_tuple);
                });
                
                // Determine sum of lineage_counts
                double total_num_lineages = accumulate(lineage_counts.begin(), lineage_counts.end(), 0.0);
                
                // Normalize lineage_counts to create a discrete probability distribution
                vector<double> probs(lineage_counts.size());
                transform(lineage_counts.begin(), lineage_counts.end(), probs.begin(), [total_num_lineages](unsigned nlineages){return (double)nlineages/total_num_lineages;});
                assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
                
                // Choose gene and species
                unsigned which = SMCGlobal::multinomialDraw(_lot, probs);
                //unsigned n                = get<0>(species_tuples[which]);
                unsigned g                = get<1>(species_tuples[which]);
                SMCGlobal::species_t spp  = get<2>(species_tuples[which]);

                GeneForest & gene_forest = _gene_forests[g];
            
                _log_weight = gene_forest.coalescentEvent(_lot, spp, compute_partial);
                
                // Adjust weight for fact that proposal differs from prior
                if (num_species_lineages > 1) {
                    _log_weight -= speciation_rate*(delta - Delta);
                    _log_weight -= log(speciation_rate);
                }
                assert(!isnan(_log_weight));
                assert(!isinf(_log_weight));
            }
            else {
                // prior post code goes here
            }

            // Only speciation events should have log_weight 0.0, unless we're
            // simulating, in which case calculate_partial will be false
            assert(_log_weight != 0.0 || !compute_partial);
        }
    
        if (is_speciation) {
            incrementSpeciations();
        }
    }

    inline void Particle::refreshHeightsInternalsPreorders() {
        _species_forest.heightsInternalsPreorders();
        for (auto & gf : _gene_forests) {
            gf.heightsInternalsPreorders();
        }
    }
    
    inline SpeciesForest & Particle::getSpeciesForest() {
        return _species_forest;
    }
    
    inline const SpeciesForest & Particle::getSpeciesForest() const {
        return _species_forest;
    }
    
    inline GeneForest & Particle::getGeneForest(unsigned gene) {
        assert(gene < SMCGlobal::_ngenes);
        return _gene_forests[gene];
    }

    inline vector<GeneForest> & Particle::getGeneForests() {
        return _gene_forests;
    }
    
    inline const vector<GeneForest> & Particle::getGeneForests() const {
        return _gene_forests;
    }
    
    inline const GeneForest & Particle::getGeneForest(unsigned gene) const {
        assert(gene < SMCGlobal::_ngenes);
        return _gene_forests[gene];
    }

    inline void Particle::copyParticleFrom(const Particle & other) {
        clear();

        // Performs a deep copy of other to this particle
        _data = other._data;
        
        // Copy gene forests
        assert(SMCGlobal::_ngenes == other._gene_forests.size());
        _gene_forests.resize(SMCGlobal::_ngenes);
        for (unsigned i = 0; i < SMCGlobal::_ngenes; i++)
            _gene_forests[i] = other._gene_forests[i];

        // Copy species forest
        _species_forest = other._species_forest;
        
        _nspeciations = other._nspeciations;
        _last_event = other._last_event;
        
        // No need to copy _log_weight
        // No need to copy _lot
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
}
