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
    };
        
    inline Particle::Particle() {
    }

    inline Particle::Particle(const Particle & other) {
        copyParticleFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
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
            gf.setGeneIndex(g++);
            gf.createTrivialForest(compute_partials);
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
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
}
