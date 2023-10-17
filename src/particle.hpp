#pragma once

#if defined(USING_MPI)
extern int my_rank;
extern int ntasks;
extern unsigned my_first_gene;
extern unsigned my_last_gene;
#endif

namespace proj {

    class SpeciesForest;

    class Particle {
    
        // A Particle can be used in estimating a particular gene tree, in which case
        //   1) its _gene_index is >= 0
        //   2) only element _gene_index of _gene_forests is used
        //   3) _species_forest is identical for all particles
        // A Particle can also be used in estimating the species tree, in which case
        //   1) its _gene_index is -1
        //   2) all elements _gene_forests are used and contain fully-resolved gene trees
        //   3) _gene_forests is identical for all particles

        friend class SpeciesForest;
        
        public:
        
            Particle();
            Particle(const Particle & other);
            ~Particle();
            
            void setGeneIndex(int g);
            int  getGeneIndex() const {return _gene_index;}

            void clear();
            void setData(Data::SharedPtr data);
            void initSpeciesForest(string newick);
            void initGeneForest(string newick);
            void initGeneForests(vector<string> & newicks);
            void resetSpeciesForest();
            void resetGeneForest();
            void digestSpeciesTree(bool append);
            void digestGeneTrees(bool append);
            void sortEpochs();
            void reconcileEpochs();
            void copyEpochsToSpeciesForest();
            double advanceGeneForest(unsigned particle, unsigned step);
            double advanceSpeciesForest(unsigned particle, unsigned step);
            SpeciesForest & getSpeciesForest() {return _species_forest;}
            GeneForest & getGeneForest();
            const GeneForest & getGeneForest() const;
            GeneForest & getGeneForest(unsigned gene);
            const GeneForest & getGeneForest(unsigned gene) const;
            const epoch_list_t & getEpochs() const {return _epochs;}
            epoch_list_t       & getEpochs()       {return _epochs;}
            
            double calcLogCoalLikeGivenTheta(double theta);
            double calcLogSpeciesTreeDensityGivenLambda(double lambda);
            double speciesForestHeight() const;
            
            double calcLogLikelihoodForGene(unsigned gene);
            
            double debugSaveTreesAsJavascript(string fnprefix) const;

            static void outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric);
            
            void copyFrom(const Particle & other);
            void operator=(const Particle & other) {copyFrom(other);}
                                
        protected:
                
            Data::SharedPtr _data;
            int _gene_index;    // -2 means just-constructed, -1 means "species" particle
            vector<GeneForest> _gene_forests;
            SpeciesForest _species_forest;
            epoch_list_t _epochs;
    };
    
    typedef vector<Particle> vect_particle_t;
    
    inline Particle::Particle() {
        _gene_index = -2;   // -2 flags a particle as just-constructed
    }

    inline Particle::Particle(const Particle & other) {
        copyFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _gene_index = -2;   // -2 flags a particle as just-constructed
        _gene_forests.clear();
        _species_forest.clear();
        _epochs.clear();
    }
    
    inline void Particle::setGeneIndex(int g) {
        // Should be first thing done to a particle (so _gene_index should equal -2)
        assert(_gene_index < -1);
        assert(g > -2);
        assert(g < (int)Forest::_ngenes);
        _gene_index = g;
    }
    
    inline void Particle::setData(Data::SharedPtr data) {
        _data = data;
    }
            
    inline void Particle::initSpeciesForest(string newick) {
        assert(newick.size() > 0);
        _species_forest.buildFromNewick(newick);
    }
    
    inline void Particle::initGeneForest(string newick) {
        assert(_gene_index > -2 && _gene_index < Forest::_ngenes);
        assert(Forest::_ntaxa > 0);
        assert(Forest::_ngenes > 0);
        _gene_forests.clear();
        _gene_forests.resize(Forest::_ngenes);
        _gene_forests[_gene_index].setData(_data);
        _gene_forests[_gene_index].setGeneIndex(_gene_index);
        _gene_forests[_gene_index].buildFromNewick(newick);
    }
    
    inline void Particle::initGeneForests(vector<string> & newicks) {
        assert(Forest::_ntaxa > 0);
        assert(Forest::_ngenes > 0);
        assert(Forest::_ngenes == newicks.size());
        _gene_forests.clear();
        _gene_forests.resize(Forest::_ngenes);
#if defined(USING_MPI)
        for (unsigned g = ::my_first_gene; g < ::my_last_gene; ++g) {
#else
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
#endif
            _gene_forests[g].setData(_data);
            _gene_forests[g].setGeneIndex(g);

            _gene_forests[g].buildFromNewick(newicks[g]);
        }
    }
    
    inline void Particle::resetSpeciesForest() {
        _species_forest.createTrivialForest();
    }

    inline void Particle::resetGeneForest() {
        assert(Forest::_ntaxa > 0);
        assert(Forest::_ngenes > 0);
        assert(_gene_index >= 0 && _gene_index < Forest::_ngenes);
        _gene_forests.clear();
        _gene_forests.resize(Forest::_ngenes);
        _gene_forests[_gene_index].setData(_data);
        _gene_forests[_gene_index].setGeneIndex(_gene_index);

        //temporary!
        //cerr << "~~> Particle::resetGeneForest" << endl;
        //cerr << str(format("~~>   _gene_index  = %d") % _gene_index) << endl;
        //cerr << str(format("~~>   rank         = %d") % ::my_rank) << endl;
        //cerr << str(format("~~>   first        = %d") % ::my_first_gene) << endl;
        //cerr << str(format("~~>   last         = %d") % ::my_last_gene) << endl;

        _gene_forests[_gene_index].createTrivialForest();
    }

    inline void Particle::digestSpeciesTree(bool append) {
        auto & e = _species_forest.digest();
        if (append) {
            move(e.begin(), e.end(), inserter(_epochs, _epochs.end()));
        }
        else {
            _epochs = _species_forest._epochs;
        }
    }
    
    inline void Particle::digestGeneTrees(bool append) {
        assert(_gene_index == -1);   // should only be called by a "species" particle
        
        // Populate _epochs based on all gene trees'
        if (!append)
            _epochs.clear();
        
#if defined(USING_MPI)
        for (unsigned g = ::my_first_gene; g < ::my_last_gene; ++g) {
#else
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
#endif
            auto & e = _gene_forests[g].digest();
            move(e.begin(), e.end(), inserter(_epochs, _epochs.end()));
        }
    }
    
    inline void Particle::sortEpochs() {
        // Sort all gene tree epochs by height
        _epochs.sort(epochLess);
    }
    
    inline void Particle::reconcileEpochs() {
        for (auto it = _epochs.begin(); it != _epochs.end(); ++it) {
            if (it->isSpeciationEpoch()) {
                Node::species_t species1    = it->_left_species;
                Node::species_t species2    = it->_right_species;
                Node::species_t new_species = it->_anc_species;
                
                // Replace species1 and species2 with new_species for all coalescent
                // epochs encountered before the next speciation epoch or _epochs.end()
                auto start_it = it;
                ++start_it;
                
                _species_forest.reconcileEpochsStartingAt(_epochs, start_it, species1, species2, new_species);
            }
        }
    }

    inline void Particle::copyEpochsToSpeciesForest() {
        _species_forest.copyEpochsFrom(_epochs);
    }
   
    inline double Particle::advanceGeneForest(unsigned particle, unsigned step) {
        assert(_gene_index < _gene_forests.size());
        GeneForest & gf = _gene_forests[_gene_index];
        if (step == 0) {
            // Each gene forest gets a copy of the species tree epochs
            gf.copyEpochsFrom(_epochs);
            gf.createInitEpoch(/*ancestral_species_only*/_epochs.size() == 0 ? true : false);
            resetAllEpochs(gf._epochs);
        }
        bool coalescent_event = false;
        double log_weight = 0.0;
        while (!coalescent_event) {
            auto result = gf.advanceGeneForest(step, particle, false);
            coalescent_event = result.first;
            log_weight = result.second;
        }
        return log_weight;
    }
    
    inline double Particle::advanceSpeciesForest(unsigned particle, unsigned step) {
        double log_weight = 0.0;
        if (step == 0) {
            // The species forest gets a copy of the combined epochs from all gene trees
            _species_forest.copyEpochsFrom(_epochs);
        }
        log_weight = _species_forest.advanceSpeciesForest(particle, step);
        return log_weight;
    }
    
    inline GeneForest & Particle::getGeneForest() {
        // Should be called only for a "gene" particle
        assert(_gene_index > -1 && _gene_index < Forest::_ngenes);
        
        return _gene_forests[_gene_index];
    }

    inline const GeneForest & Particle::getGeneForest() const {
        // Should be called only for a "gene" particle
        assert(_gene_index > -1 && _gene_index < Forest::_ngenes);
        
        return _gene_forests[_gene_index];
    }

    inline GeneForest & Particle::getGeneForest(unsigned gene) {
        // Should be called for a "species" particle or a "gene" particle only if gene is correct
        assert(_gene_index == -1 || _gene_index == gene);
        
        assert(gene < Forest::_ngenes);
        return _gene_forests[gene];
    }

    inline const GeneForest & Particle::getGeneForest(unsigned gene) const {
        // Should be called for a "species" particle or a "gene" particle only if gene is correct
        assert(_gene_index == -1 || _gene_index == gene);
        
        assert(gene < Forest::_ngenes);
        return _gene_forests[gene];
    }

    inline void Particle::copyFrom(const Particle & other) {
        clear();

        // Performs a deep copy of other to this particle
        _data = other._data;
        _gene_index = other._gene_index;
        
        // Copy gene forests
        _gene_forests.resize(other._gene_forests.size());
        if (_gene_index >= 0) {
            // "gene" particle
            _gene_forests[_gene_index] = other._gene_forests[_gene_index];
        }
        else {
            // "species" particle
            for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                _gene_forests[g] = other._gene_forests[g];
            }
        }

        // Copy species forest
        _species_forest = other._species_forest;

        // Copy _epochs;
        _epochs = other._epochs;
    }
    
    inline void Particle::outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric) {
        ofstream jsf(fn);
        jsf << "let species_translate = {\n";
        unsigned i = 1;
        for (string nm : Forest::_species_names) {
            string comma = (i < Forest::_nspecies ? "," : "");
            jsf << str(format("  %d: \"%s\"%s\n") % i % nm % comma);
            ++i;
        }
        jsf << "};\n\n";
        
        jsf << str(format("let species_newick = \"%s\";\n\n") % newick_species_tree_numeric);
        
        jsf << "let gene_translate = {\n";
        i = 1;
        for (string nm : Forest::_taxon_names) {
            string comma = (i < Forest::_ntaxa ? "," : "");
            jsf << str(format("  %d: \"%s\"%s\n") % i % nm % comma);
            ++i;
        }
        jsf << "};\n\n";

        jsf << "let gene_newicks = [\n";
        for (unsigned g = 0; g < newick_gene_trees_numeric.size(); ++g) {
            jsf << str(format("  {name:\"gene-%d\",  relrate:1.0, newick:\"%s\"},\n") % (g+1) % newick_gene_trees_numeric[g]);
        }
        jsf << "];\n";
        
        jsf.close();
    
    }

    //TODO: needs modification for MPI
    inline double Particle::debugSaveTreesAsJavascript(string fnprefix) const {
        // Should only be called on "species" particles
        assert(_gene_index == -1);
        
        string newick_species_tree_numeric = _species_forest.makeNewick(/*precision*/9, /*use names*/false, /*coal edge lengths*/false);
        vector<string> newick_gene_trees_numeric;
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
            newick_gene_trees_numeric.push_back(_gene_forests[g].makeNewick(/*precision*/9, /*use names*/false, /*coal edge lengths*/false));
        }
        string fn = str(format("%s.js") % fnprefix);
        outputJavascriptTreefile(fn, newick_species_tree_numeric, newick_gene_trees_numeric);

        double log_coalescent_likelihood = 0.0;
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            log_coalescent_likelihood += _species_forest.calcLogCoalescentLikelihood(_species_forest._epochs, g);
        }

        return log_coalescent_likelihood;
    }
    
   inline double Particle::calcLogCoalLikeGivenTheta(double theta) {
        resetAllEpochs(_species_forest._epochs);

        _species_forest.debugShowEpochs(_species_forest._epochs);
        
        double prev_theta = Forest::_theta;
        Forest::_theta = theta;
        
        double log_coal_like = 0.0;
        
        //TODO: needs modification for MPI
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
            log_coal_like += _species_forest.calcLogCoalescentLikelihood(_species_forest._epochs, g);
        }
        
        Forest::_theta = prev_theta;
        return log_coal_like;
   }
       
   inline double Particle::speciesForestHeight() const {
        return _species_forest.getHeight();
   }
   
   inline double Particle::calcLogSpeciesTreeDensityGivenLambda(double lambda) {
        double prev_lambda = Forest::_lambda;
        Forest::_lambda = lambda;
        
        double log_density = _species_forest.calcLogSpeciesTreeDensity(lambda);
        
        Forest::_lambda = prev_lambda;
        return log_density;
   }
   
   inline double Particle::calcLogLikelihoodForGene(unsigned gene) {
        // Should be called on "species" particles or, if on a "gene" particle, gene should be correct
        assert(_gene_index == -1 || _gene_index == gene);
        assert(gene < Forest::_ngenes);

        // Calculate partials
        _gene_forests[gene].computeAllPartials();
        
        // Now compute the log-likelihood at the root
        return _gene_forests[gene].calcLogLikelihood();
   }
   
}
