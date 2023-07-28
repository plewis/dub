#pragma once

namespace proj {

    class SpeciesForest;

    class Particle {

        friend class SpeciesForest;
        
        public:
        
            Particle();
            Particle(const Particle & other);
            ~Particle();

            void clear();
            void setData(Data::SharedPtr data);
            void initSpeciesForest(string newick);
            void initGeneForests(vector<string> & newicks);
            void resetSpeciesForest();
            void resetGeneForests();
            void digestSpeciesTree();
            void digestGeneTrees();
            double advanceGeneForest(unsigned particle, unsigned gene, unsigned step);
            double advanceSpeciesForest(unsigned particle, unsigned step);
            SpeciesForest & getSpeciesForest() {return _species_forest;}
            GeneForest & getGeneForest(unsigned gene);
            const GeneForest & getGeneForest(unsigned gene) const;
            const epoch_list_t & getEpochs() const {return _epochs;}
            epoch_list_t       & getEpochs()       {return _epochs;}
            
            double debugSaveTreesAsJavascript(string fnprefix) const;

            static void outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric);
            
            void copyFrom(const Particle & other);
            void operator=(const Particle & other) {copyFrom(other);}
                                
        protected:
                
            Data::SharedPtr _data;
            vector<GeneForest> _gene_forests;
            SpeciesForest _species_forest;
            
            epoch_list_t _epochs;
    };
    
    inline Particle::Particle() {
        //cerr << "Particle constructor" << endl;
    }

    inline Particle::Particle(const Particle & other) {
        //cerr << "Particle copy constructor" << endl;
        copyFrom(other);
    }

    inline Particle::~Particle() {
        //cerr << "Particle destructor" << endl;
        clear();
    }

    inline void Particle::clear() {
    }
    
    inline void Particle::setData(Data::SharedPtr data) {
        _data = data;
    }
            
    inline void Particle::initSpeciesForest(string newick) {
        assert(newick.size() > 0);
        _species_forest.buildFromNewick(newick);
    }
    
    inline void Particle::initGeneForests(vector<string> & newicks) {
        assert(Forest::_ntaxa > 0);
        assert(Forest::_ngenes > 0);
        assert(Forest::_ngenes == newicks.size());
        _gene_forests.clear();
        _gene_forests.resize(Forest::_ngenes);
        for (unsigned i = 0; i < Forest::_ngenes; i++) {
            _gene_forests[i].setData(_data);
            _gene_forests[i].setGeneIndex(i);
            _gene_forests[i].buildFromNewick(newicks[i]);
        }
    }
    
    inline void Particle::resetSpeciesForest() {
        _species_forest.createTrivialForest();
    }

    inline void Particle::resetGeneForests() {
        assert(Forest::_ntaxa > 0);
        assert(Forest::_ngenes > 0);
        _gene_forests.clear();
        _gene_forests.resize(Forest::_ngenes);
        for (unsigned i = 0; i < Forest::_ngenes; i++) {
            _gene_forests[i].setData(_data);
            _gene_forests[i].setGeneIndex(i);
            _gene_forests[i].createTrivialForest();
        }
    }

    inline void Particle::digestSpeciesTree() {
        _species_forest.digest();
        _epochs = _species_forest._epochs;
    }
    
    inline void Particle::digestGeneTrees() {
        // Populate _epochs based on all gene trees
        _epochs.clear();
        
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            auto & e = _gene_forests[g].digest();
            move(e.begin(), e.end(), inserter(_epochs, _epochs.end()));
        }

        // Sort all gene tree epochs by height (note: epochs
        // is actually _species_forest._epochs)
        _epochs.sort(epochLess);
    }

    inline double Particle::advanceGeneForest(unsigned particle, unsigned gene, unsigned step) {
        assert(gene < _gene_forests.size());
        GeneForest & gf = _gene_forests[gene];
        if (step == 0) {
            // Each gene forest gets a copy of the species tree epochs
            gf.createTrivialForest();
            gf.copyEpochsFrom(_epochs);
            gf.createInitEpoch();
            resetAllEpochs(gf._epochs);
        }
        bool coalescent_event = false;
        double log_weight = 0.0;
        while (!coalescent_event) {
            auto result = gf.advanceGeneForest(step, particle, gene, false);
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
    
    inline GeneForest & Particle::getGeneForest(unsigned gene) {
        assert(gene < Forest::_ngenes);
        return _gene_forests[gene];
    }

    inline const GeneForest & Particle::getGeneForest(unsigned gene) const {
        assert(gene < Forest::_ngenes);
        return _gene_forests[gene];
    }

    inline void Particle::copyFrom(const Particle & other) {
        //cerr << "in Particle::copyFrom" << endl;

        // Performs a deep copy of other to this particle
        _data = other._data;
        
        // Copy gene forests
        _gene_forests.resize(other._gene_forests.size());
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
            _gene_forests[g] = other._gene_forests[g];
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
    
    inline double Particle::debugSaveTreesAsJavascript(string fnprefix) const {
        string newick_species_tree_numeric = _species_forest.makeNewick(/*precision*/9, /*use names*/false, /*coal edge lengths*/false);
        vector<string> newick_gene_trees_numeric;
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
            newick_gene_trees_numeric.push_back(_gene_forests[g].makeNewick(/*precision*/9, /*use names*/false, /*coal edge lengths*/false));
        }
        string fn = str(format("%s.js") % fnprefix);
        outputJavascriptTreefile(fn, newick_species_tree_numeric, newick_gene_trees_numeric);

        double log_coalescent_likelihood = 0.0;
        for (unsigned g = 0; g < Forest::_ngenes; g++)
            _species_forest.calcLogCoalescentLikelihood(_species_forest._epochs, g);

        return log_coalescent_likelihood;
    }
    
}
