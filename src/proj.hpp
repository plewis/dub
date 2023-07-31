#pragma once

using boost::filesystem::current_path;
using boost::algorithm::join;
using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::store;
using boost::program_options::parse_command_line;
using boost::program_options::parsed_options;
using boost::program_options::parse_config_file;
using boost::program_options::reading_file;
using boost::program_options::notify;

extern proj::PartialStore ps;
extern proj::Lot rng;

namespace proj {

    class Proj {
        public:
                                       Proj();
                                       ~Proj();

            void                       clear();
            void                       processCommandLineOptions(int argc, const char * argv[]);
            void                       run();

        private:
            
            void                       debugShowStringVector(string title, const vector<string> & svect) const;
            void                       debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const;

            static string              inventName(unsigned k, bool lower_case);
            
            void                       simulateData();
            void                       readData();
            unsigned                   buildSpeciesMap();
            void                       drawStartingSpeciesTree();
            void                       readStartingSpeciesTree();
            void                       readStartingGeneTrees();
            void                       clearParticles();
            void                       initializeParticles(bool start_with_species_tree);
            double                     filterParticles(const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like);
            void                       saveGeneTreesUsingNames(string filename, Particle & p);
            void                       saveGeneTreesUsingNumbers(string filename, Particle & p);
            void                       propagateSampledParticle(bool growing_gene_forests);
            void                       debugCheckGeneTrees() const;

            void                       growGeneTrees(unsigned iter);
            void                       growSpeciesTrees(unsigned iter);
            void                       updateTheta(Particle & p, unsigned ntries, double delta);

            void                       saveSpeciesTreeUsingNames(string fn, Particle & p);
            void                       saveSpeciesTreeUsingNumbers(string fn, Particle & p);
            void                       saveUniqueSpeciesTrees(string fn, const vector<unsigned> & counts);
            void                       showSettings() const;
            void                       outputGeneTreesToFile(string fn, const vector<string> & newicks) const;
            void                       outputNexusTreefile(string fn, const vector<string> & newicks) const;
            void                       outputAnnotatedNexusTreefile(string fn, const vector<string> & newicks, const vector<string> & treenames, const vector<string> & annotations) const;
            
            string                     _data_file_name;
            string                     _species_tree_file_name;
            string                     _gene_tree_file_name;
            string                     _start_mode;
            unsigned                   _niter;
            Partition::SharedPtr       _partition;
            Data::SharedPtr            _data;
            
            //TreeSummary                _starting_gene_trees;
            //vector<Tree::SharedPtr>    _gene_trees;
            string                      _starting_species_newick;
            vector<string>              _starting_gene_newicks;
            
            bool                       _use_gpu;
            bool                       _ambig_missing;
            bool                       _verbose;
            bool                       _simulate;
            unsigned                   _nsimspecies;
            vector<unsigned>           _nsimtaxaperspecies;
            vector<unsigned>           _nsites_per_gene;
            unsigned                   _nparticles;
            int                        _track_split;
            unsigned                   _rnseed;
            bool                       _sort_forests;
            double                     _visualization_cutoff;
                        
            vector<Particle>            _particles;
#if defined(POLTMP)
            Particle                    _template_particle;
#endif

            static string               _program_name;
            static unsigned             _major_version;
            static unsigned             _minor_version;
    };

    inline Proj::Proj() {
        clear();
    }

    inline Proj::~Proj() {
    }
    
    inline void Proj::clear() {
        _use_gpu                = true;
        _ambig_missing          = true;
        _sort_forests           = false;
        _verbose                = false;
        _visualization_cutoff   = 0.99;

        // data related
        _data                   = nullptr;
        _data_file_name         = "";
        _species_tree_file_name = "stree.tre";
        _gene_tree_file_name    = "gtrees.tre";
        _start_mode             = "random"; // "random" (random species tree), "species" (species tree from file), "gene" (gene trees from file)
        _niter                  = 1;
        _partition.reset(new Partition());

        // simulation related
        _simulate = false;
        _nsimspecies = 5;
        _nsimtaxaperspecies = {2,2,2,2,2};
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> partition_subsets;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile",  value(&_data_file_name)->required(), "name of a data file in NEXUS format")
        ("speciestreefile",  value(&_species_tree_file_name), "name of a tree file in NEXUS format containing an ultrametric starting species tree (also used to store simulated species tree if simulate is true)")
        ("genetreefile",  value(&_gene_tree_file_name), "name of a tree file for: (1) obtaining starting gene trees (startspeciestree false); (2) storing gene trees (startspeciestree true); or (3) storing true gene trees for each gene (if simulate true)")
        ("startmode", value(&_start_mode), "if 'random', start with species tree drawn from prior; if 'species', start with first species tree defined in speciestreefile; if 'gene', start with gene trees defined in genetreefile")
        ("niter", value(&_niter), "number of iterations, where one iteration involves SMC of gene trees give species tree combined with an SMC of species tree given gene trees")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("gpu",           value(&_use_gpu)->default_value(true), "use GPU if available")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("verbose",  value(&_verbose)->default_value(false), "if true, saveAllForests called after each generation")
        ("simulate",  value(&_simulate)->default_value(false), "if true, data set will be simulated (note: file specified by datafile setting will be overwritten if it exists)")
        ("nspecies",  value(&_nsimspecies)->default_value(1), "number of species (only used if simulate specified)")
        ("ntaxaperspecies",  value(&_nsimtaxaperspecies), "number of taxa sampled per species (only used if simulate specified); should be _nimspecies of these entries, one for each species simulated")
        ("nparticles",  value(&_nparticles)->default_value(10), "number of particles in the population")
        ("theta",  value(&Forest::_theta)->default_value(0.05), "coalescent parameter assumed for gene trees")
        ("lambda",  value(&Forest::_lambda)->default_value(10.9), "per lineage speciation rate assumed for the species tree")
        ("rnseed",  value(&_rnseed)->default_value(13579), "pseudorandom number seed")
        ("sortforests",  value(&_sort_forests)->default_value(false), "sort forests by weight when saving to file")
        //("priorpost",  value(&Forest::_prior_post)->default_value(true), "if yes, use prior-post to choose gene tree pairings; if no, use prior-prior (species tree always uses prior-prior)")
        ("visualizationcutoff", value(&_visualization_cutoff)->default_value(0.99), "particles sorted from highest to lowest weight will be saved for visualization if cumulative weight is less than this value")
        ;
        
        store(parse_command_line(argc, argv, desc), vm);
        try {
            const parsed_options & parsed = parse_config_file< char >("proj.conf", desc, false);
            store(parsed, vm);
        }
        catch(reading_file & x) {
            cout << "Note: configuration file (proj.conf) not found" << endl;
        }
        notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            cout << desc << "\n";
            exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            cout << str(format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << endl;
            exit(1);
        }
        
        // If user specified --subset on command line, break specified partition subset
        // definition into name and character set string and add to _partition
        if (vm.count("subset") > 0) {
            _partition.reset(new Partition());
            for (auto s : partition_subsets) {
                _partition->parseSubsetDefinition(s);
            }
        }
    }
    
    inline void Proj::readData() {
        cout << "\nReading and storing the data in the file " << _data_file_name << endl;
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);
        
        Forest::_gene_names.clear();

        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        cout << "\nNumber of taxa: " << _data->getNumTaxa() << endl;
        cout << "Number of partition subsets: " << nsubsets << endl;
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(nsubsets);
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Set length of partials for gene g
            ps.setNElements(Forest::_nstates*_data->getNumPatternsInSubset(subset), subset);
            
            DataType dt = _partition->getDataTypeForSubset(subset);
            Forest::_gene_names.push_back(_data->getSubsetName(subset));
            cout << "  Subset " << (subset+1) << " (" << _data->getSubsetName(subset) << ")" << endl;
            cout << "    data type: " << dt.getDataTypeAsString() << endl;
            cout << "    sites:     " << _data->calcSeqLenInSubset(subset) << endl;
            cout << "    patterns:  " << _data->getNumPatternsInSubset(subset) << endl;
            }
    }
    
    inline void Proj::debugShowStringVector(string title, const vector<string> & svect) const {
#if defined(DEBUGGING)
        cout << "\n" << title << ":" << endl;
        for (auto s : svect) {
            cout << "  " << s << endl;
         }
#endif
    }
    
    inline void Proj::debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const {
#if defined(DEBUGGING)
        cout << "\n" << title << ":" << endl;
        for (auto su : sumap) {
            cout << "  " << su.first << ": " << su.second << endl;
         }
#endif
    }
    
    inline unsigned Proj::buildSpeciesMap() {
        // Populates Forest::_species_names and Forest::_taxon_to_species
        unsigned nspecies = 0;
        map<string, unsigned> species_name_to_index;
        Forest::_species_names.clear();
        Forest::_taxon_to_species.clear();
        
        cout << "\nMapping taxa to species" << endl;
        const Data::taxon_names_t & tnames = _data->getTaxonNames();
        unsigned ntax = (unsigned)tnames.size();
        for (auto & tname : tnames) {
            string species_name = Node::taxonNameToSpeciesName(tname);
            cout << str(format("  %s --> %s\n") % tname % species_name);
            unsigned species_index = ntax;
            if (species_name_to_index.find(species_name) == species_name_to_index.end()) {
                // species_name not found
                species_index = nspecies;
                Forest::_species_names.push_back(species_name);
                species_name_to_index[species_name] = nspecies++;
            } else {
                // species_name found
                species_index = species_name_to_index[species_name];
            }
            Forest::_taxon_to_species[tname] = species_index;
        }
        return (unsigned)Forest::_species_names.size();
    }
    
    inline void Proj::readStartingSpeciesTree() {
        cout << "\nReading and storing the species tree in the file " << _species_tree_file_name << endl;
        vector<string> tree_names;
        vector<string> newicks;
        // _nexus_taxon_map[k] provides 0-based index (into _species_names vector) of
        // the taxon having 1-based index k in the tree description
        Forest::readTreefile(   _species_tree_file_name,
                                0,
                                Forest::_species_names,
                                Forest::_nexus_taxon_map,
                                tree_names,
                                newicks);
        _starting_species_newick = newicks[0];
    }
    
    inline void Proj::readStartingGeneTrees() {
        cout << "\nReading and storing the gene trees in the file " << _gene_tree_file_name << endl;
        vector<string> tree_names;
        // _nexus_taxon_map[k] provides 0-based index (into _taxon_names vector) of
        // the taxon having 1-based index k in the tree description
        Forest::readTreefile(   _gene_tree_file_name,       // tree file name
                                0,                          // skip
                                Forest::_taxon_names,       // leaf names
                                Forest::_nexus_taxon_map,   // taxon map
                                tree_names,                 // container for tree names
                                _starting_gene_newicks);    // container for tree descriptions
    
        unsigned ntrees = (unsigned)_starting_gene_newicks.size();
        cout << "  no. trees: " << ntrees << endl;
        for (unsigned i = 0; i < ntrees; i++) {
            string tree_name = tree_names[i];
            string gene_name = Forest::_gene_names[i];
            if (gene_name != tree_name) {
                throw XProj(format("Expecting name of %dth gene tree to be \"%s\" but it was instead \"%s\"") % i % gene_name % tree_name);
            }
            cout << "  tree " << (i+1) << " has name " << tree_name << endl;
        }
    }
    
    inline void Proj::showSettings() const {
        cout << "Speciation rate (lambda): " << Forest::_lambda << endl;
        cout << "Coalescent parameter (theta): " << Forest::_theta << endl;
        cout << "Number of particles: " << _nparticles << endl;
    }
    
    inline void Proj::initializeParticles(bool start_with_species_tree) {
        assert(_data);
        
        cout << str(format("\nInitializing %d particles...") % _nparticles) << endl;

        clearParticles();
        _particles.resize(_nparticles);

        // Set up a "template" particle the hard way
#if defined(POLTMP)
        _template_particle.clear();
#else
        Particle _template_particle;
#endif
        _template_particle.setData(_data);
        if (start_with_species_tree) {
            cout << "  Particles all have starting species tree but trivial gene forests\n";
            _template_particle.initSpeciesForest(_starting_species_newick);
            _template_particle.digestSpeciesTree();
            _template_particle.resetGeneForests();
        }
        else {
            cout << "  Particles all have starting gene trees but trivial species forests\n";
            _template_particle.resetSpeciesForest();
            _template_particle.initGeneForests(_starting_gene_newicks);
            _template_particle.digestGeneTrees();
        }
        
        // Now use first as a template and copy to all other particles
        // Saves recalculating partials many times
        fill(_particles.begin(), _particles.end(), _template_particle);
    }
    
    inline void Proj::saveSpeciesTreeUsingNames(string treefname, Particle & p) {
        SpeciesForest & sf = p.getSpeciesForest();
        ofstream tmpf(treefname);
        tmpf << "#nexus\n\n";
        tmpf << "begin trees;\n";
        tmpf << str(format("  tree species_tree = [&R] %s;\n") % sf.makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
        tmpf << "end;\n";
        tmpf.close();
    }
    
    inline void Proj::saveSpeciesTreeUsingNumbers(string treefname, Particle & p) {
        SpeciesForest & sf = p.getSpeciesForest();
        ofstream tmpf(treefname);
        tmpf << "#nexus\n\n";
        tmpf << "begin trees;\n";
        tmpf << "  translate\n";
        for (unsigned s = 0; s < Forest::_nspecies; ++s) {
            tmpf << str(format("    %d %s%s\n") % (s+1) % Forest::_species_names[s] % (s == Forest::_nspecies - 1 ? ";" : ","));
        }
        tmpf << str(format("  tree species_tree = [&R] %s;\n") % sf.makeNewick(/*precision*/9, /*use names*/false, /*coalescent units*/false));
        tmpf << "end;\n";
        tmpf.close();
    }
    
    inline void Proj::saveGeneTreesUsingNames(string treefname, Particle & p) {
        vector<string> log_likes;
        ofstream tmpf(treefname);
        tmpf << "#nexus\n\n";
        tmpf << "begin trees;\n";
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            string gene_name = Forest::_gene_names[g];
            GeneForest & gf = p.getGeneForest(g);
            double lnL = gf.calcLogLikelihood();
            log_likes.push_back(str(format("%.9f") % lnL));
            tmpf << str(format("  tree %s = [lnL = %.5f] [&R] %s;\n") % gene_name % lnL % gf.makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
        }
        tmpf << "end;\n";
        tmpf.close();
    }
        
    inline void Proj::saveGeneTreesUsingNumbers(string treefname, Particle & p) {
        vector<string> log_likes;
        ofstream tmpf(treefname);
        tmpf << "#nexus\n\n";
        tmpf << "begin trees;\n";
        tmpf << "  translate\n";
        for (unsigned t = 0; t < Forest::_ntaxa; ++t) {
            tmpf << str(format("    %d %s%s\n") % (t+1) % Forest::_taxon_names[t] % (t == Forest::_ntaxa - 1 ? ";" : ","));
        }
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            string gene_name = Forest::_gene_names[g];
            GeneForest & gf = p.getGeneForest(g);
            double lnL = gf.calcLogLikelihood();
            log_likes.push_back(str(format("%.9f") % lnL));
            tmpf << str(format("  tree %s = [lnL = %.5f] [&R] %s;\n") % gene_name % lnL % gf.makeNewick(/*precision*/9, /*use names*/false, /*coalescent units*/false));
        }
        tmpf << "end;\n";
        tmpf.close();
        
#if 0
        // This part only useful if analyzing the snake example
        ofstream paupf("paup-jclikes.nex");
        paupf << "#nexus\n";
        paupf << "\n";
        paupf << "begin paup;\n";
        paupf << "  log start file=paup-jclike-log.txt replace;\n";
        paupf << "  exe data-snakes.nex;\n";
        paupf << "  outgroup ac1OUTG1^agikstrodon ac1OUTG2^agikstrodon apc1OUTG1^agikstrodon apc1OUTG2^agikstrodon;\n";
        paupf << "  set crit=like;\n";
        paupf << "  lset clock nst=1 basefreq=equal rates=equal;\n";
        paupf << str(format("  gettrees file=%s storebrlen;\n") % treefname);
        vector<unsigned> taboo({2, 6, 13, 15});
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            unsigned gg = g + 1;
            if (find(taboo.begin(), taboo.end(), gg) != taboo.end())
                continue;
            paupf << "\n";
            paupf << str(format("include %s / only;\n") % Forest::_gene_names[g]);
            paupf << str(format("lscores %d / userbrlens;\n") % gg);
            paupf << str(format("[!*** %s %.9f ***]\n") % Forest::_gene_names[g] % log_likes[g]);
        }
        paupf << "  log stop;\n";
        paupf << "  quit;\n";
        paupf << "end;\n";
        paupf.close();
#endif
    }
        
    inline void Proj::propagateSampledParticle(bool growing_gene_forests) {
        // Choose one particle at random
        assert(_nparticles == _particles.size());
        
        // Copy that particle to all particles
        unsigned which = rng.randint(0, _nparticles - 1);
        Particle & other = _particles[which];
#if defined(POLTMP)
        _template_particle = other;
#else
        Particle _template_particle(other);
#endif
        if (growing_gene_forests) {
            // Prepare for growing species forests
            _template_particle.resetSpeciesForest();
            _template_particle.digestGeneTrees();
        }
        else {
            // Just finished updating species tree, so update theta before next round
            updateTheta(_template_particle, 100, 0.5);
            
            // Prepare for growing gene forests
            _template_particle.resetGeneForests();
            _template_particle.digestSpeciesTree();
        }
        fill(_particles.begin(), _particles.end(), _template_particle);
    }
    
    inline void Proj::debugCheckGeneTrees() const {
#if defined(DEBUGGING)
        for (auto & p : _particles) {
            for (unsigned g = 0; g < Forest::_ngenes; g++) {
                const GeneForest & gf = p.getGeneForest(g);
                gf.debugCheckSpeciesSets();
            }
        }
#endif
    }
    
    inline void Proj::growGeneTrees(unsigned iter) {
        // Grows all gene trees conditional on the species tree stored in each particle
        // Assumes gene forests are trivial and species forest has been digested into epochs
        
        cout << "Growing gene trees..." << endl;
        
        // Number of joins is one fewer than the number of taxa
        unsigned nsteps = Forest::_ntaxa - 1;
                
        // Loop over steps
        double log_marg_like = 0.0;
        vector<unsigned> counts;
        for (unsigned step = 0; step < nsteps; step++) {
            cout << "  step " << (step + 1) << " of " << nsteps << ": ";
            
            vector<double> logw(_particles.size(), 0.0);

            // Loop over particles, advancing all genes one step
            unsigned p = 0;
            for (auto & particle : _particles) {
                // Loop over genes
                for (unsigned g = 0; g < Forest::_ngenes; g++) {
                    double log_weight = particle.advanceGeneForest(p, g, step);
                    logw[p] += log_weight;
                }
                ++p;
            }
            
            // Filter particles
            double ESS = filterParticles(logw, counts, log_marg_like);
            cout << str(format("ESS = %.1f%%\n") % (100.9*ESS/_nparticles));

        } // step loop
        
        cout << str(format("  log(marg. like.) = %.5f\n") % log_marg_like);
        
        // Choose one particle at random and propagate it to all other particles
        propagateSampledParticle(/*growing_gene_forests*/true);
             
        string fn = str(format("gene-trees-chosen-iter-%d.tre") % iter);
        saveGeneTreesUsingNames(fn, _particles[0]);
        //fn = str(format("gene-trees-chosen-iter-%d-numbers.tre") % iter);
        //saveGeneTreesUsingNumbers(fn, _particles[0]);
    }
    
    inline void Proj::growSpeciesTrees(unsigned iter) {
        // Grows species trees conditional on the gene trees stored in each particle

        cout << "Growing species tree..." << endl;
        
        // Number of joins is one fewer than the number of species
        unsigned nsteps = Forest::_nspecies - 1;
        
        // Compile information about sequence and timing of joins in the gene trees
        // for each particle. This info will be used in Particle::advanceSpeciesForest.
        //for (auto & particle : _particles) {
        //    particle.resetSpeciesForest();
        //    particle.digestGeneTrees();
        //}
        
        // Loop over steps
        double log_marg_like = 0.0;
        vector<unsigned> counts;
        for (unsigned step = 0; step < nsteps; step++) {
            cout << "  step " << (step + 1) << " of " << nsteps << ": ";
            vector<double> logw(_particles.size(), 0.0);

            // Loop over particles, advancing species forest one step in each
            unsigned p = 0;
            for (auto & particle : _particles) {
                double log_weight = particle.advanceSpeciesForest(p, step);
                logw[p] += log_weight;
                ++p;
            }
                        
            // Filter particles
            double ESS = filterParticles(logw, counts, log_marg_like);
            cout << str(format("ESS = %.1f%%\n") % (100.9*ESS/_nparticles));
        }

        cout << str(format("  log(marg. like.) = %.5f\n") % log_marg_like);
        
        string fn = str(format("species-trees-after-iter-%d.tre") % iter);
        saveUniqueSpeciesTrees(fn, counts);
        
        auto maxit = max_element(counts.begin(), counts.end());
        assert(maxit != counts.end());
        unsigned maxi = (unsigned)distance(counts.begin(), maxit);
        
        cout << str(format("  Best particle had index %d with count %d\n") % maxi % (*maxit));

        fn = str(format("species-tree-chosen-iter-%d.tre") % iter);
        saveSpeciesTreeUsingNames(fn, _particles[maxi]);
        
        string fnprefix = str(format("newicks-%d") % iter);
        double log_coal_like = _particles[maxi].debugSaveTreesAsJavascript(fnprefix);
        cout << str(format("  log(coal. like.) after iter %d = %.9f\n") % iter % log_coal_like);

        // Choose one particle at random and propagate it to all other particles
        propagateSampledParticle(/*growing_gene_forests*/false);
    }
    
    inline void Proj::updateTheta(Particle & p, unsigned ntries, double delta) {
        // Use multiple-try Metropolis to update theta conditional on the gene forests
        // and species forest defined in p
        double theta0 = Forest::_theta;
        
        // Sample ntries new values of theta from symmetric proposal distribution
        // Compute weights (coalescent likelihood) for each proposed_thetas value
        vector<double> proposed_thetas(ntries, 0.0);
        vector<double> logwstar(ntries, 0.0);
        for (unsigned i = 0; i < ntries; ++i) {
            double q = theta0 - delta + 2.0*delta*rng.uniform();
            if (q < 0.0)
                q = -q;
            proposed_thetas[i] = q;
            logwstar[i] = p.calcLogCoalLikeGivenTheta(q);
        }
        
        // Compute log of the sum of the weights (this sum will form the numerator of the acceptance ratio)
        double log_sum_numer_weights = Forest::calcLogSum(logwstar);

        // Normalize weights to create a discrete probability distribution
        vector<double> probs(ntries, 0.0);
        transform(logwstar.begin(), logwstar.end(), probs.begin(), [log_sum_numer_weights](double logw){
            return exp(logw - log_sum_numer_weights);
        });
        
        // Choose one theta value from the probability distribution
        unsigned which = Forest::multinomialDraw(probs);
        double theta_star = proposed_thetas[which];
        
        // Sample ntries-1 new values of theta from symmetric proposal distribution
        logwstar[0] = p.calcLogCoalLikeGivenTheta(theta0);
        for (unsigned i = 1; i < ntries; ++i) {
            double q = theta_star - delta + 2.0*delta*rng.uniform();
            if (q < 0.0)
                q = -q;
            logwstar[i] = p.calcLogCoalLikeGivenTheta(q);
        }
        
        // Compute log of the sum of the weights (this sum will form the denominator of the acceptance ratio)
        double log_sum_denom_weights = Forest::calcLogSum(logwstar);

        // Compute acceptance ratio
        double logr = log_sum_numer_weights - log_sum_denom_weights;
        bool accept = true;
        if (logr < 0.0) {
            double logu = log(rng.uniform());
            accept = logu < logr;
        }
        if (accept) {
            Forest::_theta = theta_star;
            cout << str(format("\n*** new theta: %.5f\n") % theta_star);
        }
    }
    
    inline void Proj::saveUniqueSpeciesTrees(string fn, const vector<unsigned> & counts) {
        vector<string> tree_names;
        vector<string> notes;
        vector<string> species_tree_newicks;
        unsigned p = 0;
        unsigned i = 0;
        for (auto c : counts) {
            if (c > 0) {
                unsigned nparticles = (unsigned)_particles.size();
                double pct = 100.0*c/nparticles;
                notes.push_back(str(format("This tree found in %d particles (%.1f%% of %d total particles)") % c % pct % nparticles));
                tree_names.push_back(str(format("tree%d-%d") % i % c));
                species_tree_newicks.push_back(_particles[p].getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
                ++i;
            }
            ++p;
        }
        outputAnnotatedNexusTreefile(fn, species_tree_newicks, tree_names, notes);
    }
    
    inline double Proj::filterParticles(const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like) {
        // Sanity checks
        assert(_particles.size() == _nparticles);
        assert(log_weights.size() == _nparticles);
        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = Forest::calcLogSum(log_weights);
        vector<double> probs(_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood
        log_marg_like += log_sum_weights - log(_nparticles);
        
        // Compute effective sample size
        double sum_squared_weights = 0.0;
        for (auto it = probs.begin(); it != probs.end(); it++) {
            double w = *it;
            sum_squared_weights += w*w;
        }
        double ess = 1.0/sum_squared_weights;
        
        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());
        
        // Create vector of counts storing number of darts hitting each particle
        counts.resize(_nparticles);
        counts.assign(_nparticles, 0);
        
        // Throw _nparticles darts
        for (unsigned i = 0; i < _nparticles; ++i) {
            double u = rng.uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            counts[which]++;
        }
        
        // The vector counts should represent the results of multinomial sampling, so make a copy
        // to use for filtering particles
        vector<unsigned> working_counts(counts.begin(), counts.end());
        
        // Find index of first element of working_counts with value == 0
        auto copy_to_iter = find_if(working_counts.begin(), working_counts.end(), [](unsigned c){return c == 0;});
        
        // Find index of first element of working_counts with value > 1
        auto copy_from_iter = find_if(working_counts.begin(), working_counts.end(), [](unsigned c){return c > 1;});
        
        // Copy particles with count > 1 to particles with count = 0
        while (copy_from_iter != working_counts.end()) {
            unsigned index_from = (unsigned)distance(working_counts.begin(), copy_from_iter);
            
            // n is the count for the particle
            assert(*copy_from_iter > 1);
            while (*copy_from_iter > 1) {
                // It's a zero-sum game, so there should always be a slot to copy this particle into
                assert(copy_to_iter != working_counts.end());
                
                // Find index of particle to copy to
                unsigned index_to = (unsigned)distance(working_counts.begin(), copy_to_iter);
                
                // Copy the "from" particle to the "to" partcle
                _particles[index_to] = _particles[index_from];
                
                // Adjust working_counts to reflect copying done
                (*copy_to_iter)++;
                (*copy_from_iter)--;
                
                // Advance iterator to the next "to" particle
                copy_to_iter = find_if(++copy_to_iter, working_counts.end(), [](unsigned c){return c == 0;});
            }
            copy_from_iter = find_if(++copy_from_iter, working_counts.end(), [](unsigned c){return c > 1;});
        }
        
        return ess;
    }
    
    inline void Proj::clearParticles() {
        _particles.clear();
    }
    
    inline void Proj::run() {
        cout << "Starting..." << endl;
        cout << "Current working directory: " << current_path() << endl;
        
        try {
            rng.setSeed(_rnseed);
            
            if (_simulate) {
                simulateData();
            }
            else {
                readData();
                debugShowStringVector("Gene names", Forest::_gene_names);

                Forest::_ngenes = _data->getNumSubsets();
                assert(Forest::_ngenes > 0);

                Forest::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(Forest::_taxon_names);
                
                Forest::_nspecies = buildSpeciesMap();
                debugShowStringVector("Species names", Forest::_species_names);
                
                GeneForest::computeLeafPartials(_data);
                
                if (_start_mode == "random") {
                    drawStartingSpeciesTree();
                    initializeParticles(/*start with species tree*/true);
                    for (unsigned i = 0; i < _niter; ++i) {
                        cout << str(format("\nIteration %d of %d...") % (i+1) % _niter) << endl;
                        growGeneTrees(i+1);
                        growSpeciesTrees(i+1);
                    }
                }
                else if (_start_mode == "species") {
                    readStartingSpeciesTree();
                    initializeParticles(/*start with species tree*/true);
                    for (unsigned i = 0; i < _niter; ++i) {
                        cout << str(format("\nIteration %d of %d...") % (i+1) % _niter) << endl;
                        growGeneTrees(i+1);
                        growSpeciesTrees(i+1);
                    }
                }
                else {
                    readStartingGeneTrees();
                    initializeParticles(/*start with species tree*/false);
                    for (unsigned i = 0; i < _niter; ++i) {
                        cout << str(format("\nIteration %d of %d...") % (i+1) % _niter) << endl;
                        growSpeciesTrees(i+1);
                        growGeneTrees(i+1);
                    }
                    cout << "\nFinal species tree SMC..." << endl;
                    growSpeciesTrees(_niter+1);
                }
            }
        }
        catch (XProj & x) {
            cerr << "Proj encountered a problem:\n  " << x.what() << endl;
        }

        cout << "\nFinished!" << endl;
    }
    
    inline string Proj::inventName(unsigned k, bool lower_case) {
        // If   0 <= k < 26, returns A, B, ..., Z,
        // If  26 <= k < 702, returns AA, AB, ..., ZZ,
        // If 702 <= k < 18278, returns AAA, AAB, ..., ZZZ, and so on.
        //
        // For example, k = 19009 yields ABCD:
        // ABCD 19009 = 26 + 26*26 + 26*26*26 + 0*26*26*26 + 1*26*26 + 2*26 + 3
        //              <------- base ------>   ^first       ^second   ^third ^fourth
        // base = (26^4 - 1)/25 - 1 = 18278
        //   26^1 + 26^2 + 26^3 = 26^0 + 26^1 + 26^2 + 26^3 - 1 = (q^n - 1)/(q - 1) - 1, where q = 26, n = 4
        //   n = 1 + floor(log(19009)/log(26))
        // fourth = ((19009 - 18278                           )/26^0) % 26 = 3
        // third  = ((19009 - 18278 - 3*26^0                  )/26^1) % 26 = 2
        // second = ((19009 - 18278 - 3*26^0 - 2*26^1         )/26^2) % 26 = 1
        // first  = ((19009 - 18278 - 3*26^0 - 2*26^1 - 1*26^2)/26^3) % 26 = 0
                
        // Find how long a species name string must be
        double logibase26 = log(k)/log(26);
        unsigned n = 1 + (unsigned)floor(logibase26);
        vector<char> letters;
        unsigned base = (unsigned)((pow(26,n) - 1)/25.0 - 1);
        unsigned cum = 0;
        int ordA = (unsigned)(lower_case ? 'a' : 'A');
        for (unsigned i = 0; i < n; ++i) {
            unsigned ordi = (unsigned)((k - base - cum)/pow(26,i)) % 26;
            letters.push_back(char(ordA + ordi));
            cum += (unsigned)(ordi*pow(26,i));
        }
        string species_name(letters.rbegin(), letters.rend());
        return species_name;
    }
    
    inline void Proj::outputGeneTreesToFile(string fn, const vector<string> & newicks) const {
        assert(Forest::_ngenes == newicks.size());
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto newick : newicks) {
            streef << str(format("  tree %s = [&R] %s;\n") % Forest::_gene_names[t++] % newick);
        }
        streef << "end;\n";
        streef.close();
    }
    
    inline void Proj::outputNexusTreefile(string fn, const vector<string> & newicks) const {
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto newick : newicks) {
            streef << str(format("  tree tree%d = [&R] %s;\n") % (++t) % newick);
        }
        streef << "end;\n";
        streef.close();
    }
    
    inline void Proj::outputAnnotatedNexusTreefile(string fn, const vector<string> & newicks, const vector<string> & treenames, const vector<string> & annotations) const {
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto newick : newicks) {
            streef << str(format("  tree %s = [%s] [&R] %s;\n") % treenames[t] % annotations[t] % newick);
            ++t;
        }
        streef << "end;\n";
        streef.close();
    }
    
    inline void Proj::drawStartingSpeciesTree() {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;

        cout << "\nDrawing starting species tree from Yule prior:\n";
        cout << str(format("  lambda = %.5f (per-lineage speciation rate)\n") % Forest::_lambda);
        cout << str(format("  no. species = %d\n") % Forest::_nspecies);

        // Create species tree
        SpeciesForest sf;
        epoch_list_t epochs;
        sf.simulateSpeciesTree(epochs);
        _starting_species_newick = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
        Forest::createDefaultSpeciesTreeNexusTaxonMap();

        // Output tree file containing starting species tree
        string newick_with_names = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        outputNexusTreefile("starting-species-tree.tre", {newick_with_names});
        cout << "  Species tree saved in file \"starting-species-tree.tre\"\n";
    }
    
    inline void Proj::simulateData() {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;
        
        cout << "Simulating sequence data under multispecies coalescent model:\n";
        cout << str(format("  theta  = %.5f\n") % Forest::_theta);
        cout << str(format("  lambda = %.5f\n") % Forest::_lambda);
        cout << str(format("  no. species = %d\n") % _nsimspecies);
        
        // Interrogate _partition to determine number of genes and number of sites in each gene
        Forest::_ngenes = _partition->getNumSubsets();
        vector<unsigned> nsites_per_gene(Forest::_ngenes);
        Forest::_gene_names.resize(Forest::_ngenes);
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            nsites_per_gene[g] = _partition->numSitesInSubset(g);
            Forest::_gene_names[g] = _partition->getSubsetName(g);
        }
        //unsigned total_sites = accumulate(nsites_per_gene.begin(), nsites_per_gene.end(), 0);
        
        // Set taxon/species names and numbers, and create taxpartition vector used later
        // in paup command file
        Forest::_nspecies       = _nsimspecies;
        Forest::_ntaxa          = (unsigned)accumulate(_nsimtaxaperspecies.begin(), _nsimtaxaperspecies.end(), 0);
        Forest::_species_names.resize(Forest::_nspecies);
        Forest::_taxon_names.resize(Forest::_ntaxa);
        unsigned k = 0;
        vector<string> taxpartition;
        for (unsigned i = 0; i < Forest::_nspecies; ++i) {
            string species_name = inventName(i, false);
            Forest::_species_names[i] = species_name;
            for (unsigned j = 0; j < _nsimtaxaperspecies[i]; ++j) {
                taxpartition.push_back(Forest::_species_names[i]);
                string taxon_name = str(format("%s^%s") % inventName(k, true) % Forest::_species_names[i]);
                Forest::_taxon_names[k] = taxon_name;
                Forest::_taxon_to_species[taxon_name] = i;
                ++k;
            }
        }

        // Report species and taxon names created
        cout << "  Species names:\n";
        for (unsigned i = 0; i < Forest::_nspecies; ++i) {
            cout << "    " << Forest::_species_names[i] << "\n";
        }
        cout << "  Taxon names:\n";
        for (unsigned i = 0; i < Forest::_ntaxa; ++i) {
            string taxon_name = Forest::_taxon_names[i];
            cout << "    " << taxon_name << "\n";
        }
        
        cout << "Simulating species tree from Yule prior:\n";
        cout << str(format("  lambda = %.5f (per-lineage speciation rate)\n") % Forest::_lambda);
        cout << str(format("  no. species = %d\n") % Forest::_nspecies);

        // Create species tree
        SpeciesForest sf;
        epoch_list_t epochs;
        sf.simulateSpeciesTree(epochs);
        string newick_species_tree_alpha = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        string newick_species_tree_numeric = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);

        // Output tree file containing true species tree
        outputNexusTreefile("true-species-tree.tre", {newick_species_tree_alpha});
        cout << "  True species tree saved in file \"true-species-tree.tre\"\n";
        
        // Create data object
        assert(!_data);
        _data = Data::SharedPtr(new Data());
        _data->setTaxonNames(Forest::_taxon_names);
        _data->setPartition(_partition);
        
        // Create gene trees and simulate sequence data
        unsigned starting_site = 0;
        vector<string> newick_gene_trees_alpha;
        vector<string> newick_gene_trees_numeric;
        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
            Forest::debugShow(format("gene %d of %d") % g % Forest::_ngenes);
            ps.setNElements(4*nsites_per_gene[g], g);
            GeneForest gf;
            gf.simulateGeneTree(g, epochs);
            newick_gene_trees_alpha.push_back(gf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units));
            newick_gene_trees_numeric.push_back(gf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units));
            gf.simulateData(_data, starting_site, nsites_per_gene[g]);
            starting_site += nsites_per_gene[g];
        }
                       
        // Output data to file
        _data->compressPatterns();
        _data->writeDataToFile(_data_file_name);
        cout << "  Sequence data saved to file\"" << _data_file_name << "\"\n";
         
        // Output tree file containing true gene trees
        outputGeneTreesToFile("true-gene-trees.tre", newick_gene_trees_alpha);
        cout << "  True gene trees saved in file \"true-gene-trees.tre\"\n";
        
        // Output a PAUP* command file for estimating the species tree using
        // svd quartets and qage
        ofstream paupf("svd-qage.nex");
        paupf << "#NEXUS\n\n";
        paupf << "begin paup;\n";
        paupf << "  exe simulated.nex;\n";
        paupf << "  taxpartition species (vector) = " << join(taxpartition," ") << ";\n";
        paupf << "  svd taxpartition=species;\n";
        paupf << "  roottrees;\n";
        paupf << "  qage taxpartition=species patprob=exactjc outUnits=coalescent treefile=svd.tre replace;\n";
        paupf << "  quit;\n";
        paupf << "end;\n";
        paupf.close();
        cout << "  PAUP* commands saved in file\"svd-qage.nex\"\n";

        // Output gene trees and species trees for javascript viewer
        Particle::outputJavascriptTreefile("newicks.js", newick_species_tree_numeric, newick_gene_trees_numeric);
     }

}
