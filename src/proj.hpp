#pragma once

using boost::filesystem::current_path;
using boost::algorithm::join;
using boost::is_any_of;
using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::store;
using boost::program_options::parse_command_line;
using boost::program_options::parsed_options;
using boost::program_options::parse_config_file;
using boost::program_options::reading_file;
using boost::program_options::notify;

#if defined(USING_MPI)
extern int my_rank;
extern int ntasks;
#endif

#if defined(LOG_MEMORY)
extern ofstream memfile;
#endif

extern void output(string msg);
extern void output(string msg, unsigned level);
extern proj::PartialStore ps;
extern proj::StopWatch stopwatch;
extern proj::Lot rng;

namespace proj {

    class Proj {
        public:
                                       Proj();
                                       ~Proj();

            void                       clear();
            void                       processCommandLineOptions(int argc, const char * argv[]);
            void                       run();
            
            void                       memoryReport(ofstream & memf) const;
            
            static unsigned            _verbosity;

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
            void                       clearAllParticles();
            void                       clearGeneParticles(unsigned gene);
            void                       clearSpeciesParticles();
            void                       initializeGeneParticles(unsigned gene);
            void                       initializeSpeciesParticles();
            //double                     filterParticles(int gene, const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like);
            double                     filterGeneParticles(int gene, const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like);
            double                     filterSpeciesParticles(const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like);
            void                       saveStartingGeneTrees(string fn);
            void                       saveStartingSpeciesTree(string fn);
            void                       saveGeneTreesUsingNames(string filename, Particle & p);
            void                       saveGeneTreesUsingNumbers(string filename, Particle & p);
            //void                       propagateSampledParticle(bool growing_gene_forests);
            void                       debugCheckGeneTrees() const;

            void                       growGeneTrees(unsigned iter, unsigned gene);
            Particle &                 growSpeciesTrees(unsigned iter);
            void                       updateTheta(Particle & p, unsigned ntries, double delta);
            void                       updateLambda(Particle & p, unsigned ntries, double delta);

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
            
            bool                       _starting_species_tree_from_file;
            bool                       _starting_gene_trees_from_file;
            
            string                     _starting_species_newick;
            vector<string>             _starting_gene_newicks;
            
            bool                       _use_gpu;
            bool                       _ambig_missing;
            unsigned                   _nsimspecies;
            vector<unsigned>           _nsimtaxaperspecies;
            vector<unsigned>           _nsites_per_gene;
            unsigned                   _gene_nparticles;
            unsigned                   _species_nparticles;
            int                        _track_split;
            unsigned                   _rnseed;
            bool                       _sort_forests;
            double                     _visualization_cutoff;
            
            double                      _theta;
            double                      _lambda;
                        
            double                      _theta_delta;
            double                      _lambda_delta;

            vector<vect_particle_t>     _gene_particles;
            vect_particle_t             _species_particles;
            Particle                    _template_particle;
            
#if defined(USING_MPI)
            void mpiSetSchedule();
            vector<unsigned>            _mpi_first_gene;
            vector<unsigned>            _mpi_last_gene;
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
        _nsimspecies = 5;
        _nsimtaxaperspecies = {2,2,2,2,2};
        
        // Flags to allow us to preserve Forest::_nexus_taxon_map if
        // built from taxa block in tree file
        _starting_species_tree_from_file = false;
        _starting_gene_trees_from_file = false;
        
        _theta = 0.1;
        _lambda = 1.0;
        
        _theta_delta = 1.0;
        _lambda_delta = 20.0;
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
        ("startmode", value(&_start_mode), "if 'random', start with species tree drawn from prior; if 'species', start with first species tree defined in speciestreefile; if 'gene', start with gene trees defined in genetreefile; if 'evaluate', read both starting gene trees and starting species tree and compute likelihoods; if 'simulate', simulated gene trees, species tree, and data")
        ("niter", value(&_niter), "number of iterations, where one iteration involves SMC of gene trees give species tree combined with an SMC of species tree given gene trees")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("gpu",           value(&_use_gpu)->default_value(true), "use GPU if available")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("verbosity",  value(&_verbosity)->default_value(0), "0, 1, or 2: higher number means more output")
        ("nspecies",  value(&_nsimspecies)->default_value(1), "number of species (only used if simulate specified)")
        ("ntaxaperspecies",  value(&_nsimtaxaperspecies), "number of taxa sampled per species (only used if simulate specified); should be _nimspecies of these entries, one for each species simulated")
        ("ngeneparticles",  value(&_gene_nparticles)->default_value(100), "number of gene particles in a population")
        ("nspeciesparticles",  value(&_species_nparticles)->default_value(1000), "number of species particles in a population")
        ("theta",  value(&_theta)->default_value(0.05), "coalescent parameter assumed for gene trees")
        ("lambda",  value(&_lambda)->default_value(10.9), "per lineage speciation rate assumed for the species tree")
        ("thetapriormean",  value(&Forest::_theta_prior_mean)->default_value(0.05), "mean of exponential prior for theta")
        ("lambdapriormean",  value(&Forest::_lambda_prior_mean)->default_value(1.0), "mean of exponential prior for lambda")
        ("thetadelta",  value(&_theta_delta)->default_value(1.0), "mean of exponential prior for theta")
        ("lambdadelta",  value(&_lambda_delta)->default_value(20.0), "mean of exponential prior for lambda")
        ("updatetheta",  value(&Forest::_update_theta)->default_value(true), "if yes, update theta at the end of each iteration")
        ("updatelambda",  value(&Forest::_update_lambda)->default_value(true), "if yes, update lambda at the end of each iteration")
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
            output("Note: configuration file (proj.conf) not found\n");
        }
        notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            output(str(format("%s\n") % desc));
            exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            output(str(format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version));
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
        
        if (vm.count("theta") > 0) {
            Forest::_theta = _theta;
        }
        
        if (vm.count("lambda") > 0) {
            Forest::_lambda = _lambda;
        }
        
        if (vm.count("ntaxaperspecies") > 0) {
            if (_nsimspecies > 1) {
                unsigned ntaxaperspecies = (unsigned)_nsimtaxaperspecies.size();
                if (ntaxaperspecies == 1) {
                    unsigned n = _nsimtaxaperspecies[0];
                    for (unsigned i = 1; i < _nsimspecies; ++i) {
                        _nsimtaxaperspecies.push_back(n);
                    }
                }
                else if (ntaxaperspecies != _nsimspecies)
                    throw XProj(format("Expecting either 1 or %d ntaxaperspecies entries, but found %d") % _nsimspecies % ntaxaperspecies);
            }
        }
    }
    
    inline void Proj::readData() {
        output(str(format("\nReading and storing the data in the file %s\n") % _data_file_name));
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);
        
        Forest::_gene_names.clear();

        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        output(str(format("\nNumber of taxa: %d\n") % _data->getNumTaxa()));
        output(str(format("Number of partition subsets: %d") % nsubsets));
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(nsubsets);
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Set length of partials for gene g
            ps.setNElements(Forest::_nstates*_data->getNumPatternsInSubset(subset), subset);
            
            DataType dt = _partition->getDataTypeForSubset(subset);
            Forest::_gene_names.push_back(_data->getSubsetName(subset));
            output(str(format("  Subset %d (%s)\n") % (subset+1) % _data->getSubsetName(subset)));
            output(str(format("    data type: %s\n") % dt.getDataTypeAsString()));
            output(str(format("    sites:     %d\n") % _data->calcSeqLenInSubset(subset)));
            output(str(format("    patterns:  %d\n") % _data->getNumPatternsInSubset(subset)));
            }
    }
    
    inline void Proj::debugShowStringVector(string title, const vector<string> & svect) const {
#if defined(DEBUGGING)
        output(str(format("\n%s\n") % title));
        for (auto s : svect) {
            output(str(format("  %s\n") % s));
         }
#endif
    }
    
    inline void Proj::debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const {
#if defined(DEBUGGING)
        output(str(format("\n%s:\n") % title));
        for (auto su : sumap) {
            output(str(format("  %s: %d\n") % su.first % su.second));
         }
#endif
    }
    
    inline unsigned Proj::buildSpeciesMap() {
        // Populates Forest::_species_names and Forest::_taxon_to_species
        unsigned nspecies = 0;
        map<string, unsigned> species_name_to_index;
        Forest::_species_names.clear();
        Forest::_taxon_to_species.clear();
        
        output("\nMapping taxa to species\n");
        const Data::taxon_names_t & tnames = _data->getTaxonNames();
        unsigned ntax = (unsigned)tnames.size();
        for (auto & tname : tnames) {
            string species_name = Node::taxonNameToSpeciesName(tname);
            output(str(format("  %s --> %s\n") % tname % species_name));
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
        output(str(format("\nReading and storing the species tree in the file %s\n") % _species_tree_file_name));
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
        output(str(format("\nReading and storing the gene trees in the file %s\n") % _gene_tree_file_name));
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
        output(str(format("  no. trees: %d\n") % ntrees));
        for (unsigned i = 0; i < ntrees; i++) {
            string tree_name = tree_names[i];
            string gene_name = Forest::_gene_names[i];
            if (gene_name != tree_name) {
                throw XProj(format("Expecting name of %dth gene tree to be \"%s\" but it was instead \"%s\"") % i % gene_name % tree_name);
            }
            output(str(format("  tree %d has name %s\n") % (i+1) % tree_name));
        }
    }
    
    inline void Proj::showSettings() const {
        output(str(format("Speciation rate (lambda): %.9f\n") % Forest::_lambda));
        output(str(format("Coalescent parameter (theta): %.9f\n") % Forest::_theta));
        output(str(format("Number of gene particles: %d") % _gene_nparticles));
        output(str(format("Number of species particles: %d") % _species_nparticles));
    }
        
    inline void Proj::initializeGeneParticles(unsigned gene) {
        assert(_data);
        assert(Forest::_ngenes > 0);
        assert(gene < Forest::_ngenes);
        
        clearGeneParticles(gene);
        _gene_particles.resize(Forest::_ngenes);
        _gene_particles[gene].resize(_gene_nparticles);

        // Ensure that numbers in _starting_species_newick will be interpreted correctly
        if (!_starting_species_tree_from_file)
            Forest::createDefaultSpeciesTreeNexusTaxonMap();
            
        // //temporary!
        //cerr << "~~~> setting up _template_particle\n";

        // Set up a "template" particle the hard way
        _template_particle.clear();
        _template_particle.setGeneIndex(gene);
        _template_particle.setData(_data);
        _template_particle.initSpeciesForest(_starting_species_newick);
        _template_particle.digestSpeciesTree(/*append*/false);
        _template_particle.resetGeneForest();
        _template_particle.sortEpochs();
        
        // //temporary!
        //cerr << "~~~> copying _template_particle to " << _gene_particles[gene].size() << " particles\n";

        // Now use first as a template and copy to all other particles
        fill(_gene_particles[gene].begin(), _gene_particles[gene].end(), _template_particle);
    }
    
    inline void Proj::initializeSpeciesParticles() {
        assert(_data);
        assert(Forest::_ngenes > 0);
        output(str(format("  Initializing %d species particles...\n") % _species_nparticles));

        clearSpeciesParticles();
        _species_particles.resize(_species_nparticles);

        // Ensure that numbers in _starting_gene_newicks will be interpreted correctly
        if (!_starting_gene_trees_from_file)
            Forest::createDefaultGeneTreeNexusTaxonMap();
        
        // Set up a "template" particle the hard way
        _template_particle.clear();
        _template_particle.setGeneIndex(-1);
        _template_particle.setData(_data);
        _template_particle.resetSpeciesForest();
        _template_particle.initGeneForests(_starting_gene_newicks);
        _template_particle.digestGeneTrees(/*append*/false);
        _template_particle.sortEpochs();
        
        // Now use first as a template and copy to all other particles
        fill(_species_particles.begin(), _species_particles.end(), _template_particle);
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
    
    inline void Proj::saveStartingGeneTrees(string fn) {
        ofstream tmpf(fn);
        tmpf << "#nexus\n\n";
        tmpf << "begin trees;\n";
        tmpf << "  translate\n";
        for (unsigned t = 0; t < Forest::_ntaxa; ++t) {
            tmpf << str(format("    %d %s%s\n") % (t+1) % Forest::_taxon_names[t] % (t == Forest::_ntaxa - 1 ? ";" : ","));
        }
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            string gene_name = Forest::_gene_names[g];
            tmpf << str(format("  tree %s = [&R] %s;\n") % gene_name % _starting_gene_newicks[g]);
        }
        tmpf << "end;\n";
        tmpf.close();
    }
    
    inline void Proj::saveStartingSpeciesTree(string fn) {
        ofstream tmpf(fn);
        tmpf << "#nexus\n\n";
        tmpf << "begin trees;\n";
        tmpf << "  translate\n";
        for (unsigned t = 0; t < Forest::_nspecies; ++t) {
            tmpf << str(format("    %d %s%s\n") % (t+1) % Forest::_species_names[t] % (t == Forest::_nspecies - 1 ? ";" : ","));
        }
        tmpf << str(format("  tree species_tree = [&R] %s;\n") % _starting_species_newick);
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

#if 0
    inline void Proj::propagateSampledParticle(bool growing_gene_forests) {
        // Choose one particle at random
        assert(_nparticles == _particles.size());
        
        // Copy that particle to all particles
        unsigned which = rng.randint(0, _nparticles - 1);
        Particle & other = _particles[which];
        _template_particle = other;
        if (growing_gene_forests) {
            // Prepare for growing species forests
            _template_particle.resetSpeciesForest();
            _template_particle.digestGeneTrees(/*append*/false);
            _template_particle.sortEpochs();
        }
        else {
            // Just finished updating species tree, so update theta before next round
            updateTheta(_template_particle, 100, 0.5);
            
            // Prepare for growing gene forests
            _template_particle.resetGeneForests();
            _template_particle.digestSpeciesTree(/*append*/false);
            _template_particle.sortEpochs();
        }
        fill(_particles.begin(), _particles.end(), _template_particle);
    }
#endif
    
    inline void Proj::debugCheckGeneTrees() const {
#if defined(DEBUGGING)
        for (auto & p : _species_particles) {
            for (unsigned g = 0; g < Forest::_ngenes; g++) {
                const GeneForest & gf = p.getGeneForest(g);
                gf.debugCheckSpeciesSets();
            }
        }
#endif
    }
    
    inline void Proj::growGeneTrees(unsigned iter, unsigned gene) {
        // Grows forests for the given gene conditional on the species tree stored in each particle
        // Assumes gene forests are trivial and species forest has been digested into epochs

        //output(str(format("\nGrowing forests for gene %d (%s)...\n") % gene % Forest::_gene_names[gene]));

        // Ensure all of the particles for this gene have trivial gene forests and
        // a species tree built from _starting_species_newick
        initializeGeneParticles(gene);
        
        // Number of joins is one fewer than the number of taxa
        unsigned nsteps = Forest::_ntaxa - 1;
                
        // Loop over steps
        double log_marg_like = 0.0;
        vector<unsigned> counts;
        for (unsigned step = 0; step < nsteps; step++) {
            output(str(format("  step %d of %d: ") % (step + 1) % nsteps), 1);
            
            // Loop over particles, advancing gene forests one step
            vector<double> logw(_gene_particles[gene].size(), 0.0);
            unsigned p = 0;
            for (auto & particle : _gene_particles[gene]) {
                double log_weight = particle.advanceGeneForest(p, step);
                logw[p] += log_weight;
                ++p;
            }
            
            // Filter particles
            double ESS = filterGeneParticles(gene, logw, counts, log_marg_like);
            output(str(format("ESS = %.1f%%\n") % (100.9*ESS/_gene_nparticles)), 1);

        } // step loop
        
        // Choose one particle at random
        unsigned which = rng.randint(0, _gene_nparticles - 1);
        
#if defined(USING_MPI)
        // Send newick for selected gene tree to coordinator
        GeneForest & selected = _gene_particles[gene][which].getGeneForest();
        string newick = selected.makeNewick(/*precision*/9, /*use_names*/false, /*coalunits*/false);
        
        if (my_rank == 0) {
            // Copy selected gene tree to _starting_gene_newicks
            _starting_gene_newicks[gene] = newick;
        }
        else {
            // Send the newick string to the coordinator
            int msglen = 1 + newick.size();
            newick.resize(msglen);      // Adds \0 character to the end
            MPI_Send(&newick[0],        // void* data
                msglen,                 // int count,
                MPI_CHAR,               // MPI_Datatype datatype,
                0,                      // int destination,
                gene,                   // int tag,
                MPI_COMM_WORLD);        // MPI_Comm communicator)
        }
#else
        // Copy selected gene tree to _starting_gene_newicks
        GeneForest & selected = _gene_particles[gene][which].getGeneForest();
        _starting_gene_newicks[gene] = selected.makeNewick(/*precision*/9, /*use_names*/false, /*coalunits*/false);
        
        output(str(format("  log(marg. like.) = %.5f\n") % log_marg_like));
#endif
    }
    
    inline Particle & Proj::growSpeciesTrees(unsigned iter) {
        bool last_iter = (iter == _niter - 1);
        
        // Grows species trees conditional on the gene trees stored in each particle
        
        // Ensure all of the particles for this gene have trivial gene forests and
        // a species tree built from _starting_species_newick
        initializeSpeciesParticles();
        
        // Number of joins is one fewer than the number of species
        unsigned nsteps = Forest::_nspecies - 1;
                
        // Loop over steps
        double log_marg_like = 0.0;
        vector<unsigned> counts;
        for (unsigned step = 0; step < nsteps; step++) {
            output(str(format("  step %d of %d: ") % (step + 1) % nsteps), 1);
            vector<double> logw(_species_particles.size(), 0.0);

            // Loop over particles, advancing species forest one step in each
            unsigned p = 0;
            for (auto & particle : _species_particles) {
                double log_weight = particle.advanceSpeciesForest(p, step);
                logw[p] += log_weight;
                ++p;
            }
                        
            // Filter particles
            double ESS = filterSpeciesParticles(logw, counts, log_marg_like);
            output(str(format("ESS = %.1f%%\n") % (100.9*ESS/_species_nparticles)), 1);
        }

        output(str(format("  log(marg. like.) = %.5f\n") % log_marg_like));
        
        if (_verbosity > 1 || last_iter) {
            string fn = str(format("species-trees-after-iter-%d.tre") % iter);
            saveUniqueSpeciesTrees(fn, counts);
            
            auto maxit = max_element(counts.begin(), counts.end());
            assert(maxit != counts.end());
            unsigned maxi = (unsigned)distance(counts.begin(), maxit);
            output(str(format("  Best particle had index %d with count %d\n") % maxi % (*maxit)));
            
            string fnprefix = str(format("newicks-%d") % iter);
            double log_coal_like = _species_particles[maxi].debugSaveTreesAsJavascript(fnprefix);
            output(str(format("  log(coal. like.) after iter %d = %.9f\n") % iter % log_coal_like));
        }

        // Choose one particle at random and save its tree to _starting_species_tree
        unsigned which = rng.randint(0, _species_nparticles - 1);
        SpeciesForest & sf = _species_particles[which].getSpeciesForest();
        _starting_species_newick = sf.makeNewick(/*precision*/9, /*use_names*/false, /*coalunts*/false);
                
        return _species_particles[which];
    }
    
    inline void Proj::updateTheta(Particle & p, unsigned ntries, double delta) {
        if (!Forest::_update_theta)
            return;
            
        // Use multiple-try Metropolis to update theta conditional on the gene forests
        // and species forest defined in p. Uses the algorithm presented in
        // https://en.wikipedia.org/wiki/Multiple-try_Metropolis
        // assuming a symmetric proposal (so that w(x,y) = pi(x)).
        
        output("\nUpdating theta...\n");

        // r is the rate of the theta exponential prior
        double prior_rate = 1.0/Forest::_theta_prior_mean;
        double log_prior_rate = log(prior_rate);
        double log_prior = log_prior_rate - prior_rate*Forest::_theta;

        // theta0 is the current global theta value
        double theta0 = Forest::_theta;
        
        // Sample ntries new values of theta from symmetric proposal distribution
        // (window of width 2*delta centered on theta0). Compute weights (coalescent
        // likelihood) for each proposed_theta value.
        vector<double> proposed_thetas(ntries, 0.0);
        vector<double> logwstar(ntries, 0.0);
        for (unsigned i = 0; i < ntries; ++i) {
            double q = theta0 - delta + 2.0*delta*rng.uniform();
            if (q < 0.0)
                q = -q;
            proposed_thetas[i] = q;
            log_prior = log_prior_rate - prior_rate*q;
            logwstar[i] = p.calcLogCoalLikeGivenTheta(q) + log_prior;
        }
        
        // Compute log of the sum of the weights (this sum will form the
        // denominator of the acceptance ratio)
        double log_sum_denom_weights = Forest::calcLogSum(logwstar);

        // Normalize weights to create a discrete probability distribution
        vector<double> probs(ntries, 0.0);
        transform(logwstar.begin(), logwstar.end(), probs.begin(), [log_sum_denom_weights](double logw){
            return exp(logw - log_sum_denom_weights);
        });
        
        // Choose one theta value from the probability distribution
        unsigned which = Forest::multinomialDraw(probs);
        double theta_star = proposed_thetas[which];
        
        // Sample ntries-1 new values of theta from symmetric proposal distribution
        // (window of width 2*delta centered on theta_star)
        log_prior = log_prior_rate - prior_rate*theta0;
        logwstar[0] = p.calcLogCoalLikeGivenTheta(theta0) + log_prior;
        for (unsigned i = 1; i < ntries; ++i) {
            double q = theta_star - delta + 2.0*delta*rng.uniform();
            if (q < 0.0)
                q = -q;
            log_prior = log_prior_rate - prior_rate*q;
            logwstar[i] = p.calcLogCoalLikeGivenTheta(q) + log_prior;
        }
        
        // Compute log of the sum of the weights (this sum will form
        // the numerator of the acceptance ratio)
        double log_sum_numer_weights = Forest::calcLogSum(logwstar);

        // Compute acceptance ratio
        double logr = log_sum_numer_weights - log_sum_denom_weights;
        bool accept = true;
        if (logr < 0.0) {
            double logu = log(rng.uniform());
            accept = logu < logr;
        }
        if (accept) {
            Forest::_theta = theta_star;
            output(str(format("  New theta: %.5f\n") % Forest::_theta));
        }
        else {
            output(str(format("  Theta unchanged: %.5f\n") % Forest::_theta));
        }
    }
    
    inline void Proj::updateLambda(Particle & p, unsigned ntries, double delta) {
        if (!Forest::_update_lambda)
            return;

        // Use multiple-try Metropolis to update lambda conditional on the
        // species forest defined in p. Uses the algorithm presented in
        // https://en.wikipedia.org/wiki/Multiple-try_Metropolis
        // assuming a symmetric proposal (so that w(x,y) = pi(x)).
        
        output("\nUpdating lambda...\n");

        // r is the rate of the lambda exponential prior
        double prior_rate = 1.0/Forest::_lambda_prior_mean;
        double log_prior_rate = log(prior_rate);
        double log_prior = log_prior_rate - prior_rate*Forest::_lambda;
        double log_likelihood = 0.0;

        // lambda0 is the current global lambda value
        double lambda0 = Forest::_lambda;
        
        // Sample ntries new values of lambda from symmetric proposal distribution
        // (window of width 2*delta centered on lambda0). Compute weights (coalescent
        // likelihood) for each proposed_lambdas value.
        vector<double> proposed_lambdas(ntries, 0.0);
        vector<double> logwstar(ntries, 0.0);
        for (unsigned i = 0; i < ntries; ++i) {
            double l = lambda0 - delta + 2.0*delta*rng.uniform();
            if (l < 0.0)
                l = -l;
            proposed_lambdas[i] = l;
            log_prior = log_prior_rate - prior_rate*l;
            log_likelihood =  p.calcLogSpeciesTreeDensityGivenLambda(l);
            
            logwstar[i] = log_likelihood + log_prior;
        }
        
        // Compute log of the sum of the weights (this sum will form the
        // denominator of the acceptance ratio)
        double log_sum_denom_weights = Forest::calcLogSum(logwstar);

        // Normalize weights to create a discrete probability distribution
        vector<double> probs(ntries, 0.0);
        transform(logwstar.begin(), logwstar.end(), probs.begin(), [log_sum_denom_weights](double logw){
            return exp(logw - log_sum_denom_weights);
        });
        
        // Choose one lambda value from the probability distribution
        unsigned which = Forest::multinomialDraw(probs);
        double lambda_star = proposed_lambdas[which];
        
        // Sample ntries-1 new values of lambda from symmetric proposal distribution
        // (window of width 2*delta centered on lambda_star)
        log_prior = log_prior_rate - prior_rate*lambda0;
        log_likelihood =  p.calcLogSpeciesTreeDensityGivenLambda(lambda0);
        logwstar[0] = log_likelihood + log_prior;
        for (unsigned i = 1; i < ntries; ++i) {
            double l = lambda_star - delta + 2.0*delta*rng.uniform();
            if (l < 0.0)
                l = -l;
            log_prior = log_prior_rate - prior_rate*l;
            log_likelihood =  p.calcLogSpeciesTreeDensityGivenLambda(l);
            logwstar[i] = log_likelihood + log_prior;
        }
        
        // Compute log of the sum of the weights (this sum will form
        // the numerator of the acceptance ratio)
        double log_sum_numer_weights = Forest::calcLogSum(logwstar);

        // Compute acceptance ratio
        double logr = log_sum_numer_weights - log_sum_denom_weights;
        bool accept = true;
        if (logr < 0.0) {
            double logu = log(rng.uniform());
            accept = logu < logr;
        }
        if (accept) {
            Forest::_lambda = lambda_star;
            output(str(format("  New lambda: %.5f\n") % Forest::_lambda));
        }
        else {
            output(str(format("  Lambda unchanged: %.5f\n") % Forest::_lambda));
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
                unsigned nparticles = (unsigned)_species_particles.size();
                double pct = 100.0*c/nparticles;
                notes.push_back(str(format("This tree found in %d particles (%.1f%% of %d total particles)") % c % pct % nparticles));
                tree_names.push_back(str(format("tree%d-%d") % i % c));
                species_tree_newicks.push_back(_species_particles[p].getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false));
                ++i;
            }
            ++p;
        }
        outputAnnotatedNexusTreefile(fn, species_tree_newicks, tree_names, notes);
    }

#if 1
    inline double Proj::filterSpeciesParticles(const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like) {
        // Sanity checks
        assert(_species_particles.size() == _species_nparticles);
        assert(log_weights.size() == _species_nparticles);
        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = Forest::calcLogSum(log_weights);
        vector<double> probs(_species_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood
        log_marg_like += log_sum_weights - log(_species_nparticles);
        
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
        counts.resize(_species_nparticles);
        counts.assign(_species_nparticles, 0);
        
        // Throw _species_nparticles darts
        for (unsigned i = 0; i < _species_nparticles; ++i) {
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
                
                // Copy the "from" particle to the "to" particle
                _species_particles[index_to] = _species_particles[index_from];
                
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

    inline double Proj::filterGeneParticles(int gene, const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like) {
        assert(gene >= 0);
        
        // Sanity checks
        assert(_gene_particles[gene].size() == _gene_nparticles);
        assert(log_weights.size() == _gene_nparticles);
        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = Forest::calcLogSum(log_weights);
        vector<double> probs(_gene_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood
        log_marg_like += log_sum_weights - log(_gene_nparticles);
        
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
        counts.resize(_gene_nparticles);
        counts.assign(_gene_nparticles, 0);
        
        // Throw _gene_nparticles darts
        for (unsigned i = 0; i < _gene_nparticles; ++i) {
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
                
                // //temporary!
                //cerr << str(format("~~~> copying particle %d (count = %d) to particle %d\n") % index_from % (*copy_from_iter) % index_to);
                
                // Copy the "from" particle to the "to" particle
                _gene_particles[gene][index_to] = _gene_particles[gene][index_from];
                
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
#else
    inline double Proj::filterParticles(int gene, const vector<double> & log_weights, vector<unsigned> & counts, double & log_marg_like) {
        vect_particle_t & particles  = (gene >= 0 ? _gene_particles[gene] : _species_particles);
        unsigned          nparticles = (gene >= 0 ? _gene_nparticles      : _species_nparticles);
        
        // Sanity checks
        assert(particles.size() == nparticles);
        assert(log_weights.size() == nparticles);
        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = Forest::calcLogSum(log_weights);
        vector<double> probs(nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood
        log_marg_like += log_sum_weights - log(nparticles);
        
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
        counts.resize(nparticles);
        counts.assign(nparticles, 0);
        
        // Throw nparticles darts
        for (unsigned i = 0; i < nparticles; ++i) {
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
                particles[index_to] = particles[index_from];
                
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
#endif

    inline void Proj::clearAllParticles() {
        _species_particles.clear();
        _gene_particles.clear();
    }
    
    inline void Proj::clearGeneParticles(unsigned gene) {
        assert(gene < Forest::_ngenes);
        if (_gene_particles.size() > 0) {
            assert(_gene_particles.size() > gene);
            _gene_particles[gene].clear();
        }
    }
    
    inline void Proj::clearSpeciesParticles() {
        _species_particles.clear();
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

        output("\nDrawing starting species tree from Yule prior:\n");
        output(str(format("  lambda = %.5f (per-lineage speciation rate)\n") % Forest::_lambda));
        output(str(format("  no. species = %d\n") % Forest::_nspecies));

        // Create species tree
        SpeciesForest sf;
        epoch_list_t epochs;
        sf.simulateSpeciesTree(epochs);
        _starting_species_newick = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
        Forest::createDefaultSpeciesTreeNexusTaxonMap();

        // Output tree file containing starting species tree
        string newick_with_names = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        outputNexusTreefile("starting-species-tree.tre", {newick_with_names});
        output("  Species tree saved in file \"starting-species-tree.tre\"\n");
    }
    
    inline void Proj::simulateData() {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;
        
        output("Simulating sequence data under multispecies coalescent model:\n");
        output(str(format("  theta  = %.5f\n") % Forest::_theta));
        output(str(format("  lambda = %.5f\n") % Forest::_lambda));
        output(str(format("  no. species = %d\n") % _nsimspecies));
        
        // Interrogate _partition to determine number of genes and number of sites in each gene
        Forest::_ngenes = _partition->getNumSubsets();
        vector<unsigned> nsites_per_gene(Forest::_ngenes);
        Forest::_gene_names.resize(Forest::_ngenes);
        for (unsigned g = 0; g < Forest::_ngenes; g++) {
            nsites_per_gene[g] = _partition->numSitesInSubset(g);
            Forest::_gene_names[g] = _partition->getSubsetName(g);
        }
        
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
        output("  Species names:\n");
        for (unsigned i = 0; i < Forest::_nspecies; ++i) {
            output(str(format("    %s\n") % Forest::_species_names[i]));
        }
        cout << "  Taxon names:\n";
        for (unsigned i = 0; i < Forest::_ntaxa; ++i) {
            string taxon_name = Forest::_taxon_names[i];
            output(str(format("    %s\n") % taxon_name));
        }
        
        output("Simulating species tree from Yule prior:\n");
        output(str(format("  lambda = %.5f (per-lineage speciation rate)\n") % Forest::_lambda));
        output(str(format("  no. species = %d\n") % Forest::_nspecies));

        // Create species tree
        SpeciesForest sf;
        epoch_list_t epochs;
        sf.simulateSpeciesTree(epochs);
        string newick_species_tree_alpha = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        string newick_species_tree_numeric = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);

        // Output tree file containing true species tree
        outputNexusTreefile("true-species-tree.tre", {newick_species_tree_alpha});
        output("  True species tree saved in file \"true-species-tree.tre\"\n");
        
        // Create data object
        assert(!_data);
        _data = Data::SharedPtr(new Data());
        _data->setTaxonNames(Forest::_taxon_names);
        _data->setPartition(_partition);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(Forest::_ngenes);

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
        output("  True gene trees saved in file \"true-gene-trees.tre\"\n");
        
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
        output("  PAUP* commands saved in file\"svd-qage.nex\"\n");

        // Output gene trees and species trees for javascript viewer
        Particle::outputJavascriptTreefile("newicks.js", newick_species_tree_numeric, newick_gene_trees_numeric);
     }
     
#if defined(USING_MPI)
    inline void Proj::mpiSetSchedule() {
        // Determine which genes will be handled by this processor: e.g.
        // 20 = number of genes
        //  3 = number of processors
        //  6 = 20 / 3
        //  2 = 20 % 3
        //  rank 0 gets 6, rank 1 gets 7, rank 2 gets 7
        vector<unsigned> genes_per_task(ntasks, (unsigned)(Forest::_ngenes / ntasks));

        unsigned gene_remainder = Forest::_ngenes % ntasks;

        // Each rank > 0 gets an extra job if there is any remainder
        for (unsigned rank = 1; rank < ntasks; ++rank) {
            if (gene_remainder > 0) {
                genes_per_task[rank] += 1;
                --gene_remainder;
            }
        }

        _mpi_first_gene.resize(ntasks);
        _mpi_last_gene.resize(ntasks);
        unsigned gene_cum = 0;
        output("\nGene schedule:\n");
        output(str(format("%12s %25s %25s\n") % "rank" % "first gene" % "last gene"));
        for (unsigned rank = 0; rank < ntasks; ++rank) {
            _mpi_first_gene[rank] = gene_cum;
            _mpi_last_gene[rank]  = gene_cum + genes_per_task[rank];
            gene_cum += genes_per_task[rank];
            
            output(str(format("%12d %25d %25d\n") % rank % (_mpi_first_gene[rank] + 1) % _mpi_last_gene[rank]));
        }
    }
#endif

    inline void Proj::run() {
#if defined(USING_MPI)
        output("Starting MPI parallel version...\n");
        output(str(format("No. processors: %d\n") % ntasks));
#else
        output("Starting serial version...\n");
#endif
        output(str(format("Current working directory: %s\n") % current_path()));
        
        try {
            rng.setSeed(_rnseed);
            
            if (_start_mode == "simulate") {
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
                
                // //temporary!
                //cerr << "~~~> computeLeafPartials\n";
                
                GeneForest::computeLeafPartials(_data);
                                
                if (_start_mode == "random") {
#if defined(USING_MPI)
                    mpiSetSchedule();
#endif
                    drawStartingSpeciesTree();
                    for (unsigned iter = 0; iter < _niter; ++iter) {
                        bool last_iter = (iter == _niter - 1);
                        stopwatch.start();
                        
                        output(str(format("\nIteration %d of %d...\n") % (iter+1) % _niter));
                        
                        // After each gene forest is filtered for the last time,
                        // one is chosen at random and its newick is saved to _starting_gene_newicks
                        _starting_gene_newicks.clear();
                        _starting_gene_newicks.resize(Forest::_ngenes);
                        output("\nGrowing gene trees...\n");
#if defined(USING_MPI)
                        for (unsigned g = _mpi_first_gene[my_rank]; g < _mpi_last_gene[my_rank]; ++g) {
                            growGeneTrees(iter, g);
                        }
                        
                        if (my_rank == 0) {
                            // Make a list of genes that we haven't heard from yet
                            list<unsigned> outstanding;
                            for (unsigned rank = 1; rank < ntasks; ++rank) {
                                for (unsigned g = _mpi_first_gene[rank]; g < _mpi_last_gene[rank]; ++g) {
                                    outstanding.push_back(g);
                                }
                            }
                            
                            // Receive gene trees from genes being handled by other processors
                            // until outstanding is empty
                            while (!outstand[ing.empty()) {
                                // Probe to get message status
                                int message_length = 0;
                                MPI_Status status;
                                MPI_Probe(MPI_ANY_SOURCE,   // Source rank or MPI_ANY_SOURCE
                                    MPI_ANY_TAG,            // Message tag
                                    MPI_COMM_WORLD,         // Communicator
                                    &status                 // Status object
                                );
                            
                                // Get length of message
                                MPI_Get_count(&status, MPI_CHAR, &message_length);

                                // Get gene
                                unsigned gene = (unsigned)status.MPI_TAG;

                                // Get the message itself
                                string newick;
                                newick.resize(message_length);
                                MPI_Recv(&newick[0],    // Initial address of receive buffer
                                    message_length,     // Maximum number of elements to receive
                                    MPI_CHAR,           // Datatype of each receive buffer entry
                                    MPI_ANY_SOURCE,     // Rank of source
                                    MPI_ANY_TAG,        // Message tag
                                    MPI_COMM_WORLD,     // Communicator
                                    MPI_STATUS_IGNORE   // Status object
                                );
                                
                                // Store the newick and remove gene from outstanding
                                newick.resize(message_length - 1);  // remove '\0' at end
                                _starting_gene_newicks[gene] = newick;
                                                                
                                auto it = find(outstanding.begin(), outstanding.end(), gene);
                                assert(it != outstanding.end());
                                outstanding.erase(it);
                            }
                        }
                                                
                        string message;
                        int message_length;
                        if (my_rank == 0) {
                            if (_verbosity > 1 || last_iter) {
                                // Save the chosen gene trees to file
                                saveStartingGeneTrees(str(format("gene-trees-chosen-iter-%d-numbers.tre") % iter));
                            }
                                
                            // Estimate the species tree distribution (conditional on gene trees)
                            output("\nGrowing species tree...\n");
                            Particle & chosen_particle = growSpeciesTrees(iter);
                        
                            if (_verbosity > 1 || last_iter) {
                                // Save the chosen species tree to file
                                saveStartingSpeciesTree(str(format("species-tree-chosen-iter-%d.tre") % iter));
                            }
                            
                            // Update theta
                            updateTheta(chosen_particle, 100, _theta_delta);
                        
                            // Update lambda
                            updateLambda(chosen_particle, 100, _lambda_delta);
                        
                            if (_verbosity > 1 || last_iter) {
                                // Report log-likelihood of current parameter values
                                output("\nUsing current parameter values:\n");
                                double log_coal_like = chosen_particle.calcLogCoalLikeGivenTheta(Forest::_theta);
                                output(str(format("    log(coalescent likelihood) = %.5f\n") % log_coal_like));
                                double total_log_like = 0.0;
                                for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                                    double log_like = chosen_particle.calcLogLikelihoodForGene(g);
                                    total_log_like += log_like;
                                }
                                output(str(format("    log(likelihood) = %.5f\n") % total_log_like));
                            }
                            
                            // Broadcast gene trees, species tree, theta, and lambda all in one message
                            // with values separated by | characters
                            message += _starting_species_newick;
                            message += "|";
                            for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                                message += _starting_gene_newicks[g];
                                message += "|";
                            }
                            message += str(format("%.9f|") % Forest::_theta);
                            message += str(format("%.9f") % Forest::_lambda);
                            message_length = (int)message.size() + 1;
                        }
                        
                        // Broadcast the length of the message to all processes (including coordinator process 0)
                        MPI_Bcast(
                            &message_length, // Initial address
                            1,               // Number of elements
                            MPI_INT,         // Datatype of each send buffer element
                            0,               // Rank of broadcast root
                            MPI_COMM_WORLD   // Communicator
                        );
                        
                        // Broadcast the message itself
                        message.resize(message_length);
                        MPI_Bcast(
                            &message[0],     // Initial address
                            message_length,  // Number of elements
                            MPI_CHAR,        // Datatype of each send buffer element
                            0,               // Rank of broadcast root
                            MPI_COMM_WORLD   // Communicator
                        );
                        
                        // Split message into component parts and save each part if not coordinator
                        // (who already knows all this information)
                        if (my_rank > 0) {
                            message.resize(message_length - 1);
                            vector<string> parts;
                            split(parts, message, is_any_of("|"));
                            assert(parts.size() == Forest::_ngenes + 3);
                            _starting_species_newick = parts[0];
                            assert(_starting_gene_newicks.size() == Forest::_ngenes);
                            for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                                _starting_gene_newicks[g] = parts[g+1];
                            }
                            Forest::_theta = stod(parts[Forest::_ngenes + 1]);
                            Forest::_lambda = stod(parts[Forest::_ngenes + 2]);
                        }
                        
                        double secs = stopwatch.stop();
                        output(str(format("\nTime used for iteration %d was %.5f seconds\n") % iter % secs));
#else
                        for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                            // //temporary!
                            //cerr << "~~~> growing gene tree for gene " << g << "\n";
                            
                            growGeneTrees(iter, g);
                        }
                        
                        if (_verbosity > 1 || last_iter) {
                            // Save the chosen gene trees to file
                            saveStartingGeneTrees(str(format("gene-trees-chosen-iter-%d-numbers.tre") % iter));
                        }
                        
                        // Estimate the species tree distribution (conditional on gene trees)
                        Particle & chosen_particle = growSpeciesTrees(iter);
                        
                        if (_verbosity > 1 || last_iter) {
                            // Save the chosen species tree to file
                            saveStartingSpeciesTree(str(format("species-tree-chosen-iter-%d.tre") % iter));
                        }
                        
                        // Update theta
                        updateTheta(chosen_particle, 100, 0.5);
                        
                        // Update lambda
                        updateLambda(chosen_particle, 100, 0.5);
                        
                        if (_verbosity > 1 || last_iter) {
                            // Report log-likelihood of current parameter values
                            output("\nUsing current parameter values:\n");
                            double log_coal_like = chosen_particle.calcLogCoalLikeGivenTheta(Forest::_theta);
                            output(str(format("    log(coalescent likelihood) = %.5f\n") % log_coal_like));
                            double total_log_like = 0.0;
                            for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                                double log_like = chosen_particle.calcLogLikelihoodForGene(g);
                                total_log_like += log_like;
                            }
                            output(str(format("    log(likelihood) = %.5f\n") % total_log_like));
                        }
                        
                        double secs = stopwatch.stop();
                        output(str(format("\nTime used for iteration %d was %.5f seconds\n") % iter % secs));
                        
#if defined(LOG_MEMORY)
                        memfile << "\n**** after iteration " << iter << "\n\n";
                        ps.memoryReport(memfile);
#endif
#endif
                            
                    }

#if defined(USING_MPI)
                    // Ensure no one starts on next iteration until coordinator is ready
                    MPI_Barrier(MPI_COMM_WORLD);
#endif
                }
                else if (_start_mode == "species") {
                    readStartingSpeciesTree();
                    _starting_species_tree_from_file = true;
                    for (unsigned iter = 0; iter < _niter; ++iter) {
                        output(str(format("\nIteration %d of %d...\n") % (iter+1) % _niter));
                        for (unsigned g = 0; g < Forest::_ngenes; ++g)
                            growGeneTrees(iter, g);
                        growSpeciesTrees(iter);
                    }
                }
                else if (_start_mode == "gene") {
                    readStartingGeneTrees();
                    _starting_gene_trees_from_file = true;
                    for (unsigned iter = 0; iter < _niter; ++iter) {
                        output(str(format("\nIteration %d of %d...\n") % (iter+1) % _niter));
                        growSpeciesTrees(iter);
                        for (unsigned g = 0; g < Forest::_ngenes; ++g)
                            growGeneTrees(iter, g);
                    }
                    output("\nFinal species tree SMC...\n");
                    growSpeciesTrees(_niter);
                }
                else if (_start_mode == "evaluate") {
                    Particle p;
                    p.setGeneIndex(-1);
                    p.setData(_data);

                    readStartingSpeciesTree();
                    p.initSpeciesForest(_starting_species_newick);
                    p.digestSpeciesTree(/*append*/false);

                    readStartingGeneTrees();
                    p.initGeneForests(_starting_gene_newicks);
                    p.digestGeneTrees(/*append*/true);
                    
                    p.sortEpochs();
                    p.reconcileEpochs();
                    p.copyEpochsToSpeciesForest();
                    
                    double theta = Forest::_theta;
                    double log_coal_like = p.calcLogCoalLikeGivenTheta(theta);
                    output(str(format("\nlog(coalescent likelihood) = %.9f\n") % log_coal_like));
                    
                    double total_log_like = 0.0;
                    for (unsigned g = 0; g < Forest::_ngenes; ++g) {
                        double log_like = p.calcLogLikelihoodForGene(g);
                        total_log_like += log_like;
                        output(str(format("log(likelihood for %s) = %.9f\n") % Forest::_gene_names[g] % log_like));
                    }
                    output(str(format("Total log(likelihood) = %.9f\n") % total_log_like));
                }
                else if (_start_mode == "lambdatest") {
                    Forest::_nspecies = _nsimspecies;
                    Forest::_species_names.resize(_nsimspecies);
                    for (unsigned i = 0; i < Forest::_nspecies; ++i) {
                        string species_name = inventName(i, false);
                        Forest::_species_names[i] = species_name;
                    }
                    
                    output(str(format("Starting lambda   = %.1f\n") % Forest::_lambda));
                    output(str(format("Mean lambda prior = %.1f\n") % Forest::_lambda_prior_mean));

                    drawStartingSpeciesTree();

                    Particle p;
                    p.setGeneIndex(-1);
                    p.initSpeciesForest(_starting_species_newick);
                    
                    double invseries = 0.0;
                    for (unsigned z = 2; z <= Forest::_nspecies; z++)
                        invseries += 1.0/z;
                    output(str(format("%12.5f = 1/2 + 1/3 + ... + 1/%d\n") % invseries % Forest::_nspecies));
                    double species_tree_height = p.speciesForestHeight();
                    output(str(format("%12.5f = species tree height\n") % species_tree_height));
                    double expected_species_tree_height = invseries/_lambda;
                    output(str(format("%12.5f = expected species tree height (given lambda = %.5f)\n") % expected_species_tree_height % _lambda));
                    double lambda_crude = invseries/species_tree_height;
                    output(str(format("%12.5f = lambda crude estimate\n") % lambda_crude));
                    
                    // Place starting lambda far away
                    Forest::_lambda += 20.0;
                    
                    double logLtrue =  p.calcLogSpeciesTreeDensityGivenLambda(_lambda);
                    output(str(format("%12.5f = Log likelihood given true lambda (%.5f)\n") % logLtrue % _lambda));

                    double logLcrude =  p.calcLogSpeciesTreeDensityGivenLambda(lambda_crude);
                    output(str(format("%12.5f = Log likelihood given crude estimate of lambda (%.5f)\n") % logLcrude % lambda_crude));

                    double logLstart =  p.calcLogSpeciesTreeDensityGivenLambda(Forest::_lambda);
                    output(str(format("%12.5f = Log likelihood given starting lambda (%.5f)\n") % logLstart % Forest::_lambda));
                    
                    for (unsigned iter = 0; iter < _niter; ++iter) {
                        updateLambda(p, 100, 1.0);
                    }
                }
                else {
                    throw XProj(format("Unknown start mode (\"%s\")") % _start_mode);
                }
            }
        }
        catch (XProj & x) {
            output(str(format("Proj encountered a problem:\n  %s\n") % x.what()));
        }

        output("\nFinished!\n");
    }
    
    void Proj::memoryReport(ofstream & memf) const {
        memf << "\nProj memory report:\n\n";
        memf << str(format("  Number of gene particles: %d\n") % _gene_nparticles);
        memf << str(format("  Number of species particles: %d\n") % _species_nparticles);
        _data->memoryReport(memf);
    }
    
}
