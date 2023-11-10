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

#if defined(LOG_MEMORY)
extern char dummy_char;
extern ofstream memfile;
#endif

#if defined(USING_SIGNPOSTS)
extern os_log_t log_handle;
extern os_signpost_id_t signpost_id;
#endif

#if defined(USING_MPI)
extern int my_rank;
extern int ntasks;
extern unsigned my_first_particle;
extern unsigned my_last_particle;
#endif

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
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
            
        private:
                    
            void updateCountMap(map<string, unsigned> & count_map, string & newick,    unsigned count);
            unsigned countSpeciations(const vector<Particle> & particles);
            void debugCheckCountMap(const map<string, unsigned> & count_map, unsigned target);
            void parseRelRateDefinition(string & s);
            
            void debugShowSpecies(double speciation_incr, vector<SMCGlobal::species_tuple_t> & species_tuples, unsigned which_gene, SMCGlobal::species_t which_spp);
            void                       debugShowStringVector(string title, const vector<string> & svect) const;
            void                       debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const;

            static string              inventName(unsigned k, bool lower_case);
            void                       outputNexusTreefile(string fn, const vector<string> & newicks) const;
            
            void                       simulateData();
            void                       simulateTrees(Particle & particle);
            void                       outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect);
            void                       outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric);
            void                       outputAnnotatedNexusTreefile(string fn, const vector<tuple<unsigned, double, string, string, string> > & treeinfo) const;
            void                       saveAllSpeciesTrees(string fn, const vector<Particle> particles, const map<string, unsigned> & count_map);
            void                       saveAllGeneTrees(unsigned gene, string fn, const vector<Particle> particles);
            void                       saveUniqueSpeciesTrees(string fn, const vector<Particle> particles, const vector<unsigned> & counts);

#if defined(USING_MULTITHREADING)
            void                       advanceParticleRange(unsigned step, unsigned first, unsigned last);
#endif
//            double                     advanceParticle(unsigned step, unsigned pindex, Particle & particle, bool compute_partials);
            double                     filterParticles(unsigned step, vector<Particle> & particles, vector<double> & log_weights, vector<unsigned> & counts, map<string, unsigned> & count_map);

            void                       readData();
            void                       setRelativeRates();
            unsigned                   buildSpeciesMap();
            void                       outputGeneTreesToFile(string fn,
                                            const vector<string> & newicks) const;
            void                       showSettings() const;
            
            string                     _data_file_name;
            string                     _start_mode;
            unsigned                   _niter;
            Partition::SharedPtr       _partition;
            Data::SharedPtr            _data;
            
            bool                       _use_gpu;
            bool                       _ambig_missing;
            unsigned                   _nsimspecies;
            vector<unsigned>           _nsimtaxaperspecies;
            vector<unsigned>           _nsites_per_gene;
            unsigned                   _nparticles;
            int                        _track_split;
            unsigned                   _rnseed;
            bool                       _sort_forests;
            double                     _visualization_cutoff;
            
            vector<Particle>            _particles;
            vector<double>              _log_weights;
            
            map<string, double>         _relrate_map;
            
            double                      _theta;
            double                      _lambda;
                        
            double                      _theta_delta;
            double                      _lambda_delta;
            
            double                      _ntries_theta;
            double                      _ntries_lambda;

            //vector<vect_particle_t>     _gene_particles;
            //vect_particle_t             _species_particles;
            //Particle                    _template_particle;
            
#if defined(USING_MULTITHREADING)
            void threadSetSchedule();
#endif

#if defined(USING_MPI)
            void mpiSetSchedule();
            vector<unsigned>            _mpi_first_particle;
            vector<unsigned>            _mpi_last_particle;
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
        _start_mode             = "smc";
        _niter                  = 1;
        _partition.reset(new Partition());
        
        // simulation related
        _nsimspecies = 5;
        _nsimtaxaperspecies = {2,2,2,2,2};
                
        _theta = 0.1;
        _lambda = 1.0;
        
        _theta_delta = 1.0;
        _lambda_delta = 20.0;
        
        _ntries_theta = 100;
        _ntries_lambda = 100;
        
        _nparticles = 1000;
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> partition_subsets;
        vector<string> partition_relrates;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile",  value(&_data_file_name)->required(), "name of a data file in NEXUS format")
        ("startmode", value(&_start_mode), "if 'simulate', simulate gene trees, species tree, and data; if 'smc', estimate from supplied datafile")
        ("niter", value(&_niter), "number of iterations, where one iteration involves SMC of gene trees give species tree combined with an SMC of species tree given gene trees")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("relrate",  value(&partition_relrates), "a relative rate for a previously-defined subset; format first:3.1")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("verbosity",  value(&SMCGlobal::_verbosity)->default_value(0), "0, 1, or 2: higher number means more output")
        ("nspecies",  value(&_nsimspecies)->default_value(1), "number of species (only used if simulate specified)")
        ("ntaxaperspecies",  value(&_nsimtaxaperspecies), "number of taxa sampled per species (only used if simulate specified); should be _nimspecies of these entries, one for each species simulated")
        ("nparticles",  value(&_nparticles)->default_value(1000), "number of particles in a population")
        ("nthreads",  value(&SMCGlobal::_nthreads)->default_value(3), "number of threads (each thread will handle nparticles/nthreads particle updates)")
        ("priorprior", value(&SMCGlobal::_prior_prior)->default_value(true), "use prior-prior approach to choose coalescence joins")
        ("theta",  value(&_theta)->default_value(0.05), "coalescent parameter assumed for gene trees")
        ("lambda",  value(&_lambda)->default_value(10.9), "per lineage speciation rate assumed for the species tree")
        ("thetapriormean",  value(&SMCGlobal::_theta_prior_mean)->default_value(0.05), "mean of exponential prior for theta")
        ("lambdapriormean",  value(&SMCGlobal::_lambda_prior_mean)->default_value(1.0), "mean of exponential prior for lambda")
        ("ntriestheta",  value(&_ntries_theta)->default_value(100), "number of multiple-try Metropolis trials when updating theta")
        ("ntrieslambda",  value(&_ntries_lambda)->default_value(100), "number of multiple-try Metropolis trials when updating lambda")
        ("thetadelta",  value(&_theta_delta)->default_value(1.0), "mean of exponential prior for theta")
        ("lambdadelta",  value(&_lambda_delta)->default_value(20.0), "mean of exponential prior for lambda")
        ("updatetheta",  value(&SMCGlobal::_update_theta)->default_value(true), "if yes, update theta at the end of each iteration")
        ("updatelambda",  value(&SMCGlobal::_update_lambda)->default_value(true), "if yes, update lambda at the end of each iteration")
        ("rnseed",  value(&_rnseed)->default_value(13579), "pseudorandom number seed")
        ("sortforests",  value(&_sort_forests)->default_value(false), "sort forests by weight when saving to file")
        ("visualizationcutoff", value(&_visualization_cutoff)->default_value(0.99), "particles sorted from highest to lowest weight will be saved for visualization if cumulative weight is less than this value")
        ;
        
        store(parse_command_line(argc, argv, desc), vm);
        try {
            const parsed_options & parsed = parse_config_file< char >("proj.conf", desc, false);
            store(parsed, vm);
        }
        catch(reading_file & x) {
            throw XProj("Configuration file (proj.conf) not found\n");
        }
        notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            output(format("%s\n") % desc, 2);
            exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            output(format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version, 2);
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
        
        // If user specified --relrate on command line, break specified relrate
        // definition into name and rate
        if (vm.count("relrate") > 0) {
            for (auto relrate_definition : partition_relrates) {
                parseRelRateDefinition(relrate_definition);
            }
        }
        
        if (vm.count("theta") > 0) {
            SMCGlobal::_theta = _theta;
        }
        
        if (vm.count("lambda") > 0) {
            SMCGlobal::_lambda = _lambda;
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
        
#if defined(USING_MULTITHREADING)
        // nothing to do
#else
        if (vm.count("nthreads") > 0) {
            if (SMCGlobal::_nthreads != 1) {
                output(format("\nWARNING: You specified %d threads but this non-multithreading version only allows 1 thread\nProceeding with a single thread.\n\n")  % SMCGlobal::_nthreads,1);
                SMCGlobal::_nthreads = 1;
            }
        }
#endif
    }
    
    inline void Proj::parseRelRateDefinition(string & s) {
        vector<string> v;
        
        // First separate part before colon (stored in v[0])
        // from part after colon (stored in v[1])
        split(v, s, boost::is_any_of(":"));
        if (v.size() != 2)
            throw XProj("Expecting exactly one colon in relrate definition");

        string gene_name = v[0];
        string relrate_string = v[1];
        double relrate;
        try {
            relrate = stod(relrate_string);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert \"%s\" to a relative rate for gene \"%s\"") % relrate_string % gene_name);
        }
        _relrate_map[gene_name] = relrate;
    }
    
    inline void Proj::readData() {
        output(format("\nReading and storing the data in the file %s\n") % _data_file_name, 2);
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);
        
        SMCGlobal::_gene_names.clear();

        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        output(format("\nNumber of taxa: %d\n") % _data->getNumTaxa(), 2);
        output(format("Number of partition subsets: %d") % nsubsets, 2);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(nsubsets);
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Set length of partials for gene g
            ps.setNElements(SMCGlobal::_nstates*_data->getNumPatternsInSubset(subset), subset);
            
            DataType dt = _partition->getDataTypeForSubset(subset);
            SMCGlobal::_gene_names.push_back(_data->getSubsetName(subset));
            output(format("  Subset %d (%s)\n") % (subset+1) % _data->getSubsetName(subset), 2);
            output(format("    data type: %s\n") % dt.getDataTypeAsString(), 2);
            output(format("    sites:     %d\n") % _data->calcSeqLenInSubset(subset), 2);
            output(format("    patterns:  %d\n") % _data->getNumPatternsInSubset(subset), 2);
        }
    }
    
    inline void Proj::setRelativeRates() {
        bool relrates_specified = !_relrate_map.empty();
        if (relrates_specified) {
            output("\nRelative rates specified for each gene:\n",2);
            output(format("%12s %12s\n") % "gene" % "rate",2);
            output(format("%12s %12s\n") % "-----------" % "-----------",2);
            
            SMCGlobal::_relrate_for_gene.clear();
            unsigned total_nsites = _partition->getNumSites();
            double mean_rate = 0.0;
            for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                unsigned gnsites = _partition->numSitesInSubset(g);
                string gname = _partition->getSubsetName(g);
                double r = _relrate_map.at(gname);
                SMCGlobal::_relrate_for_gene[g] = r;
                output(format("%12s %12.5f\n") % gname % r,2);
                mean_rate += r*gnsites/total_nsites;
            }
            output(format("%12s %12s\n") % "-----------" % "-----------",2);
            output(format("%12s %12.5f\n") % "mean" % mean_rate,2);
            if (fabs(mean_rate - 1.0) > 0.001) {
                XProj("The mean rate is more than 0.001 away from 1.0");
            }
        }
        else {
            output("\nRelative rates not specified (assuming no rate heterogeneity across genes)",2);
            for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                SMCGlobal::_relrate_for_gene[g] = 1.0;
            }
        }
    }

    inline unsigned Proj::buildSpeciesMap() {
        // Populates Forest::_species_names and Forest::_taxon_to_species
        unsigned nspecies = 0;
        map<string, unsigned> species_name_to_index;
        SMCGlobal::_species_names.clear();
        SMCGlobal::_taxon_to_species.clear();
        
        output("\nMapping taxon names to species\n", 2);
        const Data::taxon_names_t & tnames = _data->getTaxonNames();
        unsigned ntax = (unsigned)tnames.size();
        for (auto & tname : tnames) {
            string species_name = Node::taxonNameToSpeciesName(tname);
            output(format("  %s --> %s\n") % tname % species_name, 2);
            unsigned species_index = ntax;
            if (species_name_to_index.find(species_name) == species_name_to_index.end()) {
                // species_name not found
                species_index = nspecies;
                SMCGlobal::_species_names.push_back(species_name);
                species_name_to_index[species_name] = nspecies++;
            } else {
                // species_name found
                species_index = species_name_to_index[species_name];
            }
            SMCGlobal::_taxon_to_species[tname] = species_index;
        }

        output("\nMapping species names to species index\n", 2);
        for (auto & sname : SMCGlobal::_species_names) {
            unsigned species_index = species_name_to_index.at(sname);
            output(format("  %s --> %d\n") % sname % species_index, 2);
            SMCGlobal::_taxon_to_species[sname] = species_index;
        }

        return (unsigned)SMCGlobal::_species_names.size();
    }
    
    inline void Proj::debugShowStringVector(string title, const vector<string> & svect) const {
        output(format("\n%s\n") % title, 2);
        for (auto s : svect) {
            output(format("  %s\n") % s, 2);
         }
    }
    
    inline void Proj::debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const {
        output(format("\n%s:\n") % title, 2);
        for (auto su : sumap) {
            output(format("  %s: %d\n") % su.first % su.second, 2);
         }
    }
    
    inline void Proj::outputGeneTreesToFile(string fn, const vector<string> & newicks) const {
        assert(SMCGlobal::_ngenes == newicks.size());
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto newick : newicks) {
            streef << str(format("  tree %s = [&R] %s;\n") % SMCGlobal::_gene_names[t++] % newick);
        }
        streef << "end;\n";
        streef.close();
    }
    
    inline void Proj::showSettings() const {
        output(format("Speciation rate (lambda): %.9f\n") % SMCGlobal::_lambda, 2);
        output(format("Coalescent parameter (theta): %.9f\n") % SMCGlobal::_theta, 2);
        output(format("Number of particles: %d") % _nparticles, 2);
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
    
    inline void Proj::debugShowSpecies(double speciation_incr, vector<SMCGlobal::species_tuple_t> & species_tuples, unsigned which_gene, SMCGlobal::species_t which_spp) {
        output(format("\ndebugShowSpecies (speciation increment = %.5f):\n") % speciation_incr, 1);

        unsigned ntuples = (unsigned)species_tuples.size();
        output(format("%12s %12s %s\n")  % "nlineages" % "gene" % "species", 1);
        
        for (unsigned i = 0; i < ntuples; i++) {
            SMCGlobal::species_tuple_t & species_tuple = species_tuples[i];
            unsigned              n = get<0>(species_tuple);
            unsigned              g = get<1>(species_tuple);
            SMCGlobal::species_t  s = get<2>(species_tuple);
            string               ss = SMCGlobal::speciesStringRepresentation(s);
            assert(g > 0);
            output(format("%12d %12d %s\n")  % n % g % ss, 1);
        }
    }
    
    inline void Proj::simulateTrees(Particle & particle) {
        unsigned nspecies = SMCGlobal::_nspecies;
        unsigned ngenes   = SMCGlobal::_ngenes;
        unsigned ntaxa    = SMCGlobal::_ntaxa;
        
        // Start with trivial species forest
        particle.resetSpeciesForest();

        // Start with trivial gene forests
        particle.resetGeneForests(/*compute_partials*/false);
        
        assert(particle.getGeneForests().size() == ngenes);
        
        // Determine total number of steps required to build all gene trees and the species tree
        unsigned nsteps = (nspecies - 1) + ngenes*(ntaxa - 1);
        
        // Add coalescent events and speciation events until all trees are built
        for (unsigned step = 0; step < nsteps; step++) {
            particle.advance(step, 0, /*calculate_partial*/false);
        }
        
        particle.refreshHeightsInternalsPreorders();
    }
    
    inline void Proj::outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect) {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;
        
        // Output tree file containing true species tree
        string newick_species_tree_alpha = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        string newick_species_tree_numeric = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
        if (SMCGlobal::_nspecies > 1) {
            outputNexusTreefile("true-species-tree.tre", {newick_species_tree_alpha});
            output(format("  True species tree height = %.9f\n") % sf.getHeight(),2);
            output("  True species tree saved in file \"true-species-tree.tre\"\n", 1);
        }
        else {
            output("  True species tree not saved because it is just a single node!\n", 1);
        }

        // Output tree file containing true gene trees
        vector<string> newick_gene_trees_alpha;
        vector<string> newick_gene_trees_numeric;
        output("  True gene tree height (height/theta):\n",2);
        for (auto & gf : gfvect) {
            string newick_alpha = gf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
            newick_gene_trees_alpha.push_back(newick_alpha);
            
            string newick_numeric = gf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
            newick_gene_trees_numeric.push_back(newick_numeric);
            output(format("  %12.9f (%.9f)\n") % gf.getHeight() % (gf.getHeight()/SMCGlobal::_theta), 2);
        }
        outputGeneTreesToFile("true-gene-trees.tre", newick_gene_trees_alpha);
        output("  True gene trees saved in file \"true-gene-trees.tre\"\n", 1);

        // Output gene trees and species trees for javascript viewer
        outputJavascriptTreefile("newicks.js", newick_species_tree_numeric, newick_gene_trees_numeric);
    }
    
    inline void Proj::outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric) {
        ofstream jsf(fn);
        jsf << "let species_translate = {\n";
        unsigned i = 1;
        for (string nm : SMCGlobal::_species_names) {
            string comma = (i < SMCGlobal::_nspecies ? "," : "");
            jsf << str(format("  %d: \"%s\"%s\n") % i % nm % comma);
            ++i;
        }
        jsf << "};\n\n";
        
        jsf << str(format("let species_newick = \"%s\";\n\n") % newick_species_tree_numeric);
        
        jsf << "let gene_translate = {\n";
        i = 1;
        for (string nm : SMCGlobal::_taxon_names) {
            string comma = (i < SMCGlobal::_ntaxa ? "," : "");
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
    
    inline void Proj::simulateData() {
        SMCGlobal::_nspecies = _nsimspecies;
        SMCGlobal::_ntaxa = (unsigned)accumulate(_nsimtaxaperspecies.begin(), _nsimtaxaperspecies.end(), 0);

        output("Simulating sequence data under multispecies coalescent model:\n", 2);
        output(format("  theta  = %.5f\n") % SMCGlobal::_theta, 2);
        output(format("  lambda = %.5f\n") % SMCGlobal::_lambda, 2);
        output(format("  Number of species = %d\n") % _nsimspecies, 2);
        
        // Expected height of species tree is (1/2 + 1/3 + ... + 1/nspecies)/lambda
        // Approximation to (1/2 + 1/3 + ... + 1/nspecies) is
        //   ln(nspecies) + 0.58 (Euler's constant)
        //   see http://www.dimostriamogoldbach.it/en/inverses-integers-sum/
        // 0.25 = Exact (nspecies = 5, lambda = 5)
        // 0.44 = Approximate
        double expected_species_tree_height = 0.0;
        for (unsigned i = 2; i <= SMCGlobal::_nspecies; i++) {
            expected_species_tree_height += 1.0/i;
        }
        expected_species_tree_height /= SMCGlobal::_lambda;
        output(format("  Expected species tree height = %.5f\n") % expected_species_tree_height, 2);
        
        // Expected gene tree height without species tree constraints
        // is theta/n(n-1) + theta/(n-1)(n-2) + ... + theta/(2)(1)
        double expected_gene_tree_height = 0.0;
        for (unsigned i = 2; i <= SMCGlobal::_ntaxa; i++) {
            expected_gene_tree_height += 1.0/(i*(i-1));
        }
        expected_gene_tree_height *= SMCGlobal::_theta;
        output(format("  Expected gene tree height (unconstrained) = %.5f\n") % expected_gene_tree_height, 2);
        
        // Interrogate _partition to determine number of genes, gene names, and
        // number of sites in each gene
        SMCGlobal::_ngenes = _partition->getNumSubsets();
        SMCGlobal::_nsites_per_gene.resize(SMCGlobal::_ngenes);
        SMCGlobal::_gene_names.resize(SMCGlobal::_ngenes);
        unsigned total_nsites = 0;
        for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
            unsigned gnsites = _partition->numSitesInSubset(g);
            total_nsites += gnsites;
            string gname = _partition->getSubsetName(g);
            SMCGlobal::_nsites_per_gene[g] = gnsites;
            SMCGlobal::_gene_names[g] = gname;
        }

        setRelativeRates();

        // Invent taxon/species names and numbers, and create
        // taxpartition vector used later in paup command file
        SMCGlobal::_species_names.resize(SMCGlobal::_nspecies);
        SMCGlobal::_taxon_names.resize(SMCGlobal::_ntaxa);
        unsigned k = 0;
        vector<string> taxpartition;
        for (unsigned i = 0; i < SMCGlobal::_nspecies; ++i) {
            string species_name = inventName(i, false);
            SMCGlobal::_species_names[i] = species_name;
            for (unsigned j = 0; j < _nsimtaxaperspecies[i]; ++j) {
                taxpartition.push_back(SMCGlobal::_species_names[i]);
                string taxon_name = str(format("%s^%s") % inventName(k, true) % SMCGlobal::_species_names[i]);
                SMCGlobal::_taxon_names[k] = taxon_name;
                SMCGlobal::_taxon_to_species[taxon_name] = i;
                ++k;
            }
        }

        // Report species names created
        output("  Species names:\n", 2);
        for (unsigned i = 0; i < SMCGlobal::_nspecies; ++i) {
            output(format("    %s\n") % SMCGlobal::_species_names[i], 2);
        }

        // Report taxon names created
        cout << "  Taxon names:\n";
        for (unsigned i = 0; i < SMCGlobal::_ntaxa; ++i) {
            string taxon_name = SMCGlobal::_taxon_names[i];
            output(format("    %s\n") % taxon_name, 2);
        }
        
        // Create data object
        assert(!_data);
        _data = Data::SharedPtr(new Data());
        _data->setTaxonNames(SMCGlobal::_taxon_names);
        _data->setPartition(_partition);
        
        // Build species tree and gene trees jointly
        Particle particle;
        particle.setData(_data);
        simulateTrees(particle);
        
        SpeciesForest      & species_forest = particle.getSpeciesForest();
        vector<GeneForest> & gene_forests   = particle.getGeneForests();
        outputTrees(species_forest, gene_forests);

        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(SMCGlobal::_ngenes);

        // Simulate sequence data
        unsigned starting_site = 0;
        for (unsigned g = 0; g < SMCGlobal::_ngenes; ++g) {
            ps.setNElements(4*SMCGlobal::_nsites_per_gene[g], g);
            gene_forests[g].simulateData(_data, starting_site, SMCGlobal::_nsites_per_gene[g]);
            starting_site += SMCGlobal::_nsites_per_gene[g];
        }
                       
        // Output data to file
        _data->compressPatterns();
        _data->writeDataToFile(_data_file_name);
        cout << "  Sequence data saved in file \"" << _data_file_name << "\"\n";
                 
        // Output a PAUP* command file for estimating the species tree using
        // svd quartets and qage
        output("  PAUP* commands saved in file \"svd-qage.nex\"\n", 1);
        ofstream paupf("svd-qage.nex");
        paupf << "#NEXUS\n\n";
        paupf << "begin paup;\n";
        paupf << "  log start file=svdout.txt replace;\n";
        paupf << "  exe simulated.nex;\n";
        paupf << "  taxpartition species (vector) = " << join(taxpartition," ") << ";\n";
        paupf << "  svd taxpartition=species;\n";
        paupf << "  roottrees;\n";
        paupf << "  qage taxpartition=species patprob=exactjc outUnits=substitutions treefile=svd.tre replace;\n";
        paupf << "  log stop;\n";
        paupf << "  quit;\n";
        paupf << "end;\n";
        paupf.close();
     }
     
#if 0
//    inline double Proj::advanceParticle(unsigned step, unsigned pindex, Particle & particle, bool compute_partials) {
//        double log_weight = 0.0;
//
//        // Get pointer to the pseudorandom number generator to use for this update
//        Lot::SharedPtr lot = particle.getLot();
//
//        SpeciesForest & species_forest = particle.getSpeciesForest();
//        vector<GeneForest> & gene_forests = particle.getGeneForests();
//
//        // Create species_tuples vector. Each 3-tuple entry stores:
//        //  1. number of lineages
//        //  2. gene index
//        //  3. species within gene (ignored if species tree)
//        vector<SMCGlobal::species_tuple_t> species_tuples;
//
//        // Draw a speciation increment Delta. Note: speciation_increment will
//        // equal "infinity" if species tree is complete.
//        auto incr_rate = species_forest.drawIncrement(lot);
//        double Delta = incr_rate.first;
//        double speciation_rate = incr_rate.second;
//
//        // Visit each species within each locus, storing the number of lineages in
//        // species_tuples and computing the total coalescence rate (total_rate).
//        // Next increment will be drawn from an Exponential distribution with rate
//        // total_rate.
//        double total_rate = 0.0;
//        for (auto & gene_forest : gene_forests) {
//            total_rate += gene_forest.calcTotalRate(species_tuples, Delta);
//        }
//
//        // Draw a coalescence increment delta
//        double delta = -log(1.0 - lot->uniform())/total_rate;
//
//        // If delta > Delta, then a speciation event happened; otherwise a coalescent event happened.
//        bool is_speciation = (delta > Delta);
//
//        // Advance all forests
//        double dt = is_speciation ? Delta : delta;
//
//        unsigned num_species_lineages = 1;
//        if (speciation_rate > 0.0) {
//            num_species_lineages = species_forest.advanceAllLineagesBy(dt);
//            if (fabs(SMCGlobal::_lambda*num_species_lineages - speciation_rate) > 0.0001) {
//                throw XProj(format("%d | %.5f | %.5f | %.5f | %.5f > 0.0001") % num_species_lineages % SMCGlobal::_lambda % (SMCGlobal::_lambda*num_species_lineages) % speciation_rate % fabs(SMCGlobal::_lambda*num_species_lineages - speciation_rate));
//            }
//        }
//        for (auto & gene_forest : gene_forests) {
//            gene_forest.advanceAllLineagesBy(dt);
//        }
//
//        if (is_speciation) {
//            particle.setLastEvent(Particle::LAST_EVENT_SPECIATION);
//
//            // Create speciation event
//            SMCGlobal::species_t left;
//            SMCGlobal::species_t right;
//            SMCGlobal::species_t anc;
//            species_forest.speciationEvent(lot, left, right, anc);
//
//            // Advise all gene trees of the change in the species tree
//            for (auto & gene_forest : gene_forests) {
//                gene_forest.mergeSpecies(left, right, anc);
//            }
//        }
//        else {
//            particle.setLastEvent(Particle::LAST_EVENT_COALESCENCE);
//
//            assert(SMCGlobal::_prior_prior);    // prior-post not yet ready
//            if (SMCGlobal::_prior_prior) {
//                // Choose which locus and species in which to have a coalescence
//
//                // Create vector of lineage counts
//                vector<double> lineage_counts(species_tuples.size());
//                transform(species_tuples.begin(), species_tuples.end(), lineage_counts.begin(), [](SMCGlobal::species_tuple_t & spp_tuple){
//                    return get<0>(spp_tuple);
//                });
//
//                // Determine sum of lineage_counts
//                double total_num_lineages = accumulate(lineage_counts.begin(), lineage_counts.end(), 0.0);
//
//                // Normalize lineage_counts to create a discrete probability distribution
//                vector<double> probs(lineage_counts.size());
//                transform(lineage_counts.begin(), lineage_counts.end(), probs.begin(), [total_num_lineages](unsigned nlineages){return (double)nlineages/total_num_lineages;});
//                assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
//
//                // Choose gene and species
//                unsigned which = SMCGlobal::multinomialDraw(lot, probs);
//                //unsigned n                = get<0>(species_tuples[which]);
//                unsigned g                = get<1>(species_tuples[which]);
//                SMCGlobal::species_t spp  = get<2>(species_tuples[which]);
//
//                GeneForest & gene_forest = gene_forests[g];
//
//                log_weight = gene_forest.coalescentEvent(lot, spp, compute_partials);
//
//                // Adjust weight for fact that proposal differs from prior
//                if (num_species_lineages > 1) {
//                    log_weight -= speciation_rate*(delta - Delta);
//                    log_weight -= log(speciation_rate);
//                }
//                assert(!isnan(log_weight));
//                assert(!isinf(log_weight));
//            }
//            else {
//#if 0   // old prior post code
////                // Store log weight for all possible coalescence events in all loci
////                // at the current height
////                vector<GeneForest::coal_tuple_t> coals;
////                for (auto & gene_forest : gene_forests) {
////                    gene_forest.possibleCoalescentEvents(coals);
////                }
////
////                // Compute log(sum of weights)
////                vector<double> log_weights(coals.size());
////                transform(coals.begin(), coals.end(), log_weights.begin(), [](GeneForest::coal_tuple_t & coal_tuple){
////                    return get<0>(coal_tuple);
////                });
////
////                // Compute sum of weights on log scale
////                double log_sum_weights = SMCGlobal::calcLogSum(log_weights);
////
////                // Normalize log weights to create a discrete probability distribution
////                vector<double> probs(log_weights.size());
////                transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
////                assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
////
////                // Choose one pair to join
////                unsigned which = SMCGlobal::multinomialDraw(probs);
////                unsigned g                = get<1>(coals[which]);
////                SMCGlobal::species_t spp  = get<2>(coals[which]);
////                unsigned i                = get<3>(coals[which]);
////                unsigned j                = get<4>(coals[which]);
////
////                // Carry out the chosen join
////                GeneForest & gf = gene_forests[g];
////                gf.coalescePair(spp, i, j);
////
////                // Log sum of weights is the log weight
////                log_weight = log_sum_weights;
//#endif
//            }
//
//            // Only speciation events should have log_weight 0.0, unless we're
//            // simulating, in which case compute_partials will be false
//            assert(log_weight != 0.0 || !compute_partials);
//
//            // If any log_weight is not a number, now's the time to find out about it
//            assert(!isnan(log_weight));
//        }
//
//        if (is_speciation) {
//            particle.incrementSpeciations();
//        }
//        return log_weight;
//    }
#endif

    inline void Proj::updateCountMap(map<string, unsigned> & count_map, string & newick, unsigned count) {
        // If newick is not already in count_map, insert it and set value to count.
        // If it does exist, increase its current value by count.
        // (see item 24, p. 110, in Meyers' Efficient STL for more info on the technique used here)
        auto lowb = count_map.lower_bound(newick);
        if (lowb != count_map.end() && !(count_map.key_comp()(newick, lowb->first))) {
            // this pattern has already been seen
            lowb->second += count;
        }
        else
            {
            // this pattern has not yet been seen
            count_map.insert(lowb, map<string, unsigned>::value_type(newick, count));
        }
    }

    inline double Proj::filterParticles(unsigned step, vector<Particle> & particles, vector<double> & log_weights, vector<unsigned> & counts, map<string, unsigned> & count_map) {
        // Sanity checks
        assert(particles.size() == _nparticles);
        assert(log_weights.size() == _nparticles);
                
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = SMCGlobal::calcLogSum(log_weights);
        vector<double> probs(_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood
        //log_marg_like += log_sum_weights - log(nparticles);
        
        // Compute effective sample size
        double sum_squared_weights = 0.0;
        for (auto it = probs.begin(); it != probs.end(); it++) {
            double w = *it;
            sum_squared_weights += w*w;
        }
        double ess = 1.0/sum_squared_weights;
        
        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());
        
        // Zero vector of counts storing number of darts hitting each particle
        counts.assign(_nparticles, 0);
        
        // Clear count_map from previous filtering
        count_map.clear();
        
        // Throw _nparticles darts
        for (unsigned i = 0; i < _nparticles; ++i) {
            double u = rng.uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            counts[which]++;
        }
        
        // Build count_map
        unsigned count_check = 0;
        for (unsigned i = 0; i < _nparticles; ++i) {
            if (counts[i] > 0) {
                unsigned c = counts[i];
                count_check += c;
                string newick = particles[i].getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false);
                updateCountMap(count_map, newick, c);
            }
        }
        assert(count_check == _nparticles);
        
        // The vector counts represents the results of multinomial sampling, and
        // will need to be intact after this function is called, so make a copy
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

    inline void Proj::outputAnnotatedNexusTreefile(string fn, const vector<tuple<unsigned, double, string, string, string> > & treeinfo) const {
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto tinfo : treeinfo) {
            //unsigned  count = get<0>(tinfo);
            //double      pct = get<1>(tinfo);
            string     note = get<2>(tinfo);
            string treename = get<3>(tinfo);
            string   newick = get<4>(tinfo);
            streef << str(format("  tree %s = [%s] [&R] %s;\n") % treename % note % newick);
            ++t;
        }
        streef << "end;\n";
        streef.close();
    }

    inline void Proj::saveAllGeneTrees(unsigned gene_index, string fn, const vector<Particle> particles) {
    
        // First pass builds map of newick (key) and count (value)
        map<string, unsigned> count_map;
        for (auto & p : particles) {
            assert(gene_index < p.getGeneForests().size());
            const GeneForest & gf = p.getGeneForests()[gene_index];
            string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false);
            updateCountMap(count_map, newick, 1);
        }
        
        // Now create treeinfo tuples
        vector<tuple<unsigned, double, string, string, string> > treeinfo;
        //unsigned i = 0;
        for (auto & kv : count_map) {
            double pct = 100.0*kv.second/_nparticles;
            string note = str(format("freq = %d") % kv.second);
            string treename = str(format("tree%d") % kv.second);
            treeinfo.push_back(make_tuple(kv.second, pct, note, treename, kv.first));

            // Compute log-likelihood and store tree description of
            // gene tree from locus gene_index in this particle
            //assert(gene_index < p.getGeneForests().size());
            //const GeneForest & gf = p.getGeneForests()[gene_index];
            //double log_likelihood = gf.calcLogLikelihood();
            //string note = str(format("lnL = %.5f") % log_likelihood);
            //string treename = str(format("tree%d") % i);
            //string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false);
            //treeinfo.push_back(make_tuple(c, pct, note, treename, newick));
            //++i;
        }
        outputAnnotatedNexusTreefile(fn, treeinfo);
    }
    
    inline void Proj::saveAllSpeciesTrees(string fn, const vector<Particle> particles, const map<string, unsigned> & count_map) {
        vector<tuple<unsigned, double, string, string, string> > treeinfo;
        unsigned i = 0;
        for (auto & kv : count_map) {
            unsigned c = kv.second;
            double pct = 100.0*c/_nparticles;
            string note = str(format("freq = %d") % c);
            string treename = str(format("'tree%d-freq%d'") % i % c);
            string newick = kv.first;
            treeinfo.push_back(make_tuple(c, pct, note, treename, newick));
            ++i;
        }
        outputAnnotatedNexusTreefile(fn, treeinfo);
    }
    
    inline void Proj::saveUniqueSpeciesTrees(string fn, const vector<Particle> particles, const vector<unsigned> & counts) {
        vector<tuple<unsigned, double, string, string, string> > treeinfo;
        unsigned p = 0;
        unsigned i = 0;
        for (auto c : counts) {
            if (c > 0) {
                double pct = 100.0*c/_nparticles;
                string note = str(format("This tree found in %d particles (%.1f%% of %d total particles)") % c % pct % _nparticles);
                string treename = str(format("tree%d-%d") % i % c);
                string newick = particles[p].getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false);
                treeinfo.push_back(make_tuple(c, pct, note, treename, newick));
                ++i;
            }
            ++p;
        }
        sort(treeinfo.begin(), treeinfo.end());
        reverse(treeinfo.begin(), treeinfo.end());
        outputAnnotatedNexusTreefile(fn, treeinfo);
    }
    
    inline unsigned Proj::countSpeciations(const vector<Particle> & particles) {
        unsigned nspeciations = 0;
        for_each(particles.begin(), particles.end(), [&nspeciations](const Particle & p){
            nspeciations += (p.getSpeciations() ? 1 : 0);
        });
        return nspeciations;
    }
    
    inline void Proj::debugCheckCountMap(const map<string, unsigned> & count_map, unsigned target) {
        unsigned total_count = 0;
        for (auto & kv : count_map) {
            total_count += kv.second;
        }
        assert(total_count == target);
    }
    
#if defined(USING_MULTITHREADING)
    inline void Proj::advanceParticleRange(unsigned step, unsigned first, unsigned last) {
        for (unsigned k = first; k < last; ++k) {
            Particle & p = _particles[k];
            p.advance(step, k, /*calculate_partial*/true);
            while (p.lastEventSpeciation()) {
                p.advance(step, k, /*calculate_partial*/true);
            }
        }
    }
#endif
    
#if defined(USING_MPI)
    inline void Proj::run() {
    
        output("Starting...\n", 2);
        output(format("Current working directory: %s\n") % current_path(), 2);
        
        try {
            rng.setSeed(_rnseed);
            
            if (_start_mode == "simulate") {
                throw XProj("Cannot simulate with MPI version");
            }
            else {
                if (_start_mode != "smc")
                    throw XProj("startmode must be either \"simulate\" or \"smc\"");
                    
                readData();
                debugShowStringVector("Gene names", SMCGlobal::_gene_names);
                
                SMCGlobal::_ngenes = _data->getNumSubsets();
                assert(SMCGlobal::_ngenes > 0);

                // Set relative rates if specified
                setRelativeRates();

                SMCGlobal::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(SMCGlobal::_taxon_names);
                
                SMCGlobal::_nspecies = buildSpeciesMap();
                Node::setSpeciesMask(SMCGlobal::_species_mask, SMCGlobal::_nspecies);
                debugShowStringVector("Species names", SMCGlobal::_species_names);

                // Compute leaf partials
                for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                    GeneForest::computeLeafPartials(g, _data);
                }
                
                // Create counts vector used in filtering
                vector<unsigned> counts(_nparticles);
                
                // Create count map that keeps track of unique trees during filtering
                // Key is number of times the particle appeared in the sample
                // Value is the newick description of the species tree.
                map<string, unsigned> count_map;
                
                // Create log_weights vector that stores log weight of each particle
                vector<double> log_weights(_nparticles);
                
                // Create particles vector
                vector<Particle> particles(_nparticles);
                
                // Initialize particles
                for (auto & p : particles) {
                    p.setData(_data);
                    p.resetSpeciesForest();
                    p.resetGeneForests(/*compute_partials*/true);
                }
                
                // Determine total number of steps required to build all
                // gene trees and the species tree
                unsigned nsteps = SMCGlobal::_ngenes*(SMCGlobal::_ntaxa - 1);
                //nsteps += (SMCGlobal::_nspecies - 1);
                
                output(format("\nSMC will require %d steps.\n") % nsteps, 2);

                output(format("\n%12s %12d %12d %12s %12s %12s %12s %12s %12s\n") % "Step" % "before" % "after" % "ESS" % "minlogwt" % "maxlogwt" % "partials" % "secs" % "wait", 2);
                
                double cum_secs = 0.0;
                for (unsigned step = 0; step < nsteps; ++step) {
                    stopwatch.start();
                    
                    log_weights.assign(_nparticles, 0.0);
                    
                    // Clear the number of speciations stored by each particle
                    for_each(particles.begin(), particles.end(), [](Particle & p){p.clearSpeciations();});
                    
                    // Advance each particle by one coalescent event
                    for (unsigned k = ::my_first_particle; k < ::my_last_particle; ++k) {
                        Particle & p = particles[k];
                        double logw = advanceParticle(step, k, p, /*compute_partials*/true);
                        log_weights[k] = logw;
                        while (p.lastEventSpeciation()) {
                            logw = advanceParticle(step, k, p, /*compute_partials*/true);
                            log_weights[k] = logw;
                        }
                    }
                    
                    double minlogw = *min_element(log_weights.begin(), log_weights.end());
                    double maxlogw = *max_element(log_weights.begin(), log_weights.end());

                    // Filter particles using normalized weights and multinomial sampling
                    double ess = _nparticles;

                    unsigned nspeciations_before = countSpeciations(particles);
                    ess = filterParticles(step, particles, log_weights, counts, count_map);
                    unsigned nspeciations_after = countSpeciations(particles);
                    
                    debugCheckCountMap(count_map, _nparticles);
                    double secs = stopwatch.stop();
                    
                    // Calculate waiting time until done
                    cum_secs += secs;
                    double avg_per_step = cum_secs/(step + 1);
                    unsigned steps_to_go = nsteps - (step + 1);
                    double wait = avg_per_step*steps_to_go;
                    
                    unsigned npartials = 0;
#if defined(LOG_MEMORY)
                    npartials = ps.getNumberConstructed();
#endif
                    
                    output(format("%12d %12d %12d %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n") % (step+1) % nspeciations_before % nspeciations_after % ess % minlogw % maxlogw % npartials % secs % wait, 2);
                }
                
                output("\nSpecies trees saved to file \"final-species-trees.tre\"\n", 1);
                saveAllSpeciesTrees("final-species-trees.tre", particles, count_map);
                
                for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                    output(format("Gene trees for locus %d saved to file \"final-gene%d-trees.tre\"\n") % (g+1) % (g+1), 1);
                    saveAllGeneTrees(g, str(format("final-gene%d-trees.tre") % (g+1)), particles);
                }
            }
        }
        catch (XProj & x) {
            output(format("Proj encountered a problem:\n  %s\n") % x.what(), 2);
        }
        
        output("\nFinished!\n", 2);
    }
#else // *not* USING_MPI
    inline void Proj::run() {
    
        output("Starting...\n", 2);
        output(format("Current working directory: %s\n") % current_path(), 2);
        
        try {
            rng.setSeed(_rnseed);
            
            if (_start_mode == "simulate") {
                simulateData();
            }
            else {
                if (_start_mode != "smc")
                    throw XProj("startmode must be either \"simulate\" or \"smc\"");

#if defined(USING_MULTITHREADING)
                threadSetSchedule();
#endif
                    
                readData();
                debugShowStringVector("Gene names", SMCGlobal::_gene_names);
                
                SMCGlobal::_ngenes = _data->getNumSubsets();
                assert(SMCGlobal::_ngenes > 0);

                // Set relative rates if specified
                setRelativeRates();

                SMCGlobal::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(SMCGlobal::_taxon_names);
                
                SMCGlobal::_nspecies = buildSpeciesMap();
                Node::setSpeciesMask(SMCGlobal::_species_mask, SMCGlobal::_nspecies);
                debugShowStringVector("Species names", SMCGlobal::_species_names);

                // Compute leaf partials
                for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                    GeneForest::computeLeafPartials(g, _data);
                }
                
                // Create counts vector used in filtering
                vector<unsigned> counts(_nparticles);
                
                // Create count map that keeps track of unique trees during filtering
                // Key is number of times the particle appeared in the sample
                // Value is the newick description of the species tree.
                map<string, unsigned> count_map;
                
                // Create log_weights vector that stores log weight of each particle
                _log_weights.resize(_nparticles);
                
                // Create particles vector
                _particles.resize(_nparticles);
                
                // Initialize particles
                for (auto & p : _particles) {
                    p.setData(_data);
                    p.resetSpeciesForest();
                    p.resetGeneForests(/*compute_partials*/true);
                }
                
                // Determine total number of steps required to build all
                // gene trees and the species tree
                unsigned nsteps = SMCGlobal::_ngenes*(SMCGlobal::_ntaxa - 1);
                //nsteps += (SMCGlobal::_nspecies - 1);
                
                output(format("\nSMC will require %d steps.\n") % nsteps, 2);

                output(format("\n%12s %12d %12d %12s %12s %12s %12s %12s %12s\n") % "Step" % "before" % "after" % "ESS" % "minlogwt" % "maxlogwt" % "partials" % "secs" % "wait", 2);
                
                double cum_secs = 0.0;
                for (unsigned step = 0; step < nsteps; ++step) {
                    stopwatch.start();
                    
                    _log_weights.assign(_nparticles, 0.0);
                    
                    // Clear the number of speciations stored by each particle
                    for_each(_particles.begin(), _particles.end(), [](Particle & p){p.clearSpeciations();});
                    
                    // Assign each particle a random number seed to use for its update
                    // Multi- and single-threading versions should thus give the same
                    // results, and, in multithreaded version, threads all use
                    // different Lot objects so no mutex needed.
                    unsigned psuffix = 1;
                    for (auto & p : _particles) {
                        p.setSeed(rng.randint(1,9999) + psuffix);
                        psuffix += 2;    // pure superstition (I always use odd seeds)
                    }
                    
                    // Advance each particle by one coalescent event
#if defined(USING_MULTITHREADING)
                    // Advance each particle by one coalescent event
                    vector<thread> threads;
                    for (unsigned i = 0; i < SMCGlobal::_nthreads; i++) {
                        threads.push_back(thread(&Proj::advanceParticleRange,
                            this,
                            step,
                            SMCGlobal::_thread_first_particle[i],
                            SMCGlobal::_thread_last_particle[i]));
                    }
                    
                    // The join function causes this loop to pause until the ith thread finishes
                    for (unsigned i = 0; i < threads.size(); i++) {
                        threads[i].join();
                    }

                    // Copy log weights
                    _log_weights.resize(_nparticles);
                    unsigned pindex = 0;
                    for (auto & p : _particles) {
                        _log_weights[pindex] = p.getLogWeight();
                        ++pindex;
                    }
#else
                    unsigned i = 0;
                    for (auto & p : _particles) {
                        p.advance(step, i, partial);
                        while (p.lastEventSpeciation()) {
                            p.advance(step, i, partial);
                        }
                        i++;
                    }
#endif
                    
                    double minlogw = *min_element(_log_weights.begin(), _log_weights.end());
                    double maxlogw = *max_element(_log_weights.begin(), _log_weights.end());

                    // Filter particles using normalized weights and multinomial sampling
                    double ess = _nparticles;

                    unsigned nspeciations_before = countSpeciations(_particles);
                    ess = filterParticles(step, _particles, _log_weights, counts, count_map);
                    unsigned nspeciations_after = countSpeciations(_particles);
                    
                    debugCheckCountMap(count_map, _nparticles);
                    double secs = stopwatch.stop();
                    
                    // Calculate waiting time until done
                    cum_secs += secs;
                    double avg_per_step = cum_secs/(step + 1);
                    unsigned steps_to_go = nsteps - (step + 1);
                    double wait = avg_per_step*steps_to_go;
                    
                    unsigned npartials = 0;
#if defined(LOG_MEMORY)
                    npartials = ps.getNumberConstructed();
#endif
                    
                    output(format("%12d %12d %12d %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n") % (step+1) % nspeciations_before % nspeciations_after % ess % minlogw % maxlogw % npartials % secs % wait, 2);
                    
                    //temporary!
                    if (nspeciations_after == _nparticles) {
                        output(format("%s\n") % _particles[0].getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false),3);
                    }
                }
                
                output("\nSpecies trees saved to file \"final-species-trees.tre\"\n", 1);
                saveAllSpeciesTrees("final-species-trees.tre", _particles, count_map);
                
#if 0
                for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                    output(format("Gene trees for locus %d saved to file \"final-gene%d-trees.tre\"\n") % (g+1) % (g+1), 1);
                    saveAllGeneTrees(g, str(format("final-gene%d-trees.tre") % (g+1)), _particles);
                }
#endif
            }
        }
        catch (XProj & x) {
            output(format("Proj encountered a problem:\n  %s\n") % x.what(), 2);
        }
        
        output("\nFinished!\n", 2);
    }
#endif // USING_MPI
    
    void Proj::memoryReport(ofstream & memf) const {
        memf << "\nProj memory report:\n\n";
        memf << str(format("  Size of int:           %d\n") % sizeof(int));
        memf << str(format("  Size of char:          %d\n") % sizeof(char));
        memf << str(format("  Size of double:        %d\n") % sizeof(double));
        memf << str(format("  Size of unsigned:      %d\n") % sizeof(unsigned));
        memf << str(format("  Size of unsigned long: %d\n") % sizeof(unsigned long));
        memf << str(format("  Size of Node *:        %d\n") % sizeof(Node *));
        memf << str(format("  Number of particles: %d\n") % _nparticles);
        _data->memoryReport(memf);
    }
    
#if defined(USING_MULTITHREADING)
    inline void Proj::threadSetSchedule() {
        // Determine which particles will be handled by which thread: e.g.
        // 1000 = number of particles
        //  3   = number of threads
        //  333 = 1000 / 3
        //  1   = 1000 % 3
        //  thread 0 gets 333, thread 1 gets 334, thread 2 gets 333
        vector<unsigned> particles_per_thread(SMCGlobal::_nthreads, (unsigned)(_nparticles / SMCGlobal::_nthreads));

        unsigned remainder = _nparticles % SMCGlobal::_nthreads;

        // Each thread > 0 gets an extra job if there is any remainder
        for (unsigned thread = 1; thread < SMCGlobal::_nthreads; ++thread) {
            if (remainder > 0) {
                particles_per_thread[thread] += 1;
                --remainder;
            }
        }

        SMCGlobal::_thread_first_particle.resize(SMCGlobal::_nthreads);
        SMCGlobal::_thread_last_particle.resize(SMCGlobal::_nthreads);
        unsigned cum = 0;
        output("\nParticle schedule:\n", 1);
        output(format("%12s %25s %25s\n") % "thread" % "first particle" % "last particle", 1);
        for (unsigned thread = 0; thread < SMCGlobal::_nthreads; ++thread) {
            SMCGlobal::_thread_first_particle[thread] = cum;
            SMCGlobal::_thread_last_particle[thread]  = cum + particles_per_thread[thread];
            cum += particles_per_thread[thread];
            
            output(format("%12d %25d %25d\n") % thread % (SMCGlobal::_thread_first_particle[thread] + 1) % SMCGlobal::_thread_last_particle[thread], 1);
        }
    }
#endif

#if defined(USING_MPI)
    inline void Proj::mpiSetSchedule() {
        // Determine which particles will be handled by this processor: e.g.
        // 1000 = number of particles
        //  3   = number of processors
        //  333 = 1000 / 3
        //  1   = 1000 % 3
        //  rank 0 gets 333, rank 1 gets 334, rank 2 gets 333
        vector<unsigned> particles_per_task(ntasks, (unsigned)(_nparticles / ntasks));

        unsigned remainder = _nparticles % ntasks;

        // Each rank > 0 gets an extra job if there is any remainder
        for (unsigned rank = 1; rank < ntasks; ++rank) {
            if (remainder > 0) {
                particles_per_task[rank] += 1;
                --remainder;
            }
        }

        _mpi_first_particle.resize(ntasks);
        _mpi_last_particle.resize(ntasks);
        unsigned cum = 0;
        output("\nParticle schedule:\n", 1);
        output(format("%12s %25s %25s\n") % "rank" % "first particle" % "last particle", 1);
        for (unsigned rank = 0; rank < ntasks; ++rank) {
            _mpi_first_particle[rank] = cum;
            _mpi_last_particle[rank]  = cum + particles_per_task[rank];
            cum += particles_per_task[rank];
            
            output(format("%12d %25d %25d\n") % rank % (_mpi_first_particle[rank] + 1) % _mpi_last_particle[rank], 1);
        }
        ::my_first_particle = _mpi_first_particle[::my_rank];
        ::my_last_particle  = _mpi_last_particle[::my_rank];

        // //temporary!
        // cerr << "~~> Proj::mpiSetSchedule\n";
        // cerr << str(format("~~> my_rank           = %d\n") % ::my_rank );
        // cerr << str(format("~~> my_first_particle = %d\n") % ::my_first_particle );
        // cerr << str(format("~~> my_last_particle  = %d\n") % ::my_last_particle );
    }
#endif

}
