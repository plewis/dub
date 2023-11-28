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

            void                         clear();
            void                         processCommandLineOptions(int argc, const char * argv[]);

#if defined(USING_MPI)
            void                         particleLoopMPI(unsigned step, const vector<unsigned> & update_seeds);
#elif defined(USING_MULTITHREADING)
            void                         particleLoopMT(unsigned step, const vector<pair<unsigned, unsigned> > & thread_schedule, const vector<unsigned> & update_seeds);
#else
            void                         particleLoopStd(unsigned step, const vector<unsigned> & update_seeds);
#endif
            void                         run();
            
            void                         generateUpdateSeeds(vector<unsigned> & seeds) const;
            void                         memoryReport(ofstream & memf) const;
            
        private:
                    
            void                         updateCountMap(map<string, unsigned> & count_map, string & newick,    unsigned count);
            //unsigned                    countSpeciations(const vector<Particle> & particles);
            //void debugCheckCountMap(const map<string, unsigned> & count_map, unsigned target);
            void                         parseRelRateDefinition(string & s);
            
            void                         debugShowSpecies(double speciation_incr, vector<Node::species_tuple_t> & species_tuples, unsigned which_gene, SMCGlobal::species_t which_spp);
            void                         debugShowStringVector(string title, const vector<string> & svect) const;
            void                         debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const;

            static string                inventName(unsigned k, bool lower_case);
            void                         outputNexusTreefile(string fn, const vector<string> & newicks) const;
            
            void                         simulateData();
            void                         simulateTrees(Particle & particle);
            void                         outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect);
            void                         outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric);
            void                         outputAnnotatedNexusTreefile(string fn, const vector<tuple<unsigned, double, string, string, string> > & treeinfo) const;
            void                         saveAllSpeciesTrees(string fn, const list<Particle> particle_list);
            void                         saveAllGeneTrees(unsigned gene, string fn, const vector<Particle> particles);
            void                         saveUniqueSpeciesTrees(string fn, const vector<Particle> particles, const vector<unsigned> & counts);

#if defined(USING_MULTITHREADING)
            void                         debugCheckThreadSchedule(const vector<pair<unsigned, unsigned> > & thread_schedule) const;
            void                         debugShowThreadSchedule(const vector<pair<unsigned, unsigned> > & thread_schedule, double percent_of_max_entropy) const;
            void                         divideLargestParticle();
            double                       buildThreadSchedule(vector<pair<unsigned, unsigned> > & thread_schedule);
            void                         balanceThreads(vector<pair<unsigned, unsigned> > & thread_schedule);
            void                         advanceParticleRange(unsigned step, unsigned first_particle, unsigned last_particle, const vector<unsigned> & update_seeds);
#endif
            double                       computeEffectiveSampleSize(const vector<double> & probs) const;
            void                         pruneParticles(list<Particle> & particle_list);
            void                         findNonZeroCountsInRange(vector<unsigned> & nonzeros, const vector<unsigned> & counts, unsigned begin_index, unsigned end_index) const;
            double                       filterParticles(unsigned step, list<Particle> & particle_list, vector<double> & log_weights, vector<unsigned> & counts, vector<unsigned> & rnseeds);

            void                         readData();
            void                         setRelativeRates();
            unsigned                     buildSpeciesMap();
            void                         outputGeneTreesToFile(string fn,
                                            const vector<string> & newicks) const;
            void                         showSettings() const;
            
            string                       _data_file_name;
            string                       _start_mode;
            unsigned                     _niter;
            Partition::SharedPtr         _partition;
            Data::SharedPtr              _data;
            
            double                       _entropy_percent_cutoff;
            bool                         _use_gpu;
            bool                         _ambig_missing;
            unsigned                     _nsimspecies;
            vector<unsigned>             _nsimtaxaperspecies;
            vector<unsigned>             _nsites_per_gene;
            unsigned                     _nparticles;
            int                          _track_split;
            unsigned                     _rnseed;
            bool                         _sort_forests;
            double                       _visualization_cutoff;
            
            list<Particle>               _particle_list;
            vector<double>               _log_weights;
            vector<unsigned>             _proposed_gene;
            vector<SMCGlobal::species_t> _proposed_spp;
            vector<unsigned>             _proposed_first;
            vector<unsigned>             _proposed_second;

            double                       _log_marg_like;
            
            map<string, double>          _relrate_map;
            
            double                       _theta;
            double                       _lambda;
                        
            double                       _theta_delta;
            double                       _lambda_delta;
            
            double                       _ntries_theta;
            double                       _ntries_lambda;
            
#if defined(USING_MPI)
            void mpiSetSchedule();
            vector<unsigned>             _mpi_first_particle;
            vector<unsigned>             _mpi_last_particle;
            vector<unsigned>             _mpi_num_particles;
#endif

            static string                _program_name;
            static unsigned              _major_version;
            static unsigned              _minor_version;
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
        _entropy_percent_cutoff = 99.7;

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
        _log_marg_like = 0.0;
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
        ("entcut", value(&_entropy_percent_cutoff)->default_value(99.7), "when multithreading, determines evenness of distribution of particles to threads (should be 90-100, default 99.7)")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("verbosity",  value(&SMCGlobal::_verbosity)->default_value(0), "0, 1, or 2: higher number means more output")
        ("nspecies",  value(&_nsimspecies)->default_value(1), "number of species (only used if simulate specified)")
        ("ntaxaperspecies",  value(&_nsimtaxaperspecies), "number of taxa sampled per species (only used if simulate specified); should be _nimspecies of these entries, one for each species simulated")
        ("nparticles",  value(&_nparticles)->default_value(1000), "number of particles in a population")
        ("nthreads",  value(&SMCGlobal::_nthreads)->default_value(3), "number of threads (each thread will handle nparticles/nthreads particle updates)")
        ("priorpost", value(&SMCGlobal::_prior_post)->default_value(false), "use prior-post approach to choose coalescence joins (default is prior-prior)")
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
        
        // If user specified --theta on command line...
        if (vm.count("theta") > 0) {
            SMCGlobal::_theta = _theta;
        }
        
        // If user specified --lambda on command line...
        if (vm.count("lambda") > 0) {
            SMCGlobal::_lambda = _lambda;
        }
        
        // If user specified --priorpost on command line, check to ensure
        // that the user is not trying to simulate with priorpost set
        // (that won't work because priorpost requires data).
        if (vm.count("priorpost") > 0) {
            if (SMCGlobal::_prior_post && _start_mode == "simulate") {
                throw XProj("Cannot choose simulate for startmode and yes for prior_post at the same time");
            }
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
                double r = 0.0;
                try {
                    r = _relrate_map.at(gname);
                } catch(const out_of_range &) {
                    throw XProj(format("Proj::setRelativeRates failed because key \"%s\" does not exist in _relrate_map") % gname);
                }
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
            unsigned species_index = 0;
            try {
                species_index = species_name_to_index.at(sname);
            } catch(const out_of_range & ) {
                throw XProj(format("Proj::buildSpeciesMap failed because key \"%s\" does not exist in species_name_to_index map") % sname);
            }
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
    
    inline void Proj::debugShowSpecies(double speciation_incr, vector<Node::species_tuple_t> & species_tuples, unsigned which_gene, SMCGlobal::species_t which_spp) {
        output(format("\ndebugShowSpecies (speciation increment = %.5f):\n") % speciation_incr, 1);

        unsigned ntuples = (unsigned)species_tuples.size();
        output(format("%12s %12s %s\n")  % "nlineages" % "gene" % "species", 1);
        
        for (unsigned i = 0; i < ntuples; i++) {
            Node::species_tuple_t & species_tuple = species_tuples[i];
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
            particle.advance(step, 0, /*calculate_partial*/false, /*make_permanent*/true);
            //TODO: deal with removal of make_permanent
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
        //output("  True gene tree height (height/theta):\n",2);
        for (auto & gf : gfvect) {
            string newick_alpha = gf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
            newick_gene_trees_alpha.push_back(newick_alpha);
            
            string newick_numeric = gf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
            newick_gene_trees_numeric.push_back(newick_numeric);
            //output(format("  %12.9f (%.9f)\n") % gf.getHeight() % (gf.getHeight()/SMCGlobal::_theta), 2);
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
        particle.setSeed(rng.randint(1,9999));
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
            gene_forests[g].simulateData(particle.getLot(), _data, starting_site, SMCGlobal::_nsites_per_gene[g]);
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

    inline double Proj::computeEffectiveSampleSize(const vector<double> & probs) const {
        double ss = 0.0;
        for_each(probs.begin(), probs.end(), [&ss](double w){ss += w*w;});
        double ess = 1.0/ss;
        return ess;
    }

    inline void Proj::pruneParticles(list<Particle> & particle_list) {
        list<Particle>::iterator it = particle_list.begin();
        while (it != particle_list.end()) {
            if (it->getCount() == 0) {
                // Iterator post-increment returns previous iterator,
                // which can now be erased since the it has moved on
                particle_list.erase(it++);
            }
            else {
                // Leave this Particle intact because it has non-zero count
                ++it;
            }
        }
    }

    inline void Proj::findNonZeroCountsInRange(vector<unsigned> & nonzeros, const vector<unsigned> & counts, unsigned begin_index, unsigned end_index) const {
        nonzeros.clear();
        for (unsigned k = begin_index; k < end_index; k++) {
            if (counts[k] > 0) {
                nonzeros.push_back(k);
            }
        }
    }

    inline double Proj::filterParticles(unsigned step, list<Particle> & particle_list, vector<double> & log_weights, vector<unsigned> & counts, vector<unsigned> & rnseeds) {
        //TODO: Proj::filterParticles
        // Sanity checks
        assert(counts.size() == _nparticles);
        assert(log_weights.size() == _nparticles);
                
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = SMCGlobal::calcLogSum(log_weights);
        vector<double> probs(_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood
        _log_marg_like += log_sum_weights - log(_nparticles);
        
        // Compute effective sample size
        double ess = computeEffectiveSampleSize(probs);
        
        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());
        
        // Zero vector of counts storing number of darts hitting each particle
        counts.assign(_nparticles, 0);
                
        // Throw _nparticles darts
        for (unsigned i = 0; i < _nparticles; ++i) {
            double u = rng.uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            counts[which]++;
        }
        
        // Suppose _nparticles = 20 and counts looks like this:
        //                   <----------- 1st particle -----------><-- 2nd particle --->
        //         counts = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 12}
        //                   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18  19
        //                              nonzeros = [3,6,11]            nonzeros = [19]
        //
        // Here, _particle_list originally contains 2 particles:
        //
        //   Particle  Count
        //          0     12
        //          1      8
        //
        // After filtering, _particle_list will contain the following particle counts:
        //
        //                    Copied   Using
        //   Particle  Count    from    seed
        //          0      1       0       3
        //          1      2       0       6
        //          2      5       0      11
        //          3     12       1      19
        unsigned num_parent_particles = (unsigned)particle_list.size();
        unsigned i = 0;
        for (unsigned j = 0; j < num_parent_particles; j++) {
            auto pit = particle_list.begin();
            advance(pit, j);
            Particle & p = *pit;
            
            // Next n elements of counts belong to particle p
            unsigned n = p.getCount();
            
            // Find all non-zero counts associated with particle p
            vector<unsigned> nonzeros;
            findNonZeroCountsInRange(nonzeros, counts, i, i + n);
                        
            for (auto k : nonzeros) {
                // Make a copy of p
                particle_list.push_back(p);
                Particle & plast = *(particle_list.rbegin());
                                    
                // Advance to same coalescence event created when log weight was calculated
                // but this time keep it by setting compute_partial to false. Fills coal_proposal
                // struct with information about the coalescence event as well as any
                // speciation events created beforehand.
                
                plast.setLastProposedGene(_proposed_gene[k]);
                plast.setLastProposedSpecies(_proposed_spp[k]);
                plast.setLastProposedFirstIndex(_proposed_first[k]);
                plast.setLastProposedSecondIndex(_proposed_second[k]);

#if defined(MINIMIZE_PARTIALS)
                plast.proposeCoalescence(rnseeds[k], step, k,  /*compute_partial*/false, /*make_permanent*/true);
#else
                plast.proposeCoalescence(rnseeds[k], step, k,  /*compute_partial*/true, /*make_permanent*/true);
#endif
                                
                // Set count for new particle
                plast.setCount(counts[k]);

                //temporary! Useful for debugging filtering
                //output(format("~~~| Copying %d (%d) -> %d (%d)\n") % j % p.getCount() % (particle_list.size()-1) % plast.getCount(), 2);
            }
                            
            // Flag original particle for deletion
            p.setCount(0);
            
            i += n;
        }
        
        // Eliminate all particles with a count of 0
        pruneParticles(particle_list);

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
    
    inline void Proj::saveAllSpeciesTrees(string fn, const list<Particle> particle_list) {
        map<string,unsigned> tree_map;
        auto it = tree_map.begin();
        for (const Particle & p : particle_list) {
            unsigned c = p.getCount();
            string newick = p.getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
            try {
                unsigned curr = tree_map.at(newick);
                tree_map[newick] += c;
            }
            catch(const out_of_range &) {
                tree_map[newick] = c;
            }
        }
        
        vector<tuple<unsigned, double, string, string, string> > treeinfo;
        unsigned i = 0;
        for (auto it = tree_map.begin(); it != tree_map.end(); ++it) {
            const string & newick = it->first;
            unsigned c = it->second;
            double pct = 100.0*c/_nparticles;
            string note = str(format("freq = %d") % c);
            string treename = str(format("'tree%d-freq%d'") % i % c);
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
    
    //inline unsigned Proj::countSpeciations(const vector<Particle> & particles) {
    //    unsigned nspeciations = 0;
    //    for_each(particles.begin(), particles.end(), [&nspeciations](const Particle & p){
    //        nspeciations += (p.getSpeciations() ? 1 : 0);
    //    });
    //    return nspeciations;
    //}
    
    //inline void Proj::debugCheckCountMap(const map<string, unsigned> & count_map, unsigned target) {
    //    unsigned total_count = 0;
    //    for (auto & kv : count_map) {
    //        total_count += kv.second;
    //    }
    //    assert(total_count == target);
    //}
    
#if defined(USING_MULTITHREADING)
    inline void Proj::advanceParticleRange(unsigned step, unsigned first_particle, unsigned last_particle, const vector<unsigned> & update_seeds) {
        // Uncomment line below to ensure each thread gets through this entire function
        // uninterrupted. This is good for debugging but of course non-sensical for a
        // multithreading application.
        //lock_guard<mutex> guard(SMCGlobal::_debug_mutex);
        
        // Get iterator to start of range
        auto it_start = _particle_list.begin();
        advance(it_start, first_particle);
        
        int first_index = it_start->getBeginIndex();
        assert(first_index >= 0 && first_index < _log_weights.size());

        // Get iterator to (one past) end of range
        auto it_end = _particle_list.begin();
        advance(it_end, last_particle);
        
        // Loop over particles and process each
        unsigned i = first_index;
        for (auto it = it_start; it != it_end; ++it) {
            Particle & p = *it;
            
            // Second element of pair is number of copies of particle
            unsigned n = p.getCount();
        
#if defined(MINIMIZE_PARTIALS)
            // Recompute all partials
            p.computeAllPartials();
#endif
        
            while (n > 0) {
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugStoreForests();
#endif
                
                // Propose a coalescence event (which may involve also proposing
                // one or more intervening speciation events)
                _log_weights[i] = p.proposeCoalescence(update_seeds[i], step, i,  /*compute_partial*/true, /*make_permanent*/false);
                _proposed_gene[i] = p.getLastProposedGene();
                _proposed_spp[i] = p.getLastProposedSpecies();
                _proposed_first[i] = p.getLastProposedFirstIndex();
                _proposed_second[i] = p.getLastProposedSecondIndex();
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugCheckForests();
#endif
                
                n--;
                i++;
            }
            
#if defined(MINIMIZE_PARTIALS)
            p.stowAllPartials();
#endif
        }
    }
#endif
    
    inline void Proj::generateUpdateSeeds(vector<unsigned> & seeds) const {
        unsigned psuffix = 1;
        for (auto & s : seeds) {
            s = rng.randint(1,9999) + psuffix;
            psuffix += 2;    // pure superstition (I always use odd seeds)
        }
    }
    
#if defined(USING_MPI)
    inline void Proj::particleLoopMPI(unsigned step, const vector<unsigned> & update_seeds) {
        // Advance each particle by one coalescent event
        for (unsigned k = ::my_first_particle; k < ::my_last_particle; ++k) {
            Particle & p = _particles[k];
            p.advance(step, k, /*calculate_partial*/true);
            while (p.lastEventSpeciation()) {
                p.advance(step, k, /*calculate_partial*/true);
            }
        }
        
        // Erect barrier causing this loop to pause
        // until the last process finishes
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (my_rank == 0) {
            // The outstanding variable equals the number of
            // processors we haven't heard from yet
            unsigned outstanding = ntasks;
            
            // Receive particle weights from other processors
            while (outstanding) {
                // Probe to get message status
                int message_length = 0;
                MPI_Status status;
                MPI_Probe(MPI_ANY_SOURCE,   // Source rank or MPI_ANY_SOURCE
                    MPI_ANY_TAG,            // Message tag
                    MPI_COMM_WORLD,         // Communicator
                    &status                 // Status object
                );
            
                // Get length of message
                MPI_Get_count(&status, MPI_DOUBLE, &message_length);

                // Get processor sending message
                unsigned which = (unsigned)status.MPI_TAG;

                // Get the message itself
                vector<double> log_weights(_mpi_num_particles[which], 0.0);
                MPI_Recv(&log_weights[0],    // Initial address of receive buffer
                    message_length,     // Maximum number of elements to receive
                    MPI_DOUBLE,           // Datatype of each receive buffer entry
                    status.MPI_SOURCE,     // Rank of source
                    status.MPI_TAG,        // Message tag
                    MPI_COMM_WORLD,     // Communicator
                    MPI_STATUS_IGNORE   // Status object
                );
                
                // Copy log_weights to appropriate section of
                // the _log_weights vector
                unsigned kfrom = 0;
                for (unsigned kto = ::my_first_particle; kto < ::my_last_particle; ++kto) {
                    assert(kfrom < log_weights.size());
                    _log_weights[kto] = log_weights[kfrom];
                }
                
                --outstanding;
            }
        }
    }
#elif defined(USING_MULTITHREADING)
    inline void Proj::particleLoopMT(unsigned step, const vector<pair<unsigned, unsigned> > & thread_schedule, const vector<unsigned> & update_seeds) {
        //TODO: multithreading particle loop

        vector<thread> threads;
        for (unsigned i = 0; i < SMCGlobal::_nthreads; i++) {
            threads.push_back(thread(&Proj::advanceParticleRange,
                this,
                step,
                thread_schedule[i].first,
                thread_schedule[i].second,
                update_seeds)
            );
        }

        // The join function causes this loop to pause until the ith thread finishes
        for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
        }
    }
#else   // not USING_MPI or USING_MULTITHREADING
    inline void Proj::particleLoopStd(unsigned step, const vector<unsigned> & update_seeds) {
        //TODO: standard particle loop
        unsigned i = 0;
        for (auto & p : _particle_list) {
            // Second element of pair is number of copies of particle
            unsigned n = p.getCount();

#if defined(MINIMIZE_PARTIALS)
            // Recompute all partials
            p.computeAllPartials();
#endif
            
            while (n > 0) {
            
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugStoreForests();
#endif
                
                // Propose a coalescence event (which may involve also proposing
                // one or more intervening speciation events)
                _log_weights[i] = p.proposeCoalescence(update_seeds[i], step, i, /*compute_partial*/true, /*make_permanent*/false);
                _proposed_gene[i] = p.getLastProposedGene();
                _proposed_spp[i] = p.getLastProposedSpecies();
                _proposed_first[i] = p.getLastProposedFirstIndex();
                _proposed_second[i] = p.getLastProposedSecondIndex();
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugCheckForests();
#endif

                n--;
                i++;
            }
            
#if defined(MINIMIZE_PARTIALS)
            p.stowAllPartials();
#endif
        }
        assert(i == _nparticles);
    }
#endif

#if defined(USING_MULTITHREADING)
    inline void Proj::debugCheckThreadSchedule(const vector<pair<unsigned, unsigned> > & thread_schedule) const {
    
        //auto it = _particle_list.begin();
        //for (auto & p : thread_schedule) {
        //    unsigned begin_particle  = p.first;
        //    unsigned end_particle = p.second;
        //    unsigned begin_index = it->getBeginIndex();
        //}
        // begin again here
    }
    
    inline void Proj::debugShowThreadSchedule(const vector<pair<unsigned, unsigned> > & thread_schedule, double percent_of_max_entropy) const {
        unsigned i = 0;
        auto it = _particle_list.begin();
        unsigned sum_pcount = 0;
        output(format("\n%12s %12s %12s %12s %12s %12s\n") % "index" % "begin" % "end" % "begin-index" % "end-index" % "count", 2);
        for (auto & p : thread_schedule) {
            unsigned begin_particle  = p.first;
            unsigned end_particle = p.second;
            unsigned begin_index     = it->getBeginIndex();
            unsigned n = end_particle - begin_particle;
            unsigned pcount = 0;
            for (unsigned j = 0; j < n; j++) {
                pcount += it->getCount();
                ++it;
            }
            sum_pcount += pcount;
            unsigned end_index = begin_index + pcount;
            output(format("%12d %12d %12d %12d %12d %12d\n") % i % begin_particle % end_particle % begin_index % end_index % pcount, 2);
            i++;
        }
        assert(it == _particle_list.end());
        output(format("Total actual particles:  %d\n") % _particle_list.size(), 2);
        output(format("Total virtual particles: %d\n") % sum_pcount, 2);
        output(format("Percentage of maximum entropy: %.1f\n") % percent_of_max_entropy, 2);
        cerr << endl;
    }
#endif

#if defined(USING_MULTITHREADING)
    inline void Proj::divideLargestParticle() {
        auto it = max_element(_particle_list.begin(), _particle_list.end(), [](const Particle & first, const Particle & second){return first.getCount() < second.getCount() ? true : false;});

        unsigned begin_index = it->getBeginIndex();
        unsigned count = it->getCount();
        unsigned first_half = count/2;
        unsigned second_half = count - first_half;

        auto it0 = _particle_list.insert(it, *it);
                
        it0->setCount(first_half);
        it0->setXtra(0);
        it0->setBeginIndex(begin_index);
        
        it->setCount(second_half);
        it->setXtra(0);
        it->setBeginIndex(begin_index + first_half);
    }
#endif
    
#if defined(USING_MULTITHREADING)
    inline double Proj::buildThreadSchedule(vector<pair<unsigned, unsigned> > & thread_schedule) {
        // Calculate thread_schedule and its associated entropy and return the percentage
        // entropy relative to the maximum possible entropy

        // Create thread schedule using the prefix-sum algorithm
        thread_schedule.clear();
        unsigned current_thread = 0;
        pair<unsigned, unsigned> begin_end = make_pair(0,0);
        unsigned prefix_sum = 0;
        unsigned max_count = (unsigned)floor(1.0*_nparticles/SMCGlobal::_nthreads);
        
        vector<unsigned> freqs(SMCGlobal::_nthreads, 0);
        for (auto & p : _particle_list) {
            // Add particle count to prefix sum
            unsigned count = p.getCount();
            if (count > max_count) {
                thread_schedule.clear();
                return 0.0;
            }
            
            p.setBeginIndex(prefix_sum);
            prefix_sum += count;
                
            // Calculate thread index
            unsigned thread_index = (unsigned)floor(1.0*SMCGlobal::_nthreads*(prefix_sum - 1)/_nparticles);
            
            if (thread_index > current_thread) {
                // Add thread to the schedule
                thread_schedule.push_back(begin_end);
                current_thread++;
                
                // Start work on next thread
                begin_end.first  = begin_end.second;
                begin_end.second = begin_end.first + 1;
            }
            else {
                begin_end.second += 1;
            }
            freqs[current_thread] += count;
        }
        thread_schedule.push_back(begin_end);

        // Calculate entropy, max_entropy, and percentage
        double max_entropy = log(SMCGlobal::_nthreads);
        double entropy = 0.0;
        assert(_nparticles == accumulate(freqs.begin(), freqs.end(), 0));
        for_each(freqs.begin(), freqs.end(), [&entropy](unsigned f){entropy -= 1.0*f*log(f);});
        entropy /= _nparticles;
        entropy += log(_nparticles);
        double percentage = 100.0*entropy/max_entropy;
        assert(!isnan(percentage));
        
        return percentage;
    }
#endif
    
#if defined(USING_MULTITHREADING)
    inline void Proj::balanceThreads(vector<pair<unsigned, unsigned> > & thread_schedule) {
        // Make copies of particles as necessary to achieve a more even distribution of work.
        // Example: 2 threads, 100 (virtual) particles, 4 (actual) particles
        //
        // Original particles:
        //  index  count  cumulative   prefix-sum   calculation  thread     count
        // ------ ------ ------------ ----------- ------------- ------- ---------
        //      0     11          11           11 2*10/100=0.20       0
        //      1     25          25           36 2*35/100=0.70       0  11+25=36
        //      2     48          48           84 2*83/100=1.66       1
        //      3     16          16          100 2*99/100=1.98       1  48+16=64
        // ------ ------ ------------ ----------- ------------- ------- ---------
        // entropy      = 0.653 = -(0.36*ln(0.36) + 0.64*ln(0.64))
        // max(entropy) = 0.693 = log(2)
        // percentage   = 94.2
        //
        // Continue dividing largest particle until % max entropy > target (e.g. 99.6)
        //   99.971 = -100.0*(0.49*ln(0.49) + 0.51*ln(0.51))/ln(2)
        //   99.885 = -100.0*(0.48*ln(0.48) + 0.52*ln(0.52))/ln(2)
        //   99.740 = -100.0*(0.47*ln(0.47) + 0.53*ln(0.53))/ln(2) *
        //   99.538 = -100.0*(0.46*ln(0.46) + 0.54*ln(0.54))/ln(2)
        //   99.277 = -100.0*(0.45*ln(0.45) + 0.55*ln(0.55))/ln(2)
        assert(SMCGlobal::_nthreads > 0);
        double pct = 0.0;
        if (SMCGlobal::_nthreads == 1) {
            thread_schedule.clear();
            thread_schedule.push_back(make_pair(0,(unsigned)_particle_list.size()));
        }
        else {
            while (_particle_list.size() < SMCGlobal::_nthreads)
                divideLargestParticle();
            pct = buildThreadSchedule(thread_schedule);
            while (pct < _entropy_percent_cutoff) {
                divideLargestParticle();
                pct = buildThreadSchedule(thread_schedule);
            }
        }
        
        //debugShowThreadSchedule(thread_schedule, pct);
        //cerr << endl;
    }
#endif

//    inline void Proj::balanceThreads(vector<pair<unsigned, unsigned> > & thread_schedule) {
//        // Make copies of particles as necessary to achieve a more even distribution of work.
//        //    6 nthreads
//        //  100 total
//        // 16.0 threshold
//        //
//        // Original particles:
//        //           id        count   cumulative
//        // ------------ ------------ ------------
//        //            a           67           67
//        //            b           21           88
//        //            c            7           95
//        //            d            3           98
//        //            e            2          100
//        //
//        // Assigments:
//        //   id   count   extra   total   cumulative   calc   thread
//        // ---- ------- ------- ------- ------------ ------ --------
//        //    b      16       0      16           16   0.90        0
//        //    a      16       0      16           32   1.86        1
//        //    a      16       1      17           49   2.82        2
//        //    a      16       1      17           66   3.78        3
//        //    a      16       1      17           83   4.74        4
//        //    c       7       0       7           90   5.16        5
//        //    b       5       0       5           95   5.46        5
//        //    d       3       0       3           98   5.64        5
//        //    e       2       0       2          100   5.76        5
//        //
//        // Number of particles per thread:
//        //       Thread    Particles
//        // ------------ ------------
//        //            0           16
//        //            1           16
//        //            2           17
//        //            3           17
//        //            4           17
//        //            5           17
//        // ------------ ------------
//        //        total          100
//        //TODO: Proj::balanceThreads
//        assert(SMCGlobal::_nthreads > 0);
//        if (SMCGlobal::_nthreads == 1) {
//            thread_schedule.clear();
//            thread_schedule.push_back(make_pair(0,(unsigned)_particle_list.size()));
//        }
//        else {
//            unsigned total_xtra = 0;
//            unsigned begin_index = 0;
//            unsigned threshold = (unsigned)floor(_nparticles/SMCGlobal::_nthreads);
//
//            for (auto it = _particle_list.begin(); it != _particle_list.end(); ++it) {
//                unsigned count = it->getCount();
//                if (count > threshold) {
//                    unsigned d = count / threshold;
//                    unsigned r = count % threshold;
//                    bool augment = (r <= d ? true : false);
//                    unsigned remaining = count;
//                    while (remaining >= threshold) {
//                        remaining -= threshold;
//
//                        // Calculate xtra
//                        unsigned x = 0;
//                        if (augment && r > 0) {
//                            x = 1;
//                            remaining -= 1;
//                            r--;
//                        }
//                        total_xtra += x;
//
//                        unsigned new_begin_index = begin_index + threshold + x;
//
//                        // Make a copy of the current particle if count remaining >= threshold
//                        if (remaining >= threshold) {
//                            _particle_list.push_front(*it);
//
//                            // Adjust current particle counts
//                            it->setCount(remaining);
//                            it->setXtra(0);
//                            it->setBeginIndex(new_begin_index);
//
//                            // Set copied particle counts
//                            _particle_list.begin()->setCount(threshold);
//                            _particle_list.begin()->setXtra(x);
//                            _particle_list.begin()->setBeginIndex(begin_index);
//                        }
//                        else {
//                            // Adjust current particle counts
//                            it->setCount(count - x);
//                            it->setXtra(x);
//                            it->setBeginIndex(begin_index); // redundant
//                        }
//
//                        begin_index = new_begin_index;
//                        count = remaining;
//                    }
//                }
//                else {
//                    // count not greater than threshold
//                    it->setCount(count);    // redundant
//                    it->setXtra(0);
//                    it->setBeginIndex(begin_index);
//                    begin_index += count;
//                }
//            }
//
//            // Sort particles by increasing begin index
//            _particle_list.sort([](const Particle & first, const Particle & second){
//                return first.getBeginIndex() < second.getBeginIndex() ? true : false;
//            });
//
//            unsigned reduced_total = _nparticles - total_xtra;
//
//            // Create thread schedule
//            thread_schedule.clear();
//            unsigned current_thread = 0;
//            pair<unsigned, unsigned> begin_end = make_pair(0,0);
//            unsigned sum_counts = 0;
//            for (auto & p : _particle_list) {
//                // Move xtra into count on Particle itself
//                unsigned count = p.getCount();
//                unsigned xtra = p.getXtra();
//                p.setCount(count + xtra);
//                p.setXtra(0);
//
//                // For purposes of determining thread index, however, ignore xtra
//                sum_counts += count;
//
//                // Calculate thread index
//                unsigned thread_index = (unsigned)floor(1.0*SMCGlobal::_nthreads*(sum_counts - 1)/reduced_total);
//
//                if (thread_index > current_thread) {
//                    // Adding a new thread to the schedule
//                    thread_schedule.push_back(begin_end);
//                    unsigned b = begin_end.second;
//                    unsigned e = b + 1;
//                    begin_end.first = b;
//                    begin_end.second = e;
//                    current_thread++;
//                }
//                else {
//                    begin_end.second += 1;
//                }
//            }
//            thread_schedule.push_back(begin_end);
//        }
//
//        //debugShowThreadSchedule(thread_schedule);
//    }
                    
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

                readData();
                debugShowStringVector("Gene names", SMCGlobal::_gene_names);
                
                SMCGlobal::_ngenes = _data->getNumSubsets();
                assert(SMCGlobal::_ngenes > 0);

//#if defined(USING_MULTITHREADING)
//                threadSetSchedule();
#if defined(USING_MPI)
                mpiSetSchedule();
#endif
                    
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
                // Stores number of darts that hit each of the _nparticles in
                // multinomial sampling.
                vector<unsigned> counts(_nparticles);
                
                // Create vector of random number seeds to be reused each step.
                vector<unsigned> update_seeds(_nparticles);
                
                // Create count map that keeps track of unique trees during filtering
                // and which is used in saveAllSpeciesTrees function.
                // Key is number of times the particle appeared in the sample
                // Value is the newick description of the species tree.
                //map<string, unsigned> count_map;
                
                // Create log_weights vector that stores log weight of each particle
                _log_weights.resize(_nparticles);
                
                // Initialize first particle with trivial forests
                Particle template_particle;
                template_particle.setCount(_nparticles);
                template_particle.setData(_data);
                template_particle.resetSpeciesForest();
                template_particle.resetGeneForests(/*compute_partials*/true);
                
                // Compute initial log likelihood
                double starting_log_likelihood = template_particle.calcLogLikelihood();
                
#if defined(MINIMIZE_PARTIALS)
                stowAllPartials();
#endif
                
                // Initialize particle map with one particle having count _nparticles
                _particle_list.clear();
                _particle_list.push_back(template_particle);
                                
                // Determine total number of steps required to build all
                // gene trees and the species tree
                unsigned nsteps = SMCGlobal::_ngenes*(SMCGlobal::_ntaxa - 1);
                
                output(format("\nSMC will require %d steps.\n") % nsteps, 2);

#if defined(LOG_MEMORY)
                output(format("\n%12s %12s %12s %12s %24s %12s %12s %12s %12s\n") % "Step" % "ESS" % "minlogwt" % "maxlogwt" % "logml" % "in-use" % "stored" % "secs" % "wait", 2);
#else
                output(format("\n%12s %12s %12s %12s %24s %12s %12s\n") % "Step" % "ESS" % "minlogwt" % "maxlogwt" % "logml" % "secs" % "wait", 2);
#endif
                
                _log_marg_like = starting_log_likelihood;
                double cum_secs = 0.0;
                for (unsigned step = 0; step < nsteps; ++step) {
                    //TODO: main step loop
                    stopwatch.start();

#if defined(USING_MULTITHREADING)
                    vector<pair<unsigned, unsigned> > thread_schedule;
                    balanceThreads(thread_schedule);
#endif

                    _log_weights.assign(_nparticles, 0.0);
                    _proposed_gene.assign(_nparticles, 0);
                    _proposed_spp.assign(_nparticles, 0);
                    _proposed_first.assign(_nparticles, 0);
                    _proposed_second.assign(_nparticles, 0);
                    
                    // Assign _nparticles random number seeds to use for generating the next step
                    generateUpdateSeeds(update_seeds);
                    
                    // Advance each particle by one coalescent event
#if defined(USING_MPI)
                    particleLoopMPI(step, update_seeds);
#elif defined(USING_MULTITHREADING)
                    particleLoopMT(step, thread_schedule, update_seeds);
#else
                    particleLoopStd(step, update_seeds);
#endif
                    
                    double minlogw = *min_element(_log_weights.begin(), _log_weights.end());
                    double maxlogw = *max_element(_log_weights.begin(), _log_weights.end());

                    // Filter particles using normalized weights and multinomial sampling
                    double ess = _nparticles;
                    //unsigned nspeciations_before = countSpeciations(_particles);
                    ess = filterParticles(step, _particle_list, _log_weights, counts, update_seeds);
                    //unsigned nspeciations_after = countSpeciations(_particles);
                    
                    //debugCheckCountMap(count_map, _nparticles);
                    double secs = stopwatch.stop();
                    
                    // Calculate waiting time until done
                    cum_secs += secs;
                    double avg_per_step = cum_secs/(step + 1);
                    unsigned steps_to_go = nsteps - (step + 1);
                    double wait = avg_per_step*steps_to_go;
                    
#if defined(LOG_MEMORY)
                    //npartials = ps.getNumberConstructed();
                    npartials_inuse = ps.getInUse();
                    npartials_stored = ps.getStored();
                    
                    output(format("%12d %12.3f %12.3f %12.3f %24.6f %12d %12d %12.3f %12.3f\n") % (step+1) % ess % minlogw % maxlogw % _log_marg_like % npartials_inuse % npartials_stored % secs % wait, 2);
#else
                    output(format("%12d %12.3f %12.3f %12.3f %24.6f %12.3f %12.3f\n") % (step+1) % ess % minlogw % maxlogw % _log_marg_like % secs % wait, 2);
#endif
                    
                    // Show species tree for iterations in which a species
                    // tree join ended up being copied to all particles
                    //if (nspeciations_after == _nparticles) {
                    //    output(format("%s\n") % _particles[0].getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalescent units*/false),3);
                    //}
                }
                
                output(format("log(marginal likelihood) = %.6f\n") % _log_marg_like, 1);
                output("\nSpecies trees saved to file \"final-species-trees.tre\"\n", 1);
                saveAllSpeciesTrees("final-species-trees.tre", _particle_list);
                
                //for (unsigned g = 0; g < SMCGlobal::_ngenes; g++) {
                //    output(format("Gene trees for locus %d saved to file \"final-gene%d-trees.tre\"\n") % (g+1) % (g+1), 1);
                //    saveAllGeneTrees(g, str(format("final-gene%d-trees.tre") % (g+1)), _particles);
                //}
            }
        }
        catch (XProj & x) {
            output(format("Proj encountered a problem:\n  %s\n") % x.what(), 2);
        }
        
        output("\nFinished!\n", 2);
    }
    
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
    
//#if defined(USING_MULTITHREADING)
//    inline void Proj::threadSetSchedule() {
//        if (SMCGlobal::_ngenes < SMCGlobal::_nthreads) {
//            throw XProj(format("Aborting because the number of genes (%d) is less than the number of threads (%d) requested") % SMCGlobal::_ngenes % SMCGlobal::_nthreads);
//        }
//
//        // Determine which genes will be handled by which thread: e.g.
//        // 1000 = number of genes
//        //  3   = number of threads
//        //  333 = 1000 / 3
//        //  1   = 1000 % 3
//        //  thread 0 gets 333, thread 1 gets 334, thread 2 gets 333
//        vector<unsigned> genes_per_thread(SMCGlobal::_nthreads, (unsigned)(SMCGlobal::_ngenes / SMCGlobal::_nthreads));
//
//        unsigned remainder = SMCGlobal::_ngenes % SMCGlobal::_nthreads;
//
//        // Each thread > 0 gets an extra job if there is any remainder
//        for (unsigned thread = 1; thread < SMCGlobal::_nthreads; ++thread) {
//            if (remainder > 0) {
//                genes_per_thread[thread] += 1;
//                --remainder;
//            }
//        }
//
//        SMCGlobal::_thread_first_gene.resize(SMCGlobal::_nthreads);
//        SMCGlobal::_thread_last_gene.resize(SMCGlobal::_nthreads);
//        unsigned cum = 0;
//        output("\nGene schedule:\n", 1);
//        output(format("%12s %25s %25s\n") % "thread" % "first gene" % "last gene", 1);
//        for (unsigned thread = 0; thread < SMCGlobal::_nthreads; ++thread) {
//            SMCGlobal::_thread_first_gene[thread] = cum;
//            SMCGlobal::_thread_last_gene[thread]  = cum + genes_per_thread[thread];
//            cum += genes_per_thread[thread];
//
//            output(format("%12d %25d %25d\n") % thread % (SMCGlobal::_thread_first_gene[thread] + 1) % SMCGlobal::_thread_last_gene[thread], 1);
//        }
//    }
//#endif

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
        _mpi_num_particles.resize(ntasks);
        unsigned cum = 0;
        output("\nParticle schedule:\n", 1);
        output(format("%12s %25s %25s %25s\n") % "rank" % "first particle" % "last particle" % "no. particles", 1);
        for (unsigned rank = 0; rank < ntasks; ++rank) {
            _mpi_first_particle[rank] = cum;
            _mpi_last_particle[rank]  = cum + particles_per_task[rank];
            _mpi_num_particles[rank]  = particles_per_task[rank];
            cum += particles_per_task[rank];
            
            output(format("%12d %25d %25d %25d\n") % rank % (_mpi_first_particle[rank] + 1) % _mpi_last_particle[rank] % _mpi_num_particles[rank], 1);
        }
        ::my_first_particle = _mpi_first_particle[::my_rank];
        ::my_last_particle  = _mpi_last_particle[::my_rank];
    }
#endif

}
