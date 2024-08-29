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

#if defined(POLTMP)
// Define column_vector type used by DLib
typedef dlib::matrix<double,0,1> column_vector;
#endif

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
//POLWAS extern proj::Lot rng;
extern proj::Lot::SharedPtr rng;

namespace proj {

    class Proj {
        public:
                                         Proj();
                                         ~Proj();

            void                         clear();
            void                         processCommandLineOptions(int argc, const char * argv[]);

#if defined(USING_MPI)
            void                         particleLoopMPI(unsigned step, const vector<unsigned> & update_seeds);
#else
            void                         particleLoopStd(unsigned step, const vector<unsigned> & update_seeds);
#endif
            void                         run();
            
            void                         memoryReport(ofstream & memf) const;
            void                         debugCheckEnsemble(const SMC & ensemble) const;
            
        private:
        
            void                         updateCountMap(map<string, unsigned> & count_map, string & newick,    unsigned count);
            //unsigned                    countSpeciations(const vector<Particle> & particles);
            //void debugCheckCountMap(const map<string, unsigned> & count_map, unsigned target);
            void                         parseRelRateDefinition(string & s);
            
            void                         debugShowSpecies(double speciation_incr, vector<Node::species_tuple_t> & species_tuples, unsigned which_gene, G::species_t which_spp);
            void                         debugShowStringVector(string title, const vector<string> & svect) const;
            void                         debugShowStringUnsignedMap(string title, const map<string, unsigned> & sumap) const;

            static string                inventName(unsigned k, bool lower_case);
            void                         outputNexusTreefile(string fn, const vector<string> & newicks) const;
            
            void                         simulateData(bool chib);
            void                         simulateTrees(Particle & particle);
            void                         outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect);
            void                         outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric, const vector<string> & newick_gene_trees_numeric);
            void                         calcCoalLikeForSpecified();
            void                         testSecondLevelSMC();
            void                         chibSim(/*const Particle & test_particle*/);
            
#if defined(POLTMP)
            void                         geneTreeExperiment(bool forest, bool smctree);
#endif

            void                         selectParticlesToKeep(list<Particle> & first_level_particles, vector<unsigned> & kept);

            void                         readData();
            void                         setRelativeRates();
            unsigned                     buildSpeciesMap(bool taxa_from_data);
            void                         outputGeneTreesToFile(string fn,
                                            const vector<string> & newicks) const;
            
            string                       _data_file_name;
            string                       _start_mode;
            unsigned                     _chib_nreps;
            unsigned                     _niter;
            Partition::SharedPtr         _partition;
            Data::SharedPtr              _data;
            
            SpeciesForest::SharedPtr        _species_tree_ref;
            vector<GeneForest::SharedPtr>   _gene_tree_refs;
            
            double                       _entropy_percent_cutoff;
            bool                         _use_gpu;
            bool                         _ambig_missing;
            unsigned                     _nsimspecies;
            vector<unsigned>             _nsimtaxaperspecies;
            vector<unsigned>             _nsites_per_gene;
            int                          _track_split;
            unsigned                     _rnseed;
            bool                         _sort_forests;
            double                       _visualization_cutoff;
                        
            map<string, double>          _relrate_map;
            
            //double                       _theta;
            //double                       _lambda;
                        
            //double                       _theta_delta;
            //double                       _lambda_delta;
            
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
        _use_gpu                    = true;
        _ambig_missing              = true;
        _sort_forests               = false;
        _visualization_cutoff       = 0.99;
        _entropy_percent_cutoff     = 99.7;

        // data related
        _data                       = nullptr;
        _data_file_name             = "";
        _species_tree_ref           = nullptr;
        _gene_tree_refs.clear();
        _start_mode                 = "smc";
        _niter                      = 1;
        _chib_nreps                 = 1000;
        _partition.reset(new Partition());
        
        // simulation related
        _nsimspecies = 5;
        _nsimtaxaperspecies = {2,2,2,2,2};
                
        //_theta = 0.1;
        //_lambda = 1.0;
        
        //_theta_delta = 1.0;
        //_lambda_delta = 20.0;
        
        _ntries_theta = 100;
        _ntries_lambda = 100;
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> partition_subsets;
        vector<string> partition_relrates;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile",  value(&_data_file_name), "name of a data file in NEXUS format")
        ("speciestreeref",  value(&G::_species_tree_ref_file_name), "name of a tree file containing a single reference species tree")
        ("genetreeref",  value(&G::_gene_trees_ref_file_name), "name of a tree file containing a reference gene tree for each locus")
        ("startmode", value(&_start_mode), "if 'sim', simulate gene trees, species tree, and data; if 'smc', estimate from supplied datafile; if 'chib', computes prior probability of species species and gene tree topologies; if 'spec', computes coalescent likelihood for specified speciestreeref and genetreeref; if '2ndlevel', tests second-level SMC from gene trees supplied by genetreeref")
        ("chibreps", value(&_chib_nreps)->default_value(1000), "number of simulation replicates (default 1000) used in computing log prior probability of species and gene tree topologies (only used if startmode is \"chib\")")
        ("niter", value(&_niter), "number of iterations, where one iteration involves SMC of gene trees give species tree combined with an SMC of species tree given gene trees")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("relrate",  value(&partition_relrates), "a relative rate for a previously-defined subset; format first:3.1")
        ("entcut", value(&_entropy_percent_cutoff)->default_value(99.7), "when multithreading, determines evenness of distribution of particles to threads (should be 90-100, default 99.7)")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("verbosity",  value(&G::_verbosity)->default_value(0), "0, 1, or 2: higher number means more output")
        ("nspecies",  value(&_nsimspecies)->default_value(1), "number of species (only used if startmode sim specified)")
        ("ntaxaperspecies",  value(&_nsimtaxaperspecies), "number of taxa sampled per species (only used if simulate specified); should be _nimspecies of these entries, one for each species simulated")
        ("nparticles",  value(&G::_nparticles)->default_value(500), "number of particles in a population for joint estimation")
        ("nkept",  value(&G::_nkept)->default_value(100), "number of particles from joint smc kept for conditional smc")
        ("nspeciesparticles",  value(&G::_nparticles2)->default_value(1000), "number of particles in a population for species tree only estimation")
        ("nthreads",  value(&G::_nthreads)->default_value(3), "number of threads (each thread will handle nparticles/nthreads particle updates)")
        ("priorpost", value(&G::_prior_post)->default_value(false), "use prior-post approach to choose coalescence joins (default is prior-prior)")
        ("phi",  value(&G::_phi)->default_value(1.0), "power to which particle weight is raised")
        ("lambda",  value(&G::_lambda)->default_value(10.9), "per lineage speciation rate assumed for the species tree")
        ("compression",  value(&G::_treefile_compression)->default_value(0), "0 = no compression; 1 = intermediate compression; 2 = only unique trees saved")
#if defined(EST_THETA)
        ("invgammashape", value(&G::_invgamma_shape)->default_value(2.0), "shape parameter of inverse gamma prior distribution used for species-specific theta values (default 2.0)")
        ("freezethetamean",  value(&G::_theta_mean_frozen)->default_value(false), "if true, every species uses fixedthetamean (ignored if fixedthetamean is negative); if false, theta for any given species is drawn from InverseGamma(invgammashape, fixedthetamean) (default false)")
        ("fixedthetamean",  value(&G::_theta_mean_fixed)->default_value(-1.0), "mean theta will be fixed to this value if specified; if not specified, theta mean will vary among particles (default -1.0: negative value indicates \"not specified\")")
#else
        ("theta",  value(&G::_theta)->default_value(0.05), "coalescent parameter assumed for gene trees")
#endif
        ("thetaproposalmean",  value(&G::_theta_proposal_mean)->default_value(0.1), "mean of exponential proposal distribution for theta")
        ("thetapriormean",  value(&G::_theta_prior_mean)->default_value(1.0), "mean of exponential prior for theta")
        ("lambdapriormean",  value(&G::_lambda_prior_mean)->default_value(1.0), "mean of exponential prior for lambda")
        ("ntriestheta",  value(&_ntries_theta)->default_value(100), "number of multiple-try Metropolis trials when updating theta")
        ("ntrieslambda",  value(&_ntries_lambda)->default_value(100), "number of multiple-try Metropolis trials when updating lambda")
        //("thetadelta",  value(&_theta_delta)->default_value(1.0), "window width for updating theta")
        //("lambdadelta",  value(&_lambda_delta)->default_value(20.0), "window width for updating lambda")
        //("updatetheta",  value(&G::_update_theta)->default_value(true), "if yes, update theta at the end of each iteration")
        //("updatelambda",  value(&G::_update_lambda)->default_value(true), "if yes, update lambda at the end of each iteration")
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
            output(format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version, 1);
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
        
        // If user specified --priorpost on command line, check to ensure
        // that the user is not trying to simulate with priorpost set
        // (that won't work because priorpost requires data).
        if (vm.count("priorpost") > 0) {
            if (G::_prior_post && _start_mode == "sim") {
                throw XProj("Cannot choose simulate for startmode and yes for prior_post at the same time");
            }
        }
        
        // If user specified --nkept on command line, check to ensure
        // that nkept <= nparticles.
        if (vm.count("nkept") > 0) {
            if (G::_nkept > G::_nparticles) {
                throw XProj("nkept cannot be greater than nparticles");
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
        
        G::_gene_names.clear();

        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        output(format("\nNumber of taxa: %d\n") % _data->getNumTaxa(), 2);
        output(format("Number of partition subsets: %d") % nsubsets, 2);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(nsubsets);
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Set length of partials for gene g
            ps.setNElements(G::_nstates*_data->getNumPatternsInSubset(subset), subset);
            
            DataType dt = _partition->getDataTypeForSubset(subset);
            G::_gene_names.push_back(_data->getSubsetName(subset));
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
            
            G::_relrate_for_gene.clear();
            unsigned total_nsites = _partition->getNumSites();
            double mean_rate = 0.0;
            for (unsigned g = 0; g < G::_ngenes; g++) {
                unsigned gnsites = _partition->numSitesInSubset(g);
                string gname = _partition->getSubsetName(g);
                double r = 0.0;
                if (_relrate_map.count(gname) == 0)
                    throw XProj(format("Proj::setRelativeRates failed because key \"%s\" does not exist in _relrate_map") % gname);
                else {
                    r = _relrate_map.at(gname);
                }
                G::_relrate_for_gene[g] = r;
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
            for (unsigned g = 0; g < G::_ngenes; g++) {
                G::_relrate_for_gene[g] = 1.0;
            }
        }
    }

    inline unsigned Proj::buildSpeciesMap(bool taxa_from_data) {
        // Populates G::_species_names and G::_taxon_to_species
        // For example, given these taxon names
        //  t1^A, t2^A, t3^B, t4^B, t5^B, t6^C, t7^C, t8^D, t9^D, t10^D
        // G::_species_names = ["A", "B", "C", "D"]
        // G::_taxon_to_species["t1^A"]  = 0
        // G::_taxon_to_species["t2^A"]  = 0
        // G::_taxon_to_species["t3^B"]  = 1
        // G::_taxon_to_species["t4^B"]  = 1
        // G::_taxon_to_species["t5^B"]  = 1
        // G::_taxon_to_species["t6^C"]  = 2
        // G::_taxon_to_species["t7^C"]  = 2
        // G::_taxon_to_species["t8^D"]  = 3
        // G::_taxon_to_species["t9^D"]  = 3
        // G::_taxon_to_species["t10^D"] = 3
        // G::_taxon_to_species["A"]     = 0
        // G::_taxon_to_species["B"]     = 1
        // G::_taxon_to_species["C"]     = 2
        // G::_taxon_to_species["D"]     = 3
        unsigned nspecies = 0;
        map<string, unsigned> species_name_to_index;
        G::_taxon_to_species.clear();

        output("\nMapping taxon names to species index:\n", 2);
        unsigned ntax = (unsigned)G::_taxon_names.size();
        assert(ntax > 0);
        if (taxa_from_data) {
            // Assume taxon names are already stored in _data object and no
            // species names have yet been stored
            G::_species_names.clear();
            const Data::taxon_names_t & tnames = _data->getTaxonNames();
            assert(tnames.size() > 0);
            ntax = (unsigned)tnames.size();
            for (auto & tname : tnames) {
                string species_name = Node::taxonNameToSpeciesName(tname);
                output(format("  %s --> %s\n") % tname % species_name, 2);
                unsigned species_index = ntax;
                if (species_name_to_index.find(species_name) == species_name_to_index.end()) {
                    // species_name not found
                    species_index = nspecies;
                    G::_species_names.push_back(species_name);
                    species_name_to_index[species_name] = nspecies++;
                } else {
                    // species_name found
                    species_index = species_name_to_index[species_name];
                }
                G::_taxon_to_species[tname] = species_index;
            }
        }
        else {
            // Assume G::_species_names and G::_taxon_names are already populated
            // and _data is null because no data is needed for this analysis ("spec" mode)
            
            // First build species_name_to_index from G::_species_names
            unsigned s = 0;
            for (auto & species_name : G::_species_names) {
                species_name_to_index[species_name] = s++;
            }
            
            // Now build G::_taxon_to_species from G::_taxon_names and species_name_to_index
            for (auto & tname : G::_taxon_names) {
                string species_name = Node::taxonNameToSpeciesName(tname);
                unsigned species_index = 0;
                if (species_name_to_index.count(species_name) == 0) {
                    // species_name not found
                    throw XProj(format("Expecting species name \"%s\" to be found in global species names vector but it was not") % species_name);
                }
                else {
                    // species_name found
                    species_index = species_name_to_index[species_name];
                }
                output(format("  %s --> %s (%d)\n") % tname % species_name % species_index, 2);
                G::_taxon_to_species[tname] = species_index;
            }
        }
        
        output("\nMapping species names to species index:\n", 2);
        for (auto & sname : G::_species_names) {
            unsigned species_index = 0;
            if (species_name_to_index.count(sname) == 0)
                throw XProj(format("Proj::buildSpeciesMap failed because key \"%s\" does not exist in species_name_to_index map") % sname);
            else {
                species_index = species_name_to_index.at(sname);
            }
            output(format("  %s --> %d\n") % sname % species_index, 2);
            
            // Note: despite appearances, this next line does not
            // overwrite anything. We need to be able to map taxon
            // names in species trees to a species index as well as
            // taxon names in gene trees
            G::_taxon_to_species[sname] = species_index;
        }

        return (unsigned)G::_species_names.size();
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
        assert(G::_ngenes == newicks.size());
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto newick : newicks) {
            streef << str(format("  tree %s = [&R] %s;\n") % G::_gene_names[t++] % newick);
        }
        streef << "end;\n";
        streef.close();
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
        double logibase26 = (k > 0 ? log(k)/log(26) : 0);
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
    
    inline void Proj::debugShowSpecies(double speciation_incr, vector<Node::species_tuple_t> & species_tuples, unsigned which_gene, G::species_t which_spp) {
        output(format("\ndebugShowSpecies (speciation increment = %.5f):\n") % speciation_incr, 1);

        unsigned ntuples = (unsigned)species_tuples.size();
        output(format("%12s %12s %s\n")  % "nlineages" % "gene" % "species", 1);
        
        for (unsigned i = 0; i < ntuples; i++) {
            Node::species_tuple_t & species_tuple = species_tuples[i];
            unsigned              n = get<0>(species_tuple);
            unsigned              g = get<1>(species_tuple);
            G::species_t  s = get<2>(species_tuple);
            string               ss = G::speciesStringRepresentation(s);
            assert(g > 0);
            output(format("%12d %12d %s\n")  % n % g % ss, 1);
        }
    }

    inline void Proj::simulateTrees(Particle & particle) {
        unsigned ngenes   = G::_ngenes;
        unsigned ntaxa    = G::_ntaxa;
        
        // Start with trivial species forest
        particle.resetSpeciesForest();

#if defined(EST_THETA)
        particle.setThetas();
#endif

        // Start with trivial gene forests
        particle.resetGeneForests(/*compute_partials*/false);
        
        assert(particle.getGeneForests().size() == ngenes);
        
        // Determine total number of steps required to build all gene trees and the species tree
        unsigned nsteps = ngenes*(ntaxa - 1);
        
        vector<unsigned> update_seeds(nsteps);
        G::generateUpdateSeeds(update_seeds);
        
#if defined(ONE_LOCUS_PER_STEP)
        vector< pair<double, unsigned> > locus_ordering(G::_ngenes);
#endif

        // Add coalescent events and speciation events until all trees are built
        for (unsigned step = 0; step < nsteps; step++) {
#if defined(ONE_LOCUS_PER_STEP)
            unsigned step_modulus = step % G::_ngenes;
            if (step_modulus == 0) {
                // Time to choose a new random locus ordering
                for (unsigned i = 0; i < G::_ngenes; i++) {
                    double u = rng->uniform();
                    locus_ordering[i] = make_pair(u, i);
                }
                sort(locus_ordering.begin(), locus_ordering.end());
            }
            particle.proposeCoalescence(update_seeds[step], step, locus_ordering[step_modulus].second, /*particle index*/0, /*compute_partial*/false, /*make_permanent*/true);
#else
            particle.proposeCoalescence(update_seeds[step], step, /*particle index*/0, /*compute_partial*/false, /*make_permanent*/true);
#endif
        }
        
        particle.refreshHeightsInternalsPreorders();
    }
    
    inline void Proj::outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect) {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;
        
        // Output tree file containing true species tree
        string newick_species_tree_alpha = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        string newick_species_tree_numeric = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
        if (G::_nspecies > 1) {
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
            //output(format("  %12.9f (%.9f)\n") % gf.getHeight() % (gf.getHeight()/G::_theta), 2);
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
        for (string nm : G::_species_names) {
            string comma = (i < G::_nspecies ? "," : "");
            jsf << str(format("  %d: \"%s\"%s\n") % i % nm % comma);
            ++i;
        }
        jsf << "};\n\n";
        
        jsf << str(format("let species_newick = \"%s\";\n\n") % newick_species_tree_numeric);
        
        jsf << "let gene_translate = {\n";
        i = 1;
        for (string nm : G::_taxon_names) {
            string comma = (i < G::_ntaxa ? "," : "");
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
    
    inline void Proj::simulateData(bool chib) {
//#if defined(EST_THETA)
//        throw XProj("simulating data not yet implemented when EST_THETA is #defined");
//#endif
        G::_nspecies = _nsimspecies;
        G::_ntaxa = (unsigned)accumulate(_nsimtaxaperspecies.begin(), _nsimtaxaperspecies.end(), 0);

        output("Simulating sequence data under multispecies coalescent model:\n", 2);
        output(format("  theta  = %.5f\n") % G::_theta, 2);
        output(format("  lambda = %.5f\n") % G::_lambda, 2);
        output(format("  Number of species = %d\n") % _nsimspecies, 2);
        
        // Expected height of species tree is (1/2 + 1/3 + ... + 1/nspecies)/lambda
        // Approximation to (1/2 + 1/3 + ... + 1/nspecies) is
        //   ln(nspecies) + 0.58 (Euler's constant)
        //   see http://www.dimostriamogoldbach.it/en/inverses-integers-sum/
        // 0.25 = Exact (nspecies = 5, lambda = 5)
        // 0.44 = Approximate
        double expected_species_tree_height = 0.0;
        for (unsigned i = 2; i <= G::_nspecies; i++) {
            expected_species_tree_height += 1.0/i;
        }
        expected_species_tree_height /= G::_lambda;
        output(format("  Expected species tree height = %.5f\n") % expected_species_tree_height, 2);
        
        // Expected gene tree height without species tree constraints
        // is theta/n(n-1) + theta/(n-1)(n-2) + ... + theta/(2)(1)
        double expected_gene_tree_height = 0.0;
        for (unsigned i = 2; i <= G::_ntaxa; i++) {
            expected_gene_tree_height += 1.0/(i*(i-1));
        }
        expected_gene_tree_height *= G::_theta;
        output(format("  Expected gene tree height (unconstrained) = %.5f\n") % expected_gene_tree_height, 2);
        
        // Interrogate _partition to determine number of genes, gene names, and
        // number of sites in each gene
        G::_ngenes = _partition->getNumSubsets();
        G::_nsites_per_gene.resize(G::_ngenes);
        G::_gene_names.resize(G::_ngenes);
        unsigned total_nsites = 0;
        for (unsigned g = 0; g < G::_ngenes; g++) {
            unsigned gnsites = _partition->numSitesInSubset(g);
            total_nsites += gnsites;
            string gname = _partition->getSubsetName(g);
            G::_nsites_per_gene[g] = gnsites;
            G::_gene_names[g] = gname;
        }

        setRelativeRates();

        // Invent taxon/species names and numbers, and create
        // taxpartition vector used later in paup command file
        G::_species_names.resize(G::_nspecies);
        G::_taxon_names.resize(G::_ntaxa);
        unsigned k = 0;
        vector<string> taxpartition;
        for (unsigned i = 0; i < G::_nspecies; ++i) {
            string species_name = inventName(i, false);
            G::_species_names[i] = species_name;
            for (unsigned j = 0; j < _nsimtaxaperspecies[i]; ++j) {
                taxpartition.push_back(G::_species_names[i]);
                string taxon_name = str(format("%s^%s") % inventName(k, true) % G::_species_names[i]);
                G::_taxon_names[k] = taxon_name;
                G::_taxon_to_species[taxon_name] = i;
                ++k;
            }
        }

        // Report species names created
        output("  Species names:\n", 2);
        for (unsigned i = 0; i < G::_nspecies; ++i) {
            output(format("    %s\n") % G::_species_names[i], 2);
        }

        // Report taxon names created
        output("  Taxon names:\n", 2);
        for (unsigned i = 0; i < G::_ntaxa; ++i) {
            string taxon_name = G::_taxon_names[i];
            output(format("    %s\n") % taxon_name, 2);
        }
        
        // Create data object
        assert(!_data);
        _data = Data::SharedPtr(new Data());
        _data->setTaxonNames(G::_taxon_names);
        _data->setPartition(_partition);
        
        // Build species tree and gene trees jointly
        Particle particle;
        particle.setData(_data);
        //POLWAS particle.setSeed(rng.randint(1,9999));
        particle.setSeed(rng->randint(1,9999));
        simulateTrees(particle);
        
        SpeciesForest      & species_forest = particle.getSpeciesForest();
        vector<GeneForest> & gene_forests   = particle.getGeneForests();
        outputTrees(species_forest, gene_forests);
        
        if (chib) {
            // Make a copy of the species forest that is managed by shared pointer _species_tree_ref
            _species_tree_ref = SpeciesForest::SharedPtr(new SpeciesForest());
            *_species_tree_ref = species_forest;
            
            // Make a copy of each gene tree forest that is managed by a shared pointer
            _gene_tree_refs.resize(gene_forests.size());
            for (unsigned g = 0; g < gene_forests.size(); g++) {
                _gene_tree_refs[g] = GeneForest::SharedPtr(new GeneForest());
                *(_gene_tree_refs[g]) = gene_forests[g];
            }
        }

        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(G::_ngenes);

        // Simulate sequence data
        unsigned starting_site = 0;
        for (unsigned g = 0; g < G::_ngenes; ++g) {
            ps.setNElements(4*G::_nsites_per_gene[g], g);
            gene_forests[g].simulateData(particle.getLot(), _data, starting_site, G::_nsites_per_gene[g]);
            starting_site += G::_nsites_per_gene[g];
        }
                       
        // Output data to file
        _data->compressPatterns();
        _data->writeDataToFile(_data_file_name);
        output(format("  Sequence data saved in file \"%s\"\n") % _data_file_name, 1);

        if (!chib) {
            // Output a PAUP* command file for estimating the species tree using
            // svd quartets and qage
            output("  PAUP* commands saved in file \"svd-qage.nex\"\n", 1);
            ofstream paupf("svd-qage.nex");
            paupf << "#NEXUS\n\n";
            paupf << "begin paup;\n";
            paupf << "  log start file=svdout.txt replace;\n";
            paupf << "  exe sim.nex;\n";
            paupf << "  taxpartition species (vector) = " << join(taxpartition," ") << ";\n";
            paupf << "  svd taxpartition=species;\n";
            paupf << "  roottrees;\n";
            paupf << "  qage taxpartition=species patprob=exactjc outUnits=substitutions treefile=svd.tre replace;\n";
            paupf << "  log stop;\n";
            paupf << "  quit;\n";
            paupf << "end;\n";
            paupf.close();
        }
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

    inline void Proj::calcCoalLikeForSpecified() {
        // Read in the reference species tree
        G::_nexus_taxon_map.clear();
        map<unsigned,unsigned> species_taxon_map;
        vector<string> species_tree_names;
        vector<string> species_newicks;
        Forest::readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, species_taxon_map, species_tree_names, species_newicks);
        G::_nspecies = (unsigned)G::_species_names.size();
                
        // Read in the reference gene trees
        map<unsigned,unsigned> gene_taxon_map;
        vector<string> gene_tree_names;
        vector<string> gene_newicks;
        GeneForest::readTreefile(G::_gene_trees_ref_file_name, /*skip*/0, G::_taxon_names, gene_taxon_map, gene_tree_names, gene_newicks);
        G::_ntaxa = (unsigned)G::_taxon_names.size();
        
        // Copy taxon names to global variable _taxon_names
        buildSpeciesMap(/*taxa_from_data*/false);
        
        // Create a particle in which to store species and gene trees
        Particle particle;
        
        // Build the species tree
        SpeciesForest & sppref = particle.getSpeciesForest();
        G::_nexus_taxon_map = species_taxon_map;
        sppref.buildFromNewick(species_newicks[0]);
        
        // Build the gene trees
        G::_ngenes = (unsigned)gene_newicks.size();
        vector<GeneForest> & gfref = particle.getGeneForests();
        gfref.resize(gene_newicks.size());
        unsigned g = 0;
        G::_nexus_taxon_map = gene_taxon_map;
        for (auto newick : gene_newicks) {
            gfref[g].setParticle(&particle);
            gfref[g].setGeneIndex(g);
            gfref[g++].buildFromNewick(newick);
        }

        vector<Forest::coalinfo_t> coalinfo_vect;
        particle.calcLogCoalescentLikelihood(coalinfo_vect,
            /*integrate_out_thetas*/true, /*verbose*/true);
    }

#if defined(POLTMP)
    GeneForest * optgf = nullptr;
        
    Node * sppA = nullptr;
    Node * sppAl = nullptr;
    Node * sppAr = nullptr;
    Node * sppB = nullptr;
    Node * sppBl = nullptr;
    Node * sppBr = nullptr;
    Node * sppC = nullptr;
    Node * sppCl = nullptr;
    Node * sppCr = nullptr;
    Node * sppD = nullptr;
    Node * sppDl = nullptr;
    Node * sppDr = nullptr;
    Node * sppE = nullptr;
    Node * sppEl = nullptr;
    Node * sppEr = nullptr;
    Node * sppBC = nullptr;
    Node * sppDE = nullptr;
    Node * sppCDE = nullptr;
    Node * sppBCDE = nullptr;
    Node * sppABCDE = nullptr;
    
    double minTrueForest(const column_vector & m) {
        assert(optgf);
        GeneForest & gf = *optgf;

        // Unpack variables
        double hA = m(0);
        double hBC = m(1);
        double hDE = m(2);
        double logitB = m(3);
        double logitC = m(4);
        double logitD = m(5);
        double logitE = m(6);
        
        double hB = hBC*exp(logitB)/(1.0 + exp(logitB));
        double hC = hBC*exp(logitC)/(1.0 + exp(logitC));
        double hD = hDE*exp(logitD)/(1.0 + exp(logitD));
        double hE = hDE*exp(logitE)/(1.0 + exp(logitE));
        
        assert(hB < hBC);
        assert(hC < hBC);
        assert(hD < hDE);
        assert(hE < hDE);

        sppA->setHeight(hA);
        sppAl->setHeight(0.0);
        sppAr->setHeight(0.0);
        
        sppBC->setHeight(hBC);
        
        sppB->setHeight(hB);
        sppBl->setHeight(0.0);
        sppBr->setHeight(0.0);

        sppC->setHeight(hC);
        sppCl->setHeight(0.0);
        sppCr->setHeight(0.0);

        sppDE->setHeight(hDE);
        
        sppD->setHeight(hD);
        sppDl->setHeight(0.0);
        sppDr->setHeight(0.0);

        sppE->setHeight(hE);
        sppEl->setHeight(0.0);
        sppEr->setHeight(0.0);
        
        sppAl->setEdgeLength(hA);
        sppAr->setEdgeLength(hA);

        sppBl->setEdgeLength(hB);
        sppBr->setEdgeLength(hB);

        sppCl->setEdgeLength(hC);
        sppCr->setEdgeLength(hC);

        sppDl->setEdgeLength(hD);
        sppDr->setEdgeLength(hD);

        sppEl->setEdgeLength(hE);
        sppEr->setEdgeLength(hE);
        
        sppB->setEdgeLength(hBC - hB);
        sppC->setEdgeLength(hBC - hC);
        sppD->setEdgeLength(hDE - hD);
        sppE->setEdgeLength(hDE - hE);

        sppA->setEdgeLength(0.0);
        sppBC->setEdgeLength(0.0);
        sppDE->setEdgeLength(0.0);
        
        gf.calcPartialArray(sppA);
        gf.calcPartialArray(sppB);
        gf.calcPartialArray(sppC);
        gf.calcPartialArray(sppD);
        gf.calcPartialArray(sppE);
        gf.calcPartialArray(sppBC);
        gf.calcPartialArray(sppDE);
        return gf.calcLogLikelihood();
    }
    
    double minWrongForest(const column_vector & m) {
        assert(optgf);
        GeneForest & gf = *optgf;

        double hA       = m(0); // starting height of A, ancestor of a^A and b^A
        double hB       = m(1); // starting height of B, ancestor of c^B and d^B
        double hCDE     = m(2); // starting height of ancestor of C and D+E
        double logitpDE = m(3); // log-proportion of height of D+E ancestor
                              //   to height of ancestor of C and D+E
        double logitpC  = m(4); // log-proportion of height of e^C and f^C ancestor
                              //   to height of ancestor of C and D+E
        double logitpD  = m(5); // log-proportion of height of g^D and h^D ancestor
                              //   to height of D+E ancestor
        double logitpE  = m(6); // log-proportion of height of i^E and j^E ancestor

        sppA->setHeight(hA);
        sppAl->setHeight(0.0);
        sppAr->setHeight(0.0);
        sppB->setHeight(hB);
        sppBl->setHeight(0.0);
        sppBr->setHeight(0.0);
        sppCDE->setHeight(hCDE);
        double hC = hCDE*exp(logitpC)/(1.0 + exp(logitpC));
        sppC->setHeight(hC);
        sppCl->setHeight(0.0);
        sppCr->setHeight(0.0);
        double hDE = hCDE*exp(logitpDE)/(1.0 + exp(logitpDE));
        sppDE->setHeight(hDE);
        double hD = hDE*exp(logitpD)/(1.0 + exp(logitpD));
        sppD->setHeight(hD);
        sppDl->setHeight(0.0);
        sppDr->setHeight(0.0);
        double hE = hDE*exp(logitpE)/(1.0 + exp(logitpE));
        sppE->setHeight(hE);
        sppEl->setHeight(0.0);
        sppEr->setHeight(0.0);
        
        sppA->setEdgeLength(0.0);
        sppAl->setEdgeLength(hA);
        sppAr->setEdgeLength(hA);
        
        sppB->setEdgeLength(0.0);
        sppBl->setEdgeLength(hB);
        sppBr->setEdgeLength(hB);
        
        double edgelenC = hCDE - hC;
        assert(edgelenC >= 0.0);
        sppC->setEdgeLength(edgelenC);
        sppCl->setEdgeLength(hC);
        sppCr->setEdgeLength(hC);
        
        double edgelenD = hDE - hD;
        assert(edgelenD >= 0.0);
        sppD->setEdgeLength(edgelenD);
        sppDl->setEdgeLength(hD);
        sppDr->setEdgeLength(hD);
        
        double edgelenE = hDE - hE;
        assert(edgelenE >= 0.0);
        sppE->setEdgeLength(edgelenE);
        sppEl->setEdgeLength(hE);
        sppEr->setEdgeLength(hE);
        
        double edgelenDE = hCDE - hDE;
        assert(edgelenDE >= 0.0);
        sppDE->setEdgeLength(edgelenDE);
        
        sppCDE->setEdgeLength(0.0);
                
        gf.calcPartialArray(sppA);
        gf.calcPartialArray(sppB);
        gf.calcPartialArray(sppC);
        gf.calcPartialArray(sppD);
        gf.calcPartialArray(sppE);
        gf.calcPartialArray(sppDE);
        gf.calcPartialArray(sppCDE);
        return gf.calcLogLikelihood();
    }
    
    inline void Proj::geneTreeExperiment(bool forest, bool smctree) {
        // forest   smctree   newick
        //    yes        no   A, (B,C), (D,E)
        //    yes       yes   A, B, (C,(D,E))
        //     no        no   (A,((B,C),(D,E))
        //     no       yes   (A,(B,(C,(D,E)))
        readData();
        G::_ngenes = _data->getNumSubsets();
        assert(G::_ngenes > 0);

        // Copy taxon names to global variable _taxon_names
        G::_ntaxa = _data->getNumTaxa();
        _data->copyTaxonNames(G::_taxon_names);
                    
        for (unsigned g = 0; g < G::_ngenes; g++) {
            GeneForest::computeLeafPartials(g, _data);
        }

        // Read in the true species tree
        G::_nexus_taxon_map.clear();
        map<unsigned,unsigned> species_taxon_map;
        vector<string> species_tree_names;
        vector<string> species_newicks;
        Forest::readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, species_taxon_map, species_tree_names, species_newicks);
        G::_nspecies = (unsigned)G::_species_names.size();
                
        // Read in the true gene tree
        map<unsigned,unsigned> gene_taxon_map;
        vector<string> gene_tree_names;
        vector<string> gene_newicks;
        GeneForest::readTreefile(G::_gene_trees_ref_file_name, /*skip*/0, G::_taxon_names, gene_taxon_map, gene_tree_names, gene_newicks);
        G::_ntaxa = (unsigned)G::_taxon_names.size();
        
        // Copy taxon names to global variable _taxon_names
        buildSpeciesMap(/*taxa_from_data*/false);
        
        // Create a particle in which to store species and gene trees
        Particle particle;
        particle.setData(_data);
        
        // Build the species tree
        SpeciesForest & sppref = particle.getSpeciesForest();
        G::_nexus_taxon_map = species_taxon_map;
        sppref.buildFromNewick(species_newicks[0]);
        
        // Build the gene trees
        G::_ngenes = (unsigned)gene_newicks.size();
        vector<GeneForest> & gfref = particle.getGeneForests();
        gfref.resize(gene_newicks.size());
        
        // Only considering locus 2
        unsigned g = 1;
        const string & newick = gene_newicks[g];
        output(format("locus %d:\n") % (g+1), 2);
        
        G::_nexus_taxon_map = gene_taxon_map;
        GeneForest & gf = gfref[g];
        gf.setParticle(&particle);
        gf.setGeneIndex(g);
        output(format("  newick = \"%s\"\n") %  newick, 2);
        gf.buildFromNewick(newick);
        
        // Calculate log likelihood for locus 2 gene tree
        gf.setData(_data);
        gf.computeAllPartials();
        double logL = gf.calcLogLikelihood();
        output(format("  logL = %.5f\n") % logL, 2);
        
        // Create pointers to important nodes
        //gf.debugShowPreorder();
        
        sppA = gf.findNodeNumbered(17);
        sppAl = gf.findNodeNumbered(0);
        sppAr = gf.findNodeNumbered(1);

        sppB = gf.findNodeNumbered(13);
        sppBl = gf.findNodeNumbered(2);
        sppBr = gf.findNodeNumbered(3);

        sppC = gf.findNodeNumbered(14);
        sppCl = gf.findNodeNumbered(4);
        sppCr = gf.findNodeNumbered(5);

        sppD = gf.findNodeNumbered(11);
        sppDl = gf.findNodeNumbered(6);
        sppDr = gf.findNodeNumbered(7);

        sppE = gf.findNodeNumbered(10);
        sppEl = gf.findNodeNumbered(8);
        sppEr = gf.findNodeNumbered(9);

        sppBC = gf.findNodeNumbered(15);
        sppDE = gf.findNodeNumbered(12);
        sppBCDE = gf.findNodeNumbered(16);
        sppABCDE = gf.findNodeNumbered(18);

        double sppAheight  = sppA->getHeight();
        double sppBheight  = sppB->getHeight();
        double sppCheight  = sppC->getHeight();
        double sppDheight  = sppD->getHeight();
        double sppEheight  = sppE->getHeight();
        double sppBCheight = sppBC->getHeight();
        double sppDEheight = sppDE->getHeight();
        double sppBCDEheight = sppBCDE->getHeight();
        double sppABCDEheight = sppABCDE->getHeight();
        double hBCDE  = sppBCDEheight; //(sppABCDEheight + sppBCDEheight)/2.0;
        double hstart = forest ? sppABCDEheight : sppBCDEheight;
        double hstop  = sppDEheight;
        double hincr  = (hstart - hstop)/20.0;
        output(format("sppAheight    = %.5f\n") % sppAheight, 1);
        output(format("sppBheight    = %.5f\n") % sppBheight, 1);
        output(format("sppCheight    = %.5f\n") % sppCheight, 1);
        output(format("sppDheight    = %.5f\n") % sppDheight, 1);
        output(format("sppEheight    = %.5f\n") % sppEheight, 1);
        output(format("sppBCheight   = %.5f\n") % sppBCheight, 1);
        output(format("sppDEheight   = %.5f\n") % sppDEheight, 1);
        output(format("sppBCDEheight = %.5f\n") % sppBCDEheight, 1);
        output(format("sppaBCDEheight = %.5f\n") % sppABCDEheight, 1);
        output(format("hstart = %.5f\n") % hstart, 1);
        output(format("hstop  = %.5f\n") % hstop, 1);
        output(format("hincr  = %.5f\n") % hincr, 1);
        
        if (forest) {
            if (smctree) {
                // (A,((B,C),(D,E)) --> (A,B,(C,(D,E)))
                // Separate sppABCDE into sppA and sppBCDE
                gf.addTwoRemoveOne(gf.getLineages(), sppA, sppBCDE, sppABCDE);
                sppA->setParent(nullptr);
                sppA->setRightSib(nullptr);
                sppBCDE->setParent(nullptr);
                sppBCDE->setRightSib(nullptr);
                sppABCDE->setLeftChild(nullptr);
                //gf.refreshAllPreorders();
                //gf.debugShowPreorder();
                                
                // Separate sppBCDE into sppBC and sppDE
                gf.addTwoRemoveOne(gf.getLineages(), sppBC, sppDE, sppBCDE);
                sppBC->setParent(nullptr);
                sppBC->setRightSib(nullptr);
                sppDE->setParent(nullptr);
                sppDE->setRightSib(nullptr);
                sppBCDE->setLeftChild(nullptr);
                sppBCDE->setRightSib(nullptr);
                sppBCDE->setParent(nullptr);

                // Separate sppBC into sppB and sppC
                gf.addTwoRemoveOne(gf.getLineages(), sppB, sppC, sppBC);
                sppB->setParent(nullptr);
                sppB->setRightSib(nullptr);
                sppC->setParent(nullptr);
                sppC->setRightSib(nullptr);
                sppBC->setLeftChild(nullptr);
                sppBC->setRightSib(nullptr);
                sppBC->setParent(nullptr);
                //gf.refreshAllPreorders();
                //gf.debugShowPreorder();

                // Join sppC and sppDE to form sppCDE
                // sppBC is no longer being used, so recyle it as sppCDE
                sppCDE = sppBC;
                sppCDE->setLeftChild(sppC);
                sppCDE->setParent(nullptr);
                sppCDE->setRightSib(nullptr);
                sppCDE->setEdgeLength(0.0);

                // Finish connecting new trio of nodes
                sppC->setRightSib(sppDE);
                sppC->setParent(sppCDE);
                sppDE->setParent(sppCDE);
                
                // Add sppCDE to and remove sppC from _lineages
                gf.removeTwoAddOne(gf.getLineages(), sppC, sppDE, sppCDE);
                gf.refreshAllPreorders();
                gf.debugShowPreorder();
            }
            else /* not smctree */ {
                // Separate sppABCDE into sppA and sppBCDE
                gf.addTwoRemoveOne(gf.getLineages(), sppA, sppBCDE, sppABCDE);
                
                // Separate sppBCDE into sppBC and sppDE
                gf.addTwoRemoveOne(gf.getLineages(), sppBC, sppDE, sppBCDE);
                
                // Tree has been converted to a forest, so recalculate preorders
                gf.refreshAllPreorders();
                //gf.debugShowPreorder();
            }
        }
        else /* not forest */ {
            if (smctree) {
                // (A,((B,C),(D,E)) --> (A,(B,(C,(D,E))))
                // Separate sppABCDE into sppA and sppBCDE
                gf.addTwoRemoveOne(gf.getLineages(), sppA, sppBCDE, sppABCDE);
                sppA->setParent(nullptr);
                sppA->setRightSib(nullptr);
                sppBCDE->setParent(nullptr);
                sppBCDE->setRightSib(nullptr);
                sppABCDE->setLeftChild(nullptr);
                //gf.refreshAllPreorders();
                //gf.debugShowPreorder();
                                
                // Separate sppBCDE into sppBC and sppDE
                gf.addTwoRemoveOne(gf.getLineages(), sppBC, sppDE, sppBCDE);
                sppBC->setParent(nullptr);
                sppBC->setRightSib(nullptr);
                sppDE->setParent(nullptr);
                sppDE->setRightSib(nullptr);
                sppBCDE->setLeftChild(nullptr);
                sppBCDE->setRightSib(nullptr);
                sppBCDE->setParent(nullptr);

                // Separate sppBC into sppB and sppC
                gf.addTwoRemoveOne(gf.getLineages(), sppB, sppC, sppBC);
                sppB->setParent(nullptr);
                sppB->setRightSib(nullptr);
                sppC->setParent(nullptr);
                sppC->setRightSib(nullptr);
                sppBC->setLeftChild(nullptr);
                sppBC->setRightSib(nullptr);
                sppBC->setParent(nullptr);
                //gf.refreshAllPreorders();
                //gf.debugShowPreorder();

                // Join sppC and sppDE to form sppCDE
                // sppBC is no longer being used, so recyle it as sppCDE
                sppCDE = sppBC;
                sppCDE->setLeftChild(sppC);
                //sppCDE->setParent(sppBCDE);
                sppCDE->setRightSib(nullptr);
                sppC->setRightSib(sppDE);
                sppC->setParent(sppCDE);
                sppDE->setParent(sppCDE);
                gf.removeTwoAddOne(gf.getLineages(), sppC, sppDE, sppCDE);
                //gf.refreshAllPreorders();
                //gf.debugShowPreorder();
                
                // Join sppB and sppCDE to form sppBCDE
                sppBCDE->setLeftChild(sppB);
                sppBCDE->setParent(sppABCDE);
                sppBCDE->setRightSib(nullptr);
                sppB->setRightSib(sppCDE);
                sppB->setParent(sppBCDE);
                sppCDE->setParent(sppBCDE);
                gf.removeTwoAddOne(gf.getLineages(), sppB, sppCDE, sppBCDE);
                //gf.refreshAllPreorders();
                //gf.debugShowPreorder();

                // Join sppA and sppBCDE to form sppABCDE
                sppABCDE->setLeftChild(sppA);
                sppABCDE->setParent(nullptr);
                sppABCDE->setRightSib(nullptr);
                sppA->setRightSib(sppBCDE);
                sppA->setParent(sppABCDE);
                sppBCDE->setParent(sppABCDE);
                gf.removeTwoAddOne(gf.getLineages(), sppA, sppBCDE, sppABCDE);
                gf.refreshAllPreorders();
                //gf.debugShowPreorder();
            }
            else {
                // nothing to do because true tree already built
            }
        }
        
        bool find_MLEs = true;
        if (find_MLEs) {
            assert(forest);
            optgf = &gf;
            if (smctree) {
                double sppCDEheight = sppBCheight;
                assert(sppCDEheight > sppDEheight);
                double sppDElogitprop = log(sppDEheight) - log(sppCDEheight-sppDEheight);
                assert(sppCDEheight > sppCheight);
                double sppClogitprop  = log(sppCheight) - log(sppCDEheight-sppCheight);
                assert(sppDEheight > sppDheight);
                double sppDlogitprop  = log(sppDheight) - log(sppDEheight-sppDheight);
                assert(sppDEheight > sppEheight);
                double sppElogitprop  = log(sppEheight) - log(sppDEheight-sppEheight);
                column_vector starting_point = {
                    sppAheight,      // starting height of A, ancestor of a^A and b^A
                    sppBheight,      // starting height of B, ancestor of c^B and d^B
                    sppCDEheight,    // starting height of ancestor of C and D+E
                    sppDElogitprop,  // logit-proportion of height of D+E ancestor
                                     //   to height of ancestor of C and D+E
                    sppClogitprop,   // logit-proportion of height of e^C and f^C ancestor
                                     //   to height of ancestor of C and D+E
                    sppDlogitprop,   // logit-proportion of height of g^D and h^D ancestor
                                     //   to height of D+E ancestor
                    sppElogitprop    // logit-proportion of height of i^E and j^E ancestor
                                     //   to height of D+E ancestor
                };
                
                double maximized_log_likelihood = dlib::find_min_using_approximate_derivatives(
                    dlib::bfgs_search_strategy(),
                    dlib::objective_delta_stop_strategy(1e-7),
                    minWrongForest,
                    starting_point,
                    -1);
                
                // Unpack variables
                double hA = starting_point(0);
                double hB = starting_point(1);
                double hCDE = starting_point(2);
                double logitDE = starting_point(3);
                double hDE = hCDE*exp(logitDE)/(1.0 + exp(logitDE));
                double logitC = starting_point(4);
                double hC = hCDE*exp(logitC)/(1.0 + exp(logitC));
                double logitD = starting_point(5);
                double hD = hDE*exp(logitD)/(1.0 + exp(logitD));
                double logitE = starting_point(6);
                double hE = hDE*exp(logitE)/(1.0 + exp(logitE));
                
                string Aclade = str(format("(a^A:%.9f,b^A:%.9f)A:0.0") % hA % hA);
                string Bclade = str(format("(c^B:%.9f,d^B:%.9f)B:0.0") % hB % hB);
                string Cclade = str(format("(e^C:%.9f,f^C:%.9f)C:%.9f") % hC % hC % (hCDE-hC));
                string Dclade = str(format("(g^D:%.9f,h^D:%.9f)D:%.9f") % hD % hD % (hDE-hD));
                string Eclade = str(format("(i^E:%.9f,j^E:%.9f)E:%.9f") % hE % hE % (hDE-hE));
                string CDEclade = str(format("(%s,(%s,%s)DE:%.9f)CDE:0.0") % Cclade % Dclade % Eclade % (hCDE-hDE));
                string newick = str(format("(%s,%s,%s)\n") % Aclade % Bclade % CDEclade);
                output("Wrong forest maximizing log-likelihood:\n  ", 1);
                output(newick, 1);
                output(format("Maximized log-likelihood: %.9f\n") % maximized_log_likelihood,1);
            }
            else {
                double hA = sppAheight;
                double hB = sppBheight;
                double hC = sppCheight;
                double hD = sppDheight;
                double hE = sppEheight;
                double hBC = sppBCheight;
                double hDE = sppDEheight;

                assert(hB < hBC);
                assert(hC < hBC);
                assert(hD < hDE);
                assert(hE < hDE);
                
                double logitB = log(hB) - log(hBC-hB);
                double logitC = log(hC) - log(hBC-hC);
                double logitD = log(hD) - log(hDE-hD);
                double logitE = log(hE) - log(hDE-hE);

                column_vector starting_point = {
                    hA,     // starting height of A, ancestor of a^A and b^A
                    hBC,    // starting height of ancestor of B and C
                    hDE,    // starting height of ancestor of D and E
                    logitB, // logit of proportion hB to hBC
                    logitC, // logit of proportion hC to hBC
                    logitD, // logit of proportion hD to hDE
                    logitE  // logit of proportion hE to hDE
                };
                
                double maximized_log_likelihood = dlib::find_min_using_approximate_derivatives(
                    dlib::bfgs_search_strategy(),
                    dlib::objective_delta_stop_strategy(1e-7),
                    minTrueForest,
                    starting_point,
                    -1);
                
                // Unpack variables
                hA = starting_point(0);
                hBC = starting_point(1);
                hDE = starting_point(2);
                logitB = starting_point(3);
                logitC = starting_point(4);
                logitD = starting_point(5);
                logitE = starting_point(6);
                
                hB = hBC*exp(logitB)/(1.0 + exp(logitB));
                hC = hBC*exp(logitC)/(1.0 + exp(logitC));
                hD = hDE*exp(logitD)/(1.0 + exp(logitD));
                hE = hDE*exp(logitE)/(1.0 + exp(logitE));
                
                assert(hB < hBC);
                assert(hC < hBC);
                assert(hD < hDE);
                assert(hE < hDE);
                
                string Aclade = str(format("(a^A:%.9f,b^A:%.9f)A:0.0") % hA % hA);
                string Bclade = str(format("(c^B:%.9f,d^B:%.9f)B:%.9f") % hB % hB % (hBC - hB));
                string Cclade = str(format("(e^C:%.9f,f^C:%.9f)C:%.9f") % hC % hC % (hBC-hC));
                string Dclade = str(format("(g^D:%.9f,h^D:%.9f)D:%.9f") % hD % hD % (hDE-hD));
                string Eclade = str(format("(i^E:%.9f,j^E:%.9f)E:%.9f") % hE % hE % (hDE-hE));
                string BCclade = str(format("(%s,%s)BC:0.0") % Bclade % Cclade);
                string DEclade = str(format("(%s,%s)DE:0.0") % Dclade % Eclade);
                string newick = str(format("(%s,%s,%s)\n") % Aclade % BCclade % DEclade);
                output("True forest maximizing log-likelihood:\n  ", 1);
                output(newick, 1);
                output(format("Maximized log-likelihood: %.9f\n") % maximized_log_likelihood,1);
            }
        }
        else {
            // Change height of sppBC and record log-likelihood for each increment
            if (forest) {
                if (smctree) {
                    output("Adjusting height of C+(D+E) node in gene forest:\n", 1);
                }
                else {
                    output("Adjusting height of B+C node in gene forest:\n", 1);
                }
            }
            else {
                if (smctree) {
                    output("Adjusting height of C+(D+E) node in full gene tree:\n", 1);
                }
                else {
                    output("Adjusting height of B+C node in full gene tree:\n", 1);
                }
            }
            
            vector< pair<double, double> > wrong_forest;
            vector< pair<double, double> > true_forest;
            vector< pair<double, double> > wrong_tree;
            vector< pair<double, double> > true_tree;

            output(format("%20s %12s\n") % "logL" % "height", 1);
            for (double h = hstart; h >= hstop; h -= hincr) {
                if (forest) {
                    if (smctree) {
                        sppCDE->setEdgeLength(0.0);
                        sppCDE->setHeight(h);
                        sppC->setEdgeLength(h - sppC->getHeight());
                        sppDE->setEdgeLength(h - sppDEheight);
                        gf.calcPartialArray(sppCDE);
                        logL = gf.calcLogLikelihood();
                        wrong_forest.push_back(make_pair(h, logL));
                    }
                    else {
                        sppBC->setEdgeLength(0.0);
                        sppBC->setHeight(h);
                        sppB->setEdgeLength(h - sppB->getHeight());
                        sppC->setEdgeLength(h - sppC->getHeight());
                        gf.calcPartialArray(sppBC);
                        gf.calcPartialArray(sppDE);
                        logL = gf.calcLogLikelihood();
                        true_forest.push_back(make_pair(h, logL));
                    }
                }
                else {
                    if (smctree) {
                        sppA->setEdgeLength(sppABCDE->getHeight() - sppA->getHeight());
                        sppB->setEdgeLength(hBCDE - sppB->getHeight());
                        sppC->setEdgeLength(h - sppC->getHeight());
                        sppDE->setEdgeLength(h - sppDE->getHeight());
                        sppCDE->setHeight(h);
                        sppCDE->setEdgeLength(hBCDE - h);
                        sppBCDE->setHeight(hBCDE);
                        sppBCDE->setEdgeLength(sppABCDE->getHeight() - hBCDE);
                        gf.calcPartialArray(sppCDE);
                        gf.calcPartialArray(sppBCDE);
                        gf.calcPartialArray(sppABCDE);
                        logL = gf.calcLogLikelihood();
                        wrong_tree.push_back(make_pair(h, logL));
                    }
                    else {
                        sppBC->setHeight(h);
                        sppBC->setEdgeLength(sppBCDE->getHeight() - h);
                        sppB->setEdgeLength(h - sppB->getHeight());
                        sppC->setEdgeLength(h - sppC->getHeight());
                        gf.calcPartialArray(sppBC);
                        gf.calcPartialArray(sppBCDE);
                        gf.calcPartialArray(sppABCDE);
                        logL = gf.calcLogLikelihood();
                        true_tree.push_back(make_pair(h, logL));
                    }
                }
                string newick = gf.makeNewick(9, true, false);
                output(format("%20.9f %12.5f %s\n") % logL % h % newick, 1);
            }
            
            // Create commands for plotting in R
            double ymin = G::_infinity;
            double ymax = G::_negative_infinity;
            double best_logL   = G::_negative_infinity;
            double best_height = 0.0;
            if (forest) {
                if (smctree) {
                    bool first = true;
                    output("height_wrong_forest <- c(", 1);
                    for (auto & p : wrong_forest) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        output(format("%.5f") % p.first, 1);
                    }
                    output(")\n", 1);
                    first = true;
                    output("logL_wrong_forest <- c(", 1);
                    for (auto & p : wrong_forest) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        if (p.second < ymin)
                            ymin = p.second;
                        if (p.second > ymax)
                            ymax = p.second;
                        if (p.second > best_logL) {
                            best_height = p.first;
                            best_logL = p.second;
                        }
                        output(format("%.5f") % p.second, 1);
                    }
                    output(")\n", 1);
                    output(format("ymin <- %.5f\n") % ymin, 1);
                    output(format("ymax <- %.5f\n") % ymax, 1);
                    output(format("mle_wrong_forest  <- %.5f\n") % best_height, 1);
                    output(format("mle_logL_wrong_forest  <- %.5f\n") % best_logL, 1);
                    output(format("plot(height_wrong_forest, logL_wrong_forest, type=\"l\", lwd=2, col=\"navy\", ylim=c(%.5f,%.5f))\n") % ymin % ymax, 1);
                    output("#lines(height_wrong_forest, logL_wrong_forest, lwd=2, col=\"red\")\n", 1);
                    output("abline(v=mle_wrong_forest, lwd=2, lty=\"dotted\", col=\"gray\")\n", 1);
                }
                else {
                    bool first = true;
                    output("height_true_forest <- c(", 1);
                    for (auto & p : true_forest) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        output(format("%.5f") % p.first, 1);
                    }
                    output(")\n", 1);
                    first = true;
                    double ymin = G::_infinity;
                    double ymax = G::_negative_infinity;
                    output("logL_true_forest <- c(", 1);
                    for (auto & p : true_forest) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        if (p.second < ymin)
                            ymin = p.second;
                        if (p.second > ymax)
                            ymax = p.second;
                        if (p.second > best_logL) {
                            best_height = p.first;
                            best_logL = p.second;
                        }
                        output(format("%.5f") % p.second, 1);
                    }
                    output(")\n", 1);
                    output(format("ymin <- %.5f\n") % ymin, 1);
                    output(format("ymax <- %.5f\n") % ymax, 1);
                    output(format("mle_true_forest  <- %.5f\n") % best_height, 1);
                    output(format("mle_logL_true_forest  <- %.5f\n") % best_logL, 1);
                    output(format("plot(height_true_forest, logL_true_forest, type=\"l\", lwd=2, col=\"navy\", ylim=c(%.5f,%.5f))\n") % ymin % ymax, 1);
                    output("#lines(height_true_forest, logL_true_forest, lwd=2, col=\"red\")\n", 1);
                    output("abline(v=mle_true_forest, lwd=2, lty=\"dotted\", col=\"gray\")\n", 1);
                }
            }
            else /* not forest */ {
                if (smctree) {
                    bool first = true;
                    output("height_wrong_tree <- c(", 1);
                    for (auto & p : wrong_tree) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        output(format("%.5f") % p.first, 1);
                    }
                    output(")\n", 1);
                    first = true;
                    output("logL_wrong_tree <- c(", 1);
                    for (auto & p : wrong_tree) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        if (p.second < ymin)
                            ymin = p.second;
                        if (p.second > ymax)
                            ymax = p.second;
                            if (p.second > best_logL) {
                                best_height = p.first;
                                best_logL = p.second;
                            }
                        output(format("%.5f") % p.second, 1);
                    }
                    output(")\n", 1);
                    output(format("ymin <- %.5f\n") % ymin, 1);
                    output(format("ymax <- %.5f\n") % ymax, 1);
                    output(format("mle_wrong_tree  <- %.5f\n") % best_height, 1);
                    output(format("mle_logL_wrong_tree  <- %.5f\n") % best_logL, 1);
                    output(format("plot(height_wrong_tree, logL_wrong_tree, type=\"l\", lwd=2, col=\"navy\", ylim=c(%.5f,%.5f))\n") % ymin % ymax, 1);
                    output("#lines(height_wrong_tree, logL_wrong_tree, lwd=2, col=\"red\")\n", 1);
                    output("abline(v=mle_wrong_tree, lwd=2, lty=\"dotted\", col=\"gray\")\n", 1);
                }
                else /* not smctree */ {
                    bool first = true;
                    output("height_true_tree <- c(", 1);
                    for (auto & p : true_tree) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        output(format("%.5f") % p.first, 1);
                    }
                    output(")\n", 1);
                    first = true;
                    output("logL_true_tree <- c(", 1);
                    for (auto & p : true_tree) {
                        if (first)
                            first = false;
                        else
                            output(",", 1);
                        if (p.second < ymin)
                            ymin = p.second;
                        if (p.second > ymax)
                            ymax = p.second;
                            if (p.second > best_logL) {
                                best_height = p.first;
                                best_logL = p.second;
                            }
                        output(format("%.5f") % p.second, 1);
                    }
                    output(")\n", 1);
                    output(format("ymin <- %.5f\n") % ymin, 1);
                    output(format("ymax <- %.5f\n") % ymax, 1);
                    output(format("mle_true_tree  <- %.5f\n") % best_height, 1);
                    output(format("mle_logL_true_tree  <- %.5f\n") % best_logL, 1);
                    output(format("plot(height_true_tree, logL_true_tree, type=\"l\", lwd=2, col=\"navy\", ylim=c(%.5f,%.5f))\n") % ymin % ymax, 1);
                    output("#lines(height_true_tree, logL_true_tree, lwd=2, col=\"red\")\n", 1);
                    output("abline(v=mle_true_tree, lwd=2, lty=\"dotted\", col=\"gray\")\n", 1);
                }
            }
        }   // end if find_MLEs else ...
    }
#endif

    inline void Proj::testSecondLevelSMC() {
        // Read data (only taxon names are used)
        readData();
        G::_ngenes = _data->getNumSubsets();
        assert(G::_ngenes > 0);

        // Copy taxon names to global variable _taxon_names
        G::_ntaxa = _data->getNumTaxa();
        _data->copyTaxonNames(G::_taxon_names);
        
        // Populate G::_species_names using taxon names from data file
        G::_species_names.clear();
        for (auto tname : G::_taxon_names) {
            string species_name = Node::taxonNameToSpeciesName(tname);
            if (find(G::_species_names.begin(), G::_species_names.end(), species_name) == G::_species_names.end()) {
                // Species_name not found in G::_species_names
                G::_species_names.push_back(species_name);
            }
        }

        // Read in the reference species tree
        G::_nexus_taxon_map.clear();
        map<unsigned,unsigned> species_taxon_map;
        vector<string> species_tree_names;
        vector<string> species_newicks;
        // if G::_species_names is empty, readTreeFile will populate it using the taxa block
        Forest::readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, species_taxon_map, species_tree_names, species_newicks);
        G::_nspecies = (unsigned)G::_species_names.size();
                
        // Read in the reference gene trees
        map<unsigned,unsigned> gene_taxon_map;
        vector<string> gene_tree_names;
        vector<string> gene_newicks;
        // if G::_taxon_names is empty, readTreeFile will populate it using the taxa block
        GeneForest::readTreefile(G::_gene_trees_ref_file_name, /*skip*/0, G::_taxon_names, gene_taxon_map, gene_tree_names, gene_newicks);
        G::_ntaxa = (unsigned)G::_taxon_names.size();
        
        // Populates G::_taxon_to_species
        // Assumes G::_taxon_names and G::_species_names are already populated
        buildSpeciesMap(/*taxa_from_data*/false);
        
        // Create a particle in which to store gene trees
        Particle particle;
                
        // Build the species tree
        SpeciesForest & sppref = particle.getSpeciesForest();
        G::_nexus_taxon_map = species_taxon_map;
        sppref.buildFromNewick(species_newicks[0]);

        // Build the gene trees
        G::_ngenes = (unsigned)gene_newicks.size();
        vector<GeneForest> & gfref = particle.getGeneForests();
        gfref.resize(gene_newicks.size());
        unsigned g = 0;
        G::_nexus_taxon_map = gene_taxon_map;
        for (auto newick : gene_newicks) {
            gfref[g].setParticle(&particle);
            gfref[g].setGeneIndex(g);
            gfref[g++].buildFromNewick(newick);
        }
        
        // Build coal info vectors
        sppref.buildCoalInfoVect();
        //vector<GeneForest> & gtvect = particle.getGeneForests();
        for (auto & gt : gfref) {
            gt.buildCoalInfoVect();
        }
        
        particle.setThetas();
        
        //temporary!
        // Calculate log-coalescent-likelihood for true species tree
        vector<Forest::coalinfo_t> coalinfo_vect;
        particle.recordAllForests(coalinfo_vect);
        double log_coallike = particle.calcLogCoalescentLikelihood(coalinfo_vect,  /*integrate_out_thetas*/false, /*verbose*/true);
        output(format("log-coallike = %.9f (true species tree)\n") % log_coallike, 2);
        
        // Calculate Jones 2017) log-coalescent-likelihood for true species tree
        //coalinfo_vect.clear();
        //particle.recordAllForests(coalinfo_vect);
        //double log_coallike_jones = particle.calcLogCoalescentLikelihood(coalinfo_vect,  /*integrate_out_thetas*/true, /*verbose*/true);
        //output(format("log-coallike-jones = %.9f (true species tree)\n") % log_coallike_jones, 2);
        
        //temporary!
        G::_speclog.clear();

        SMC smc2;
        smc2.setMode(SMC::SPECIES_GIVEN_GENE);
        smc2.setNParticles(G::_nparticles2);
        smc2.initFromParticle(particle);
        smc2.run();
        smc2.summarize();
        
        //temporary!
        // Save proposal distributions at each step for plotting in R
        unsigned nsteps = (unsigned)G::_speclog.size();
        vector<string> tmp;
        vector<string> ablines;
        vector<string> joincolors;
        vector<double> logweights;
        vector<double> increments;
        vector<double> maxheights;
        vector<double> forestheights;
        vector<double> upperbounds;
        typedef tuple<double, double, double, double, double, string> list_t;
        vector<list_t> list;
        map<set<G::species_t>, unsigned> precombos;
        map<set<G::species_t>, unsigned> postcombos;
        for (unsigned step = 0; step < nsteps; step++) {
            logweights.clear();
            increments.clear();
            maxheights.clear();
            forestheights.clear();
            upperbounds.clear();
            joincolors.clear();
            ablines.clear();
            precombos.clear();
            postcombos.clear();
            
            ofstream tmpf(str(format("step-%d.R") % step));
            
            tmpf << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
            tmpf << "setwd(cwd)\n";
            tmpf << str(format("pdf(\"plot-%d.pdf\")\n") % step);
            
            vector<G::SpecLog> & v = G::_speclog.at(step);
                        
            unsigned prefilter_count = 0;
            unsigned postfilter_count = 0;
            for (auto x : v) {
                if (x._filtered) {
                    postfilter_count += x._freq;
                }
                else {
                    prefilter_count += 1;
                }
            }
            tmpf << "# step = " << step << "\n";
            tmpf << "# pre-filtering  = " << prefilter_count << "\n";
            tmpf << "# post-filtering = " << postfilter_count << "\n";
            
            // Save log weights
            tmp.clear();
            for (auto x : v) {
                if (!x._filtered) {
                    tmp.push_back(str(format("%.9f") % x._logw));
                    logweights.push_back(x._logw);
                }
            }
            tmpf << "logw <- c(" << boost::algorithm::join(tmp, ",") << ")\n";
            
            // Save log increments
            tmp.clear();
            for (auto x : v) {
                if (!x._filtered) {
                    tmp.push_back(str(format("%.9f") % x._incr));
                    increments.push_back(x._incr);
                    joincolors.push_back(x.calcColor());
                    set<G::species_t> s = {x._left, x._right};
                    precombos[s] += 1;
                }
                else {
                    assert(x._freq > 0);
                    set<G::species_t> s = {x._left, x._right};
                    postcombos[s] += x._freq;
                    ablines.push_back(str(format("abline(v=%.9f, col=%s)\n") % x._incr % x.calcColor()));
                }
            }
            tmpf << "incr <- c(" << boost::algorithm::join(tmp, ",") << ")\n";
            tmpf << "cols <- c(" << boost::algorithm::join(joincolors, ",") << ")\n";
            
            // Save maxh
            tmp.clear();
            for (auto x : v) {
                if (!x._filtered) {
                    tmp.push_back(str(format("%.9f") % x._maxh));
                    maxheights.push_back(x._maxh);
                    forestheights.push_back(x._height);
                    upperbounds.push_back(x._maxh - x._height);
                }
            }
            tmpf << "maxh <- c(" << boost::algorithm::join(tmp, ",") << ")\n";
            
            tmpf << "plot(incr, logw, type=\"p\", pch=19, col=cols)\n";
            for (auto abline : ablines) {
                tmpf << abline;
            }
            
            tmpf << "\n# pre-combos:\n";
            for (auto c : precombos) {
                auto s = c.first;
                unsigned count = c.second;
                tmpf << "# " << count << ": ";
                for (auto z : s) {
                    tmpf << z << " ";
                }
                tmpf << "\n";
            }
            
            tmpf << "\n# post-combos:\n";
            for (auto c : postcombos) {
                auto s = c.first;
                unsigned count = c.second;
                tmpf << "# " << count << ": ";
                for (auto z : s) {
                    tmpf << z << " ";
                }
                tmpf << "\n";
            }
            
            tmpf << "\n# (logw, incr, maxh, height, bound, color) from highest to lowest log-weight:\n";
            unsigned n = (unsigned)logweights.size();
            assert(n == increments.size());
            assert(n == maxheights.size());
            assert(n == forestheights.size());
            assert(n == upperbounds.size());
            assert(n == joincolors.size());
            list.clear();
            for (unsigned i = 0; i < n; i++) {
                list.push_back(make_tuple(logweights[i], increments[i], maxheights[i], forestheights[i], upperbounds[i], joincolors[i]));
            }
            sort(list.begin(), list.end(), greater<list_t>());
            for (auto l : list) {
                tmpf << "# " << get<0>(l) << ", " << get<1>(l) << ", " << get<2>(l) << ", " << get<3>(l) << ", " << get<4>(l) << ", " << get<5>(l) << ", " <<"\n";
            }

            tmpf << "dev.off()\n";
            tmpf.close();
        }
    }

    inline void Proj::chibSim(/*const Particle & test_particle*/) {
        // Begin by simulating data, species trees, gene trees
        simulateData(/*chib*/true);

        unsigned nreps = _chib_nreps;
        unsigned nsame = 0;
        Particle particle;
        for (unsigned i = 0; i < nreps; i++) {
            // if nreps = 1000,
            // i = 0: seed between 1001 and 2000
            // i = 1: seed between 2001 and 3000
            // ...
            // i = 999: seed between 1000001 and 1001000
            //POLWAS particle.setSeed((i+1)*nreps + rng.randint(1,nreps));
            particle.setSeed((i+1)*nreps + rng->randint(1,nreps));
            simulateTrees(particle);
            SpeciesForest      & sf = particle.getSpeciesForest();
            vector<GeneForest> & gfvect   = particle.getGeneForests();
            double rfsum = Forest::calcTreeDistances(*_species_tree_ref, sf).second;
            assert(gfvect.size() == _gene_tree_refs.size());
            for (unsigned g = 0; g < gfvect.size(); g++) {
                GeneForest::SharedPtr gfref = _gene_tree_refs[g];
                rfsum += Forest::calcTreeDistances(*gfref, gfvect[g]).second;
            }
            if (rfsum == 0.0) {
                nsame++;
            }
        }

        double logp = G::_negative_infinity;
        assert(nreps > 0);
        if (nsame > 0) {
            logp = log(nsame) - log(nreps);
        }
        output(format("Log(prob. topol.lo) = %.9f") % logp,1);
    }
    
    inline void Proj::selectParticlesToKeep(list<Particle> & first_level_particles, vector<unsigned> & kept) {
        assert(kept.size() == first_level_particles.size());
        assert(G::_nkept <= G::_nparticles);
        if (G::_nkept == G::_nparticles) {
            // No need to select, just keep every particle
            unsigned which = 0;
            for (auto & p : first_level_particles) {
                unsigned n = p.getCount();
                kept[which++] = n;
            }
        }
        else {
            double n = (double)G::_nparticles;
            
            // Choose _nkept particles
            vector<double> prob;
            for (auto & p : first_level_particles) {
                prob.push_back(1.0*p.getCount()/n);
            }
            double sum_probs = accumulate(prob.begin(), prob.end(), 0.0);
            assert(fabs(sum_probs - 1.0) < 0.0001);
            
            for (unsigned k = 0; k < G::_nkept; ++k) {
                unsigned which = G::multinomialDraw(rng, prob);
                kept[which]++;
            }
        }
    }

    inline void Proj::run() {
    
        output("Starting...\n", 2);

        if (G::_nthreads > 1) {
            output(format("Multithreaded version with %d threads\n") % G::_nthreads, 2);
        }
        else {
            output("Singlethreaded version\n", 2);
        }

        output(format("Current working directory: %s\n") % current_path(), 2);
        
        try {
            //POLWAS rng.setSeed(_rnseed);
            rng->setSeed(_rnseed);
            
#if defined(POLTMP)
            geneTreeExperiment(/*forest*/true, /*smctree*/false);
#else
            if (_start_mode == "sim") {
                simulateData(/*chib*/false);
            }
            else if (_start_mode == "chib") {
                // Estimate prior probability of choosing specified species
                // and gene tree topologies for purposes of using Chib's method
                chibSim();
            }
            else if (_start_mode == "spec") {
                // Estimate coalescent log likelihood for specified species tree
                // and gene trees
                calcCoalLikeForSpecified();
            }
            else if (_start_mode == "2ndlevel") {
                if (G::_theta_mean_fixed < 0.0) {
                    throw XProj("Must specify fixedthetamean if startmode is 2ndlevel");
                }
        
                // Test 2nd-level SMC using specified gene trees
                testSecondLevelSMC();
            }
            else {
                if (_start_mode != "smc")
                    throw XProj("startmode must be either \"sim\", \"smc\", or \"chib\"");

                readData();
                debugShowStringVector("Gene names", G::_gene_names);
                
                G::_ngenes = _data->getNumSubsets();
                assert(G::_ngenes > 0);

#if defined(USING_MPI)
                mpiSetSchedule();
#endif
                    
                // Set relative rates if specified
                setRelativeRates();

                // Copy taxon names to global variable _taxon_names
                G::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(G::_taxon_names);
                
                // Save species names to global variable _species_names
                // and create global _taxon_to_species map that provides
                // the species index for each taxon name
                G::_nspecies = buildSpeciesMap(/*taxa_from_data*/true);
                
                // Create global _species_mask that has "on" bits only for
                // the least significant _nspecies bits
                Node::setSpeciesMask(G::_species_mask, G::_nspecies);
                
                // Show species names if debugging
                debugShowStringVector("Species names", G::_species_names);

                // Precompute leaf partials (never have to be computed again)
                for (unsigned g = 0; g < G::_ngenes; g++) {
                    GeneForest::computeLeafPartials(g, _data);
                }
                
                // First-level particle filtering
                SMC smc;
                smc.setMode(SMC::SPECIES_AND_GENE);
                smc.setNParticles(G::_nparticles);
                smc.setData(_data);
                smc.init();
                smc.run();
                smc.summarize();
                
                if (G::_nparticles2 > 0) {
                    output(format("Performing 2nd-level SMC on %d 1st-level particles...\n") % G::_nkept, 2);
                    
                    list<Particle> & first_level_particles = smc.getParticles();

                    // Choose G::_nkept 1st-level particles for use in 2nd level
                    vector<unsigned> kept(first_level_particles.size(), 0);
                    selectParticlesToKeep(first_level_particles, kept);
                    
                    // Second-level particle filtering
                    SMC ensemble;
                    ensemble.setMode(SMC::SPECIES_GIVEN_GENE);
                    unsigned which_first_level = 0;
                    unsigned total_first_level = 0;
                    unsigned which_second_level = 0;
                    for (auto & p : first_level_particles) {
                        unsigned n = kept[which_first_level]; //p.getCount();
                                                
                        if (n > 0) {
                            // Rebuild coal info vectors, stripping effect of previous species tree
                            vector<GeneForest> & gtvect = p.getGeneForests();
                            for (auto & gt : gtvect) {
                                gt.buildCoalInfoVect();
                            }
                                                    
                            for (unsigned j = 0; j < n; j++) {
                            
                                ++total_first_level;
                                if (total_first_level % 10 == 0)
                                    output(format("  working on 1st-level %d...\n") % total_first_level, 2);
                                    
                                SMC smc2;
                                smc2.setMode(SMC::SPECIES_GIVEN_GENE);
                                smc2.setNParticles(G::_nparticles2);
                                smc2.initFromParticle(p);
                                smc2.run();
                                smc2.dumpParticles(ensemble);
                                //debugCheckEnsemble(ensemble);
                                which_second_level++;
                            }
                        }
                        which_first_level++;
                    }
                    //debugCheckEnsemble(ensemble);
                    ensemble.summarize();
                }
            }
#endif      //POLTMP
        }
        catch (XProj & x) {
            output(format("Proj encountered a problem:\n  %s\n") % x.what(), 2);
        }
        
        output("\nFinished!\n", 2);
    }
    
    void Proj::debugCheckEnsemble(const SMC & ensemble) const {
        // Check to make sure all particles in ensemble have gene trees
        output("Checking ensemble:\n",2);
        const list<Particle> & plist = ensemble.getParticlesConst();
        unsigned which = 0;
        for (auto & p : plist) {
            which++;
            unsigned c = p.getCount();
            const vector<GeneForest> * gf = p.getGeneTreesConst();
            if (gf) {
                unsigned n = (unsigned)(gf->size());
                output(format("  particle %d (count %d) has %d gene trees\n") % which % c % n,2);
            }
            else {
                output(format("  particle %d (count = %d) has NO gene trees\n") % which % c,2);
            }
        }
    }
    
    void Proj::memoryReport(ofstream & memf) const {
        memf << "\nProj memory report:\n\n";
        memf << str(format("  Size of int:           %d\n") % sizeof(int));
        memf << str(format("  Size of char:          %d\n") % sizeof(char));
        memf << str(format("  Size of double:        %d\n") % sizeof(double));
        memf << str(format("  Size of unsigned:      %d\n") % sizeof(unsigned));
        memf << str(format("  Size of unsigned long: %d\n") % sizeof(unsigned long));
        memf << str(format("  Size of Node *:        %d\n") % sizeof(Node *));
        memf << str(format("  Number of particles: %d\n") % G::_nparticles);
        _data->memoryReport(memf);
    }
    
#if defined(USING_MPI)
    inline void Proj::mpiSetSchedule() {
        // Determine which particles will be handled by this processor: e.g.
        // 1000 = number of particles
        //  3   = number of processors
        //  333 = 1000 / 3
        //  1   = 1000 % 3
        //  rank 0 gets 333, rank 1 gets 334, rank 2 gets 333
        vector<unsigned> particles_per_task(ntasks, (unsigned)(G::_nparticles / ntasks));

        unsigned remainder = G::_nparticles % ntasks;

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
