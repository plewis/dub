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

extern proj::PartialStore ps;
extern proj::StopWatch stopwatch;
extern proj::Lot::SharedPtr rng;

namespace proj {

    class Proj {
        public:
                                            Proj();
                                            ~Proj();
   
            void                            clear();
            void                            processCommandLineOptions(int argc, const char * argv[]);
            void                            secondLevelRange(vector<Particle> & first_level_particles, SMC & ensemble, unsigned nparticles, unsigned nkept, unsigned first, unsigned last);
            void                            secondLevelConditionedOn(unsigned rnseed, const Particle & p, unsigned nparticles, unsigned nkept, SMC & ensemble);
        void                                buildEnsembleCoalInfo(vector<Particle> & first_level_particles);
            void                            stripPartials(vector<Particle> & first_level_particles);
            void                            run();
            void                            readData();
            void                            setRelativeRates();
            void                            parseRelRateDefinition(string & s);
            void                            simulateData();
            void                            outputGeneTreesToFile(string fn, const vector<string> & newicks) const;
            void                            outputNexusTreefile(string fn, const vector<string> & newicks) const;
            void                            outputJavascriptTreefile(string fn, const string & newick_species_tree_numeric,
                                                const vector<string> & newick_gene_trees_numeric);
#if defined(LAZY_COPYING)
            void                            outputTrees(SpeciesForest & sf, vector<GeneForest::SharedPtr> & gfpvect);
#else
            void                            outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect);
#endif
            void                            simulateTrees(Particle & particle);
            void                            reportDeepCoalescences(Particle & particle);
            
#if defined(FOSSILS)
            Fossil                          parseFossilDefinition(string & fossil_def);
            TaxSet                          parseTaxsetDefinition(string & taxset_def);
#endif
    
            void                            selectParticlesToKeep(list<Particle> & first_level_particles, vector<unsigned> & kept);
               
            string                          _data_file_prefix;
            string                          _start_mode;
            Partition::SharedPtr            _partition;
            Data::SharedPtr                 _data;
            map<string, double>             _relrate_map;
            
            SpeciesForest::SharedPtr        _species_tree_ref;
            vector<GeneForest::SharedPtr>   _gene_tree_refs;
            
            unsigned                        _nsimspecies;
            vector<unsigned>                _nsimtaxaperspecies;

            bool                            _use_gpu;
            bool                            _ambig_missing;
            
            static string                   _program_name;
            static unsigned                 _major_version;
            static unsigned                 _minor_version;
            
        protected:
            
            unsigned                        buildSpeciesMap(bool taxa_from_data);
    };

    inline Proj::Proj() {
        clear();
    }

    inline Proj::~Proj() {
    }
    
    inline void Proj::clear() {
        _use_gpu                    = true;
        _ambig_missing              = true;
        _data                       = nullptr;
        _data_file_prefix           = "";
        _species_tree_ref           = nullptr;
        _start_mode                 = "smc";
        _gene_tree_refs.clear();
        _partition.reset(new Partition());
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> log_include;
        vector<string> partition_subsets;
#if defined(SPECIES_IN_CONF)
        vector<string> species_definitions;
#endif
        vector<string> partition_relrates;
#if defined(FOSSILS)
        vector<string> taxsets;
        vector<string> fossils;
#endif
        bool do_compress;
        //bool dummy_bool;
        int dummy_int;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("log", value(&log_include), "categories of output to include: INFO = include normal output; VERBOSE = include debugging output; SECONDLEVEL = show progress in second level SMC analyses; DEBUGGING = show debugging output")
        ("datafile",  value(&_data_file_prefix), "prefix to be used in creating data file names (extension, e.g. '.nex', will be appended)")
        ("speciestreeref",  value(&G::_species_tree_ref_file_name), "name of a tree file containing a single reference species tree")
        ("genetreeref",  value(&G::_gene_trees_ref_file_name), "name of a tree file containing a reference gene tree for each locus")
        ("startmode", value(&_start_mode), "if 'sim', simulate gene trees, species tree, and data; if 'smc', estimate from supplied datafile; if 'chib', computes prior probability of species species and gene tree topologies; if 'spec', computes coalescent likelihood for specified speciestreeref and genetreeref; if '2ndlevel', tests second-level SMC from gene trees supplied by genetreeref")
#if defined(SPECIES_IN_CONF)
        ("species", value(&species_definitions), "a string defining a species, e.g. 'A:x,y,z' says that taxa x, y, and z are in species A")
#endif
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
#if defined(FOSSILS)
        ("fossil",  value(&fossils), "a string defining a fossil, e.g. 'Ursus_abstrusus         1.8–5.3 4.3' (4.3 is time, 1.8-5.3 is prior range)")
        ("taxset",  value(&taxsets), "a string defining a taxon set, e.g. 'Ursinae: Helarctos_malayanus Melursus_ursinus Ursus_abstrusus Ursus_americanus Ursus_arctos Ursus_maritimus Ursus_spelaeus Ursus_thibetanus'")
#endif
        ("relrate",  value(&partition_relrates), "a relative rate for a previously-defined subset; format first:3.1")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("nparticles",  value(&G::_nparticles)->default_value(500), "number of particles in a population for joint estimation (1st level SMC)")
        ("nsubpops",  value(&G::_nsubpops)->default_value(1), "number of subpopulations in 1st level SMC (must divide evenly into nparticles)")
        ("nkept",  value(&G::_nkept)->default_value(0), "number of particles from joint smc kept for conditional (2nd level) smc")
        ("nspeciesparticles",  value(&G::_nparticles2)->default_value(0), "number of particles in a population for (2nd level) species tree only estimation")
        ("nspecieskept",  value(&G::_nkept2)->default_value(0), "number of particles from conditional smc kept")
        ("treefilecompression", value(&do_compress)->default_value(0), "no=no compression, yes=maximum compression")
        ("lambda",  value(&G::_lambda)->default_value(10.0), "per lineage speciation rate assumed for the species tree")
        ("simnspecies",  value(&_nsimspecies)->default_value(1), "number of species (only used if startmode is 'sim')")
        ("simntaxaperspecies",  value(&_nsimtaxaperspecies), "number of taxa sampled per species (only used if startmode is 'sim'); should be _nimspecies of these entries, one for each species simulated")
        ("simedgeratevar",  value(&G::_edge_rate_variance)->default_value(0.0), "variance of lognormal relative rate distribution across edges in gene trees (used only if startmode is 'sim')")
        ("simasrvshape",  value(&G::_asrv_shape)->default_value(G::_infinity), "Shape of gamma among-site rate heterogeneity within a locus (used only if startmode is 'sim')")
        ("simoccupancy",  value(&G::_occupancy)->default_value(1.0), "probability that any given taxon will have data for any given locus; 1-_occupancy is prob. all missing data for a taxon (used only if startmode is 'sim')")
        ("simcomphet",  value(&G::_comphet)->default_value(G::_infinity), "Dirichlet parameter governing compositional heterogeneity (default value results in compositional homogeneity (used only if startmode is 'sim')")
        ("rnseed",  value(&G::_rnseed)->default_value(13579), "pseudorandom number seed")
        ("theta",  value(&G::_theta)->default_value(0.05), "coalescent parameter assumed for gene trees")
#if defined(USE_HEATING)
        ("heatingpower",  value(&G::_heating_power)->default_value(1.0), "power to which weights are raised (default 1.0 means no heating)")
#endif
#if defined(USING_MULTITHREADING)
        ("nthreads",  value(&G::_nthreads)->default_value(1), "specify number of threads")
#endif
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

        if (vm.count("help") > 0) {
            output(format("%s\n") % desc, G::LogCateg::ALWAYS);
            exit(1);
        }

        if (vm.count("version") > 0) {
            output(format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version);
            exit(1);
        }
        
#if defined(USING_MULTITHREADING)
        if (vm.count("nthreads") > 0) {
            if (G::_nthreads < 1) {
                output(format("Number of threads specified cannot be less than 1 (you specified %d)\n") % G::_nthreads);
                exit(1);
            }
        }
#endif

        if (vm.count("nspeciesparticles") > 0 || vm.count("nspecieskept") > 0) {
            if (G::_nparticles2 > 0 && G::_nkept2 == 0) {
                output(format("You specified nspeciesparticles = %d but nspecieskept is still set to default value 0\n") % G::_nparticles2);
                exit(1);
            }
            else if (G::_nparticles2 == 0 && G::_nkept2 > 0) {
                output(format("You specified nspecieskept = %d but nspeciesparticles is still set to default value 0\n") % G::_nkept2);
                exit(1);
            }
        }
        
        // If user specified --relrate on command line, break specified relrate
        // definition into name and rate
        if (vm.count("relrate") > 0) {
            for (auto relrate_definition : partition_relrates) {
                parseRelRateDefinition(relrate_definition);
            }
        }
        
        // If user specified --treefilecompression on command line, set G::_treefile_compression to 2
        if (vm.count("treefilecompression") > 0) {
            G::_treefile_compression = (do_compress ? 2 : 0);
        }
        
        // If user specified --nsubpops on command line, check to
        // ensure that nparticles is a multiple of nsubpops
        if (vm.count("nsubpops") > 0) {
            if (G::_nparticles % G::_nsubpops != 0)
                throw XProj(str(format("nsubpops (%d) must divide evenly into nparticles (%d)") % G::_nsubpops % G::_nparticles));
        }
        
#if defined(USE_HEATING)
        // If user specified --heatingpower on command line, check to
        // ensure that heatingpower is between 0.0 and 1.0 (inclusive)
        if (vm.count("heatingpower") > 0) {
            if (G::_heating_power < 0.0 || G::_heating_power > 1.0)
                throw XProj(str(format("heatingpower should be in range [0.0, 1.0] but you specified (%g)") % G::_heating_power));
        }
#endif
        
        // If user specified --log on command line, set G::_log_includes
        // according to the values specified
        if (vm.count("log") > 0) {
            G::_log_include = G::LogCateg::NONE;
            for (auto s : log_include) {
                if (s == "INFO" || s == "info" || s == "Info")
                    G::_log_include |= G::LogCateg::INFO;
                else if (s == "VERBOSE" || s == "verbose" || s == "Verbose")
                    G::_log_include |= G::LogCateg::VERBOSE;
                else if (s == "SECONDLEVEL" || s == "secondlevel" || s == "SecondLevel")
                    G::_log_include |= G::LogCateg::SECONDLEVEL;
                else if (s == "DEBUGGING" || s == "debugging" || s == "Debugging")
                    G::_log_include |= G::LogCateg::DEBUGGING;
                else if (s == "CONDITIONALS" || s == "conditionals" || s == "Conditionals")
                    G::_log_include |= G::LogCateg::CONDITIONALS;
                else
                    throw XProj(str(format("unrecognized log category \"%s\"") % s));
            }
        }
        
        // If user specified --subset on command line, break specified partition subset
        // definition into name and character set string and add to _partition
        if (vm.count("subset") > 0) {
            _partition.reset(new Partition());
            for (auto s : partition_subsets) {
                _partition->parseSubsetDefinition(s);
            }
        }
        
#if defined(SPECIES_IN_CONF)
        if (vm.count("species") > 0) {
            for (auto s : species_definitions) {
                G::parseSpeciesDefinition(s);
            }
        }
#endif
        
#if defined(FOSSILS)
        // If user specified --fossil on command line, break specified
        // fossil definition into species name, taxset name, and age
        if (vm.count("fossil") > 0) {
            G::_fossils.clear();
            for (auto fdef : fossils) {
                G::_fossils.push_back(parseFossilDefinition(fdef));
            }
        }
        
        // If user specified --taxset on command line, break specified
        // taxset definition into name and species included
        if (vm.count("taxset") > 0) {
            G::_taxsets.clear();
            for (auto tdef : taxsets) {
                G::_taxsets.push_back(parseTaxsetDefinition(tdef));
            }
        }
#endif
        
        // If user specified --nkept on command line, check to ensure
        // that nkept <= nparticles.
        if (vm.count("nkept") > 0) {
            if (G::_nkept > G::_nparticles) {
                throw XProj("nkept cannot be greater than nparticles");
            }
        }
        
        // If user specified --nspecieskept on command line, check to ensure
        // that nkept2 <= nspeciesparticles.
        if (vm.count("nspecieskept") > 0) {
            if (G::_nkept2 > G::_nparticles2) {
                throw XProj("nspecieskept cannot be greater than nspeciesparticles");
            }
        }
    }

#if defined(FOSSILS)
    // Note: This code is unduly complex because of the fact that there are
    // many different kinds of dashes that can be used in a utf8-encoded text file.
    // One solution is to eliminate the dash separating the lower and upper
    // bounds for the age range, but I felt that using a dash makes it clear that
    // it is an age range and I didn't want to put the onus on the user to use one
    // particular kind of dash (especially since I could not figure out how to
    // use the right kind of dash myself in a text file created using BBEdit).
    
    // The ws_to_utf8 and utf8_to_ws functions below are
    // slightly modified from Galik's answer at
    // stackoverflow.com/questions/43302279
    //   /any-good-solutions-for-c-string-code-point-and-code-unit
    //   /43302460#43302460
    string ws_to_utf8(wstring const & s) {
        wstring_convert<codecvt_utf8<wchar_t>, wchar_t> cnv;
        string utf8 = cnv.to_bytes(s);
        if(cnv.converted() < s.size())
            throw XProj("incomplete conversion to utf8");
        return utf8;
    }

    wstring utf8_to_ws(string const & utf8) {
        wstring_convert<codecvt_utf8<wchar_t>, wchar_t> cnv;
        wstring s = cnv.from_bytes(utf8);
        if(cnv.converted() < utf8.size())
            throw XProj("incomplete conversion to wstring");
        return s;
    }

    inline Fossil Proj::parseFossilDefinition(string & fossil_def) {
        // Examples showing range of possible inputs:
        //   fossil_def = "Ursus_abstrusus 1.8–5.3 4.3"
        //   fossil_def = "Parictis_montanus      33.90  – 37.20  36.6"

        // Vector to hold strings obtained by parsing fossil_def
        vector<string> v;
        
        // Separate fossil_def into 4 strings:
        //  v[0] = fossil species name
        //  v[1] = fossil age lower bound
        //  v[2] = fossil age upper bound
        //  v[3] = fossil age

        // regex_pattern specifies string of characters that do not
        // represent a tab (\u0009), space (\u0020), or dash of any kind:
        //  \u002D Hyphen-minus
        //  \u058A Armenian Hyphen
        //  \u05BE Hebrew Punctuation Maqaf
        //  \u2010 Hyphen
        //  \u2011 Non-Breaking Hyphen
        //  \u2012 Figure Dash
        //  \u2013 En dash
        //  \u2014 Em dash
        //  \u2015 Horizontal bar
        //  \u2E3A Two-Em Dash
        //  \u2E3B Three-Em Dash
        //  \uFE58 Small Em Dash
        //  \uFE63 Small Hyphen-Minus
        //  \uFF0D Fullwidth Hyphen-Minus
        wregex regex_pattern(utf8_to_ws("[^\\u0020\\u0009\\u002D\\u058A\\u05BE\\u2010\\u2011\\u2012\\u2013\\u2014\\u2015\\u2E3A\\u2E3B\\uFE58\\uFE63\\uFF0D]+"));

        // 0 means keep whole match:
        //     1 would mean keep first match only
        //     2 would mean keep second match only
        //     {1,2} would mean keep first and secod matches only
        //    -1 means use regex_pattern as delimiter
        wstring ws_fossil_def = utf8_to_ws(fossil_def);
        wsregex_token_iterator regex_iter(ws_fossil_def.begin(), ws_fossil_def.end(), regex_pattern, 0);

        // regex_end is, by default, signifies the end of the input
        wsregex_token_iterator regex_end;

        // iterate to get matches
        for ( ; regex_iter != regex_end; ++regex_iter) {
            wstring s = *regex_iter;
            v.push_back(ws_to_utf8(s));
        }

        string fossil_species_name = v[0];
        string fossil_age_lower_str  = v[1];
        string fossil_age_upper_str  = v[2];
        string fossil_age_str   = v[3];
        
        double fossil_age_lower;
        try {
            fossil_age_lower = stof(fossil_age_lower_str);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert fossil age lower bound \"%s\" to a floating point value") % fossil_age_lower_str);
        }

        double fossil_age_upper;
        try {
            fossil_age_upper = stof(fossil_age_upper_str);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert fossil age upper bound \"%s\" to a floating point value") % fossil_age_upper_str);
        }
        
        double fossil_age;
        try {
            fossil_age = stof(fossil_age_str);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert fossil age \"%s\" to a floating point value") % fossil_age_str);
        }
        
        return Fossil(fossil_species_name, fossil_age_lower, fossil_age_upper, fossil_age);
    }
    
    inline TaxSet Proj::parseTaxsetDefinition(string & taxset_def) {
        // Example:
        //   taxset_def = "Ursinae = Helarctos_malayanus Melursus_ursinus Ursus_abstrusus Ursus_americanus Ursus_arctos Ursus_maritimus Ursus_spelaeus Ursus_thibetanus;"
        vector<string> v;
        
        // Separate taxset_def into 2 strings at equal sign:
        //  v[0] = taxset name
        //  v[1] = list of species in taxset
        split(v, taxset_def, boost::is_any_of(":"));
        if (v.size() != 2)
            throw XProj(format("Expecting exactly 2 items separated by colon (:) in taxset definition but instead found %d") % v.size());
        string taxset_name = v[0];
        string taxset_species_list  = v[1];
        
        // Trim whitespace from both ends
        trim(taxset_name);
        trim(taxset_species_list);
        
        // Separate taxset_species_list into strings separated by spaces
        
        // regex_pattern specifies string of characters not a hyphen or whitespace
        regex regex_pattern("[^\\s]+");
        
        // 0 means keep whole match:
        //     1 would mean keep first match only
        //     2 would mean keep second match only
        //     {1,2} would mean keep first and secod matches only
        //    -1 means use regex_pattern as delimiter
        sregex_token_iterator regex_iter(taxset_species_list.begin(), taxset_species_list.end(), regex_pattern, 0);
        
        // regex_end is, by default, signifies the end of the input
        sregex_token_iterator regex_end;
        
        v.clear();
        for ( ; regex_iter != regex_end; ++regex_iter) {
            v.push_back(*regex_iter);
        }
        return TaxSet(taxset_name, v);
    }
    
#endif

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
        output(format("\nReading and storing the data in the file %s\n") % (_data_file_prefix + ".nex"));
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_prefix + ".nex");
        
        G::_gene_names.clear();

        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        output(format("\nNumber of taxa: %d\n") % _data->getNumTaxa());
        output(format("Number of partition subsets: %d") % nsubsets);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNLoci(nsubsets);
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Set length of partials for gene g
            ps.setNElements(G::_nstates*_data->getNumPatternsInSubset(subset), subset);
            
            DataType dt = _partition->getDataTypeForSubset(subset);
            G::_gene_names.push_back(_data->getSubsetName(subset));
            output(format("  Subset %d (%s)\n") % (subset+1) % _data->getSubsetName(subset));
            output(format("    data type: %s\n") % dt.getDataTypeAsString());
            output(format("    sites:     %d\n") % _data->calcSeqLenInSubset(subset));
            output(format("    patterns:  %d\n") % _data->getNumPatternsInSubset(subset));
        }
    }
             
    inline void Proj::setRelativeRates() {
        bool relrates_specified = !_relrate_map.empty();
        if (relrates_specified) {
            output("\n  Relative rates specified for each locus:\n");
            output(format("  %12s %12s\n") % "locus" % "rate");
            output(format("  %12s %12s\n") % "-----------" % "-----------");
            
            G::_relrate_for_gene.clear();
            unsigned total_nsites = _partition->getNumSites();
            double mean_rate = 0.0;
            for (unsigned g = 0; g < G::_nloci; g++) {
                unsigned gnsites = _partition->numSitesInSubset(g);
                string gname = _partition->getSubsetName(g);
                double r = 0.0;
                if (_relrate_map.count(gname) == 0)
                    throw XProj(format("Proj::setRelativeRates failed because key \"%s\" does not exist in _relrate_map") % gname);
                else {
                    r = _relrate_map.at(gname);
                }
                G::_relrate_for_gene[g] = r;
                output(format("  %12s %12.5f\n") % gname % r);
                mean_rate += r*gnsites/total_nsites;
            }
            output(format("  %12s %12s\n") % "-----------" % "-----------");
            output(format("  %12s %12.5f\n") % "mean" % mean_rate);
            if (fabs(mean_rate - 1.0) > 0.001) {
                XProj("The mean rate is more than 0.001 away from 1.0");
            }
        }
        else {
            output("\n  Relative rates not specified:\n    assuming rate homogeneity across loci\n");
            for (unsigned g = 0; g < G::_nloci; g++) {
                G::_relrate_for_gene[g] = 1.0;
            }
        }
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
    
    inline void Proj::outputGeneTreesToFile(string fn, const vector<string> & newicks) const {
        assert(G::_nloci == newicks.size());
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
    
#if defined(LAZY_COPYING)
    inline void Proj::outputTrees(SpeciesForest & sf, vector<GeneForest::SharedPtr> & gfpvect) {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;
        
        // Output tree file containing true species tree
        string newick_species_tree_alpha = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        string newick_species_tree_numeric = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
        if (G::_nspecies > 1) {
            outputNexusTreefile("true-species-tree.tre", {newick_species_tree_alpha});
            output(format("  True species tree height = %.9f\n") % sf.getHeight());
            output("  True species tree saved in file \"true-species-tree.tre\"\n");
        }
        else {
            output("  True species tree not saved because it is just a single node!\n");
        }

        // Output tree file containing true gene trees
        vector<string> newick_gene_trees_alpha;
        vector<string> newick_gene_trees_numeric;
        //output("  True gene tree height (height/theta):\n");
        for (auto gfp : gfpvect) {
            string newick_alpha = gfp->makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
            newick_gene_trees_alpha.push_back(newick_alpha);
            
            string newick_numeric = gfp->makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
            newick_gene_trees_numeric.push_back(newick_numeric);
            //output(format("  %12.9f (%.9f)\n") % gf.getHeight() % (gfp->getHeight()/G::_theta));
        }
        outputGeneTreesToFile("true-gene-trees.tre", newick_gene_trees_alpha);
        output("  True gene trees saved in file \"true-gene-trees.tre\"\n");

        // Output gene trees and species trees for javascript viewer
        outputJavascriptTreefile("newicks.js", newick_species_tree_numeric, newick_gene_trees_numeric);
    }
#else
    inline void Proj::outputTrees(SpeciesForest & sf, vector<GeneForest> & gfvect) {
        // This should be a setting
        bool edgelens_in_coalescent_units = false;
        
        // Output tree file containing true species tree
        string newick_species_tree_alpha = sf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
        string newick_species_tree_numeric = sf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
        if (G::_nspecies > 1) {
            outputNexusTreefile("true-species-tree.tre", {newick_species_tree_alpha});
            output(format("  True species tree height = %.9f\n") % sf.getHeight());
            output("  True species tree saved in file \"true-species-tree.tre\"\n");
        }
        else {
            output("  True species tree not saved because it is just a single node!\n");
        }

        // Output tree file containing true gene trees
        vector<string> newick_gene_trees_alpha;
        vector<string> newick_gene_trees_numeric;
        //output("  True gene tree height (height/theta):\n");
        for (auto & gf : gfvect) {
            string newick_alpha = gf.makeNewick(/*precision*/9, /*use names*/true, edgelens_in_coalescent_units);
            newick_gene_trees_alpha.push_back(newick_alpha);
            
            string newick_numeric = gf.makeNewick(/*precision*/9, /*use names*/false, edgelens_in_coalescent_units);
            newick_gene_trees_numeric.push_back(newick_numeric);
            //output(format("  %12.9f (%.9f)\n") % gf.getHeight() % (gf.getHeight()/G::_theta));
        }
        outputGeneTreesToFile("true-gene-trees.tre", newick_gene_trees_alpha);
        output("  True gene trees saved in file \"true-gene-trees.tre\"\n");

        // Output gene trees and species trees for javascript viewer
        outputJavascriptTreefile("newicks.js", newick_species_tree_numeric, newick_gene_trees_numeric);
    }
#endif
    
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
    
    inline void Proj::simulateTrees(Particle & particle) {
        unsigned ngenes   = G::_nloci;
        unsigned ntaxa    = G::_ntaxa;
        
        // Start with trivial species forest
        particle.resetSpeciesForest();

#if defined(EST_THETA)
        particle.setThetas();
#endif

        // Start with trivial gene forests
        particle.resetGeneForests(/*compute_partials*/false);
        
#if defined(LAZY_COPYING)
        assert(particle.getGeneForestPtrs().size() == ngenes);
#else
        assert(particle.getGeneForests().size() == ngenes);
#endif
        
        // Determine total number of steps required to build all gene trees and the species tree
        unsigned nsteps = ngenes*(ntaxa - 1);
        
        G::generateUpdateSeeds(nsteps);
        
        vector< pair<double, unsigned> > locus_ordering(G::_nloci);

        // Add coalescent events and speciation events until all trees are built
        for (unsigned step = 0; step < nsteps; step++) {
            unsigned step_modulus = step % G::_nloci;
            if (step_modulus == 0) {
#if defined(RANDOM_LOCUS_ORDERING)
                // Time to choose a new random locus ordering
                for (unsigned i = 0; i < G::_nloci; i++) {
                    double u = rng->uniform();
                    locus_ordering[i] = make_pair(u, i);
                }
#else
                // Order loci sequentially
                for (unsigned i = 0; i < G::_nloci; i++) {
                    locus_ordering[i] = make_pair(i, i);
                }
#endif
                sort(locus_ordering.begin(), locus_ordering.end());
            }
            unsigned locus = locus_ordering[step_modulus].second;
            particle.proposeCoalescence(step, locus, G::_seed_bank[step], /*rebuild_species_tree*/locus == locus_ordering[0].second);
        }
        
        particle.refreshHeightsInternalsPreorders();
    }
    
    inline void Proj::reportDeepCoalescences(Particle & particle) {
        vector<Forest::coalinfo_t> coalinfo_vect;
        particle.rebuildCoalInfo();
        particle.recordAllForests(coalinfo_vect);
        
#if defined(OUTPUT_FOR_DEBUGGING)
        output("\nContents of coalinfo_vect:\n");
        output(format("%12s %6s %s\n") % "height" % "locus" % "species");
        for (auto & cinfo : coalinfo_vect) {
            double      height = get<0>(cinfo);
            unsigned     locus = get<1>(cinfo);
            auto & sppvect = get<2>(cinfo);
            vector<string> s;
            for (auto x : sppvect) {
                s.push_back(to_string(x));
            }
            output(format("%12.9f %6d %s\n") % height % locus % boost::algorithm::join(s, ","));
        }
#endif

        // Count starting number of lineages in each species
        map<G::species_t, unsigned> lineages_per_species;
        for (auto t : G::_taxon_names) {
            unsigned i = G::_taxon_to_species[t];
            G::species_t s = (G::species_t)1 << i;
            lineages_per_species[s] += G::_nloci;
        }
        
        // Make a copy that will not be affected by coalescences
        // that is used to count the maximum number of possible deep
        // coalescences
        map<G::species_t, unsigned> lineages_per_species0 = lineages_per_species;
        
        //output("lineages_per_species map:\n");
        //for (auto & p : lineages_per_species) {
        //    output(format("%12d lineages in species %d\n") % p.second % p.first);
        //}

        // Count deep coalscences at each locus
        unsigned ndeep = 0;
        unsigned maxdeep = 0;
        for (auto & cinfo : coalinfo_vect) {
            unsigned     locus = get<1>(cinfo);
            auto & sppvect = get<2>(cinfo);
            if (locus > 0) {
                // Coalescence
                assert(sppvect.size() == 2);
                assert(sppvect[0] == sppvect[1]);
                G::species_t s = sppvect[0];
                lineages_per_species[s]--;
                
                //output(format("~~> Coalescence in species %d at locus %d\n") % s % locus);
            }
            else {
                // Speciation
                assert(sppvect.size() == 2);
                
                // Count actual number of deep coalescences
                unsigned nleft = lineages_per_species.at(sppvect[0]);
                unsigned nright = lineages_per_species.at(sppvect[1]);
                unsigned deep = nleft + nright - 2;
                ndeep += deep;
                
                // Count maximum number of deep coalescences
                unsigned nleft0 = lineages_per_species0.at(sppvect[0]);
                unsigned nright0 = lineages_per_species0.at(sppvect[1]);
                unsigned deep0 = nleft0 + nright0 - 2;
                maxdeep += deep0;
                
                // Merge lineages from left and right species
                G::species_t ancspp = sppvect[0] | sppvect[1];

                lineages_per_species[sppvect[0]] = 0;
                lineages_per_species[sppvect[1]] = 0;
                lineages_per_species[ancspp] = nleft + nright;

                lineages_per_species0[sppvect[0]] = 0;
                lineages_per_species0[sppvect[1]] = 0;
                lineages_per_species0[ancspp] = nleft0 + nright0;
                
                //output(format("~~> Speciation event combines species %d and %d\n") % sppvect[0] % sppvect[1]);
                //output(format("~~> No. deep coalescences here is %d (total is %d)\n") % deep % ndeep);
                //output(format("~~> Max. deep coalescences here is %d (total is %d)\n") % deep0 % maxdeep);
            }
        }
        
        output(format("\n  Number of deep coalescences = %d\n") % ndeep);
        output(format("  Maximum number of deep coalescences = %d\n") % maxdeep);
    }
            
    inline void Proj::simulateData() {
        if (_nsimtaxaperspecies.size() == 0) {
            throw XProj("You must define nspecies and ntaxaperspecies in order to simulate data");
        }
        G::_simulating = true;
        G::_nspecies = _nsimspecies;
        G::_ntaxa = (unsigned)accumulate(_nsimtaxaperspecies.begin(), _nsimtaxaperspecies.end(), 0);

        output("Simulating sequence data under multispecies coalescent model:\n");
        output(format("  random number seed = %.5f\n") % G::_rnseed);
        output(format("  theta  = %.5f\n") % G::_theta);
        output(format("  lambda = %.5f\n") % G::_lambda);
        output(format("  Number of species = %d\n") % _nsimspecies);
        
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
        output(format("  Expected species tree height = %.5f\n") % expected_species_tree_height);
        
        // Expected gene tree height without species tree constraints
        // is theta/n(n-1) + theta/(n-1)(n-2) + ... + theta/(2)(1)
        double expected_gene_tree_height = 0.0;
        for (unsigned i = 2; i <= G::_ntaxa; i++) {
            expected_gene_tree_height += 1.0/(i*(i-1));
        }
        expected_gene_tree_height *= G::_theta;
        output(format("  Expected gene tree height (unconstrained) = %.5f\n") % expected_gene_tree_height);
        
        // Interrogate _partition to determine number of genes, gene names, and
        // number of sites in each gene
        G::_nloci = _partition->getNumSubsets();
        G::_nsites_per_gene.resize(G::_nloci);
        G::_gene_names.resize(G::_nloci);
        unsigned total_nsites = 0;
        for (unsigned g = 0; g < G::_nloci; g++) {
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
            string species_name = G::inventName(i, false);
            G::_species_names[i] = species_name;
            G::_taxon_to_species[species_name] = i;
            for (unsigned j = 0; j < _nsimtaxaperspecies[i]; ++j) {
                taxpartition.push_back(G::_species_names[i]);
                string taxon_name = str(format("%s^%s") % G::inventName(k, true) % G::_species_names[i]);
                G::_taxon_names[k] = taxon_name;
                G::_taxon_to_species[taxon_name] = i;
                ++k;
            }
        }

        // Report species names created
#if defined(OUTPUT_FOR_DEBUGGING)
        output("  Species info (index is into G::_species_names):\n");
        output(format("    %12s %12s %s\n") % "index" % "species" % "name");
        for (unsigned i = 0; i < G::_nspecies; ++i) {
            G::species_t spp = (G::species_t)1 << i;
            output(format("    %12d %12d %s\n") % i % spp % G::_species_names[i]);
        }
#else
        output("\n  Species names:\n");
        for (unsigned i = 0; i < G::_nspecies; ++i) {
            output(format("    %s\n") % G::_species_names[i]);
        }
#endif

        // Report taxon names created
        output("\n  Taxon names:\n");
        for (unsigned i = 0; i < G::_ntaxa; ++i) {
            string taxon_name = G::_taxon_names[i];
            output(format("    %s\n") % taxon_name);
        }
        
        // Create data object
        assert(!_data);
        _data = Data::SharedPtr(new Data());
        _data->setTaxonNames(G::_taxon_names);
        _data->setPartition(_partition);
        
        // Build species tree and gene trees jointly
        Particle particle;
        particle.setData(_data);
        //particle.setSeed(rng->randint(1,9999));
        simulateTrees(particle);
        
#if defined(OUTPUT_FOR_DEBUGGING)
        output(format("Simulated species tree:\n  %s\n") % particle.getSpeciesForestConst().makeNewick(9, true, false));
        for (unsigned g = 0; g < G::_nloci; g++) {
            output(format("Simulated gene tree for locus %d:\n  %s\n") % (g+1) % particle.getGeneForestConst(g).makeNewick(9, true, false));
        }
#endif
        
        reportDeepCoalescences(particle);
        
        SpeciesForest      & species_forest = particle.getSpeciesForest();
#if defined(LAZY_COPYING)
        vector<GeneForest::SharedPtr> & gene_forest_ptrs   = particle.getGeneForestPtrs();
        outputTrees(species_forest, gene_forest_ptrs);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNLoci(G::_nloci);

        // Simulate sequence data
        unsigned starting_site = 0;
        for (unsigned g = 0; g < G::_nloci; ++g) {
            ps.setNElements(4*G::_nsites_per_gene[g], g);
            gene_forest_ptrs[g]->simulateData(::rng, _data, starting_site, G::_nsites_per_gene[g]);
            starting_site += G::_nsites_per_gene[g];
        }
#else
        vector<GeneForest> & gene_forests   = particle.getGeneForests();
        outputTrees(species_forest, gene_forests);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNLoci(G::_nloci);

        // Simulate sequence data
        unsigned starting_site = 0;
        for (unsigned g = 0; g < G::_nloci; ++g) {
            ps.setNElements(4*G::_nsites_per_gene[g], g);
            gene_forests[g].simulateData(::rng, _data, starting_site, G::_nsites_per_gene[g]);
            starting_site += G::_nsites_per_gene[g];
        }
#endif
                       
        // Output data for each locus to a separate file
        
        
        // Output data to file
        _data->compressPatterns();
        _data->writeDataToFile(_data_file_prefix + ".nex");
        output(format("  Sequence data saved in file \"%s.nex\"\n") % _data_file_prefix, 1);

        // Output data for each locus to a separate file for use with beast
        //_data->writeLocusSpecificDataFiles(_data_file_prefix);
        
        // Output a PAUP* command file for estimating the species tree using
        // svd quartets and qage
        output("  PAUP* commands saved in file \"svd-qage.nex\"\n");
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
     
    inline void Proj::secondLevelRange(vector<Particle> & first_level_particles, SMC & ensemble, unsigned nparticles, unsigned nkept, unsigned first, unsigned last) {
        for (unsigned i = first; i < last; i++) {
            unsigned rnseed = 1;
            Particle * p = nullptr;
            {
#if defined(USING_MULTITHREADING)
                lock_guard<mutex> guard(G::_mutex);
#endif
                p = &first_level_particles[i];
                rnseed = G::_seed_bank[i];
                //output(format("  kept %d\n") % (i+1), G::LogCateg::INFO);
            }
            secondLevelConditionedOn(rnseed, *p, nparticles, nkept, ensemble);
        }
    }
    
    inline void Proj::secondLevelConditionedOn(unsigned rnseed, const Particle & p, unsigned nparticles, unsigned nkept, SMC & ensemble) {
        
        SMC smc2;
        smc2.setNParticles(G::_nparticles2, 1); //TODO: let 2nd level have subpops?
        smc2.setMode(SMC::SPECIES_GIVEN_GENE);
        smc2.initFromParticle(p, _data);
        smc2.setRandomNumberSeed(rnseed);
        smc2.run();
        
        // Choose (discrete-uniform-randomly) G::_nkept2 2nd-level particles
        vector<unsigned> kept2;
        smc2.keepSecondLevel(kept2);

        {
#if defined(USING_MULTITHREADING)
            lock_guard<mutex> guard(G::_mutex);
#endif
            smc2.dumpParticles(ensemble, kept2);
        }
    }
    
    inline void Proj::buildEnsembleCoalInfo(vector<Particle> & first_level_particles) {
        for (auto & p : first_level_particles) {
            // Gather gene tree info into particle's _ensemble_coalinfo vector
            p.rebuildCoalInfo();
            p.buildEnsembleCoalInfo();
        }
    }
    
    inline void Proj::stripPartials(vector<Particle> & first_level_particles) {
        // Strip all partials from gene trees (not needed in second level

#if 0
        //TODO: needs work
        // First create map whose keys store unique gene forest pointers
        // and whose values store the particle and locus of that gene forest
        map<const void *, vector<pair<unsigned, unsigned> > > unique_gene_forests;
        map<unsigned, unsigned> locus_counts;
        map<unsigned, unsigned> particle_counts;
        unsigned particle_index = 0;
        for (auto & p : first_level_particles) {
            unsigned locus_index = 0;
            vector<GeneForest::SharedPtr> & gene_forest_ptrs = p.getGeneForestPtrs();
            for (auto & gfp : gene_forest_ptrs) {
                const void * ptr = (const void *)gfp.get();
                if (unique_gene_forests.count(ptr) > 0) {
                    unique_gene_forests[ptr].push_back(make_pair(particle_index, locus_index));
                }
                else {
                    unique_gene_forests[ptr] = {make_pair(particle_index, locus_index)};
                }
                locus_index++;
            }
            particle_index++;
        }
        
        //temporary!
        for (auto & m : unique_gene_forests) {
            //output(format("%s\n") % G::memoryAddressAsString(m.first));
            for (auto & v : m.second) {
                //output(format("  %d, %d\n") % v.first % v.second);
                particle_counts[v.first]++;
                locus_counts[v.second]++;
            }
        }
        //temporary!
        output("\nNumber of unique gene forests per particle:\n");
        for (auto & m : particle_counts) {
            output(format("  %d: %d\n") % m.first % m.second);
        }
        //temporary!
        output("\nNumber of unique gene forests per locus:\n");
        for (auto & m : locus_counts) {
            output(format("  %d: %d\n") % m.first % m.second);
        }
        cerr << endl;
#endif
    }
     
    inline void Proj::run() {
    
        ::rng->setSeed(G::_rnseed);
        
        if (_start_mode == "sim") {
            simulateData();
        }
        else if (_start_mode == "smc") {
            output("Starting...\n");

            output(format("Current working directory: %s\n") % current_path());
        
            readData();
            
            G::_nloci = _data->getNumSubsets();
            assert(G::_nloci > 0);
            
            setRelativeRates();
            
#if defined(SPECIES_IN_CONF)
            if (G::_nspecies > 0) {
                // Species specified in the conf file
                // Check that taxon names are the same as those
                // in the data file
                _data->checkTaxonNames(G::_taxon_names);
            }
            else {
                // Copy taxon names to global variable _taxon_names
                G::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(G::_taxon_names);
                
                // Save species names to global variable _species_names
                // and create global _taxon_to_species map that provides
                // the species index for each taxon name
                G::_nspecies = buildSpeciesMap(/*taxa_from_data*/true);
            }
#else
                // Copy taxon names to global variable _taxon_names
                G::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(G::_taxon_names);
                
                // Save species names to global variable _species_names
                // and create global _taxon_to_species map that provides
                // the species index for each taxon name
                G::_nspecies = buildSpeciesMap(/*taxa_from_data*/true);
#endif
            
            // Create global _species_mask that has "on" bits only for
            // the least significant _nspecies bits
            Node::setSpeciesMask(G::_species_mask, G::_nspecies);
            
            // Precompute leaf partials (never have to be computed again)
            for (unsigned g = 0; g < G::_nloci; g++) {
                GeneForest::computeLeafPartials(g, _data);
            }
            
            output(format("Performing 1st-level SMC using %d particles...\n") % G::_nparticles, G::LogCateg::INFO);

            // First-level particle filtering
            SMC smc;
            smc.setNParticles(G::_nparticles, G::_nsubpops);
            smc.setData(_data);
            smc.init();
            smc.run();
            smc.summarize();
            
            bool second_level = (G::_nparticles2 > 0 && G::_nkept2 > 0);
            if (second_level) {
                output(format("\nPerforming 2nd-level SMC on %d 1st-level particles...\n") % G::_nkept, G::LogCateg::INFO);
                
                vector<Particle> & first_level_particles = smc.getParticles();
                assert(first_level_particles.size() == G::_nparticles);
                
                // No longer need partials stored in gene trees
                //stripPartials(first_level_particles);
                
                // Calculate ensemble coal info vectors for each first-level
                // particle to avoid accessing gene tree pointers in multithreaded
                // version
                buildEnsembleCoalInfo(first_level_particles);

                // Choose (discrete-uniform-randomly) G::_nkept
                // 1st-level particles for use in 2nd level
                // Save indices in kept vector
                vector<unsigned> kept;
                kept.reserve(first_level_particles.size());
                for (unsigned k = 0; k < G::_nkept; k++) {
                    unsigned i = rng->randint(0, G::_nparticles-1);
                    kept.push_back(i);
                }

                // Second-level particle filtering
                SMC ensemble;
                ensemble.setMode(SMC::SPECIES_GIVEN_GENE);
                
                // Replace seeds in seed bank
                G::generateUpdateSeeds(G::_nkept);

                unsigned nparticles = G::_nparticles2;
                unsigned nkept = G::_nkept2;

#if defined(USING_MULTITHREADING)
                G::buildThreadSchedule(G::_nkept, "kept particle");
                
                vector<thread> threads;
                for (unsigned i = 0; i < G::_nthreads; i++) {
                    threads.push_back(thread(&Proj::secondLevelRange,
                        this,
                        ref(first_level_particles),
                        ref(ensemble),
                        nparticles,
                        nkept,
                        G::_thread_sched[i].first,
                        G::_thread_sched[i].second)
                    );
                }

                // The join function causes this loop to pause until
                // the ith thread finishes
                for (unsigned i = 0; i < threads.size(); i++) {
                    threads[i].join();
                }
#else
                secondLevelRange(first_level_particles, ensemble, nparticles, nkept, 0, G::_nkept);
#endif
                ensemble.summarize();
            }
            
            ps.debugReport();
        }
        else {
            throw XProj("This program currently only accepts \"sim\" and \"smc\" as start mode");
        }
        
        output("\nFinished!\n");
    }
    
    inline unsigned Proj::buildSpeciesMap(bool taxa_from_data) {
        // Assumes G::_taxon_names is already filled.
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

        output("\nMapping taxon names to species index:\n");
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
                output(format("  %s --> %s\n") % tname % species_name);
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
        
        output("\nMapping species names to species index:\n");
        for (auto & sname : G::_species_names) {
            unsigned species_index = 0;
            if (species_name_to_index.count(sname) == 0)
                throw XProj(format("Proj::buildSpeciesMap failed because key \"%s\" does not exist in species_name_to_index map") % sname);
            else {
                species_index = species_name_to_index.at(sname);
            }
            output(format("  %s --> %d\n") % sname % species_index);
            
            // Note: despite appearances, this next line does not
            // overwrite anything. We need to be able to map taxon
            // names in species trees to a species index as well as
            // taxon names in gene trees
            G::_taxon_to_species[sname] = species_index;
        }

        return (unsigned)G::_species_names.size();
    }
    
}

