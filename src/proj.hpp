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

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
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
            void                            readData();
            void                            run();
               
            string                          _data_file_name;
            string                          _start_mode;
            Partition::SharedPtr            _partition;
            Data::SharedPtr                 _data;
            
            SpeciesForest::SharedPtr        _species_tree_ref;
            vector<GeneForest::SharedPtr>   _gene_tree_refs;
            
            bool                            _use_gpu;
            bool                            _ambig_missing;
            unsigned                        _rnseed;
                                       
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
        _data_file_name             = "";
        _species_tree_ref           = nullptr;
        _start_mode                 = "smc";
        _gene_tree_refs.clear();
        _partition.reset(new Partition());
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> partition_subsets;
        vector<string> partition_relrates;
        bool dummy_bool;
        int dummy_int;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile",  value(&_data_file_name), "name of a data file in NEXUS format")
        ("speciestreeref",  value(&G::_species_tree_ref_file_name), "name of a tree file containing a single reference species tree")
        ("genetreeref",  value(&G::_gene_trees_ref_file_name), "name of a tree file containing a reference gene tree for each locus")
        ("startmode", value(&_start_mode), "if 'sim', simulate gene trees, species tree, and data; if 'smc', estimate from supplied datafile; if 'chib', computes prior probability of species species and gene tree topologies; if 'spec', computes coalescent likelihood for specified speciestreeref and genetreeref; if '2ndlevel', tests second-level SMC from gene trees supplied by genetreeref")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("verbosity",  value(&G::_verbosity)->default_value(0), "0, 1, or 2: higher number means more output")
        ("nparticles",  value(&G::_nparticles)->default_value(500), "number of particles in a population for joint estimation")
        ("nspeciesparticles",  value(&G::_nparticles2)->default_value(1000), "number of particles in a population for species tree only estimation")
        ("lambda",  value(&G::_lambda)->default_value(10.9), "per lineage speciation rate assumed for the species tree")
        ("theta",  value(&G::_theta)->default_value(0.05), "coalescent parameter assumed for gene trees")
        ("rnseed",  value(&_rnseed)->default_value(13579), "pseudorandom number seed")
        ("freezethetamean",  value(&dummy_bool)->default_value(true), "this option is not used in this version of the program")
        ("fixedthetamean",  value(&G::_theta)->default_value(0.05), "synonym of theta in this version of the program")
        ("nthreads",  value(&dummy_int)->default_value(3), "this option is not used in this version of the program")
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
             
    inline void Proj::run() {
    
        if (_start_mode != "smc") {
            throw XProj("This program currently only accepts \"smc\" as start mode");
        }
    
        output("Starting...\n", 2);

#if defined(UPGMA_WEIGHTS)
        output("Using UPGMA completion\n", 2);
#endif

        output(format("Current working directory: %s\n") % current_path(), 2);
        
        try {
            ::rng->setSeed(_rnseed);
            
            readData();
            
            G::_nloci = _data->getNumSubsets();
            assert(G::_nloci > 0);

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
            
            // Precompute leaf partials (never have to be computed again)
            for (unsigned g = 0; g < G::_nloci; g++) {
                GeneForest::computeLeafPartials(g, _data);
            }
            
            // First-level particle filtering
            SMC smc;
            smc.setData(_data);
            smc.init();
            smc.run();
            smc.summarize();
        }
        catch (XProj & x) {
            output(format("Proj encountered a problem:\n  %s\n") % x.what(), 2);
        }
        
        output("\nFinished!\n", 2);
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
    
}
