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

extern void output(string msg, proj::G::verbosity_t verb);
extern void output(format & fmt, proj::G::verbosity_t level);
extern proj::PartialStore         ps;
extern proj::StopWatch            stopwatch;
extern proj::Lot::SharedPtr       rng;
extern proj::Partition::SharedPtr partition;
extern proj::Data::SharedPtr      data;
extern vector<string>             rcolors;

namespace proj {

    class Proj {
        public:
                                         Proj();
                                         ~Proj();

            void                         clear();
            void                         processCommandLineOptions(int argc, const char * argv[]);
            void                         run();
                        
#if defined(LOG_MEMORY)
            void                         memoryReport(ofstream & memf) const;
#endif
            
        private:
        
            static string                inventName(unsigned k, bool lower_case);


            void                         readTreefile(  const string                filename,
                                                        unsigned                    skip,
                                                        vector<string> &            leaf_names,
                                                        map<unsigned,unsigned> &    taxon_map,
                                                        vector<string> &            tree_names,
                                                        vector<string> &            newicks) const;
            unsigned                     buildSpeciesMap();
            void                         debugShowStringVector(string title, const vector<string> & svect, G::verbosity_t verb) const;
            void                         compareWithReferenceTrees(const SMC & smc) const;
            void                         readData();

            string                       _data_file_name;
            string                       _start_mode;
                        
            bool                         _use_gpu;
            bool                         _ambig_missing;
            unsigned                     _rnseed;
                        
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
        _data_file_name             = "";
        _start_mode                 = "smc";
        ::partition.reset(new Partition());
        ::data                      = nullptr;
        _use_gpu                    = true;
        _ambig_missing              = true;
        _rnseed                     = 1;
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        bool debugging = false;
        bool temp = false;
        vector<string> partition_subsets;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile",  value(&_data_file_name), "name of a data file in NEXUS format")
        ("startmode", value(&_start_mode), "if 'sim', simulate gene trees, species tree, and data; if 'smc', estimate from supplied datafile; if 'chib', computes prior probability of species species and gene tree topologies; if 'spec', computes coalescent likelihood for specified speciestreeref and genetreeref; if '2ndlevel', tests second-level SMC from gene trees supplied by genetreeref")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("speciestreeref",  value(&G::_species_tree_ref_file_name), "name of a tree file containing a single reference species tree")
        ("genetreeref",  value(&G::_gene_trees_ref_file_name), "name of a tree file containing a reference gene tree for each locus")
        ("nthreads",  value(&G::_nthreads)->default_value(1), "number of threads (each thread will handle nparticles/nthreads particle updates)")
        ("nbundles", value(&G::_nbundles)->default_value(100), "number of bundles")
        ("nparticles", value(&G::_nparticles)->default_value(100), "number of gene-tree particles per locus per bundle")
        ("lambda", value(&G::_lambda)->default_value(1.0), "speciation rate")
        ("theta", value(&G::_theta)->default_value(0.1), "mutation-scaled population size parameter")
        ("ambigmissing",  value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("debugging",  value(&debugging)->default_value(false), "if yes, shows debugging output")
        ("temp",  value(&temp)->default_value(false), "if yes, shows output from temporary debugging code")
        ("rnseed",  value(&_rnseed)->default_value(13579), "pseudorandom number seed");
        
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
            output(format("%s\n") % desc, G::VSTANDARD);
            exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            output(format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version, G::VSTANDARD);
            exit(1);
        }
        
        // If user specified --debugging yes on command line, set G::VDEBUG bit in G::_verbosity
        if (debugging) {
            G::_verbosity = (G::verbosity_t)(G::_verbosity | G::VDEBUG);
        }
        
        // If user specified --temp yes on command line, set G::VTEMP bit in G::_verbosity
        if (temp) {
            G::_verbosity = (G::verbosity_t)(G::_verbosity | G::VTEMP);
        }
        
        // If user specified --subset on command line, break specified partition subset
        // definition into name and character set string and add to G::_partition
        if (vm.count("subset") > 0) {
            ::partition.reset(new Partition());
            for (auto s : partition_subsets) {
                ::partition->parseSubsetDefinition(s);
            }
        }

#if defined(USING_MULTITHREADING)
        // nothing to do
#else
        if (vm.count("nthreads") > 0) {
            if (G::_nthreads != 1) {
                output(format("\nWARNING: You specified %d threads but this non-multithreading version only allows 1 thread\nProceeding with a single thread.\n\n")  % G::_nthreads,G::VSTANDARD);
                G::_nthreads = 1;
            }
        }
#endif
    }
    
    inline void Proj::readData() {
        output(format("\nReading and storing the data in the file %s\n") % _data_file_name, G::VSTANDARD);
        ::data = Data::SharedPtr(new Data());
        ::data->setPartition(::partition);
        ::data->getDataFromFile(_data_file_name);
        
        // Report information about data partition subsets
        unsigned nsubsets = ::data->getNumSubsets();
        output(format("\nNumber of taxa: %d\n") % ::data->getNumTaxa(), G::VSTANDARD);
        output(format("Number of partition subsets: %d") % nsubsets, G::VSTANDARD);
        
        // Inform PartialStore of number of genes so that it can allocate
        // its _nelements and _storage vectors
        ps.setNGenes(nsubsets);
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Set length of partials for gene g
            ps.setNElements(G::_nstates*::data->getNumPatternsInSubset(subset), subset);
            
            DataType dt = ::partition->getDataTypeForSubset(subset);
            G::_locus_names.push_back(::data->getSubsetName(subset));
            output(format("  Subset %d (%s)\n") % (subset+1) % ::data->getSubsetName(subset), G::VSTANDARD);
            output(format("    data type: %s\n") % dt.getDataTypeAsString(), G::VSTANDARD);
            output(format("    sites:     %d\n") % ::data->calcSeqLenInSubset(subset), G::VSTANDARD);
            output(format("    patterns:  %d\n") % ::data->getNumPatternsInSubset(subset), G::VSTANDARD);
        }

        debugShowStringVector("Locus names", G::_locus_names, G::VSTANDARD);
        
        G::_nloci = ::data->getNumSubsets();
        assert(G::_nloci > 0);

        // Populate G::_taxon_names
        G::_ntaxa = ::data->getNumTaxa();
        ::data->copyTaxonNames(G::_taxon_names);
        
        // Populate G::_species_names and G::_taxon_to_species
        G::_nspecies = buildSpeciesMap();
        
        // Calculate partials for leaves and store in G::_leaf_partials[locus] for all loci
        GParticle::computeLeafPartials();
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

    inline void Proj::debugShowStringVector(string title, const vector<string> & svect, G::verbosity_t verb) const {
        output(format("\n%s\n") % title, verb);
        for (auto s : svect) {
            output(format("  %s\n") % s, verb);
         }
    }

    inline unsigned Proj::buildSpeciesMap() {
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

        output("\nMapping taxon names to species index:\n", G::VSTANDARD);
        unsigned ntax = (unsigned)G::_taxon_names.size();
        assert(ntax > 0);
        
        // Assume taxon names are already stored in _data object and no
        // species names have yet been stored
        G::_species_names.clear();
        const Data::taxon_names_t & tnames = ::data->getTaxonNames();
        assert(tnames.size() > 0);
        ntax = (unsigned)tnames.size();
        for (auto & tname : tnames) {
            string species_name = Node::taxonNameToSpeciesName(tname);
            output(format("  %s --> %s\n") % tname % species_name, G::VSTANDARD);
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
        
        output("\nMapping species names to species index:\n", G::VSTANDARD);
        for (auto & sname : G::_species_names) {
            unsigned species_index = 0;
            if (species_name_to_index.count(sname) == 0)
                throw XProj(format("Proj::buildSpeciesMap failed because key \"%s\" does not exist in species_name_to_index map") % sname);
            else {
                species_index = species_name_to_index.at(sname);
            }
            output(format("  %s --> %d\n") % sname % species_index, G::VSTANDARD);
            
            // Note: despite appearances, this next line does not
            // overwrite anything. We need to be able to map taxon
            // names in species trees to a species index as well as
            // taxon names in gene trees
            G::_taxon_to_species[sname] = species_index;
        }

        return (unsigned)G::_species_names.size();
    }
    
    inline void Proj::readTreefile(const string filename,
                                     unsigned skip,
                                     vector<string> & leaf_names,
                                     map<unsigned,unsigned> & taxon_map,
                                     vector<string> & tree_names,
                                     vector<string> & newicks) const {
                                     
        // Reads file, skipping the first skip trees in each trees block and
        // storing newick tree descriptions (in newicks) and tree names (in tree_names).
        // User-supplied taxon_map will be filled with entries such that taxon_map[k]
        // provides the index (into leaf_names) of the taxon whose number in the taxa block
        // (and thus also in the newick tree description) equals k. (If the tree file does
        // not supply a taxa block, one will be created in which the ordering of taxa is
        // arbitrary (perhaps the order in which taxa are encountered in the first tree?).

        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

        // Read and process the file
        //MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
        MultiFormatReader nexus_reader(-1, NxsReader::IGNORE_WARNINGS);

        // Both of these needed to suppress "storing read block" messages
        // see NxsReader::statusMessage in nxsreader.cpp
        nexus_reader.SetAlwaysReportStatusMessages(false);
        nexus_reader.SetWarningOutputLevel(NxsReader::SUPPRESS_WARNINGS_LEVEL);
        
        try {
            nexus_reader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...) {
            nexus_reader.DeleteBlocksFromFactories();
            throw;
        }

        // Get number of taxa blocks
        int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
        
        // Process the trees blocks associated with each taxa block
        for (int i = 0; i < num_taxa_blocks; ++i) {
            NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(i);
            string taxa_block_title = taxa_block->GetTitle();
            
            // Ensure that the number of leaves accords with expectations
            unsigned num_leaves = (unsigned)taxa_block->GetNumTaxonLabels();
            if (num_leaves != leaf_names.size()) {
                if (leaf_names.empty()) {
                    // Populate leaf_names using taxa block
                    leaf_names.resize(num_leaves);
                    for (unsigned t = 0; t < num_leaves; t++) {
                        leaf_names[t] = taxa_block->GetTaxonLabel(t);
                    }
                    output(format("Assuming %d taxon names are those found in taxa block") % num_leaves, G::VSTANDARD);
                }
                else {
                    throw XProj(format("Found %d taxa in taxa block %d but expecting %d") % num_leaves % (i+1) % leaf_names.size());
                }
            }
            
            // It is possible (probable even) that the order of taxa in taxa_block will differ from
            // the order in leaf_names. Therefore, taxon_map is constructed the keys of which are the
            // taxon numbers in the newick tree description (corresponding to taxa_block); the values
            // are the index of the taxon in leaf_names. For example (branch lengths omitted):
            //
            // #NEXUS
            // begin trees;
            //   tree tree1 = [&R] (A,(((D,B),C),E));
            // end;
            //
            // Supplied to readTreeFile function:
            //   leaf_names = {"A","B","C","D","E"}
            //                  0   1   2   3   4
            //
            // Populated by readTreeFile function:
            //   tree_names = {"tree1"}
            //   newicks = {"(1,(((2,3),4),5))"}
            //   taxon_map: (newick)    (leaf_names)
            //                key          value
            //                 1    -->      0 (i.e. "A")
            //                 2    -->      3 (i.e. "D")
            //                 3    -->      1 (i.e. "B")
            //                 4    -->      2 (i.e. "C")
            //                 5    -->      4 (i.e. "E")
            
            taxon_map.clear();
            for (unsigned t = 0; t < num_leaves; t++) {
                string taxon_label = taxa_block->GetTaxonLabel(t);
                auto it = find(leaf_names.begin(), leaf_names.end(), taxon_label);
                if (it == leaf_names.end()) {
                    throw XProj(format("Taxon label \"%s\" was not expected ") % taxon_label);
                }
                unsigned d = (unsigned)(distance(leaf_names.begin(), it));
                taxon_map[t+1] = d;
            }
            
            // Get number of trees blocks for this taxa block
            const unsigned num_trees_blocks = nexus_reader.GetNumTreesBlocks(taxa_block);
            
            // Get all trees for each trees block, skipping the first skip trees
            for (unsigned j = 0; j < num_trees_blocks; ++j) {
                // GetTreeName is inexplicably non-const, necessitating
                // making trees_block non-const
                NxsTreesBlock * trees_block = nexus_reader.GetTreesBlock(taxa_block, j);
                unsigned num_trees = trees_block->GetNumTrees();
                if (skip < num_trees) {
                    for (unsigned t = skip; t < num_trees; ++t) {
                        const NxsFullTreeDescription & d = trees_block->GetFullTreeDescription(t);

                        if (!d.IsRooted()) {
                            throw XProj(format("Tree %d in trees block %d and taxa block %d is unrooted; expecting only rooted trees") % (t+1) % (j+1) % (i+1));
                        }
                        
                        // Store the tree name
                        string tree_name = trees_block->GetTreeName(t);
                        tree_names.push_back(tree_name);
                        
                        // Store the newick tree description
                        string newick = d.GetNewick();
                        newicks.push_back(newick);
                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

        // No longer any need to store raw data from nexus file
        nexus_reader.DeleteBlocksFromFactories();
    }
            
    inline void Proj::compareWithReferenceTrees(const SMC & smc) const {
        assert(G::_nspecies == (unsigned)G::_species_names.size());
        assert(G::_ntaxa == (unsigned)G::_taxon_names.size());
        
        if (G::_species_tree_ref_file_name.length() > 0) {
            // Read in the reference species tree
            G::_nexus_taxon_map.clear();
            vector<string> species_tree_names;
            vector<string> species_newicks;
            readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, G::_nexus_taxon_map, species_tree_names, species_newicks);
            
            // Create SParticle with reference species tree (assumed to be the first tree in the file)
            SParticle species_ref;
            species_ref.buildFromNewick(species_newicks[0], G::_nexus_taxon_map, G::_species_names);

            // Save log marginal likelihood for every bundle
            vector<double> log_marg_like;
            smc.saveLogMargLike(log_marg_like);
            unsigned nbundles = (unsigned)log_marg_like.size();

            // Create R file that plots RF distances (y-axis) against
            // log marginal likelihood (x-axis)
            ofstream Rfile("species-tree-scatterplot.R");
            Rfile << "rf <- c(";
            for (unsigned b = 0; b < nbundles; b++) {
                auto kf_rf = Particle::calcTreeDistances(species_ref, smc.getSpeciesTreeConst(b));
                Rfile << (b > 0 ? "," : "") << kf_rf.second;
            }
            Rfile << ")\n";
            Rfile << "lnl <- c(";
            for (unsigned i = 0; i < nbundles; i++) {
                Rfile << (i > 0 ? "," : "") << str(format("%.5f") % log_marg_like[i]);
            }
            Rfile << ")\n";
            Rfile << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
            Rfile << "setwd(cwd)\n";
            Rfile << "pdf(\"species-tree-scatterplot.pdf\")\n";
            Rfile << "plot(rf, lnl, type=\"p\")\n";
            Rfile << "dev.off()\n";
            Rfile.close();
            
            // Create nexus file with species trees along with
            // PAUP commands for calculating distances and likelihoods
            ofstream nexf("species-tree-distances.nex");
            nexf << "#nexus\n\n";
            nexf << "begin paup;\n";
            nexf << "  log start file=pauplog.txt replace;\n";
            nexf << "  set maxtrees=1000 increase=auto autoinc=1000;\n";
            nexf << "end;\n\n";
            nexf << "begin trees;\n";
            nexf << str(format("  tree strue = [&R] %s;\n") % species_ref.makeNewick(5, true));

            // Save newick tree descriptions for every bundle's species tree
            vector<string> newicks;
            smc.saveSpeciesTrees(newicks, /*compress*/false);
            assert(nbundles == newicks.size());

            for (unsigned b = 0; b < nbundles; b++) {
                nexf << str(format("  tree s%d [log marg. like. = %.9f] = [&R] %s;\n") % (b+1) % log_marg_like[b] % newicks[b]);
            }
            nexf << "end;\n\n";
            nexf << "begin paup;\n";
            nexf << "  treedist reftree=1 measure=rfsymdiff;\n";
            nexf << "  treedist reftree=1 measure=pathlength;\n";
            nexf << "  log stop;\n";
            nexf << "  quit;\n";
            nexf << "end;\n\n";
            nexf.close();
        }
        
        if (G::_gene_trees_ref_file_name.length() > 0) {
            // Read in the reference gene trees
            G::_nexus_taxon_map.clear();
            vector<string> gene_tree_names;
            vector<string> gene_newicks;
            readTreefile(G::_gene_trees_ref_file_name, /*skip*/0, G::_taxon_names, G::_nexus_taxon_map, gene_tree_names, gene_newicks);
        
            // Create GParticle for every reference gene tree
            vector<GParticle> gene_ref_vect(G::_nloci);
            for (unsigned g = 0; g < G::_nloci; g++) {
                gene_ref_vect[g].buildFromNewick(gene_newicks[g], G::_nexus_taxon_map, G::_taxon_names);
            }

#if 1
            // Create R file for each locus that plots RF distance (y-axis)
            // between each gene tree (from every bundle) and the true gene tree
            // against gene tree log likelihood (x-axis)
            unsigned nrcolors = (unsigned)::rcolors.size();
            unsigned begin_site = 1;
            unsigned end_site = 0;
            for (unsigned g = 0; g < G::_nloci; g++) {
                end_site = begin_site + ::partition->numSitesInSubset(g) - 1;

                // Number of points is number of bundles times number of particles per bundle for locus g
                unsigned npoints = G::_nbundles*G::_nparticles;
                
                // log-likelihoods for every gene tree in every bundle for locus g
                vector<double> lnL;
                lnL.reserve(npoints);
                
                // RF distances between every gene tree in every bundle to true gene tree for locus g
                vector<double> dRF;
                dRF.reserve(npoints);
                
                // Use a different random color for every bundle
                vector<string> rcolor;
                rcolor.reserve(npoints);
                for (unsigned b = 0; b < G::_nbundles; b++) {
                    string random_rcolor = ::rcolors[rng->randint(0, nrcolors-1)];
                    fill_n(back_inserter(rcolor), G::_nparticles, random_rcolor);
                }
                
                for (unsigned b = 0; b < G::_nbundles; b++) {
                    // Save vector of log likelihoods for bundle b, locus g
                    vector<double> log_likes;
                    smc.saveLogLikesForLocus(b, g, log_likes);
                    
                    // Dump these G::_nparticles log likelihoods into lnL
                    lnL.insert(lnL.end(), log_likes.begin(), log_likes.end());
                    
                    for (unsigned i = 0; i < G::_nparticles; i++) {
                        auto kf_rf = Particle::calcTreeDistances(gene_ref_vect[g], smc.getGeneTreeConst(b,g,i));
                        dRF.push_back(kf_rf.second);
                    }
                } // bundle loop
                
                assert(npoints == lnL.size());
                assert(npoints == dRF.size());
                assert(npoints == rcolor.size());
                
                // Create R file for this locus
                ofstream Rfile(str(format("locus%d-scatterplot.R") % g));
                Rfile << "rf <- c(";
                for (unsigned i = 0; i < npoints; i++) {
                    Rfile << (i > 0 ? "," : "") << dRF[i];
                }
                Rfile << ")\n";
                Rfile << "lnl <- c(";
                for (unsigned i = 0; i < npoints; i++) {
                    Rfile << (i > 0 ? "," : "") << str(format("%.5f") % lnL[i]);
                }
                Rfile << ")\n";
                Rfile << "pointcolor <- c(";
                for (unsigned i = 0; i < npoints; i++) {
                    Rfile << (i > 0 ? "," : "") << "\"" << rcolor[i] << "\"";
                }
                Rfile << ")\n";
                Rfile << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
                Rfile << "setwd(cwd)\n";
                Rfile << str(format("pdf(\"locus%d-scatterplot.pdf\")\n") % g);
                Rfile << "plot(rf, lnl, type=\"p\", col=pointcolor)\n";
                Rfile << "dev.off()\n";
                Rfile.close();

                // Create nexus file with gene trees for locus g along with
                // PAUP commands for calculating distances and likelihoods
                ofstream nexf(str(format("locus%d-distances.nex") % g));
                nexf << "#nexus\n\n";
                nexf << "begin paup;\n";
                nexf << str(format("  log start file=locus%d-pauplog.txt replace;\n") % g);
                nexf << "  exe ../simulate/simulated.nex;\n";
                nexf << "  set crit=like maxtrees=1000 increase=auto autoinc=1000;\n";
                nexf << "  lset nst=1 basefreq=equal rates=equal pinvar=0 clock;\n";
                nexf << "end;\n\n";
                nexf << "begin trees;\n";
                nexf << str(format("  tree gtrue = [&R] %s;\n") % gene_ref_vect[g].makeNewick(9, true));

                for (unsigned b = 0; b < G::_nbundles; b++) {
                    // Save newick tree descriptions for every bundle's gene trees for locus g
                    vector<string> newicks;
                    smc.saveGeneTreeBundleLocus(b, g, newicks, /*compress*/false);
                    assert(G::_nparticles == newicks.size());
                    
                    vector<double> log_likes;
                    smc.saveLogLikesForLocus(b, g, log_likes);
                    
                    for (unsigned i = 0; i < G::_nparticles; i++) {
                        nexf << str(format("  tree t%d [lnL = %.9f] = [&R] %s;\n") % (i+1) % log_likes[i] % newicks[i]);
                    }
                }
                nexf << "end;\n\n";
                nexf << "begin paup;\n";
                nexf << "  exclude all;\n";
                nexf << str(format("  include %d-%d;\n") % begin_site % end_site);
                nexf << str(format("  lscores all / userbrlens scorefile=locus%d-lscores.txt replace;\n") % g);
                if (npoints > 65536) {
                    // PAUP imposes limit of 65536
                    nexf << str(format("  treedist 1-65536 / reftree=1 measure=rfsymdiff file=locus%d-rfdist.txt replace;\n") % g);
                    nexf << str(format("  treedist 1-65536 / reftree=1 measure=pathlength file=locus%d-pldist.txt replace;\n") % g);
                }
                else {
                    nexf << str(format("  treedist all / reftree=1 measure=rfsymdiff file=locus%d-rfdist.txt replace;\n") % g);
                    nexf << str(format("  treedist all / reftree=1 measure=pathlength file=locus%d-pldist.txt replace;\n") % g);
                }
                nexf << "  log stop;\n";
                nexf << "  quit;\n";
                nexf << "end;\n\n";
                nexf.close();

                begin_site = end_site + 1;
            } // locus loop
#else
            unsigned nbundles = G::_nbundles;
            for (unsigned b = 0; b < nbundles; b++) {
                // Save vector of log likelihoods for every locus
                vector< vector<double> > log_likes;
                smc.saveLogLikes(b, log_likes);

                // Create R file for each locus that plots RF distance (y-axis)
                // between gene tree and the true gene tree against
                // gene tree log likelihood (x-axis)
                unsigned begin_site = 1;
                unsigned end_site = 0;
                for (unsigned g = 0; g < G::_nloci; g++) {
                    ofstream Rfile(str(format("bundle%d-locus%d-scatterplot.R") % b % g));
                    Rfile << "rf <- c(";
                    for (unsigned i = 0; i < G::_nparticles; i++) {
                        auto kf_rf = Particle::calcTreeDistances(gene_ref_vect[g], smc.getGeneTreeConst(b,g,i));
                        Rfile << (i > 0 ? "," : "") << kf_rf.second;
                    }
                    Rfile << ")\n";
                    Rfile << "lnl <- c(";
                    for (unsigned i = 0; i < G::_nparticles; i++) {
                        Rfile << (i > 0 ? "," : "") << str(format("%.5f") % log_likes[g][i]);
                    }
                    Rfile << ")\n";
                    Rfile << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
                    Rfile << "setwd(cwd)\n";
                    Rfile << str(format("pdf(\"bundle%d-locus%d-scatterplot.pdf\")\n") % b % g);
                    Rfile << "plot(rf, lnl, type=\"p\")\n";
                    Rfile << "dev.off()\n";
                    Rfile.close();
            
                    // Create nexus file with gene trees along with
                    // PAUP commands for calculating distances and likelihoods
                    ofstream nexf(str(format("bundle%d-locus%d-distances.nex") % b % g));
                    nexf << "#nexus\n\n";
                    nexf << "begin paup;\n";
                    nexf << "  log start file=pauplog.txt replace;\n";
                    nexf << "  exe ../simulate/simulated.nex;\n";
                    nexf << "  set crit=like maxtrees=1000 increase=auto autoinc=1000;\n";
                    nexf << "  lset nst=1 basefreq=equal rates=equal pinvar=0 clock;\n";
                    nexf << "end;\n\n";
                    nexf << "begin trees;\n";
                    nexf << str(format("  tree gtrue = [&R] %s;\n") % gene_ref_vect[g].makeNewick(9, true));

                    // Save newick tree descriptions for every bundle's gene trees
                    // for locus g
                    vector<string> newicks;
                    smc.saveGeneTrees(b, g, newicks);
                    assert(G::_nparticles == newicks.size());

                    for (unsigned i = 0; i < G::_nparticles; i++) {
                        nexf << str(format("  tree h%d [lnL = %.9f] = [&R] %s;\n") % (i+1) % log_likes[g][i] % newicks[i]);
                    }
                    nexf << "end;\n\n";
                    nexf << "begin paup;\n";
                    nexf << "exclude all;\n";
                    end_site = begin_site + ::partition->numSitesInSubset(g) - 1;
                    nexf << str(format("include %d-%d;\n") % begin_site % end_site);
                    nexf << "lscores all / userbrlens;\n";
                    nexf << "  treedist reftree=1 measure=rfsymdiff;\n";
                    nexf << "  treedist reftree=1 measure=pathlength;\n";
                    nexf << "  log stop;\n";
                    nexf << "  quit;\n";
                    nexf << "end;\n\n";
                    nexf.close();
                    
                    begin_site = end_site + 1;
                } // locus loop
            } // bundle loop
#endif
        } // gene trees
    }
    
    inline void Proj::run() {
        output("Starting...\n", G::VSTANDARD);

        output(format("Current working directory: %s\n") % current_path(), G::VSTANDARD);
        
        try {
            rng->setSeed(_rnseed);
            
            if (_start_mode != "smc")
                throw XProj("\"sim\" is the only startmode currently");

            readData();
            
            SMC smc;
            smc.run();
            
            smc.saveSpeciesTreesToFile("species-tree.tre", /*compress*/true);
            
            for (unsigned l = 0; l < G::_nloci; l++) {
                smc.saveGeneTreeLocusToFile(str(format("gene-trees-locus-%d.tre") % l), l, /*compress*/true);
            }
            
            compareWithReferenceTrees(smc);
        }
        catch (XProj & x) {
            output(format("Proj encountered a problem:\n  %s\n") % x.what(), G::VSTANDARD);
        }
        
        output("\nFinished!\n", G::VSTANDARD);
    }

#if defined(LOG_MEMORY)
    inline void Proj::memoryReport(ofstream & memf) const {
        unsigned sizeofNodePtr = sizeof(Node *);
        unsigned sizeofNode = sizeof(Node);
        unsigned sizeofSParticle = sizeof(SParticle);
        unsigned sizeofGParticle = sizeof(GParticle);
        double Smem = sizeofSParticle*(sizeofNode*(2*G::_nspecies-1) + sizeofNodePtr*G::_nspecies);
        double Gmem = sizeofGParticle*(sizeofNode*(2*G::_ntaxa-1) + sizeofNodePtr*G::_ntaxa);
        memf << "\nMemory report:\n\n";
        memf << str(format("  Size of int:                %d\n") % sizeof(int));
        memf << str(format("  Size of char:               %d\n") % sizeof(char));
        memf << str(format("  Size of double:             %d\n") % sizeof(double));
        memf << str(format("  Size of unsigned:           %d\n") % sizeof(unsigned));
        memf << str(format("  Size of unsigned long:      %d\n") % sizeof(unsigned long));
        memf << str(format("  Size of Node:               %d\n") % sizeofNode);
        memf << str(format("  Size of SParticle:          %d\n") % sizeofSParticle);
        memf << str(format("  Size of GParticle:          %d\n") % sizeofGParticle);
        memf << str(format("  Size of Node *:             %d\n") % sizeof(Node *));
        memf << str(format("  No. species:                %d\n") % G::_nspecies);
        memf << str(format("  No. taxa:                   %d\n") % G::_ntaxa);
        memf << str(format("  No. loci:                   %d\n") % G::_nloci);
        memf << str(format("  No. species tree particles: %d\n") % G::_nbundles);
        memf << str(format("  No. gene tree particles:    %d\n") % (G::_nbundles*G::_nloci*G::_nparticles));
        memf << str(format("  SParticle memory used:      %.1f MB\n") % (Smem*G::_nbundles/1048576));
        memf << str(format("  GParticle memory used:      %.1f GB\n") % (Gmem*G::_nbundles*G::_nloci*G::_nparticles/1073741824));
        ::data->memoryReport(memf);
    }
#endif

}
