#pragma once

extern proj::Lot::SharedPtr rng;

namespace proj {

    struct G {
        typedef unsigned long                           species_t;
        typedef pair<unsigned, unsigned>                uint_pair_t;
        typedef tuple<double,G::species_t,G::species_t> merge_t; // <speciation height, spp 1, spp 2>
        
        //TODO: Xcode' clang allows vector<const merge_t> but not other compilers
        typedef vector<merge_t>                         merge_vect_t;
        
        enum LogCateg : unsigned {
            NONE           = 0,
            INFO           = 1 << 1,
            VERBOSE        = 1 << 2,
            DEBUGGING      = 1 << 3,
            SECONDLEVEL    = 1 << 4,
            CONDITIONALS   = 1 << 5,
            ALWAYS         = 1 << 6
        };

        //https://stackoverflow.com/questions/24297992/is-it-possible-to-make-a-scoped-enumeration-enum-class-contextually-converti
        //explicit operator bool() const {
        //    return *this != LogCateg::NONE;
        //}

        static unsigned                 _log_include;
        
        static unsigned long            _npartials_calculated;
        
#if defined(USE_HEATING)
        static double                   _heating_power;
#endif
        
#if defined(HACK_FOR_SNAKE_ATP)
        static int                      _hack_atp_index;
#endif
                
#if defined(FOSSILS)
        static vector<Fossil>           _fossils;
        static vector<TaxSet>           _taxsets;
#endif

        static vector<unsigned>         _seed_bank;
                
        static string                   _species_tree_ref_file_name;
        static string                   _gene_trees_ref_file_name;

        static bool                     _save_gene_trees;
        static bool                     _save_species_trees;
        static unsigned                 _treefile_compression;
        
        static unsigned                 _rnseed;

        static bool                     _simulating;
        static bool                     _debugging;
        
#if defined(USING_MULTITHREADING)
        static unsigned                  _nthreads;
        static mutex                     _mutex;
        static vector<pair<unsigned,unsigned> > _thread_sched;
#endif
        
        static unsigned                 _nstates;
        
        static unsigned                 _ntaxa;
        static vector<string>           _taxon_names;
        static map<string, unsigned>    _taxon_to_species;

        static unsigned                 _nspecies;
        static species_t                _species_mask;
        static vector<string>           _species_names;
        static map<unsigned,unsigned>   _nexus_taxon_map;

        static unsigned                 _nloci;
        static vector<string>           _gene_names;
        static vector<unsigned>         _nsites_per_gene;
        static map<unsigned, double>    _relrate_for_gene;
        
        static double                   _phi;
        static double                   _theta;
        static double                   _lambda;
        
        // These only used for simulating data
        static double                   _edge_rate_variance;
        static double                   _asrv_shape;
        static double                   _occupancy;
        static double                   _comphet;
                        
        static double                   _small_enough;
        static double                   _infinity;
        static double                   _negative_infinity;
        
        static unsigned                 _nparticles;
        static unsigned                 _nparticles2;
        static unsigned                 _nkept;
        static unsigned                 _nkept2;
        
        static vector<vector<double> >  _second_level_log_likes;
        
        static unsigned                 _nsubpops;
        
        static void     saveSecondLevelReport(string fn);
                
#if defined(USING_MULTITHREADING)
        static void     buildThreadSchedule(unsigned n, string entity_name);
#endif

        static string   inventName(unsigned k, bool lower_case);
        static void     showSettings();
        static double   inverseGammaVariate(double shape, double rate, Lot::SharedPtr lot);
        static void     getAllParamNames(vector<string> & names);
        static void     generateUpdateSeeds(unsigned num_seeds_needed);
#if defined(SPECIES_IN_CONF)
        static void     debugShowSpeciesDefinitions();
        static void     parseSpeciesDefinition(string s);
#endif
        static string   unsignedVectToString(const vector<unsigned> & v);
        static void     createDefaultGeneTreeNexusTaxonMap();
        static void     createDefaultSpeciesTreeNexusTaxonMap();
        static string   memoryAddressAsString(const void * ptr);
        static void     normalizeRates(const vector<double> rates, vector<double> & probs);
        static void     normalizeCounts(const vector<unsigned> counts, vector<double> & probs);
        static unsigned multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);
        static string   speciesStringRepresentation(G::species_t species);
        static set<unsigned> speciesToUnsignedSet(G::species_t species);
    };
    
    void output(format & fmt, unsigned level = G::LogCateg::INFO) {
        if ((G::_log_include & level) || (level & G::LogCateg::ALWAYS))
            cout << str(fmt);
    }

    void output(string msg, unsigned level = G::LogCateg::INFO) {
        if ((G::_log_include & level) || (level & G::LogCateg::ALWAYS))
            cout << msg;
    }

    G::LogCateg operator|(G::LogCateg a, G::LogCateg b) {
        return static_cast<G::LogCateg>(static_cast<unsigned>(a) | static_cast<unsigned>(b));
    }

    void operator|=(G::LogCateg a, G::LogCateg b) {
        a = static_cast<G::LogCateg>(static_cast<unsigned>(a) | static_cast<unsigned>(b));
    }
    
    G::LogCateg operator&(G::LogCateg a, G::LogCateg b) {
        return static_cast<G::LogCateg>(static_cast<unsigned>(a) & static_cast<unsigned>(b));
    }
    
    inline string G::inventName(unsigned k, bool lower_case) {
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

    inline void G::showSettings() {
        output(format("Speciation rate (lambda): %.9f\n") % G::_lambda);
        output(format("Coalescent parameter (theta): %.9f\n") % G::_theta);
        output(format("Number of 1st-level particles: %d\n") % G::_nparticles);
        output(format("Number of 1st-level particles kept: %d\n") % G::_nkept);
        output(format("Number of 2nd-level particles: %d\n") % G::_nparticles2);
    }

    inline double G::inverseGammaVariate(double shape, double rate, Lot::SharedPtr lot) {
        double gamma_variate = lot->gamma(shape, 1.0/rate);
        double invgamma_variate = 1.0/gamma_variate;
        return invgamma_variate;
    }
    
    inline void G::getAllParamNames(vector<string> & names) {
        // Get species tree increment names
        for (unsigned i = 0; i < G::_nspecies - 1; i++) {
            names.push_back("sincr." + to_string(i));
        }
        
        // Get gene tree increment names
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned i = 0; i < G::_ntaxa - 1; i++) {
                names.push_back(G::_gene_names[g] + "." + to_string(i));
            }
        }
    }

    inline void G::generateUpdateSeeds(unsigned num_seeds_needed) {
        _seed_bank.resize(num_seeds_needed);
        unsigned maxseed = num_seeds_needed*100 - 1;
        for_each(_seed_bank.begin(), _seed_bank.end(), [maxseed](unsigned & x){x = ::rng->randint(1,maxseed);});
    }

#if defined(SPECIES_IN_CONF)
    inline void G::debugShowSpeciesDefinitions() {
        for (unsigned i = 0; i < G::_nspecies; i++) {
            string s = G::_species_names[i];
            output(format("Species \"%s\":\n") % s);
            for (auto ts : G::_taxon_to_species) {
                if (ts.second == i) {
                    output(format("  \"%s\"\n") % ts.first);
                }
            }
        }
    }
#endif

#if defined(SPECIES_IN_CONF)
    inline void G::parseSpeciesDefinition(string s) {
        // Given these definitions in the conf file:
        //   species = A: a^A, b^A, c^A
        //   species = B: d^B, e^B, f^B
        //   species = C: g^C, h^C, i^C
        //   species = D: j^D, k^D, l^D
        //   species = E: m^E, n^E, o^E
        // This would be the result:
        //   G::_nspecies = 5
        //   G::_species_names = ["A", "B", "C", "D", "E"]
        //   G::_ntaxa = 15
        //   G::_taxon_names = [
        //      "a^A", "b^A", "c^A",
        //      "d^B", "e^B", "f^B",
        //      "g^C", "h^C", "i^C",
        //      "j^D", "k^D", "l^D",
        //      "m^E", "n^E", "o^E"]
        //   G::_taxon_to_species["a^A"] = 0
        //   G::_taxon_to_species["b^A"] = 0
        //   G::_taxon_to_species["c^A"] = 0
        //   G::_taxon_to_species["d^B"] = 1
        //   G::_taxon_to_species["e^B"] = 1
        //   G::_taxon_to_species["f^B"] = 1
        //   G::_taxon_to_species["g^C"] = 2
        //   G::_taxon_to_species["h^C"] = 2
        //   G::_taxon_to_species["i^C"] = 2
        //   G::_taxon_to_species["j^D"] = 3
        //   G::_taxon_to_species["k^D"] = 3
        //   G::_taxon_to_species["l^D"] = 3
        //   G::_taxon_to_species["m^E"] = 4
        //   G::_taxon_to_species["n^E"] = 4
        //   G::_taxon_to_species["o^E"] = 4
        //   G::_taxon_to_species["A"]   = 0
        //   G::_taxon_to_species["B"]   = 1
        //   G::_taxon_to_species["C"]   = 2
        //   G::_taxon_to_species["D"]   = 3
        //   G::_taxon_to_species["E"]   = 4
        
        vector<string> v;
        
        // First separate part before colon (stored in v[0])
        // from the part after colon (stored in v[1])
        split(v, s, boost::is_any_of(":"));
        if (v.size() != 2)
            throw XProj("Expecting exactly one colon in species definition");

        string species_name = v[0];
        string taxon_list = v[1];

        // Now separate the part after the colon at commas to
        // yield the taxa that are in that species
        boost::trim(taxon_list);
        split(v, taxon_list, boost::is_any_of(","));
        for_each(v.begin(), v.end(), [](string & s){boost::trim(s);});
        
        //output(format("Species \"%s\":\n") % species_name, LogCateg::ALWAYS);
        //for (string s : v) {
        //    output(format("   Taxon \"%s\"\n") % s, LogCateg::ALWAYS);
        //}
        
        G::_species_names.push_back(species_name);
        G::_nspecies = (unsigned)G::_species_names.size();
        unsigned i = (unsigned)(G::_nspecies - 1);
        G::_taxon_to_species[species_name] = i;
        for (auto t : v) {
            G::_taxon_names.push_back(t);
            G::_taxon_to_species[t] = i;
            G::_ntaxa++;
        }
    }
#endif

    inline string G::unsignedVectToString(const vector<unsigned> & v) {
        ostringstream oss;
        
        // Create string with values separated by commas
        copy(v.begin(), v.end(), ostream_iterator<unsigned>(oss, ","));

        // Remove the trailing comma
        oss.seekp(-1, oss.end);
        oss.put('\0');
        
        // Return string
        return oss.str();
    }
    
    inline void G::createDefaultGeneTreeNexusTaxonMap() {
        // Build default _nexus_taxon_map used by buildFromNewick
        // that assumes taxon indices in the newick tree description are in
        // the same order as _taxon_names
        _nexus_taxon_map.clear();
        for (unsigned i = 0; i < G::_ntaxa; ++i) {
            _nexus_taxon_map[i+1] = i;
        }
    }
    
    inline void G::createDefaultSpeciesTreeNexusTaxonMap() {
        // Build default _nexus_taxon_map used by buildFromNewick
        // that assumes taxon indices in the newick tree description are in
        // the same order as _species_names
        _nexus_taxon_map.clear();
        for (unsigned i = 0; i < G::_nspecies; ++i) {
            _nexus_taxon_map[i+1] = i;
        }
    }

    inline string G::memoryAddressAsString(const void * ptr) {
        ostringstream memory_address;
        memory_address << ptr;
        return memory_address.str();
    }
    
    inline void G::normalizeCounts(const vector<unsigned> counts, vector<double> & probs) {
        // Determine sum of counts
        double total_count = accumulate(counts.begin(), counts.end(), 0.0);
        
        // Normalize counts to create a discrete probability distribution
        transform(counts.begin(), counts.end(), probs.begin(), [total_count](unsigned count){return (double)count/total_count;});
        assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
    }
    
    inline void G::normalizeRates(const vector<double> rates, vector<double> & probs) {
        // Determine sum of rates
        double total_rate = accumulate(rates.begin(), rates.end(), 0.0);
        
        // Normalize rates to create a discrete probability distribution
        transform(rates.begin(), rates.end(), probs.begin(), [total_rate](unsigned rate){return (double)rate/total_rate;});
        assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
    }
    
    inline unsigned G::multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs) {
        // Compute cumulative probababilities
        vector<double> cumprobs(probs.size());
        partial_sum(probs.begin(), probs.end(), cumprobs.begin());
        assert(fabs(*(cumprobs.rbegin()) - 1.0) < 0.0001);

        // Draw a Uniform(0,1) random deviate
        double u = lot->uniform();

        // Find first element in cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6 <-- u
        //
        // cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(cumprobs.begin(), cumprobs.end(), [u](double cumpr){return cumpr > u;});
        if (it == cumprobs.end()) {
            double last_cumprob = *(cumprobs.rbegin());
            throw XProj(format("G::multinomialDraw failed: u = %.9f, last cumprob = %.9f") % u % last_cumprob);
        }

        auto d = distance(cumprobs.begin(), it);
        assert(d >= 0);
        assert(d < probs.size());
        return (unsigned)d;
    }
    
    inline string G::speciesStringRepresentation(G::species_t species) {
        species_t species_copy = species;
        unsigned bits_avail = (unsigned)sizeof(species_t);
        string s;
        for (unsigned i = 0; i < bits_avail; ++i) {
            species_t bitmask = ((species_t)1 << i);
            bool bit_is_set = ((species_copy & bitmask) > (species_t)0);
            if (bit_is_set) {
                // Add species i to the string
                s += to_string(i);
                
                // Zero that bit so we know when we are done
                species_copy &= ~bitmask;
            }
            if (!species_copy) {
                // If species_copy is zero, there are no more bits set
                break;
            }
        }
        return s;
    }
    
    inline set<unsigned> G::speciesToUnsignedSet(G::species_t species) {
        species_t species_copy = species;
        unsigned bits_avail = (unsigned)sizeof(species_t);
        set<unsigned> s;
        for (unsigned i = 0; i < bits_avail; ++i) {
            species_t bitmask = ((species_t)1 << i);
            bool bit_is_set = ((species_copy & bitmask) > (species_t)0);
            if (bit_is_set) {
                // Add species i to the string
                s.insert(i);
                
                // Zero that bit so we know when we are done
                species_copy &= ~bitmask;
            }
            if (!species_copy) {
                // If species_copy is zero, there are no more bits set
                break;
            }
        }
        return s;
    }
    
    inline void G::saveSecondLevelReport(string fn) {
        unsigned n = (unsigned)G::_second_level_log_likes.size();
        assert(n > 0);
        assert(n == G::_nkept);
        
        // This vector will hold p(D|G) p(G) for each 1st-level kept particle
        vector<double> log_values(n);

        ofstream outf(fn);
        
        // Save headers
        outf << "log-marg-like\t";
        for (unsigned i = 0; i < G::_nloci; i++) {
            outf << str(format("locus-%d\t") % (i+1));
        }
        outf << "\n";
        
        // Save likelihoods
        unsigned i = 0;
        for (auto & v : G::_second_level_log_likes) {
            double sum_log_likes = 0.0;
            for (auto logL : v) {
                outf << setprecision(9) << logL << "\t";
                sum_log_likes += logL;
            }
            log_values[i++] = sum_log_likes;
            outf << "\n";
        }
        
        outf.close();
        
        // Calculate log sum of elements in log_values
        double max_logv = *max_element(log_values.begin(), log_values.end());
        
        double factored_sum = 0.0;
        for (auto & logv : log_values) {
            factored_sum += exp(logv - max_logv);
        }
        double log_sum_values = max_logv + log(factored_sum);
        
        double log_marg_like = log_sum_values - log(n);
        
        output(format("%20.9f Log marginal likelihood (2nd-level estimate)") % log_marg_like, LogCateg::INFO);
    }
    
#if defined(USING_MULTITHREADING)
    inline void G::buildThreadSchedule(unsigned n, string entity_name) {
        // Populate _thread_sched for the specified number of entities
        assert(n > 0);
        G::_thread_sched.clear();
#if 1
        // Example: n = 11, _nthreads = 3
        // entities_per_thread = 3, remainder = 2
        // thread 0: begin = 0, end = 4  (= 0 + 3 + 1) <-- 0, 1, 2, 3
        // thread 1: begin = 4, end = 8  (= 4 + 3 + 1) <-- 4, 5, 6, 7
        // thread 2: begin = 8, end = 11 (= 8 + 3 + 0) <-- 8, 9, 10
        unsigned entities_per_thread = (unsigned)floor(1.0*n/G::_nthreads);
        unsigned remainder = n - G::_nthreads*entities_per_thread;
        unsigned begin = 0;
        for (unsigned i = 0; i < G::_nthreads; i++) {
            unsigned end = begin + entities_per_thread;
            if (remainder > 0) {
                end++;
                remainder--;
            }
            G::_thread_sched.push_back(make_pair(begin,end));
            begin = end;
        }
#else
//        // Old version: used when there was a list of particles each of which
//        // had a count associated with it of identical particles.
//        //
//        // Suppose we are using 3 threads and there are 10 particles with these counts:
//        // _nparticles = 100 = 20 + 10 + 5 + 1 + 4 + 20 + 23 + 12 + 1 + 4 = 100
//        // prefix_sum        = 20   30  35  36  40   60   83   95  96  100
//        // thread_index      =  0    0   1   1   1    1    2    2   2    2
//        // max_count = 100/3 = 33
//
//        // Create thread schedule using the prefix-sum algorithm
//        thread_schedule.clear();
//        unsigned current_thread = 0;
//        pair<unsigned, unsigned> begin_end = make_pair(0,0);
//        unsigned prefix_sum = 0;
//        unsigned max_count = (unsigned)floor(1.0*_nparticles/G::_nthreads);
//        
//        vector<unsigned> freqs(SMCGlobal::_nthreads, 0);
//        for (auto & p : _particle_list) {
//            // Add particle count to prefix sum
//            unsigned count = p.getCount();
//            
//            // If count greater than max_count, we cannot create
//            // a schedule because the schedule cannot split up
//            // the count of a single particle
//            if (count > max_count) {
//                thread_schedule.clear();
//                return 0.0;
//            }
//            
//            p.setBeginIndex(prefix_sum);
//            prefix_sum += count;
//                
//            // Calculate thread index
//            unsigned thread_index = (unsigned)floor(1.0*G::_nthreads*(prefix_sum - 1)/_nparticles);
//            
//            if (thread_index > current_thread) {
//                // Add thread to the schedule
//                thread_schedule.push_back(begin_end);
//                current_thread++;
//                
//                // Start work on next thread
//                begin_end.first  = begin_end.second;
//                begin_end.second = begin_end.first + 1;
//            }
//            else {
//                begin_end.second += 1;
//            }
//            freqs[current_thread] += count;
//        }
//        thread_schedule.push_back(begin_end);
//
//        // Calculate entropy, max_entropy, and percentage
//        double max_entropy = log(SMCGlobal::_nthreads);
//        double entropy = 0.0;
//        assert(_nparticles == accumulate(freqs.begin(), freqs.end(), 0));
//        for_each(freqs.begin(), freqs.end(), [&entropy](unsigned f){entropy -= 1.0*f*log(f);});
//        entropy /= _nparticles;
//        entropy += log(_nparticles);
//        double percentage = 100.0*entropy/max_entropy;
//        assert(!isnan(percentage));
//        
//        return percentage;
#endif

        output(format("\nThread schedule (n = %d %s%s spread over %d thread%s):\n")
            % n
            % entity_name
            % (n > 1 ? "s" : "")
            % _thread_sched.size()
            % (G::_nthreads > 1 ? "s" : "")
        );
        unsigned t = 0;
        output(format("%24s %12s %12s\n") % entity_name % "begin" % "end", G::LogCateg::INFO);
        for (auto be : _thread_sched) {
            output(format("%24d %12d %12d\n") % (t++) % be.first % (be.second - 1), G::LogCateg::INFO);
        }
    }
#endif
}

