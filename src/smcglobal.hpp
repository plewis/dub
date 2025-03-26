#pragma once

extern proj::Lot::SharedPtr rng;

namespace proj {

    struct G {
        typedef unsigned long           species_t;
        
        enum LogCateg : unsigned {
            NONE           = 0,
            INFO           = 1 << 1,
            VERBOSE        = 1 << 2,
            DEBUGGING      = 1 << 3,
            SECONDLEVEL    = 1 << 4
        };

        //https://stackoverflow.com/questions/24297992/is-it-possible-to-make-a-scoped-enumeration-enum-class-contextually-converti
        //explicit operator bool() const {
        //    return *this != LogCateg::NONE;
        //}

        static unsigned                 _log_include;
        
        static unsigned long            _npartials_calculated;
        
#if defined(HACK_FOR_SNAKE_ATP)
        static int                      _hack_atp_index;
#endif
                
#if defined(FOSSILS)
        static vector<Fossil>           _fossils;
        static vector<TaxSet>           _taxsets;
#endif
                
        static string                   _species_tree_ref_file_name;
        static string                   _gene_trees_ref_file_name;
        
        static unsigned                 _treefile_compression;
        
        static unsigned                 _rnseed;

        static bool                     _simulating;
        static bool                     _debugging;
        
        static unsigned                 _nthreads;
        
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
        
        static unsigned                 _nsubpops;

#if defined(UPGMA_WEIGHTS)
        static vector<vector<double> > _dmatrix;
        static vector<Split>           _dmatrix_rows;
#endif

        static string   inventName(unsigned k, bool lower_case);
        static void     showSettings();
        static double   inverseGammaVariate(double shape, double rate, Lot::SharedPtr lot);
        static void     getAllParamNames(vector<string> & names);
        static void     generateUpdateSeeds(vector<unsigned> & seeds);
        static double   calcLogSum(const vector<double> & log_values);
        static string   unsignedVectToString(const vector<unsigned> & v);
        static void     createDefaultGeneTreeNexusTaxonMap();
        static void     createDefaultSpeciesTreeNexusTaxonMap();
        static string   memoryAddressAsString(const void * ptr);
        static void     normalizeRates(const vector<double> rates, vector<double> & probs);
        static void     normalizeCounts(const vector<unsigned> counts, vector<double> & probs);
        static unsigned multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);
        static string   speciesStringRepresentation(G::species_t species);
        static set<unsigned> speciesToUnsignedSet(G::species_t species);
        
#if defined(UPGMA_WEIGHTS)
        static void     mergeDMatrixPair(vector<Split> & dmatrows, vector<double> & dmatrix, Split & s1, Split & s2);
        //static void     copyDMatrixTo(vector<Split> & dmatrows, vector<vector<double> > & dmatrix);
        static void     debugShowDistanceMatrix(const vector<Split> rows, const vector<double> & d, unsigned locus);
#endif
    };
    
    void output(format & fmt, unsigned level = G::LogCateg::INFO) {
        if (G::_log_include & level)
            cout << str(fmt);
    }

    void output(string msg, unsigned level = G::LogCateg::INFO) {
        if (G::_log_include & level)
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

    inline void G::generateUpdateSeeds(vector<unsigned> & seeds) {
        unsigned psuffix = 1;
        for (auto & s : seeds) {
            s = ::rng->randint(1,9999) + psuffix;
            psuffix += 2;    // pure superstition (I always use odd seeds)
        }
    }

    inline double G::calcLogSum(const vector<double> & log_values) {
        double max_logv = *max_element(log_values.begin(), log_values.end());
        
        double factored_sum = 0.0;
        for (auto & logv : log_values) {
            factored_sum += exp(logv - max_logv);
        }
        double log_sum_values = max_logv + log(factored_sum);
        return log_sum_values;
    }

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
    
#if defined(UPGMA_WEIGHTS)
    void G::mergeDMatrixPair(vector<Split> & dmatrows, vector<double> & dmatrix, Split & s1, Split & s2) {
        // dmatrix = [d10, d20, d21, d31, d31, d32]
        //   index =   0    1    2    3    4    5
        //
        // where the actual 2-dimensional matrix looks like this:
        //
        //   d00  d01  d02  d03
        //   d10  d11  d12  d13
        //   d20  d21  d22  d23
        //   d30  d31  d32  d33
        //
        // Only the lower triangle (not including diagonals) is used,
        // and elements are storwd by row.
        //
        // The index k of the (i,j)th element, where i > j, can be
        // obtained as follows:
        //
        // k = i*(i-1)/2 + j
        //
        // For example, i = 3, j = 1, k = 3*2/2 + 1 = 4
        //
        
        unsigned n = (unsigned)dmatrows.size();
        
        // Find index of s1 in dmatrows
        auto it1 = find(dmatrows.begin(), dmatrows.end(), s1);
        assert(it1 != dmatrows.end());
        unsigned i = (unsigned)distance(dmatrows.begin(), it1);
        
        // Find index of s2 in dmatrows
        auto it2 = find(dmatrows.begin(), dmatrows.end(), s2);
        assert(it2 != dmatrows.end());
        unsigned j = (unsigned)distance(dmatrows.begin(), it2);
        
        // Update distance matrix
        unsigned ij = (i > j) ? (i*(i-1)/2 + j) : (j*(j-1)/2 + i);
        dmatrix[ij] = G::_infinity;
        for (unsigned k = 0; k < n; k++) {
            if (k != i && k != j) {
                unsigned ik = (i > k) ? (i*(i-1)/2 + k) : (k*(k-1)/2 + i);
                unsigned jk = (j > k) ? (j*(j-1)/2 + k) : (k*(k-1)/2 + j);
                double a = dmatrix[ik];
                double b = dmatrix[jk];
                
                // Put average in cell with the smaller index
                if (i < j) {
                    dmatrix[ik] = 0.5*(a + b);
                    dmatrix[jk] = G::_infinity;
                }
                else {
                    dmatrix[ik] = G::_infinity;
                    dmatrix[jk] = 0.5*(a + b);
                }
            }
        }

        // Example (5x5 distance matrix):
        //   dmatrix [d10, d20, d21, d30, d31, d32, d40, d41, d42, d43]
        //   number of elements: 5*4/2 = 10
        //
        // Original (5x5):
        //
        //              0       1       2       3       4
        //            ----*   ---*-   --*--   -*---   *----
        //          +-------+-------+-------+-------+-------+
        // 0 ----*  |       |       |       |       |       |
        //          +-------+-------+-------+-------+-------+
        // 1 ---*-  |  d10  |       |       |       |       |
        //          +-------+-------+-------+-------+-------+
        // 2 --*--  |  d20  |  d21  |       |       |       |
        //          +-------+-------+-------+-------+-------+
        // 3 -*---  |  d30  |  d31  |  d32  |       |       |
        //          +-------+-------+-------+-------+-------+
        // 4 *----  |  d40  |  d41  |  d42  |  d43  |       |
        //          +-------+-------+-------+-------+-------+
        //
        // Intermediate during merging of 0 and 2 (4x4)
        //
        // Set d20 = inf and make these changes as well
        //
        //      k     ik     jk  dmatrix[ik]   dmatrix[jk]
        //  -----  -----  -----  -----------  ------------
        //  i = 0
        //      1    d10    d21  (d10+d21)/2           inf
        //  j = 2
        //      3    d30    d32  (d30+d32)/2           inf
        //      4    d40    d42  (d40+d42)/2           inf
        //
        //                  0           1       2       3       4
        //                ----*       ---*-   --*--   -*---   *----
        //          +---------------+-------+-------+-------+-------+
        // 0 ----*  |               |       |       |       |       | <-- kept
        //          +---------------+-------+-------+-------+-------+
        // 1 ---*-  |  (d10+d21)/2  |       |       |       |       |
        //          +---------------+-------+-------+-------+-------+
        // 2 --*--  |  inf          |  inf  |       |       |       | <-- elim
        //          +---------------+-------+-------+-------+-------+
        // 3 -*---  |  (d30+d32)/2  |  d31  |  inf  |       |       |
        //          +---------------+-------+-------+-------+-------+
        // 4 *----  |  (d40+d42)/2  |  d41  |  inf  |  d43  |       |
        //          +---------------+-------+-------+-------+-------+
        //                 kept                elim
        //
        // New matrix after eliminating row 2 and column 2 (3x3):
        //
        //                  0           1       2       3
        //                --*-*       ---*-   -*---   *----
        //          +---------------+-------+-------+-------+
        // 0 --*-*  |               |       |       |       |
        //          +---------------+-------+-------+-------+
        // 1 ---*-  |  (d10+d21)/2  |       |       |       |
        //          +---------------+-------+-------+-------+
        // 2 -*---  |  (d30+d32)/2  |  d31  |       |       |
        //          +---------------+-------+-------+-------+
        // 3 *----  |  (d40+d42)/2  |  d41  |  d43  |       |
        //          +---------------+-------+-------+-------+

        // Build new distance matrix
        unsigned n2 = n - 1;
        unsigned dim2 = n2*(n2-1)/2;
        vector<double> dmat2(dim2, G::_infinity);
        vector<Split> drows2(n2);
        
        // Eliminated row and column is the larger of i and j
        unsigned elim = i > j ? i : j;
        unsigned kept = i > j ? j : i;
        unsigned newrow = 0;
        unsigned newcol = 0;

        for (unsigned oldrow = 0; oldrow < n; oldrow++) {
            if (oldrow != elim) {
                if (oldrow == kept) {
                    // Set new row split to union of splits for rows i and j
                    drows2[newrow] = dmatrows[i] + dmatrows[j];
                }
                else {
                    drows2[newrow] = dmatrows[oldrow];
                }
                
                // Copy matrix row
                for (unsigned oldcol = 0; oldcol < oldrow; oldcol++) {
                    if (oldcol != elim) {
                        unsigned oldk = oldrow*(oldrow-1)/2 + oldcol;
                        unsigned newk = newrow*(newrow-1)/2 + newcol;
                        dmat2[newk] = dmatrix[oldk];
                        newcol++;
                    }
                }
                newrow++;
                newcol = 0;
            }
        }

        // Replace distance matrix with new one
        dmatrix = dmat2;
        dmatrows = drows2;
    }
#endif

//#if defined(UPGMA_WEIGHTS)
//    void G::copyDMatrixTo(vector<Split> & dmatrows, vector<vector<double> > & dmatrix) {
//        dmatrows.resize(G::_dmatrix_rows.size());
//        dmatrix.resize(G::_dmatrix.size());
//        for (unsigned i = 0; i < G::_dmatrix_rows.size(); i++) {
//            dmatrows[i] = G::_dmatrix_rows[i];
//        }
//        for (unsigned i = 0; i < G::_dmatrix.size(); i++) {
//            dmatrix[i] = G::_dmatrix[i];
//        }
//    }
//#endif
 
#if defined(UPGMA_WEIGHTS)
    inline void G::debugShowDistanceMatrix(const vector<Split> rows, const vector<double> & d, unsigned locus) {
        // d is a 1-dimensional vector that stores the lower triangle of a square matrix
        // (not including diagonals) in row order
        //
        // For example, for a 4x4 matrix (- means non-applicable):
        //
        //       1  2  4  8  <-- rows and columns are Splits representing taxa 1, 2, 3, and 4
        //     +-----------
        //  1  | -  -  -  -
        //  2  | 0  -  -  -
        //  4  | 1  2  -  -
        //  8  | 3  4  5  -
        //
        // For this example, d = {0, 1, 2, 3, 4 ,5}
        //
        // See this explanation for how to index d:
        //   https://math.stackexchange.com/questions/646117/how-to-find-a-function-mapping-matrix-indices
        //
        // In short, d[k] is the (i,j)th element, where k = i(i-1)/2 + j
        //       i   j   k = i*(i-1)/2 + j
        //       1   0   0 = 1*0/2 + 0
        //       2   0   1 = 2*1/2 + 0
        //       2   1   2 = 2*1/2 + 1
        //       3   0   3 = 3*2/2 + 0
        //       3   1   4 = 3*2/2 + 1
        //       3   2   5 = 3*2/2 + 2
        //
        // Number of elements in d is n(n-1)/2
        // Solving for n, and letting x = d.size(),
        //  x = n(n-1)/2
        //  2x = n^2 - n
        //  0 = a n^2 + b n + c, where a = 1, b = -1, c = -2x
        //  n = (-b += sqrt(b^2 - 4ac))/(2a)
        //    = (1 + sqrt(1 + 8x))/2
        double x = (double)d.size();
        double dbln = (1.0 + sqrt(1.0 + 8.0*x))/2.0;
        unsigned n = (unsigned)dbln;
        
        output(format("\nJC distance matrix (%d x %d) for locus %d") % n % n % locus, G::LogCateg::DEBUGGING);

        // Column headers
        output(format("%12d") % " ", G::LogCateg::DEBUGGING);
        for (unsigned j = 0; j < n; j++) {
            string s = rows[j].createPatternRepresentation();
            output(format("%12s") % s, G::LogCateg::DEBUGGING);
        }
        output("\n", G::LogCateg::DEBUGGING);
        
        unsigned k = 0;
        for (unsigned i = 0; i < n; i++) {
            string s = rows[i].createPatternRepresentation();
            output(format("%12s") % s, G::LogCateg::DEBUGGING);
            for (unsigned j = 0; j < n; j++) {
                if (j < i) {
                    double v = d[k++];
                    if (v == G::_infinity)
                        output("         inf", G::LogCateg::DEBUGGING);
                    else
                        output(format("%12.5f") % v, G::LogCateg::DEBUGGING);
                }
                else {
                    output("         inf", G::LogCateg::DEBUGGING);
                }
            }
            output("\n", G::LogCateg::DEBUGGING);
        }
        output("\n", G::LogCateg::DEBUGGING);
    }
#endif

}

