#pragma once

extern proj::Lot rng;

namespace proj {

    struct SMCGlobal {
        typedef unsigned long                       species_t;
        typedef tuple<unsigned, unsigned, species_t>  species_tuple_t;

        static unsigned                 _nthreads;
#if defined(USING_MULTITHREADING)
        static mutex                    _mutex;
        static mutex                    _debug_mutex;
        //static vector<unsigned>         _thread_first_gene;
        //static vector<unsigned>         _thread_last_gene;
#endif
        
        static unsigned                 _verbosity;
        static bool                     _debugging;

        static unsigned                 _nstates;
        
        static unsigned                 _ntaxa;
        static vector<string>           _taxon_names;
        static map<string, unsigned>    _taxon_to_species;

        static unsigned                 _nspecies;
        static species_t                _species_mask;
        static vector<string>           _species_names;
        static map<unsigned,unsigned>   _nexus_taxon_map;

        static unsigned                 _ngenes;
        static vector<string>           _gene_names;
        static vector<unsigned>         _nsites_per_gene;
        static map<unsigned, double>    _relrate_for_gene;
        
        static double                   _theta;
        static double                   _lambda;
        
        static double                   _theta_prior_mean;
        static double                   _lambda_prior_mean;
        
        static bool                     _update_theta;
        static bool                     _update_lambda;
        
        static bool                     _prior_prior;
        
        static double                   _small_enough;
        static double                   _infinity;
        static double                   _negative_infinity;
                
        static vector<double>           _cumprobs;  // workspace used by multinomialDraw
        
        static double   calcLogSum(const vector<double> & log_values);
        static string   unsignedVectToString(const vector<unsigned> & v);
        static void     createDefaultGeneTreeNexusTaxonMap();
        static void     createDefaultSpeciesTreeNexusTaxonMap();
        static string   memoryAddressAsString(const void * ptr);
        static void     normalizeRates(const vector<double> rates, vector<double> & probs);
        static void     normalizeCounts(const vector<unsigned> counts, vector<double> & probs);
        static unsigned multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);
        static string   speciesStringRepresentation(SMCGlobal::species_t species);

    };

    inline double SMCGlobal::calcLogSum(const vector<double> & log_values) {
        double max_logv = *max_element(log_values.begin(), log_values.end());
        
        double factored_sum = 0.0;
        for (auto & logv : log_values) {
            factored_sum += exp(logv - max_logv);
        }
        double log_sum_values = max_logv + log(factored_sum);
        return log_sum_values;
    }
    
    inline string SMCGlobal::unsignedVectToString(const vector<unsigned> & v) {
        ostringstream oss;
        
        // Create string with values separated by commas
        copy(v.begin(), v.end(), ostream_iterator<unsigned>(oss, ","));

        // Remove the trailing comma
        oss.seekp(-1, oss.end);
        oss.put('\0');
        
        // Return string
        return oss.str();
    }
    
    inline void SMCGlobal::createDefaultGeneTreeNexusTaxonMap() {
        // Build default _nexus_taxon_map used by buildFromNewick
        // that assumes taxon indices in the newick tree description are in
        // the same order as _taxon_names
        _nexus_taxon_map.clear();
        for (unsigned i = 0; i < SMCGlobal::_ntaxa; ++i) {
            _nexus_taxon_map[i+1] = i;
        }
    }
    
    inline void SMCGlobal::createDefaultSpeciesTreeNexusTaxonMap() {
        // Build default _nexus_taxon_map used by buildFromNewick
        // that assumes taxon indices in the newick tree description are in
        // the same order as _species_names
        _nexus_taxon_map.clear();
        for (unsigned i = 0; i < SMCGlobal::_nspecies; ++i) {
            _nexus_taxon_map[i+1] = i;
        }
    }

    inline string SMCGlobal::memoryAddressAsString(const void * ptr) {
        ostringstream memory_address;
        memory_address << ptr;
        return memory_address.str();
    }
    
    inline void SMCGlobal::normalizeCounts(const vector<unsigned> counts, vector<double> & probs) {
        // Determine sum of counts
        double total_count = accumulate(counts.begin(), counts.end(), 0.0);
        
        // Normalize counts to create a discrete probability distribution
        transform(counts.begin(), counts.end(), probs.begin(), [total_count](unsigned count){return (double)count/total_count;});
        assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
    }
    
    inline void SMCGlobal::normalizeRates(const vector<double> rates, vector<double> & probs) {
        // Determine sum of rates
        double total_rate = accumulate(rates.begin(), rates.end(), 0.0);
        
        // Normalize rates to create a discrete probability distribution
        transform(rates.begin(), rates.end(), probs.begin(), [total_rate](unsigned rate){return (double)rate/total_rate;});
        assert(fabs(accumulate(probs.begin(), probs.end(), 0.0) - 1.0) < 0.0001);
    }
    
    inline unsigned SMCGlobal::multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs) {
//#if defined(USING_MULTITHREADING)
//        lock_guard<mutex> guard(SMCGlobal::_mutex);
//#endif
        
        // Compute cumulative probababilities
        _cumprobs.resize(probs.size());
        partial_sum(probs.begin(), probs.end(), _cumprobs.begin());

        // Draw a Uniform(0,1) random deviate
        double u = lot->uniform();

        // Find first element in _cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6                <-- u
        //
        // _cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(_cumprobs.begin(), _cumprobs.end(), [u](double cumpr){return cumpr > u;});
        if (it == _cumprobs.end()) {
            throw XProj(format("SMCGlobal::multinomialDraw failed: u = %.9f, *_cumprobs.rbegin()") % u % *(_cumprobs.rbegin()));
        }

        return (unsigned)distance(_cumprobs.begin(), it);
    }
    
    inline string SMCGlobal::speciesStringRepresentation(SMCGlobal::species_t species) {
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
    
}

