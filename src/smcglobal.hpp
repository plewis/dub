#pragma once

extern proj::Lot rng;

namespace proj {

    struct G {
        typedef unsigned long           species_t;
        
        static string                   _species_tree_ref_file_name;
        static string                   _gene_trees_ref_file_name;
        
        static bool                     _debugging;

        static unsigned                 _nthreads;
        
        static unsigned                 _verbosity;

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
        
        static double                   _phi;
        static double                   _theta;
        static double                   _lambda;
        
        static double                   _theta_prior_mean;
        static double                   _lambda_prior_mean;
        
        static bool                     _update_theta;
        static bool                     _update_lambda;
        
        static bool                     _prior_post;
        
        static double                   _small_enough;
        static double                   _infinity;
        static double                   _negative_infinity;
        
        static unsigned                 _nparticles;
        static unsigned                 _nparticles2;
                        
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
    };
    
    inline void G::getAllParamNames(vector<string> & names) {
        // Get species tree increment names
        for (unsigned i = 0; i < G::_nspecies - 1; i++) {
            names.push_back("sincr." + to_string(i));
        }
        
        // Get gene tree increment names
        for (unsigned g = 0; g < G::_ngenes; g++) {
            for (unsigned i = 0; i < G::_ntaxa - 1; i++) {
                names.push_back(G::_gene_names[g] + "." + to_string(i));
            }
        }
    }

    inline void G::generateUpdateSeeds(vector<unsigned> & seeds) {
        unsigned psuffix = 1;
        for (auto & s : seeds) {
            s = rng.randint(1,9999) + psuffix;
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
    
}

