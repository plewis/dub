#pragma once

extern proj::Lot::SharedPtr  rng;

namespace proj {

    struct G {
        // Species are represented by bits so that unions can represent ancestral species
        typedef unsigned long           species_t;
                
        // Tuple: node height, left child species, right child species
        typedef tuple<double, bool, species_t, species_t>   join_info_t;
                        
        // Verbosity is a bitset enum corresponding to the kinds of output desired
        enum verbosity_t {
            VSTANDARD = 0x01,    // show standard informational and progress output
            VDEBUG    = 0x02,    // show debugging output
            VTEMP     = 0x04     // show temporary debugging output
        };
        static verbosity_t              _verbosity;
        
        static string                   _species_tree_ref_file_name;
        static string                   _gene_trees_ref_file_name;
        static map<unsigned,unsigned>   _nexus_taxon_map;

        static double                   _log_marg_like;

        static unsigned                 _step;
        static unsigned                 _bundle;
        static unsigned                 _locus;
        static unsigned                 _particle;
        
        static species_t                _species_zero;

        static unsigned                 _nstates;
        static unsigned                 _ntaxa;
        static unsigned                 _nspecies;
        static vector<string>           _taxon_names;
        static map<string, unsigned>    _taxon_to_species;
        static vector<string>           _species_names;

        static unsigned                 _nloci;
        static vector<string>           _locus_names;
        
        static unsigned                 _nsparticles;
        static unsigned                 _ngparticles;
        
        static double                   _theta;
        static double                   _lambda;
        
        static double                   _epsilon;
        static double                   _small_enough;
        static double                   _infinity;
        static double                   _negative_infinity;
        
        static void     showSettings();
        static double   extractEdgeLen(string edge_length_string);
        static int      extractNodeNumberFromName(string node_name, set<unsigned> & used);
        static double   inverseGammaVariate(double shape, double rate, Lot::SharedPtr lot);
        static void     generateUpdateSeeds(vector<unsigned> & seeds);
        static double   calcLogSum(const vector<double> & log_values);
        static string   unsignedVectToString(const vector<unsigned> & v);
        static string   memoryAddressAsString(const void * ptr);
        static void     normalizeCounts(const vector<unsigned> counts, vector<double> & probs);
        static void     normalizeRates(const vector<double> rates, vector<double> & probs);
        static unsigned multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);
        static void     setSpeciesBit(species_t & to_species, unsigned i, bool init_to_zero_first);
    };
    
    inline void G::showSettings() {
    }

    inline double G::extractEdgeLen(string edge_length_string) {
        bool success = true;
        double d = 0.0;
        try {
            d = stod(edge_length_string);
        }
        catch(invalid_argument &) {
            // edge_length_string could not be converted to a double value
            success = false;
        }

        if (!success) {
            throw XProj(str(format("%s is not interpretable as an edge length") % edge_length_string));
        }
        
        // conversion succeeded
        return d;
    }

    inline int G::extractNodeNumberFromName(string node_name, set<unsigned> & used) {
        // Attempts to convert node_name to a number and throws exception if node's
        // name cannot be converted to an integer or if node number has already
        // been encountered previously
        bool success = true;
        int x;
        try {
            x = stoi(node_name);
        }
        catch(invalid_argument &) {
            // node_name could not be converted to an integer value
            success = false;
        }

        if (success) {
            // conversion succeeded
            // attempt to insert x into the set of node numbers already used
            pair<set<unsigned>::iterator, bool> insert_result = used.insert(x);
            if (insert_result.second) {
                // insertion was made, so x has NOT already been used
                return x;
            }
            else {
                // insertion was not made, so set already contained x
                throw XProj(str(format("leaf number %d used more than once") % x));
            }
        }
        else
            throw XProj(str(format("node name (%s) not interpretable as an integer") % node_name));
    }
    
    inline double G::inverseGammaVariate(double shape, double rate, Lot::SharedPtr lot) {
        double gamma_variate = lot->gamma(shape, 1.0/rate);
        double invgamma_variate = 1.0/gamma_variate;
        return invgamma_variate;
    }
    
    inline void G::generateUpdateSeeds(vector<unsigned> & seeds) {
        unsigned psuffix = 1;
        for (auto & s : seeds) {
            s = rng->randint(1,9999) + psuffix;
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
    
    inline void G::setSpeciesBit(G::species_t & to_species, unsigned i, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (G::species_t)0;
            
        // Set ith bit in to_species
        to_species |= ((G::species_t)1 << i);
    }

}

