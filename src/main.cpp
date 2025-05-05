//#define NDEBUG
#include <cassert>
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <unordered_map>
#include <memory>   // shared_ptr
#include <new>      // used by Mallocator: bad_alloc, bad_array_new_length
#include <numeric>
#include <queue>
#include <regex>
#include <set>
#include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/program_options.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "ncl/nxsmultiformat.h"

using namespace std;
using boost::format;

// Include this header to enable macros that Valgrind recognizes
// but which are ignored if Valgrind is not being used
#include "valgrind.h"

#include "conditionals.hpp"

#if defined(FOSSILS)
#   include <codecvt>
#   include "fossil.hpp"
#   include "taxset.hpp"
#endif

#include "xproj.hpp"
#include "lot.hpp"
#include "split.hpp"
#include "smcglobal.hpp"
#include "stopwatch.hpp"
#include "genetic-code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "data.hpp"
#include "partial_store.hpp"
#include "node.hpp"
#include "forest.hpp"
#include "species-forest.hpp"
#include "gene-forest.hpp"

#if defined(LAZY_COPYING)
#   include "gene-forest-extension.hpp"
#endif

#include "gene-forest-func.hpp"
#include "particle.hpp"
#include "smc.hpp"
#include "proj.hpp"

using namespace proj;

Lot::SharedPtr                      rng(new Lot());
StopWatch                           stopwatch;
PartialStore                        ps;

#if defined(HACK_FOR_SNAKE_ATP)
int                                 G::_hack_atp_index              = -1;
#endif

unsigned                            G::_log_include                 = G::LogCateg::NONE;

unsigned long                       G::_npartials_calculated        = 0;

#if defined(USE_HEATING)
double                              G::_heating_power               = 1.0;
#endif

string                              G::_species_tree_ref_file_name  = "";
string                              G::_gene_trees_ref_file_name    = "";

vector<unsigned>                    G::_seed_bank;
unsigned                            G::_rnseed                      = 1;
bool                                G::_debugging                   = false;
bool                                G::_simulating                  = false;

unsigned                            G::_nthreads                    = 1;

#if defined(USING_MULTITHREADING)
mutex                               G::_mutex;
#endif

unsigned                            G::_treefile_compression        = 0;

unsigned                            G::_nstates                     = 4;

unsigned                            G::_ntaxa                       = 0;
vector<string>                      G::_taxon_names;
map<string, unsigned>               G::_taxon_to_species;

unsigned                            G::_nspecies                    = 0;
G::species_t                        G::_species_mask                = (G::species_t)0;
vector<string>                      G::_species_names;
map<unsigned,unsigned>              G::_nexus_taxon_map;

unsigned                            G::_nloci                      = 0;
vector<string>                      G::_gene_names;
vector<unsigned>                    G::_nsites_per_gene;
map<unsigned, double>               G::_relrate_for_gene;

double                              G::_phi                         = 1.0;
double                              G::_theta                       = 0.05;
double                              G::_lambda                      = 1.0;

double                              G::_edge_rate_variance          = 0.0;
double                              G::_occupancy                   = 1.0;
double                              G::_comphet                     = numeric_limits<double>::infinity();
double                              G::_asrv_shape                  = numeric_limits<double>::infinity();

double                              G::_small_enough                = 0.00001;

unsigned                            G::_nsubpops                    = 1;

unsigned                            G::_nparticles                  = 500;
unsigned                            G::_nparticles2                 = 0;
unsigned                            G::_nkept                       = 0;
unsigned                            G::_nkept2                      = 0;

#if defined(FOSSILS)
vector<Fossil>                      G::_fossils;
vector<TaxSet>                      G::_taxsets;
#endif

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required in order to use infinity()");
double                              G::_infinity = numeric_limits<double>::infinity();
double                              G::_negative_infinity = -numeric_limits<double>::infinity();

PartialStore::leaf_partials_t       GeneForest::_leaf_partials;

const double                        Node::_smallest_edge_length = 1.0e-12;

string                              Proj::_program_name         = "smc";
unsigned                            Proj::_major_version        = 3;
unsigned                            Proj::_minor_version        = 0;

GeneticCode::genetic_code_definitions_t GeneticCode::_definitions = { 
                             // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};

int main(int argc, const char * argv[]) {

#if defined(SYSTEMATIC_FILTERING)
    output("SYSTEMATIC_FILTERING", G::LogCateg::CONDITIONALS);
#endif

#if defined(LAZY_COPYING)
    output("LAZY_COPYING", G::LogCateg::CONDITIONALS);
#endif

#if defined(REUSE_PARTIALS)
    output("REUSE_PARTIALS", G::LogCateg::CONDITIONALS);
#endif

#if defined(USE_HEATING)
    output("USE_HEATING", G::LogCateg::CONDITIONALS);
#endif

#if defined(RANDOM_LOCUS_ORDERING)
    output("RANDOM_LOCUS_ORDERING", G::LogCateg::CONDITIONALS);
#endif

    Proj proj;
    bool normal_termination = true;
    try {
        proj.processCommandLineOptions(argc, argv);
        
#if defined(FOSSILS)
        //temporary!
        output("\n", 1);
        if (G::_fossils.size() > 0) {
            output("\nFossils defined:\n");
            for (auto f : G::_fossils) {
                output(format("  \"%s\" | %.5f <-- %.5f --> %.5f\n") % f._name % f._lower % f._age % f._upper);
            }
        } else {
            output("\nNo fossils were defined\n");
        }
        
        if (G::_taxsets.size() > 0) {
            output("\nTaxon sets defined:\n");
            for (auto t : G::_taxsets) {
                output(format("\n  \"%s\":\n") % t._name);
                for (auto s : t._species_included) {
                    output(format("    \"%s\"\n") % s);
                }
            }
        } else {
            output("\nNo taxon sets were defined\n");
        }
        output("\n");
#endif
        
        output(format("sizeof(unsigned) = %d\nsizeof(double) = %d\n") % sizeof(unsigned) % sizeof(double), G::LogCateg::INFO);

        StopWatch sw;
        sw.start();
        proj.run();
        double total_seconds = sw.stop();
        output(format("\nTotal time: %.5f seconds\n") % total_seconds, G::LogCateg::ALWAYS);
    }
    catch(std::exception & x) {
        cerr << str(format("Exception: %s\n") % x.what());
        cerr << "Aborted.\n";
        normal_termination = false;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        normal_termination = false;
    }
    
    return 0;
}
