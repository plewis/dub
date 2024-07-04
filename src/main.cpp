//#define NDEBUG
#include <cassert>
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
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

// See https://www.jviotti.com/2022/02/21/emitting-signposts-to-instruments-on-macos-using-cpp.html
#if defined(USING_SIGNPOSTS)
#   include <os/log.h>
#   include <os/signpost.h>
    os_log_t log_handle;
    os_signpost_id_t signpost_id;
#endif

#include "conditionals.hpp"
#include "xproj.hpp"
#include "lot.hpp"
#include "split.hpp"
#include "g.hpp"
#include "stopwatch.hpp"
#include "genetic-code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "data.hpp"
#include "partial_store.hpp"
#include "node.hpp"
#include "particle.hpp"
#include "sparticle.hpp"
#include "gparticle.hpp"
#include "bundle.hpp"
#include "smc.hpp"
#include "proj.hpp"

using namespace proj;

void output(format & fmt, G::verbosity_t verb) {
    if (verb & G::_verbosity)
        cout << str(fmt);
}

void output(string msg, G::verbosity_t verb) {
    if (verb & G::_verbosity)
        cout << msg;
}

Lot::SharedPtr                      rng(new Lot());
StopWatch                           stopwatch;
PartialStore                        ps;
Partition::SharedPtr                partition;
Data::SharedPtr                     data;

G::verbosity_t                      G::_verbosity          = G::VSTANDARD;

string                              G::_species_tree_ref_file_name = "";
string                              G::_gene_trees_ref_file_name = "";
map<unsigned,unsigned>              G::_nexus_taxon_map;

double                              G::_log_marg_like      = 0.0;

unsigned                            G::_step               = 0;
unsigned                            G::_bundle             = 0;
unsigned                            G::_locus              = 0;
unsigned                            G::_particle           = 0;

G::species_t                        G::_species_zero       = (G::species_t)0;
unsigned                            G::_nstates            = 4;

unsigned                            G::_ntaxa              = 0;
vector<string>                      G::_taxon_names;

unsigned                            G::_nspecies           = 0;
vector<string>                      G::_species_names;

unsigned                            G::_nloci              = 0;
vector<string>                      G::_locus_names;

map<string, unsigned>               G::_taxon_to_species;

unsigned                            G::_nsparticles        = 0;
unsigned                            G::_ngparticles        = 0;

double                              G::_theta              = 0.01;
double                              G::_lambda             = 1.0;

double                              G::_epsilon            = 1e-8;
double                              G::_small_enough       = 1e-12;

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required in order to use infinity()");
double                              G::_infinity = numeric_limits<double>::infinity();
double                              G::_negative_infinity = -numeric_limits<double>::infinity();

PartialStore::leaf_partials_t       GParticle::_leaf_partials;

const double                        Node::_smallest_edge_length = 1.0e-12;

#if defined(LOG_MEMORY)
vector<unsigned>                    Partial::_nconstructed;
vector<unsigned>                    Partial::_ndestroyed;
vector<unsigned>                    Partial::_max_in_use;
vector<unsigned long>               Partial::_bytes_per_partial;
unsigned long                       Partial::_total_max_in_use  = 0;
unsigned long                       Partial::_total_max_bytes   = 0;
unsigned                            Partial::_nstates           = 4;
#endif

string                              Proj::_program_name         = "smc5";
unsigned                            Proj::_major_version        = 1;
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

vector<string> rcolors = {"aliceblue" ,"aquamarine" ,"aquamarine1" ,"aquamarine2" ,"aquamarine3" ,"aquamarine4" ,"azure" ,"azure1" ,"azure2" ,"azure3" ,"azure4" ,"beige" ,"bisque" ,"bisque1" ,"bisque2" ,"bisque3" ,"bisque4" ,"black" ,"blanchedalmond" ,"blue" ,"blue1" ,"blue2" ,"blue3" ,"blue4" ,"blueviolet" ,"brown" ,"brown1" ,"brown2" ,"brown3" ,"brown4" ,"burlywood" ,"burlywood1" ,"burlywood2" ,"burlywood3" ,"burlywood4" ,"cadetblue" ,"cadetblue1" ,"cadetblue2" ,"cadetblue3" ,"cadetblue4" ,"chartreuse" ,"chartreuse1" ,"chartreuse2" ,"chartreuse3" ,"chartreuse4" ,"chocolate" ,"chocolate1" ,"chocolate2" ,"chocolate3" ,"chocolate4" ,"coral" ,"coral1" ,"coral2" ,"coral3" ,"coral4" ,"cornflowerblue" ,"cornsilk" ,"cornsilk1" ,"cornsilk2" ,"cornsilk3" ,"cornsilk4" ,"cyan" ,"cyan1" ,"cyan2" ,"cyan3" ,"cyan4" ,"darkblue" ,"darkcyan" ,"darkgoldenrod" ,"darkgoldenrod1" ,"darkgoldenrod2" ,"darkgoldenrod3" ,"darkgoldenrod4" ,"darkgray" ,"darkgreen" ,"darkgrey" ,"darkkhaki" ,"darkmagenta" ,"darkolivegreen" ,"darkolivegreen1" ,"darkolivegreen2" ,"darkolivegreen3" ,"darkolivegreen4" ,"darkorange" ,"darkorange1" ,"darkorange2" ,"darkorange3" ,"darkorange4" ,"darkorchid" ,"darkorchid1" ,"darkorchid2" ,"darkorchid3" ,"darkorchid4" ,"darkred" ,"darksalmon" ,"darkseagreen" ,"darkseagreen1" ,"darkseagreen2" ,"darkseagreen3" ,"darkseagreen4" ,"darkslateblue" ,"darkslategray" ,"darkslategray1" ,"darkslategray2" ,"darkslategray3" ,"darkslategray4" ,"darkslategrey" ,"darkturquoise" ,"darkviolet" ,"deeppink" ,"deeppink1" ,"deeppink2" ,"deeppink3" ,"deeppink4" ,"deepskyblue" ,"deepskyblue1" ,"deepskyblue2" ,"deepskyblue3" ,"deepskyblue4" ,"dimgray" ,"dimgrey" ,"dodgerblue" ,"dodgerblue1" ,"dodgerblue2" ,"dodgerblue3" ,"dodgerblue4" ,"firebrick" ,"firebrick1" ,"firebrick2" ,"firebrick3" ,"firebrick4" ,"floralwhite" ,"forestgreen" ,"gainsboro" ,"ghostwhite" ,"gold" ,"gold1" ,"gold2" ,"gold3" ,"gold4" ,"goldenrod" ,"goldenrod1" ,"goldenrod2" ,"goldenrod3" ,"goldenrod4" ,"gray" ,"green" ,"green1" ,"green2" ,"green3" ,"green4" ,"greenyellow" ,"honeydew" ,"honeydew1" ,"honeydew2" ,"honeydew3" ,"honeydew4" ,"hotpink" ,"hotpink1" ,"hotpink2" ,"hotpink3" ,"hotpink4" ,"indianred" ,"indianred1" ,"indianred2" ,"indianred3" ,"indianred4" ,"ivory" ,"ivory1" ,"ivory2" ,"ivory3" ,"ivory4" ,"khaki" ,"khaki1" ,"khaki2" ,"khaki3" ,"khaki4" ,"lavender" ,"lavenderblush" ,"lavenderblush1" ,"lavenderblush2" ,"lavenderblush3" ,"lavenderblush4" ,"lawngreen" ,"lemonchiffon" ,"lemonchiffon1" ,"lemonchiffon2" ,"lemonchiffon3" ,"lemonchiffon4" ,"lightblue" ,"lightblue1" ,"lightblue2" ,"lightblue3" ,"lightblue4" ,"lightcoral" ,"lightcyan" ,"lightcyan1" ,"lightcyan2" ,"lightcyan3" ,"lightcyan4" ,"lightgoldenrod" ,"lightgoldenrod1" ,"lightgoldenrod2" ,"lightgoldenrod3" ,"lightgoldenrod4" ,"lightgoldenrodyellow" ,"lightgray" ,"lightgreen" ,"lightgrey" ,"lightpink" ,"lightpink1" ,"lightpink2" ,"lightpink3" ,"lightpink4" ,"lightsalmon" ,"lightsalmon1" ,"lightsalmon2" ,"lightsalmon3" ,"lightsalmon4" ,"lightseagreen" ,"lightskyblue" ,"lightskyblue1" ,"lightskyblue2" ,"lightskyblue3" ,"lightskyblue4" ,"lightslateblue" ,"lightslategray" ,"lightslategrey" ,"lightsteelblue" ,"lightsteelblue1" ,"lightsteelblue2" ,"lightsteelblue3" ,"lightsteelblue4" ,"lightyellow" ,"lightyellow1" ,"lightyellow2" ,"lightyellow3" ,"lightyellow4" ,"limegreen" ,"linen" ,"magenta" ,"magenta1" ,"magenta2" ,"magenta3" ,"magenta4" ,"maroon" ,"maroon1" ,"maroon2" ,"maroon3" ,"maroon4" ,"mediumaquamarine" ,"mediumblue" ,"mediumorchid" ,"mediumorchid1" ,"mediumorchid2" ,"mediumorchid3" ,"mediumorchid4" ,"mediumpurple" ,"mediumpurple1" ,"mediumpurple2" ,"mediumpurple3" ,"mediumpurple4" ,"mediumseagreen" ,"mediumslateblue" ,"mediumspringgreen" ,"mediumturquoise" ,"mediumvioletred" ,"midnightblue" ,"mintcream" ,"mistyrose" ,"mistyrose1" ,"mistyrose2" ,"mistyrose3" ,"mistyrose4" ,"moccasin" ,"navajowhite" ,"navajowhite1" ,"navajowhite2" ,"navajowhite3" ,"navajowhite4" ,"navy" ,"navyblue" ,"oldlace" ,"olivedrab" ,"olivedrab1" ,"olivedrab2" ,"olivedrab3" ,"olivedrab4" ,"orange" ,"orange1" ,"orange2" ,"orange3" ,"orange4" ,"orangered" ,"orangered1" ,"orangered2" ,"orangered3" ,"orangered4" ,"orchid" ,"orchid1" ,"orchid2" ,"orchid3" ,"orchid4" ,"palegoldenrod" ,"palegreen" ,"palegreen1" ,"palegreen2" ,"palegreen3" ,"palegreen4" ,"paleturquoise" ,"paleturquoise1" ,"paleturquoise2" ,"paleturquoise3" ,"paleturquoise4" ,"palevioletred" ,"palevioletred1" ,"palevioletred2" ,"palevioletred3" ,"palevioletred4" ,"papayawhip" ,"peachpuff" ,"peachpuff1" ,"peachpuff2" ,"peachpuff3" ,"peachpuff4" ,"peru" ,"pink" ,"pink1" ,"pink2" ,"pink3" ,"pink4" ,"plum" ,"plum1" ,"plum2" ,"plum3" ,"plum4" ,"powderblue" ,"purple" ,"purple1" ,"purple2" ,"purple3" ,"purple4" ,"red" ,"red1" ,"red2" ,"red3" ,"red4" ,"rosybrown" ,"rosybrown1" ,"rosybrown2" ,"rosybrown3" ,"rosybrown4" ,"royalblue" ,"royalblue1" ,"royalblue2" ,"royalblue3" ,"royalblue4" ,"saddlebrown" ,"salmon" ,"salmon1" ,"salmon2" ,"salmon3" ,"salmon4" ,"sandybrown" ,"seagreen" ,"seagreen1" ,"seagreen2" ,"seagreen3" ,"seagreen4" ,"seashell" ,"seashell1" ,"seashell2" ,"seashell3" ,"seashell4" ,"sienna" ,"sienna1" ,"sienna2" ,"sienna3" ,"sienna4" ,"skyblue" ,"skyblue1" ,"skyblue2" ,"skyblue3" ,"skyblue4" ,"slateblue" ,"slateblue1" ,"slateblue2" ,"slateblue3" ,"slateblue4" ,"slategray" ,"slategray1" ,"slategray2" ,"slategray3" ,"slategray4" ,"slategrey" ,"snow" ,"snow1" ,"snow2" ,"snow3" ,"snow4" ,"springgreen" ,"springgreen1" ,"springgreen2" ,"springgreen3" ,"springgreen4" ,"steelblue" ,"steelblue1" ,"steelblue2" ,"steelblue3" ,"steelblue4" ,"tan" ,"tan1" ,"tan2" ,"tan3" ,"tan4" ,"thistle" ,"thistle1" ,"thistle2" ,"thistle3" ,"thistle4" ,"tomato" ,"tomato1" ,"tomato2" ,"tomato3" ,"tomato4" ,"turquoise" ,"turquoise1" ,"turquoise2" ,"turquoise3" ,"turquoise4" ,"violet" ,"violetred" ,"violetred1" ,"violetred2" ,"violetred3" ,"violetred4" ,"wheat" ,"wheat1" ,"wheat2" ,"wheat3" ,"wheat4" ,"whitesmoke" ,"yellow" ,"yellow1" ,"yellow2" ,"yellow3" ,"yellow4" ,"yellowgreen"};

int main(int argc, const char * argv[]) {

    Proj proj;
    bool normal_termination = true;
    try {
        proj.processCommandLineOptions(argc, argv);
        
        StopWatch sw;
        sw.start();
        proj.run();
        double total_seconds = sw.stop();
        output(format("\nTotal time: %.5f seconds\n") % total_seconds, G::VSTANDARD);
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
    
#if defined(LOG_MEMORY)
    if (normal_termination) {
        ofstream memfile("allocs.txt");
        proj.memoryReport(memfile);
        ps.memoryReport(memfile);
        memfile.close();
    }
#endif

    return 0;
}
