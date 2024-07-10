#pragma once

//extern void output(string msg, proj::G::verbosity_t verb);
//extern void output(format & fmt, proj::G::verbosity_t level);
extern void doof(string msg, proj::G::verbosity_t verb);
extern void doof(format & fmt, proj::G::verbosity_t level);
//extern proj::PartialStore         ps;
//extern proj::StopWatch            stopwatch;
//extern proj::Lot::SharedPtr       rng;
//extern proj::Partition::SharedPtr partition;
//extern proj::Data::SharedPtr      data;

namespace proj {

    // Bundles are compound particles comprising one SParticle and
    // a vector of GParticles for every locus
    class Bundle {
        public:
            Bundle();
            ~Bundle() {}
            
            void operator=(const Bundle & other);
            
            // Getters and setters
            void setBundleIndex(unsigned i);
            unsigned getBundleIndex() const {return _bundle_index;};
            
            Lot::SharedPtr getLot() const {return _lot;}
            void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            
            unsigned getNumSpeciesTreeLineages() const {return (unsigned)_species_tree._lineages.size();};
            
            GParticle & getGParticle(unsigned locus_index, unsigned particle_index) {
                return _locus_vect[locus_index][particle_index];
            }
            const GParticle & getGParticleConst(unsigned locus_index, unsigned particle_index) {
                return _locus_vect[locus_index][particle_index];
            }

            double getLogMargLike() const;
            double getPrevLogMargLike() const;
            double getLogWeight() const;
            void getLogLikes(vector< vector<double> > & log_likes) const;
            void getLogLikesForLocus(unsigned g, vector<double> & log_likes) const;
            
            SParticle & getSpeciesTree() {return _species_tree;}
            const SParticle & getSpeciesTreeConst() const {return _species_tree;}
            
            GParticle & getGeneTree(unsigned g, unsigned i) {
                return _locus_vect[g][i];
            }
            const GParticle & getGeneTreeConst(unsigned g, unsigned i) const {
                return _locus_vect[g][i];
            }

            void initSpeciesTree();
            void saveJavascript(string fnprefix) const;
            void advanceAllGeneTrees();
            void filterAllGeneTrees(unsigned step);
            void randomizeEvaluationOrder();
            void shrinkWrapSpeciesTree();
            void debugSanityCheck() const;
            
            double report() const;
                        
        protected:
            vector< vector<GParticle> >     _locus_vect;
            SParticle                       _species_tree;
            unsigned                        _bundle_index;
            vector<double>                  _log_marg_like;
            vector<double>                  _prev_log_marg_like;
            vector< pair<unsigned, unsigned> > _eval_order;

            // Even though it is a shared pointer, _lot is a private random number
            // generator not shared with any other bundle and has nothing to
            // to do with the global Lot shared_ptr rng defined in main.cpp.
            // Note that _lot is excluded from copying.
            mutable Lot::SharedPtr  _lot;
    };
    
    inline Bundle::Bundle() {
        _lot.reset(new Lot());

        //auto spp_incr_rate = _species_tree.drawIncrement(_lot);
        //_species_tree.extendAllLineagesBy(spp_incr_rate.first);
        
        // Build vector
        _locus_vect.resize(G::_nloci);
        for (unsigned g = 0; g < G::_nloci; g++) {
            _locus_vect[g].resize(G::_nparticles);
            for (unsigned i = 0; i < G::_nparticles; i++) {
                // Each GParticle must have a pointer to the species tree
                _locus_vect[g][i].setSpeciesTree(&_species_tree);
                
                // Each GParticle needs to know to which locus it belongs
                _locus_vect[g][i].setLocusIndex(g);
                
                // Each GParticle needs to know to which bundle it belongs
                // but currently the bundle itself doean't know its index
                // so this needs to be set in setBundleIndex
                _locus_vect[g][i].setBundleIndex(-1);
                
                // Each GParticle needs to know its index
                _locus_vect[g][i].setIndex(i);
                
                // Build a trivial forest comprising only leaf nodes
                _locus_vect[g][i].createTrivialForest();
            }
        }
        
        // Compute starting log likelihood for all loci
        _log_marg_like.resize(G::_nloci, 0.0);
        _prev_log_marg_like.resize(G::_nloci, 0.0);
        double total_log_like = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            _log_marg_like[g] = _locus_vect[g][0].calcLogLikelihood();
            total_log_like = _log_marg_like[g];
            output(format("%12.5f = starting log-likelihood for locus %d\n") %  _log_marg_like[g] % g, G::VDEBUG);
        }
        output(format("%12.5f = starting total log-likelihood\n") %  total_log_like, G::VDEBUG);
    }
    
    inline void Bundle::shrinkWrapSpeciesTree() {
        // Remove deepest joins in species tree if they have height
        // greater than the deepest coalescent event in any gene
        // tree at any locus
        
        // Find deepest gene tree
        double deepest_gene_tree = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned i = 0; i < G::_nparticles; i++) {
                double h = _locus_vect[g][i]._height;
                if (h > deepest_gene_tree)
                    deepest_gene_tree = h;
            }
        }
        
        // Search for and remove any nodes in the species
        // tree older than deepest_gene_tree
        _species_tree.trimToHeight(deepest_gene_tree, _lot);
    }
    
    inline void Bundle::initSpeciesTree() {
        assert(_species_tree.getHeight() == 0.0);
        auto spp_incr_rate = _species_tree.drawIncrement(_lot);
        _species_tree.extendAllLineagesBy(spp_incr_rate.first);
    }
    
    inline void Bundle::debugSanityCheck() const {
        // Throws exception if any gene tree violates gene flow barriers determined by species tree
        vector<G::join_info_t> spec_info;
        vector<G::join_info_t> coal_info;
        _species_tree.recordJoinInfo(spec_info);
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned i = 0; i < G::_nparticles; i++) {
                output(format("\nlocus %d, particle %d:\n") % g % i, G::VDEBUG);
                coal_info.clear();
                _locus_vect[g][i].recordJoinInfo(coal_info);
                coal_info.insert(coal_info.end(), spec_info.begin(), spec_info.end());
                sort(coal_info.begin(), coal_info.end());
                
                for (auto & t : coal_info) {
                    double h = get<0>(t);
                    bool is_spp_tree = get<1>(t);
                    G::species_t sppL = get<2>(t);
                    G::species_t sppR = get<3>(t);
                    output(format("%12.7f %6s %12d %12d %12d\n") % h % (is_spp_tree ? "spp" : "coal") % (sppL|sppR) % sppL % sppR, G::VDEBUG);
                    if (sppL != sppR && !is_spp_tree) {
                        throw XProj(format("Coalescent event joined different species (%d, %d) at height %g in locus %d, particle %d") % sppL % sppR % h % g % i);
                    }
                }
            }
        }
    }

    inline void Bundle::setBundleIndex(unsigned i) {
        _bundle_index = i;
        _species_tree.setIndex(i);
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned k = 0; k < G::_nparticles; k++) {
                // Let each GParticle know which bundle it belongs to
                _locus_vect[g][k].setBundleIndex(i);
            }
        }
    }

    inline double Bundle::getLogWeight() const {
        double curr_log_marg_like = getLogMargLike();
        double prev_log_marg_like = getPrevLogMargLike();
        return curr_log_marg_like - prev_log_marg_like;
    }
    
    inline double Bundle::getLogMargLike() const {
        double total_log_marg_like = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            total_log_marg_like += _log_marg_like[g];
        }
        return total_log_marg_like;
    }
    
    inline double Bundle::getPrevLogMargLike() const {
        double total_log_marg_like = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            total_log_marg_like += _prev_log_marg_like[g];
        }
        return total_log_marg_like;
    }
    
    inline void Bundle::getLogLikes(vector< vector<double> > & log_likes) const {
        log_likes.resize(G::_nloci);
        for (unsigned g = 0; g < G::_nloci; g++) {
            log_likes[g].resize(G::_nparticles, 0.0);
            for (unsigned i = 0; i < G::_nparticles; i++) {
                log_likes[g][i] = _locus_vect[g][i].calcLogLikelihood();
            }
        }
    }
    
    void Bundle::getLogLikesForLocus(unsigned g, vector<double> & log_likes) const {
        log_likes.resize(G::_nparticles, 0.0);
        for (unsigned i = 0; i < G::_nparticles; i++) {
            log_likes[i] = _locus_vect[g][i].calcLogLikelihood();
        }
    }
    
    inline double Bundle::report() const {
        output(format("\nReport for bundle %d:\n") % _bundle_index, G::VSTANDARD);
        double total_log_marg_like = getLogMargLike();
        output(format("  log marginal likelihood: %.5f\n") % total_log_marg_like, G::VSTANDARD);
        output(format("%s\n") % _species_tree.info(), G::VSTANDARD);
        output(format("  %s\n") % _species_tree.makeNewick(5, true), G::VSTANDARD);
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned i = 0; i < G::_nparticles; i++) {
                string delim = (i == G::_nparticles - 1 ? "\n" : "");
                output(format("%12s%s") % _locus_vect[g][i].info() % delim, G::VDEBUG);
            }
        }
        return total_log_marg_like;
    }
    
    inline void Bundle::randomizeEvaluationOrder() {
        unsigned total_eval = G::_nloci*G::_nparticles;
        
        // Build a list of indices that we will choose from
        // Each time an index is chosen, that element will be removed from the list
        list<unsigned> indices(total_eval);
        unsigned i = 0;
        for (auto iter = indices.begin(); iter != indices.end(); ++iter) {
            *iter = i++;
        }
        
        _eval_order.resize(total_eval);
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned p = 0; p < G::_nparticles; p++) {
                // Choose index
                double u = _lot->uniform();
                unsigned n = (unsigned)indices.size();
                unsigned i = (unsigned)floor(u*n);
                assert(i < indices.size());
                auto iter = indices.begin();
                advance(iter, i);
                unsigned index = *iter;
                indices.erase(iter);
                _eval_order[index] = make_pair(g, p);
            }
        }
        
    }
    
    inline void Bundle::advanceAllGeneTrees() {
        randomizeEvaluationOrder();

        output(format("Advancing all gene trees for all loci in bundle %d\n") % _bundle_index, G::VDEBUG);
        for (auto evalpair : _eval_order) {
            unsigned locus = evalpair.first;
            unsigned particle = evalpair.second;
            _locus_vect[locus][particle].coalesce(_lot);
        }
    }
    
    inline void Bundle::filterAllGeneTrees(unsigned step) {
        output(format("Filtering gene trees for all loci in bundle %d\n") % _bundle_index, G::VDEBUG);
        for (unsigned g = 0; g < G::_nloci; g++) {
            output(format("    Gene %d\n") % g, G::VDEBUG);
            
            // Copy log weights for all particles in this locus to prob vector
            vector<double> probs(G::_nparticles, 0.0);
            for (unsigned p = 0; p < G::_nparticles; p++) {
                probs[p] = _locus_vect[g][p].getLogWeight();
                output(format(" %6d %12.5f\n") % p % probs[p], G::VDEBUG);
            }

            // Normalize log_weights to create discrete probability distribution
            double log_sum_weights = G::calcLogSum(probs);
            assert(!isnan(log_sum_weights));
            transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

            // Compute component of the log marginal likelihood due to this step
            _prev_log_marg_like[g] = _log_marg_like[g];
            _log_marg_like[g] += log_sum_weights - log(G::_nparticles);

            // Compute cumulative probabilities
            partial_sum(probs.begin(), probs.end(), probs.begin());
            
            // Initialize vector of counts storing number of darts hitting each particle
            vector<unsigned> counts(G::_nparticles, 0);
                    
            // Throw G::_nparticles darts
            for (unsigned i = 0; i < G::_nparticles; ++i) {
                double u = _lot->uniform();
                auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
                assert(it != probs.end());
                unsigned which = (unsigned)distance(probs.begin(), it);
                counts[which]++;
            }
                
            // Count number of cells that need copying to
            unsigned nzeros = 0;
            for (unsigned i = 0; i < G::_nparticles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }
            
            if (nzeros > 0) {
                // Copy particles
                // Locate first donor
                unsigned donor = 0;
                while (counts[donor] < 2) {
                    donor++;
                }
                
                // Locate first recipient
                unsigned recipient = 0;
                while (counts[recipient] != 0) {
                    recipient++;
                }
                
                while (nzeros > 0) {
                    assert(donor < G::_nparticles);
                    assert(recipient < G::_nparticles);
                    _locus_vect[g][recipient] = _locus_vect[g][donor];
                    counts[donor]--;
                    counts[recipient]++;
                    nzeros--;
                    if (counts[donor] == 1) {
                        // Move donor to next slot with count > 1
                        donor++;
                        while (donor < G::_nparticles && counts[donor] < 2) {
                            donor++;
                        }
                    }
                    
                    // Move recipient to next slot with count equal to 0
                    recipient++;
                    while (recipient < G::_nparticles && counts[recipient] > 0) {
                        recipient++;
                    }
                }
            }
        }
    }
    
    inline void Bundle::saveJavascript(string fnprefix) const {
        string fn = str(format("%s.js") % fnprefix);
        ofstream jsf(fn);
        
        jsf << "let species_translate = {\n";
        unsigned i = 1;
        for (string nm : G::_species_names) {
            string comma = (i < G::_nspecies ? "," : "");
            jsf << str(format("  %d: \"%s\"%s\n") % i % nm % comma);
            ++i;
        }
        jsf << "};\n\n";
        
        string newick_species_tree_numeric = _species_tree.makeNewick(5, false);
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

        for (unsigned l = 0; l < G::_nloci; l++) {
            // Map counts (values) onto gene tree newicks (keys)
            map<string, unsigned> m;
            for (unsigned p = 0; p < G::_nparticles; p++) {
                string newick = _locus_vect[l][p].makeNewick(9, false);
                m[newick] += 1;
            }
            
            // Save newick with highest count
            auto it = max_element(m.begin(), m.end(),
                [](const pair<string,unsigned> & a, const pair<string,unsigned> & b) {return b.second > a.second;});
                jsf << str(format("  {name:\"gene-%d\",  relrate:1.0, color:\"black\", newick:\"%s\"},\n") % (l+1) % it->first);
        }
        
        jsf << "];\n";
        jsf.close();
    }
    
    inline void Bundle::operator=(const Bundle & other) {
        // _bundle_index should not be copied
        // _eval_order does not need to be copied
        // _lot should NOT be copied
        
        // Copy species tree
        _species_tree = other._species_tree;
        
        // Copy log marginal likelihood
        _log_marg_like = other._log_marg_like;
        
        // Copy gene trees and _prev_log_marg_like
        for (unsigned g = 0; g < G::_nloci; g++) {
            _prev_log_marg_like[g] = other._prev_log_marg_like[g];
            for (unsigned i = 0; i < G::_nparticles; i++) {
                _locus_vect[g][i] = other._locus_vect[g][i];
                assert(_locus_vect[g][i].getSpeciesTree() == &_species_tree);
            }
        }
    }
    
}

