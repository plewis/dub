#pragma once

//extern void output(string msg, proj::G::verbosity_t verb);
//extern void output(format & fmt, proj::G::verbosity_t level);
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
            void setBundleIndex(unsigned i) {_bundle_index = i; _species_tree.setIndex(i);}
            unsigned getBundleIndex() const {return _bundle_index;};
            
            unsigned getNumSpeciesTreeLineages() const {return (unsigned)_species_tree._lineages.size();};

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
    };
    
    inline Bundle::Bundle() {
        auto spp_incr_rate = _species_tree.drawIncrement();
        _species_tree.extendAllLineagesBy(spp_incr_rate.first);
        
        // Build vector
        _locus_vect.resize(G::_nloci);
        for (unsigned g = 0; g < G::_nloci; g++) {
            _locus_vect[g].resize(G::_ngparticles);
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                // Each GParticle must have a pointer to the species tree
                _locus_vect[g][i].setSpeciesTree(&_species_tree);
                
                // Each GParticle needs to know which locus it belongs to
                _locus_vect[g][i].setLocusIndex(g);
                
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
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                double h = _locus_vect[g][i]._height;
                if (h > deepest_gene_tree)
                    deepest_gene_tree = h;
            }
        }
        
        // Search for and remove any nodes in the species
        // tree older than deepest_gene_tree
        _species_tree.trimToHeight(deepest_gene_tree);
    }
    
    inline void Bundle::debugSanityCheck() const {
        // Throws exception if any gene tree violates gene flow barriers determined by species tree
        vector<G::join_info_t> spec_info;
        vector<G::join_info_t> coal_info;
        _species_tree.recordJoinInfo(spec_info);
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                output(format("\nlocus %d, particle %d:\n") % g % i, G::VTEMP);
                coal_info.clear();
                _locus_vect[g][i].recordJoinInfo(coal_info);
                coal_info.insert(coal_info.end(), spec_info.begin(), spec_info.end());
                sort(coal_info.begin(), coal_info.end());
                
                for (auto & t : coal_info) {
                    double h = get<0>(t);
                    bool is_spp_tree = get<1>(t);
                    G::species_t sppL = get<2>(t);
                    G::species_t sppR = get<3>(t);
                    output(format("%12.7f %6s %12d %12d %12d\n") % h % (is_spp_tree ? "spp" : "coal") % (sppL|sppR) % sppL % sppR, G::VTEMP);
                    if (sppL != sppR && !is_spp_tree) {
                        throw XProj(format("Coalescent event joined different species (%d, %d) at height %g in locus %d, particle %d") % sppL % sppR % h % g % i);
                    }
                }
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
            log_likes[g].resize(G::_ngparticles, 0.0);
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                log_likes[g][i] = _locus_vect[g][i].calcLogLikelihood();
            }
        }
    }
    
    void Bundle::getLogLikesForLocus(unsigned g, vector<double> & log_likes) const {
        log_likes.resize(G::_ngparticles, 0.0);
        for (unsigned i = 0; i < G::_ngparticles; i++) {
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
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                string delim = (i == G::_ngparticles - 1 ? "\n" : "");
                output(format("%12s%s") % _locus_vect[g][i].info() % delim, G::VDEBUG);
            }
        }
        return total_log_marg_like;
    }
    
    inline void Bundle::randomizeEvaluationOrder() {
        unsigned total_eval = G::_nloci*G::_ngparticles;
        
        // Build a list of indices that we will choose from
        // Each time an index is chosen, that element will be removed from the list
        list<unsigned> indices(total_eval);
        unsigned i = 0;
        for (auto iter = indices.begin(); iter != indices.end(); ++iter) {
            *iter = i++;
        }
        
        _eval_order.resize(total_eval);
        for (G::_locus = 0; G::_locus < G::_nloci; G::_locus++) {
            for (G::_particle = 0; G::_particle < G::_ngparticles; G::_particle++) {
                // Choose index
                double u = rng->uniform();
                unsigned n = (unsigned)indices.size();
                unsigned i = (unsigned)floor(u*n);
                assert(i < indices.size());
                auto iter = indices.begin();
                advance(iter, i);
                unsigned index = *iter;
                indices.erase(iter);
                _eval_order[index] = make_pair(G::_locus, G::_particle);
            }
        }
        
    }
    
    inline void Bundle::advanceAllGeneTrees() {
        randomizeEvaluationOrder();

        output(format("Advancing all gene trees for all loci in bundle %d\n") % _bundle_index, G::VDEBUG);
        for (auto p : _eval_order) {
            G::_locus = p.first;
            G::_particle = p.second;
            output(format("\nWorking on particle %d from locus %d...") % G::_particle % G::_locus, G::VTEMP);
            _locus_vect[G::_locus][G::_particle].coalesce();
        }
    }
    
    inline void Bundle::filterAllGeneTrees(unsigned step) {
        output(format("Filtering gene trees for all loci in bundle %d\n") % _bundle_index, G::VDEBUG);
        for (G::_locus = 0; G::_locus < G::_nloci; G::_locus++) {
            output(format("    Gene %d\n") % G::_locus, G::VDEBUG);
            
            // Copy log weights for all particles in this locus to prob vector
            vector<double> probs(G::_ngparticles, 0.0);
            for (G::_particle = 0; G::_particle < G::_ngparticles; G::_particle++) {
                probs[G::_particle] = _locus_vect[G::_locus][G::_particle].getLogWeight();
                //output(format(" %6d %12.5f\n") % G::_particle % probs[G::_particle], G::VDEBUG);
            }

            // Normalize log_weights to create discrete probability distribution
            double log_sum_weights = G::calcLogSum(probs);
            transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

            // Compute component of the log marginal likelihood due to this step
            _prev_log_marg_like[G::_locus] = _log_marg_like[G::_locus];
            _log_marg_like[G::_locus] += log_sum_weights - log(G::_ngparticles);

            // Compute cumulative probabilities
            partial_sum(probs.begin(), probs.end(), probs.begin());
            
            // Initialize vector of counts storing number of darts hitting each particle
            vector<unsigned> counts(G::_ngparticles, 0);
                    
            // Throw _nparticles darts
            for (unsigned i = 0; i < G::_ngparticles; ++i) {
                double u = rng->uniform();
                auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
                assert(it != probs.end());
                unsigned which = (unsigned)distance(probs.begin(), it);
                counts[which]++;
            }
                
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
            
            // Count number of cells that need copying to
            unsigned nzeros = 0;
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }
            
            while (nzeros > 0) {
                assert(donor < G::_ngparticles);
                assert(recipient < G::_ngparticles);
                _locus_vect[G::_locus][recipient] = _locus_vect[G::_locus][donor];
                counts[donor]--;
                counts[recipient]++;
                nzeros--;
                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < G::_ngparticles && counts[donor] < 2) {
                        donor++;
                    }
                }
                
                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < G::_ngparticles && counts[recipient] > 0) {
                    recipient++;
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
            for (unsigned p = 0; p < G::_ngparticles; p++) {
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
        
        // Copy species tree
        _species_tree = other._species_tree;
        
        // Copy log marginal likelihood
        _log_marg_like = other._log_marg_like;
        
        // Copy gene trees
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                _locus_vect[g][i] = other._locus_vect[g][i];
            }
        }
    }
}

