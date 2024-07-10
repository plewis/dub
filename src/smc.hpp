#pragma once

namespace proj {

    class SMC {
        public:
            SMC() {}
            ~SMC() {}
            
            void run();
            void advanceBundleRange(unsigned step, unsigned thread, unsigned first, unsigned last);
            
            unsigned saveBestSpeciesTree() const;

            const SParticle & getSpeciesTreeConst(unsigned b) const;
            const GParticle & getGeneTreeConst(unsigned b, unsigned g, unsigned i) const;
            
            void saveSpeciesTrees(vector<string> & newicks, bool compress) const;
            void saveSpeciesTreesToFile(string fn, bool compress) const;
            
            void saveGeneTreeBundleLocus(unsigned b, unsigned g, vector<string> & newicks, bool compress) const;
            void saveGeneTreeBundleLocusToFile(string fn, unsigned b, unsigned g, bool compress) const;

            void saveGeneTreeLocus(unsigned g, vector<string> & newicks, bool compress) const;
            void saveGeneTreeLocusToFile(string fn, unsigned g, bool compress) const;
            
            void saveLogMargLike(vector<double> & log_marg_like) const;
            void saveLogLikesForLocus(unsigned b, unsigned g, vector<double> & log_likes) const;
            void saveLogLikes(unsigned b, vector< vector<double> > & log_likes_for_locus) const;
            
        private:
            void returnUnusedPartials();
            void filterBundles(unsigned step);
            void sanityCheckBundles() const;
        
            vector<unsigned>    _update_seeds;
            vector<Bundle *>    _bundle_vect;
    };
    
    inline const SParticle & SMC::getSpeciesTreeConst(unsigned b) const {
        assert(b < _bundle_vect.size());
        assert(_bundle_vect[b]);
        return _bundle_vect[b]->getSpeciesTreeConst();
    }
    
    inline const GParticle & SMC::getGeneTreeConst(unsigned b, unsigned g, unsigned i) const {
        assert(b < _bundle_vect.size());
        assert(_bundle_vect[b]);
        assert(g < G::_nloci);
        assert(i < G::_nparticles);
        return _bundle_vect[b]->getGeneTreeConst(g, i);
    }
    
    inline void SMC::saveSpeciesTrees(vector<string> & newicks, bool compress) const {
        if (compress) {
            // Save only unique newick species tree descriptions
            // May save same topology many times if edge lengths differ

            // Create map with newicks as keys and counts as values
            newicks.clear();
            map<string, unsigned> freq;
            for (unsigned i = 0; i < G::_nbundles; i++) {
                assert(_bundle_vect[i]);
                string newick = _bundle_vect[i]->getSpeciesTreeConst().makeNewick(9, true);
                freq[newick] += 1;
            }

            // Sort newicks by decreasing count
            vector< pair<unsigned, string> > count_newick(freq.size());
            transform(freq.begin(), freq.end(), count_newick.begin(), [](pair<string, unsigned> p){return make_pair(p.second, p.first);});
            sort(count_newick.begin(), count_newick.end(), greater< pair<unsigned, string> >());
            
            // Save the newicks in order of decreasing count
            for (auto & uniq : count_newick) {
                newicks.push_back(str(format("[freq=%d] %s") % uniq.first % uniq.second));
            }
        }
        else {
            // Save all species trees, even if they are all identical newick strings
            newicks.resize(G::_nbundles);
            for (unsigned i = 0; i < G::_nbundles; i++) {
                assert(_bundle_vect[i]);
                newicks[i] = _bundle_vect[i]->getSpeciesTreeConst().makeNewick(9, true);
            }
        }
    }
    
    inline void SMC::saveSpeciesTreesToFile(string fn, bool compress) const {
        vector<string> newicks;
        saveSpeciesTrees(newicks, compress);
        
        ofstream treef(fn);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        unsigned i = 1;
        for (auto newick : newicks) {
            treef << str(format("  tree t%d = %s;\n") % (i++) % newick);
        }
        treef << "end;\n";
        treef.close();
    }
    
    inline void SMC::saveGeneTreeBundleLocus(unsigned b, unsigned g, vector<string> & newicks, bool compress) const {
        assert(_bundle_vect[b]);
        if (compress) {
            // Save only unique newick species tree descriptions
            // May save same topology many times if edge lengths differ
            
            // Create map with newicks as keys and counts as values
            newicks.clear();
            map<string, unsigned> freq;
            for (unsigned i = 0; i < G::_nparticles; i++) {
                string newick = _bundle_vect[b]->getGeneTreeConst(g, i).makeNewick(9, true);
                freq[newick] += 1;
            }
            
            // Sort newicks by decreasing count
            vector< pair<unsigned, string> > count_newick(freq.size());
            transform(freq.begin(), freq.end(), count_newick.begin(), [](pair<string, unsigned> p){return make_pair(p.second, p.first);});
            sort(count_newick.begin(), count_newick.end(), greater< pair<unsigned, string> >());
            
            // Save the newicks in order of decreasing count
            for (auto & uniq : count_newick) {
                newicks.push_back(str(format("[freq=%d] %s") % uniq.first % uniq.second));
            }
        }
        else {
            // Save all species trees, even if they are all identical newick strings
            newicks.resize(G::_nparticles);
            for (unsigned i = 0; i < G::_nbundles; i++) {
                newicks[i] = _bundle_vect[b]->getGeneTreeConst(g, i).makeNewick(9, true);
            }
        }
    }
    
    inline void SMC::saveGeneTreeBundleLocusToFile(string fn, unsigned b, unsigned g, bool compress) const {
        vector<string> newicks;
        saveGeneTreeBundleLocus(b, g, newicks, compress);
        
        ofstream treef(fn);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        unsigned i = 1;
        for (auto newick : newicks) {
            treef << str(format("  tree t%d = %s;\n") % (i++) % newick);
        }
        treef << "end;\n";
        treef.close();
    }
    
    inline void SMC::saveGeneTreeLocus(unsigned g, vector<string> & newicks, bool compress) const {
        if (compress) {
            // Save only unique newick species tree descriptions
            // May save same topology many times if edge lengths differ

            // Create map with newicks as keys and counts as values
            newicks.clear();
            map<string, unsigned> freq;
            for (unsigned i = 0; i < G::_nbundles; i++) {
                assert(_bundle_vect[i]);
                for (unsigned j = 0; j < G::_nparticles; j++) {
                    string newick = _bundle_vect[i]->getGeneTreeConst(g, j).makeNewick(9, true);
                    freq[newick] += 1;
                }
            }

            // Sort newicks by decreasing count
            vector< pair<unsigned, string> > count_newick(freq.size());
            transform(freq.begin(), freq.end(), count_newick.begin(), [](pair<string, unsigned> p){return make_pair(p.second, p.first);});
            sort(count_newick.begin(), count_newick.end(), greater< pair<unsigned, string> >());
            
            // Save the newicks in order of decreasing count
            for (auto & uniq : count_newick) {
                newicks.push_back(str(format("[freq=%d] %s") % uniq.first % uniq.second));
            }
        }
        else {
            // Save all species trees, even if they are all identical newick strings
            newicks.resize(G::_nparticles);
            for (unsigned i = 0; i < G::_nbundles; i++) {
                assert(_bundle_vect[i]);
                for (unsigned j = 0; j < G::_nparticles; j++) {
                    newicks[i] = _bundle_vect[i]->getGeneTreeConst(g, j).makeNewick(9, true);
                }
            }
        }
    }
    
    inline void SMC::saveGeneTreeLocusToFile(string fn, unsigned g, bool compress) const {
        vector<string> newicks;
        saveGeneTreeLocus(g, newicks, compress);
        
        ofstream treef(fn);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        unsigned i = 1;
        for (auto newick : newicks) {
            treef << str(format("  tree t%d = %s;\n") % (i++) % newick);
        }
        treef << "end;\n";
        treef.close();
    }
    
    inline void SMC::saveLogMargLike(vector<double> & log_marg_like) const {
        log_marg_like.resize(G::_nbundles);
        for (unsigned i = 0; i < G::_nbundles; i++) {
            assert(_bundle_vect[i]);
            log_marg_like[i] = _bundle_vect[i]->getLogMargLike();
        }
    }
    
    inline void SMC::saveLogLikesForLocus(unsigned b, unsigned g, vector<double> & log_likes) const {
        assert(_bundle_vect[b]);
        _bundle_vect[b]->getLogLikesForLocus(g, log_likes);
    }

    inline void SMC::saveLogLikes(unsigned b, vector< vector<double> > & log_likes_for_locus) const {
        assert(_bundle_vect[b]);
        log_likes_for_locus.resize(G::_nloci);
        for (unsigned g = 0; g < G::_nloci; g++) {
            _bundle_vect[b]->getLogLikes(log_likes_for_locus);
        }
    }
    
    inline void SMC::filterBundles(unsigned step) {
        output("Filtering bundles...\n", G::VDEBUG);

        // Copy log weights for all bundles to prob vector
        vector<double> probs(G::_nbundles, 0.0);
        for (unsigned b = 0; b < G::_nbundles; b++) {
            assert(_bundle_vect[b]);
            probs[b] = _bundle_vect[b]->getLogWeight();
        }

        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = G::calcLogSum(probs);
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        G::_log_marg_like += log_sum_weights - log(G::_nbundles);

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts(G::_nbundles, 0);

        // Throw G::_nbundles darts
        for (unsigned i = 0; i < G::_nbundles; ++i) {
            double u = rng->uniform();
            
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            counts[which]++;
        }
        
        // Copy particles

        // Count number of cells with zero count that can serve as copy recipients
        unsigned nzeros = 0;
        for (unsigned i = 0; i < G::_nbundles; i++) {
            if (counts[i] == 0)
                nzeros++;
        }
        
        if (nzeros > 0) {
            // Locate first donor
            unsigned donor = 0;
            while (counts[donor] < 2) {
                donor++;
                assert(donor < counts.size());
            }

            // Locate first recipient
            unsigned recipient = 0;
            while (counts[recipient] != 0) {
                recipient++;
                assert(recipient < counts.size());
            }

            while (nzeros > 0) {
                assert(donor < G::_nbundles);
                assert(recipient < G::_nbundles);

                // Copy donor to recipient
                assert(_bundle_vect[donor]);
                assert(_bundle_vect[recipient]);
                *(_bundle_vect[recipient]) = *(_bundle_vect[donor]);

                counts[donor]--;
                counts[recipient]++;
                nzeros--;

                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < G::_nbundles && counts[donor] < 2) {
                        donor++;
                    }
                }

                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < G::_nbundles && counts[recipient] > 0) {
                    recipient++;
                }
            }
        }
    }
    
    inline unsigned SMC::saveBestSpeciesTree() const {
        unsigned best_bundle = 0;
        double best_log_marg_like = G::_negative_infinity;
        string best_species_tree = "";
        for (unsigned b = 0; b < G::_nbundles; b++) {
            assert(_bundle_vect[b]);
            double log_marg_like = _bundle_vect[b]->report();
            if (log_marg_like > best_log_marg_like) {
                best_bundle = b;
                best_log_marg_like = log_marg_like;
                best_species_tree = _bundle_vect[b]->getSpeciesTreeConst().makeNewick(5, true);
            }
        }
        
        output(format("\nBest log marginal likelihood: %.5f\n") % best_log_marg_like, G::VSTANDARD);
        output(format("Best species tree: %s\n") % best_species_tree, G::VSTANDARD);
        
        // Save best species tree to a treefile
        ofstream treef("best-species-tree.tre");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree best_species_tree [log marg. like. = " << best_log_marg_like << "] = [&R] " << best_species_tree << ";\n";
        treef << "end;\n";
        treef.close();
        
        return best_bundle;
    }
    
    inline void SMC::sanityCheckBundles() const {
        for (unsigned b = 0; b < G::_nbundles; b++) {
            assert(_bundle_vect[b]);
            _bundle_vect[b]->debugSanityCheck();
        }
    }
    
#if defined(STOW_UNUSED_PARTIALS)
    inline void SMC::returnUnusedPartials() {
        {   // scope ensures that unused will be deleted (and the partials it stores as keys
            // will have their use counts decremented) before actually deleting unused partials
            
            // Determine which partials can be deleted
            vector<map<PartialStore::partial_t, int> > unused(G::_nloci);
            for (unsigned g = 0; g < G::_nloci; g++) {
                for (unsigned b = 0; b < G::_nbundles; b++) {
                    assert(_bundle_vect[b]);
                    for (unsigned i = 0; i < G::_nparticles; i++) {
                        _bundle_vect[b]->getGParticleConst(g,i).enumerateUnusedPartials(unused);
                    }
                }
            }

            // If unused[locus][partial] > -1, stow partial
            for (unsigned g = 0; g < G::_nloci; g++) {
                for (auto & mpair : unused[g]) {
                    PartialStore::partial_t     p = mpair.first;
                    int                     count = mpair.second;
                    if (count > -1) {
                        if (!p->_in_storage) {
                            ps.stowPartial(g, p);
                        }
                    }
                }
            }
        }

        // Delete partials that have been stowed
        for (unsigned g = 0; g < G::_nloci; g++) {
            for (unsigned b = 0; b < G::_nbundles; b++) {
                assert(_bundle_vect[b]);
                for (unsigned i = 0; i < G::_nparticles; i++) {
                    _bundle_vect[b]->getGParticle(g,i).deleteStowedPartials();
                }
            }
        }
    }
#endif
    
    inline void SMC::run() {
        // Sanity checks
        assert(G::_nbundles > 0);
        assert(G::_nparticles > 0);
        assert(G::_nloci > 0);

        // Create vector of G::_nbundles Bundle objects
        _bundle_vect.resize(G::_nbundles);
        for (unsigned i = 0; i < G::_nbundles; i++) {
            _bundle_vect[i] = new Bundle();
            assert(_bundle_vect[i]);
            _bundle_vect[i]->setBundleIndex(i);
        }

#if defined(USING_MULTITHREADING)
        unsigned bundles_per_thread = G::_nbundles/G::_nthreads;
        output(format("\nUsing %d threads:\n") % G::_nthreads, G::VSTANDARD);
        unsigned bstart = 0;
        unsigned bend = bundles_per_thread;
        for (unsigned th = 0; th < G::_nthreads; th++) {
            output(format("  thread %d will process bundles %d to %d:\n") % (th+1) % (bstart+1) % bend, G::VSTANDARD);
            bstart = bend;
            bend += bundles_per_thread;
            if (bend > G::_nbundles)
                bend = G::_nbundles;
        }
#else
        output(format("\nSole thread will process bundles 1 to %d\n") % G::_nbundles, G::VSTANDARD);
#endif
        
        // *****************
        // *** Main loop ***
        // *****************
        
        // Each gene tree requires ntaxa-1 steps to be complete
        unsigned nsteps = (unsigned)G::_ntaxa - 1;
        
        // Initialize global log marginal likelihood
        G::_log_marg_like = 0.0;
        
        output("\n", G::VSTANDARD);
        for (unsigned step = 0; step < nsteps; step++) {
            output(format("Step %d of %d\n") % (step + 1) % nsteps, G::VSTANDARD);

            _update_seeds.resize(G::_nbundles);
            G::generateUpdateSeeds(_update_seeds);
            
#if defined(USING_MULTITHREADING)
            bstart = 0;
            bend   = bundles_per_thread;
            vector<thread> threads;
            for (unsigned th = 0; th < G::_nthreads; th++) {
                threads.push_back(
                    thread(&SMC::advanceBundleRange,
                        this,
                        step,
                        th,
                        bstart,
                        bend
                    )
                );
                bstart =  bend;
                bend   += bundles_per_thread;
                if (bend > G::_nbundles)
                    bend = G::_nbundles;
            }

            // The join function causes this loop to pause until the ith thread finishes
            for (unsigned i = 0; i < threads.size(); i++) {
                threads[i].join();
            }
#else
            advanceBundleRange(/*step*/step, /*thread*/0, /*first*/0, /*last*/G::_nbundles);
#endif

            filterBundles(step);
            
#if defined(STOW_UNUSED_PARTIALS)
            returnUnusedPartials();
#endif
        }
        
        unsigned best_bundle = saveBestSpeciesTree();
        assert(_bundle_vect[best_bundle]);
        _bundle_vect[best_bundle]->saveJavascript("newicks-best");
        output(format("log marginal likelihood = %.5f\n") % G::_log_marg_like, G::VSTANDARD);
    }

    inline void SMC::advanceBundleRange(unsigned step, unsigned thread, unsigned first, unsigned last) {
        for (unsigned b = first; b < last; b++) {
            assert(_bundle_vect[b]);
            _bundle_vect[b]->setSeed(_update_seeds[b]);
            if (step == 0) {
                _bundle_vect[b]->initSpeciesTree();
            }
            _bundle_vect[b]->advanceAllGeneTrees();
            _bundle_vect[b]->filterAllGeneTrees(step);
            _bundle_vect[b]->shrinkWrapSpeciesTree();
        }
    }

}
