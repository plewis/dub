#pragma once

namespace proj {

    class SMC {
        public:
            SMC() {}
            ~SMC() {}
            
            void run();
            
            const SParticle & getSpeciesTreeConst(unsigned b) const;
            const GParticle & getGeneTreeConst(unsigned b, unsigned g, unsigned i) const;
            
            void saveSpeciesTrees(vector<string> & newicks) const;
            void saveGeneTrees(unsigned b, unsigned g, vector<string> & newicks) const;
            void saveLogMargLike(vector<double> & log_marg_like) const;
            void saveLogLikesForLocus(unsigned b, unsigned g, vector<double> & log_likes) const;
            void saveLogLikes(unsigned b, vector< vector<double> > & log_likes_for_locus) const;
            
        private:
            void filterBundles(unsigned step);
        
            vector<Bundle>           _bundle;
    };
    
    inline const SParticle & SMC::getSpeciesTreeConst(unsigned b) const {
        assert(b < _bundle.size());
        return _bundle[b].getSpeciesTreeConst();
    }
    
    inline const GParticle & SMC::getGeneTreeConst(unsigned b, unsigned g, unsigned i) const {
        assert(b < _bundle.size());
        assert(g < G::_nloci);
        assert(i < G::_ngparticles);
        return _bundle[b].getGeneTreeConst(g, i);
    }
    
    inline void SMC::saveSpeciesTrees(vector<string> & newicks) const {
        newicks.resize(G::_nsparticles);
        for (unsigned i = 0; i < G::_nsparticles; i++) {
            newicks[i] = _bundle[i].getSpeciesTreeConst().makeNewick(9, true);
        }
    }
    
    inline void SMC::saveGeneTrees(unsigned b, unsigned g, vector<string> & newicks) const {
        newicks.resize(G::_ngparticles);
        for (unsigned i = 0; i < G::_ngparticles; i++) {
            newicks[i] = _bundle[b].getGeneTreeConst(g, i).makeNewick(9, true);
        }
    }
    
    inline void SMC::saveLogMargLike(vector<double> & log_marg_like) const {
        log_marg_like.resize(G::_nsparticles);
        for (unsigned i = 0; i < G::_nsparticles; i++) {
            log_marg_like[i] = _bundle[i].getLogMargLike();
        }
    }
    
    inline void SMC::saveLogLikesForLocus(unsigned b, unsigned g, vector<double> & log_likes) const {
        _bundle[b].getLogLikesForLocus(g, log_likes);
    }

    inline void SMC::saveLogLikes(unsigned b, vector< vector<double> > & log_likes_for_locus) const {
        log_likes_for_locus.resize(G::_nloci);
        for (unsigned g = 0; g < G::_nloci; g++) {
            _bundle[b].getLogLikes(log_likes_for_locus);
        }
    }
    
    inline void SMC::filterBundles(unsigned step) {
        output("Filtering bundles...\n", G::VDEBUG);

        // Copy log weights for all bundles to prob vector
        vector<double> probs(G::_nsparticles, 0.0);
        for (G::_bundle = 0; G::_bundle < G::_nsparticles; G::_bundle++) {
            probs[G::_bundle] = _bundle[G::_bundle].getLogWeight();
        }

        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = G::calcLogSum(probs);
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        //temporary!
        output("\nBundle probabilities:\n", G::VTEMP);
        for (auto p : probs) {
            output(format("%6.3f ") % p, G::VTEMP);
        }
        output("\n", G::VTEMP);

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts(G::_nsparticles, 0);

        // Throw _nparticles darts
        for (unsigned i = 0; i < G::_nsparticles; ++i) {
            double u = rng->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            counts[which]++;
        }

        //temporary!
        output("\nBundle counts:\n", G::VTEMP);
        for (auto c : counts) {
            output(format("%6d ") % c, G::VTEMP);
        }
        output("\n", G::VTEMP);

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

        // Count number of cells with zero count that can serve as copy recipients
        unsigned nzeros = 0;
        for (unsigned i = 0; i < G::_nsparticles; i++) {
            if (counts[i] == 0)
                nzeros++;
        }

        while (nzeros > 0) {
            assert(donor < G::_nsparticles);
            assert(recipient < G::_nsparticles);

            // Copy donor to recipient
            _bundle[recipient] = _bundle[donor];

            counts[donor]--;
            counts[recipient]++;
            nzeros--;

            if (counts[donor] == 1) {
                // Move donor to next slot with count > 1
                donor++;
                while (donor < G::_nsparticles && counts[donor] < 2) {
                    donor++;
                }
            }

            // Move recipient to next slot with count equal to 0
            recipient++;
            while (recipient < G::_nsparticles && counts[recipient] > 0) {
                recipient++;
            }
        }

    }
    
    inline void SMC::run() {
        // Sanity checks
        assert(G::_nsparticles > 0);
        assert(G::_ngparticles > 0);
        assert(G::_nloci > 0);

        // Create vector of G::_nsparticles Bundle objects
        _bundle.resize(G::_nsparticles);
        for (unsigned i = 0; i < G::_nsparticles; i++) {
            _bundle[i].setBundleIndex(i);
        }
        
        // *****************
        // *** Main loop ***
        // *****************
        
        // Each gene tree requires ntaxa-1 steps to be complete
        unsigned nsteps = (unsigned)G::_ntaxa - 1;
        
        output("\n", G::VSTANDARD);
        for (G::_step = 0; G::_step < nsteps; G::_step++) {
            output(format("Step %d of %d\n") % (G::_step + 1) % nsteps, G::VSTANDARD);
            for (G::_bundle = 0; G::_bundle < G::_nsparticles; G::_bundle++) {
                output(format("  Bundle %d\n") % G::_bundle, G::VDEBUG);
                _bundle[G::_bundle].advanceAllGeneTrees();
                _bundle[G::_bundle].filterAllGeneTrees(G::_step);
            }
            filterBundles(G::_step);
        }
        
        double best_log_marg_like = G::_negative_infinity;
        string best_species_tree = "";
        for (G::_bundle = 0; G::_bundle < G::_nsparticles; G::_bundle++) {
            //_bundle[G::_bundle].saveJavascript(str(format("newicks-%d") % G::_bundle));
            double log_marg_like = _bundle[G::_bundle].report();
            if (log_marg_like > best_log_marg_like) {
                best_log_marg_like = log_marg_like;
                best_species_tree = _bundle[G::_bundle].getSpeciesTreeConst().makeNewick(5, true);
            }
        }
        output(format("\nBest log marginal likelihood: %.5f\n") % best_log_marg_like, G::VSTANDARD);
        output(format("Best species tree: %s\n") % best_species_tree, G::VSTANDARD);
    }
}
