#pragma once

namespace proj {

    class SMC {
        public:
            SMC() {}
            ~SMC() {}
            
            void run();
            
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
            void filterBundles(unsigned step);
            void sanityCheckBundles() const;
        
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
    
    inline void SMC::saveSpeciesTrees(vector<string> & newicks, bool compress) const {
        if (compress) {
            // Save only unique newick species tree descriptions
            // May save same topology many times if edge lengths differ

            // Create map with newicks as keys and counts as values
            newicks.clear();
            map<string, unsigned> freq;
            for (unsigned i = 0; i < G::_nsparticles; i++) {
                string newick = _bundle[i].getSpeciesTreeConst().makeNewick(9, true);
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
            newicks.resize(G::_nsparticles);
            for (unsigned i = 0; i < G::_nsparticles; i++) {
                newicks[i] = _bundle[i].getSpeciesTreeConst().makeNewick(9, true);
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
        if (compress) {
            // Save only unique newick species tree descriptions
            // May save same topology many times if edge lengths differ
            
            // Create map with newicks as keys and counts as values
            newicks.clear();
            map<string, unsigned> freq;
            for (unsigned i = 0; i < G::_ngparticles; i++) {
                string newick = _bundle[b].getGeneTreeConst(g, i).makeNewick(9, true);
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
            newicks.resize(G::_ngparticles);
            for (unsigned i = 0; i < G::_nsparticles; i++) {
                newicks[i] = _bundle[b].getGeneTreeConst(g, i).makeNewick(9, true);
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
            for (unsigned i = 0; i < G::_nsparticles; i++) {
                for (unsigned j = 0; j < G::_ngparticles; j++) {
                    string newick = _bundle[i].getGeneTreeConst(g, j).makeNewick(9, true);
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
            newicks.resize(G::_ngparticles);
            for (unsigned i = 0; i < G::_nsparticles; i++) {
                for (unsigned j = 0; j < G::_ngparticles; j++) {
                    newicks[i] = _bundle[i].getGeneTreeConst(g, j).makeNewick(9, true);
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
    
//    inline void SMC::saveGeneTrees(unsigned bundle_index, unsigned locus, vector<string> & newicks, bool compress) {
//        typedef tuple<unsigned, double, string, string, string> treeinfo_t;
//        treeinfo_t treeinfo;
//        vector<treeinfo_t> treeinfo_vect;
//        if (compress) {
//            // Save only unique newick strings
//            map<string,vector<GeneTreeDetails> > tree_info;
//            for (const Particle & p : particle_list) {
//                GeneTreeDetails info;
//
//                // Get count for this particle
//                info._count = p.getCount();
//
//                // Get newick tree description for this gene tree
//                assert(gene_index < p.getGeneForestsConst().size());
//                const GeneForest & gf = p.getGeneForestsConst()[gene_index];
//                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
//
//                // Calculate log-likelihood for this gene tree
//                info._log_likelihood = gf.calcLogLikelihood();
//
//                // Record everything for this particle
//                tree_info[newick].push_back(info);
//            }
//
//            ofstream tmpf(str(format("%s.txt") % fnprefix));
//            unsigned i = 0;
//            for (auto it = tree_info.begin(); it != tree_info.end(); ++it) {
//                const string & newick = it->first;
//                vector<GeneTreeDetails> & details_vect = it->second;
//
//                // Calculate total count (i.e. freq)
//                unsigned total_c = 0;
//                for (auto v : details_vect) {
//                    total_c += v._count;
//                }
//
//                tmpf << total_c << " <-- " << newick << endl;
//                unsigned k = 0;
//                for (auto v : details_vect) {
//                    unsigned c = v._count;
//                    tmpf << "  particle: " << k << endl;
//                    tmpf << "    count: " << c << endl;
//                    tmpf << "    log_like: " << str(format("%.9f") % v._log_likelihood) << endl;
//                    tmpf << endl;
//                    k++;
//                }
//
//                double pct = 100.0*total_c/_nparticles;
//                string note = str(format("freq = %d") % total_c);
//                string treename = str(format("'tree%d-freq%d'") % i % total_c);
//                treeinfo = make_tuple(total_c, pct, note, treename, newick);
//                treeinfo_vect.push_back(treeinfo);
//                ++i;
//            }
//            tmpf.close();
//        }
//        else {
//            unsigned i = 0;
//            for (const Particle & p : particle_list) {
//                unsigned c = p.getCount();
//                assert(gene_index < p.getGeneForestsConst().size());
//                const GeneForest & gf = p.getGeneForestsConst()[gene_index];
//                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
//                if (compression_level == 1) {
//                    double pct = 100.0*c/_nparticles;
//                    string note = str(format("freq = %d") % c);
//                    string treename = str(format("'tree%d-freq%d'") % i % c);
//                    treeinfo = make_tuple(c, pct, note, treename, newick);
//                    treeinfo_vect.push_back(treeinfo);
//                    i++;
//                }
//                else {
//                    double pct = 100.0/_nparticles;
//                    string note = "freq = 1";
//                    for (unsigned j = 0; j < c; j++) {
//                        string treename = str(format("'tree%d-freq1'") % i);
//                        treeinfo = make_tuple(c, pct, note, treename, newick);
//                        treeinfo_vect.push_back(treeinfo);
//                        i++;
//                    }
//                }
//            }
//        }
//        string fn = str(format("%s.tre") % fnprefix);
//        outputAnnotatedNexusTreefile(fn, treeinfo_vect);
//    }
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
        
        G::_log_marg_like += log_sum_weights - log(G::_nsparticles);

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
    
    inline unsigned SMC::saveBestSpeciesTree() const {
        unsigned best_bundle = 0;
        double best_log_marg_like = G::_negative_infinity;
        string best_species_tree = "";
        for (G::_bundle = 0; G::_bundle < G::_nsparticles; G::_bundle++) {
            double log_marg_like = _bundle[G::_bundle].report();
            if (log_marg_like > best_log_marg_like) {
                best_bundle = G::_bundle;
                best_log_marg_like = log_marg_like;
                best_species_tree = _bundle[G::_bundle].getSpeciesTreeConst().makeNewick(5, true);
                
                // ///temporary!
                //_bundle[G::_bundle].getSpeciesTreeConst().sanityCheck(0, G::_bundle);
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
        for (G::_bundle = 0; G::_bundle < G::_nsparticles; G::_bundle++) {
            _bundle[G::_bundle].debugSanityCheck();
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
        
        // //temporary!
        // ofstream tmpf1("after-trimming.tre");
        // tmpf1 << "#nexus\n\n";
        // tmpf1 << "begin trees;\n";
        // tmpf1.close();
        //
        // //temporary!
        // ofstream tmpf2("after-trimming.tre");
        // tmpf2 << "#nexus\n\n";
        // tmpf2 << "begin trees;\n";
        // tmpf2.close();
        
        // Each gene tree requires ntaxa-1 steps to be complete
        unsigned nsteps = (unsigned)G::_ntaxa - 1;
        
        // Initialize global log marginal likelihood
        G::_log_marg_like = 0.0;
        
        output("\n", G::VSTANDARD);
        for (G::_step = 0; G::_step < nsteps; G::_step++) {
            output(format("Step %d of %d\n") % (G::_step + 1) % nsteps, G::VSTANDARD);
            
            for (G::_bundle = 0; G::_bundle < G::_nsparticles; G::_bundle++) {
                output(format("  Bundle %d\n") % G::_bundle, G::VDEBUG);
                
                // //temporary!
                // _bundle[G::_bundle].getSpeciesTreeConst().sanityCheck(G::_step, G::_bundle);
                // cerr << "before: " << _bundle[G::_bundle].getSpeciesTreeConst().makeNewick(7,true) << endl;
                
                _bundle[G::_bundle].advanceAllGeneTrees();
                
                // //temporary!
                // cerr << "after: " << _bundle[G::_bundle].getSpeciesTreeConst().makeNewick(7,true) << endl;
                
                _bundle[G::_bundle].filterAllGeneTrees(G::_step);
                _bundle[G::_bundle].shrinkWrapSpeciesTree();
            }
            filterBundles(G::_step);
            
            // //temporary!
            // sanityCheckBundles();
        }
        
        //  //temporary!
        //  ofstream tmpf3("after-trimming.tre", ios::out | ios::app);
        //  tmpf3 << "end;\n";
        //  tmpf3.close();
        
        // //temporary!
        // ofstream tmpf4("after-trimming.tre", ios::out | ios::app);
        // tmpf4 << "end;\n";
        // tmpf4.close();
        
        unsigned best_bundle = saveBestSpeciesTree();
        _bundle[best_bundle].saveJavascript("newicks-best");
        output(format("log marginal likelihood = %.5f\n") % G::_log_marg_like, G::VSTANDARD);
    }
    
}
