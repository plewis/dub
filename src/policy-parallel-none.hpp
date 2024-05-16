#pragma once

namespace proj {

    template <class T>
    class ParallelPolicyNone {
        protected:
            bool multithreading() {return false;}
            void particleLoop(unsigned step, const vector<unsigned> & update_seeds);
            void balanceThreads() {}
    };
    
    template <class T>
    inline void ParallelPolicyNone<T>::particleLoop(unsigned step, const vector<unsigned> & update_seeds) {
        // Get reference to the SMC object to which this policy applies
        T & smc = static_cast<T &>(*this);
        
        // If a particle has n identical copies, then each of these n copies need to be advanced
        // Since many of them will be later filtered out, we store only the log weight to use in
        // filtering and the instructions needed to recreate the advance (in case the particle
        // has a log weight high enough to allow it to survive filtering).
        
        unsigned i = 0; // index of proposal
        unsigned j = 0; // index of parent particle
        for (auto & p : smc._particle_list) {
            if (smc.isConditionalMode()) {
                // Compute _prev_log_likelihood
                double check_prev_log_coallike = p.getPrevLogCoalLike();
                vector<Forest::coalinfo_t> coalinfo_vect;
                p.recordAllForests(coalinfo_vect);
                double log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);

                //temporary!
                if (fabs(check_prev_log_coallike - log_coallike) > 0.0001) {
                    cerr << "***********************" << endl;
                    cerr << "* previous log(coalescent likelihoods) doesn't match" << endl;
                    cerr << "* check_prev_log_coallike: " << check_prev_log_coallike << endl;
                    cerr << "* log_coallike:            " << log_coallike << endl;
                    cerr << "***********************" << endl;
                    Forest::debugShowCoalInfo("===== coalinfo =====", coalinfo_vect);
                }
                
                p.setPrevLogCoalLike(log_coallike);
            }

            // n is the number of copies of particle p
            unsigned n = p.getCount();
            
#if defined(MINIMIZE_PARTIALS)
            // Recompute all partials
            p.computeAllPartials();
#endif
            
            while (n > 0) {
            
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugStoreForests();
#endif
                pair<double, unsigned> proposed;
                if (smc.isJointMode()) {
                    // Propose a coalescence event (which may involve also proposing
                    // one or more intervening speciation events)
                    // proposed is a pair in which first is the log_weight and
                    // second is the current number of species tree lineages
                    proposed = p.proposeCoalescence(update_seeds[i], step, i, /*compute_partial*/true, /*make_permanent*/false);
                    
                    // This is the current number of lineages in the species tree
                    // This information is used only for sanity checking (should
                    // be the same after reconstructing this coalescent event as
                    // it was the first time
                    //smc._proposed_species_tree_lineages[i] = proposed.second;

                    // Save proposed state so that it can be reconstructed if
                    // that particle is chosen during filtering
                    //smc._proposed_gene[i] = p.getLastProposedGene();
                    //smc._proposed_spp[i] = p.getLastProposedSpecies();
                    //smc._proposed_first[i] = p.getLastProposedFirstIndex();
                    //smc._proposed_second[i] = p.getLastProposedSecondIndex();
                }
                else if (smc.isConditionalMode()) {
                    //p.debugShowMarkVariables(str(format("particleLoop before proposeSpeciation: i = %d, j = %d, n = %d") % i % j % n));
                    proposed = p.proposeSpeciation(update_seeds[i], step, i, /*make_permanent*/false);
                    //p.debugShowMarkVariables("particleLoop after proposeSpeciation");
                }
                else {
                    throw XProj(format("Unknown SMC mode encountered (%d) in ParallelPolicyNone<T>::particleLoop function") % smc._mode);
                }
                
                // Store log weight (phi is heating factor, default = 1.0)
                smc._log_weights[i] = G::_phi*proposed.first;
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugCheckForests();
#endif

                n--;
                i++;
            }
            j++;    // parent particle index
            
#if defined(MINIMIZE_PARTIALS)
            p.stowAllPartials();
#endif
        }
        assert(i == smc._nparticles);
    }

}

