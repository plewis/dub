#pragma once

namespace proj {

    class ParallelPolicyMultithreading {
        private:
            vector<pair<unsigned, unsigned> > _thread_schedule;

            void debugCheckThreadSchedule() const;
            void debugShowThreadSchedule(double percent_of_max_entropy) const;
            void divideLargestParticle();
            double buildThreadSchedule();
            void advanceParticleRange(unsigned step, unsigned first_particle, unsigned last_particle, const vector<unsigned> & update_seeds);
        
        protected:
            bool multithreading() {return true;}
            void particleLoop(unsigned step, const vector<unsigned> & update_seeds);
            void balanceThreads();
            
            mutex  _mutex;
            mutex  _gene_forest_clear_mutex;
            mutex  _debug_mutex;            
    };
    
    inline void ParallelPolicyMultithreading::advanceParticleRange(unsigned step, unsigned first_particle, unsigned last_particle, const vector<unsigned> & update_seeds) {
        // Uncomment line below to ensure each thread gets through this entire function
        // uninterrupted. This is good for debugging but of course non-sensical for a
        // multithreading application.
        //lock_guard<mutex> guard(G::_debug_mutex);
        
        // Get iterator to start of range
        auto it_start = _particle_list.begin();
        advance(it_start, first_particle);
        
        int first_index = it_start->getBeginIndex();
        assert(first_index >= 0 && first_index < _log_weights.size());

        // Get iterator to (one past) end of range
        auto it_end = _particle_list.begin();
        advance(it_end, last_particle);
        
        // Loop over particles and process each
        unsigned i = first_index;
        for (auto it = it_start; it != it_end; ++it) {
            Particle & p = *it;
            
            // Second element of pair is number of copies of particle
            unsigned n = p.getCount();
        
            while (n > 0) {
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugStoreForests();
#endif
                
                // Propose a coalescence event (which may involve also proposing
                // one or more intervening speciation events)
                auto proposed = p.proposeCoalescence(update_seeds[i], step, i,  /*compute_partial*/true, /*make_permanent*/false);
                
                // Save the log weight
#if defined(UPGMA_WEIGHTS)
                throw XProj("UPGMA_WEIGHTS not yet implemented for ParallelPolicyMultithreading");
#else
                _log_weights[i] = G::_phi*proposed.first;
#endif
                
                // This is the current number of lineages in the species tree
                // This information is used only for sanity checking (should
                // be the same after reconstructing this coalescent event as
                // it was the first time
                //_proposed_species_tree_lineages[i] = proposed.second;
                
                // Save information necessary for recreating this proposed coalescence,
                // namely the gene, species within gene, and first and second node
                // indices joined within the species
                //_proposed_gene[i] = p.getLastProposedGene();
                //_proposed_spp[i] = p.getLastProposedSpecies();
                //_proposed_first[i] = p.getLastProposedFirstIndex();
                //_proposed_second[i] = p.getLastProposedSecondIndex();
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugCheckForests();
#endif
                
                n--;
                i++;
            }
        }
    }

    inline void ParallelPolicyMultithreading::particleLoop(unsigned step, const vector<unsigned> & update_seeds) {
        //MARK: Proj::particleLoopMT
        vector<thread> threads;
        for (unsigned i = 0; i < G::_nthreads; i++) {
            threads.push_back(thread(&Proj::advanceParticleRange,
                this,
                step,
                thread_schedule[i].first,
                thread_schedule[i].second,
                update_seeds)
            );
        }

        // The join function causes this loop to pause until the ith thread finishes
        for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
        }
    }
    
    inline void ParallelPolicyMultithreading::balanceThreads() {
        // Make copies of particles as necessary to achieve a more even distribution of work.
        // Example: 2 threads, 100 (virtual) particles, 4 (actual) particles
        //
        // Original particles:
        //  index  count  cumulative   prefix-sum   calculation  thread     count
        // ------ ------ ------------ ----------- ------------- ------- ---------
        //      0     11          11           11 2*10/100=0.20       0
        //      1     25          25           36 2*35/100=0.70       0  11+25=36
        //      2     48          48           84 2*83/100=1.66       1
        //      3     16          16          100 2*99/100=1.98       1  48+16=64
        // ------ ------ ------------ ----------- ------------- ------- ---------
        // entropy      = 0.653 = -(0.36*ln(0.36) + 0.64*ln(0.64))
        // max(entropy) = 0.693 = log(2)
        // percentage   = 94.2
        //
        // Continue dividing largest particle until % max entropy > target (e.g. 99.6)
        //   99.971 = -100.0*(0.49*ln(0.49) + 0.51*ln(0.51))/ln(2)
        //   99.885 = -100.0*(0.48*ln(0.48) + 0.52*ln(0.52))/ln(2)
        //   99.740 = -100.0*(0.47*ln(0.47) + 0.53*ln(0.53))/ln(2) *
        //   99.538 = -100.0*(0.46*ln(0.46) + 0.54*ln(0.54))/ln(2)
        //   99.277 = -100.0*(0.45*ln(0.45) + 0.55*ln(0.55))/ln(2)
        assert(G::_nthreads > 0);
        double pct = 0.0;
        if (G::_nthreads == 1) {
            _thread_schedule.clear();
            _thread_schedule.push_back(make_pair(0,(unsigned)_particle_list.size()));
        }
        else {
            while (_particle_list.size() < G::_nthreads)
                divideLargestParticle();
            pct = buildThreadSchedule(_thread_schedule);
            while (pct < _entropy_percent_cutoff) {
                divideLargestParticle();
                pct = buildThreadSchedule(_thread_schedule);
            }
        }
        
        //debugShowThreadSchedule(_thread_schedule, pct);
        //cerr << endl;
    }

    inline void ParallelPolicyMultithreading::debugCheckThreadSchedule(const vector<pair<unsigned, unsigned> > & thread_schedule) const {
    
        //auto it = _particle_list.begin();
        //for (auto & p : thread_schedule) {
        //    unsigned begin_particle  = p.first;
        //    unsigned end_particle = p.second;
        //    unsigned begin_index = it->getBeginIndex();
        //}
        // begin again here
    }
    
    inline void ParallelPolicyMultithreading::debugShowThreadSchedule(const vector<pair<unsigned, unsigned> > & thread_schedule, double percent_of_max_entropy) const {
        unsigned i = 0;
        auto it = _particle_list.begin();
        unsigned sum_pcount = 0;
        output(format("\n%12s %12s %12s %12s %12s %12s\n") % "index" % "begin" % "end" % "begin-index" % "end-index" % "count", 2);
        for (auto & p : thread_schedule) {
            unsigned begin_particle  = p.first;
            unsigned end_particle = p.second;
            unsigned begin_index     = it->getBeginIndex();
            unsigned n = end_particle - begin_particle;
            unsigned pcount = 0;
            for (unsigned j = 0; j < n; j++) {
                pcount += it->getCount();
                ++it;
            }
            sum_pcount += pcount;
            unsigned end_index = begin_index + pcount;
            output(format("%12d %12d %12d %12d %12d %12d\n") % i % begin_particle % end_particle % begin_index % end_index % pcount, 2);
            i++;
        }
        assert(it == _particle_list.end());
        output(format("Total actual particles:  %d\n") % _particle_list.size(), 2);
        output(format("Total virtual particles: %d\n") % sum_pcount, 2);
        output(format("Percentage of maximum entropy: %.1f\n") % percent_of_max_entropy, 2);
        cerr << endl;
    }

    inline void ParallelPolicyMultithreading::divideLargestParticle() {
        auto it = max_element(_particle_list.begin(), _particle_list.end(), [](const Particle & first, const Particle & second){return first.getCount() < second.getCount() ? true : false;});

        unsigned begin_index = it->getBeginIndex();
        unsigned count = it->getCount();
        unsigned first_half = count/2;
        unsigned second_half = count - first_half;

        auto it0 = _particle_list.insert(it, *it);
                
        it0->setCount(first_half);
        it0->setXtra(0);
        it0->setBeginIndex(begin_index);
        
        it->setCount(second_half);
        it->setXtra(0);
        it->setBeginIndex(begin_index + first_half);
    }
    
    inline double ParallelPolicyMultithreading::buildThreadSchedule(vector<pair<unsigned, unsigned> > & thread_schedule) {
        // Calculate thread_schedule and its associated entropy and return the percentage
        // entropy relative to the maximum possible entropy

        // Create thread schedule using the prefix-sum algorithm
        thread_schedule.clear();
        unsigned current_thread = 0;
        pair<unsigned, unsigned> begin_end = make_pair(0,0);
        unsigned prefix_sum = 0;
        unsigned max_count = (unsigned)floor(1.0*G::_nparticles/G::_nthreads);
        
        vector<unsigned> freqs(G::_nthreads, 0);
        for (auto & p : _particle_list) {
            // Add particle count to prefix sum
            unsigned count = p.getCount();
            if (count > max_count) {
                thread_schedule.clear();
                return 0.0;
            }
            
            p.setBeginIndex(prefix_sum);
            prefix_sum += count;
                
            // Calculate thread index
            unsigned thread_index = (unsigned)floor(1.0*G::_nthreads*(prefix_sum - 1)/G::_nparticles);
            
            if (thread_index > current_thread) {
                // Add thread to the schedule
                thread_schedule.push_back(begin_end);
                current_thread++;
                
                // Start work on next thread
                begin_end.first  = begin_end.second;
                begin_end.second = begin_end.first + 1;
            }
            else {
                begin_end.second += 1;
            }
            freqs[current_thread] += count;
        }
        thread_schedule.push_back(begin_end);

        // Calculate entropy, max_entropy, and percentage
        double max_entropy = log(G::_nthreads);
        double entropy = 0.0;
        assert(G::_nparticles == accumulate(freqs.begin(), freqs.end(), 0));
        for_each(freqs.begin(), freqs.end(), [&entropy](unsigned f){entropy -= 1.0*f*log(f);});
        entropy /= G::_nparticles;
        entropy += log(G::_nparticles);
        double percentage = 100.0*entropy/max_entropy;
        assert(!isnan(percentage));
        
        return percentage;
    }

}

