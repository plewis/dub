#pragma once

namespace proj {

    class ParallelPolicyMPI {
        protected:
            bool multithreading() {return false;}
            void particleLoop(unsigned step, const vector<unsigned> & update_seeds);
            void balanceThreads() {}
    };
    
    inline void ParallelPolicyMPI::advanceParticleRange(unsigned step, unsigned first_particle, unsigned last_particle, const vector<unsigned> & update_seeds) {
        // Uncomment line below to ensure each thread gets through this entire function
        // uninterrupted. This is good for debugging but of course non-sensical for a
        // multithreading application.
        //lock_guard<mutex> guard(SMCGlobal::_debug_mutex);
        
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
        
#if defined(MINIMIZE_PARTIALS)
            // Recompute all partials
            p.computeAllPartials();
#endif
        
            while (n > 0) {
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugStoreForests();
#endif
                
                // Propose a coalescence event (which may involve also proposing
                // one or more intervening speciation events)
                auto proposed = p.proposeCoalescence(update_seeds[i], step, i,  /*compute_partial*/true, /*make_permanent*/false);
                _log_weights[i] = SMCGlobal::_phi*proposed.first;
                _proposed_species_tree_lineages[i] = proposed.second;
                _proposed_gene[i] = p.getLastProposedGene();
                _proposed_spp[i] = p.getLastProposedSpecies();
                _proposed_first[i] = p.getLastProposedFirstIndex();
                _proposed_second[i] = p.getLastProposedSecondIndex();
                
#if defined(DEBUGGING_SANITY_CHECK)
                p.debugCheckForests();
#endif
                
                n--;
                i++;
            }
            
#if defined(MINIMIZE_PARTIALS)
            p.stowAllPartials();
#endif
        }
    }
    
    inline void ParallelPolicyMPI::particleLoop(unsigned step, const vector<unsigned> & update_seeds) {
        // Advance each particle by one coalescent event
        for (unsigned k = ::my_first_particle; k < ::my_last_particle; ++k) {
            Particle & p = _particles[k];
            p.advance(step, k, /*calculate_partial*/true);
            while (p.lastEventSpeciation()) {
                p.advance(step, k, /*calculate_partial*/true);
            }
        }
        
        // Erect barrier causing this loop to pause
        // until the last process finishes
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (my_rank == 0) {
            // The outstanding variable equals the number of
            // processors we haven't heard from yet
            unsigned outstanding = ntasks;
            
            // Receive particle weights from other processors
            while (outstanding) {
                // Probe to get message status
                int message_length = 0;
                MPI_Status status;
                MPI_Probe(MPI_ANY_SOURCE,   // Source rank or MPI_ANY_SOURCE
                    MPI_ANY_TAG,            // Message tag
                    MPI_COMM_WORLD,         // Communicator
                    &status                 // Status object
                );
            
                // Get length of message
                MPI_Get_count(&status, MPI_DOUBLE, &message_length);

                // Get processor sending message
                unsigned which = (unsigned)status.MPI_TAG;

                // Get the message itself
                vector<double> log_weights(_mpi_num_particles[which], 0.0);
                MPI_Recv(&log_weights[0],    // Initial address of receive buffer
                    message_length,     // Maximum number of elements to receive
                    MPI_DOUBLE,           // Datatype of each receive buffer entry
                    status.MPI_SOURCE,     // Rank of source
                    status.MPI_TAG,        // Message tag
                    MPI_COMM_WORLD,     // Communicator
                    MPI_STATUS_IGNORE   // Status object
                );
                
                // Copy log_weights to appropriate section of
                // the _log_weights vector
                unsigned kfrom = 0;
                for (unsigned kto = ::my_first_particle; kto < ::my_last_particle; ++kto) {
                    assert(kfrom < log_weights.size());
                    _log_weights[kto] = log_weights[kfrom];
                }
                
                --outstanding;
            }
        }
    }

    inline void ParallelPolicyMPI::balanceThreads() {
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
        assert(SMCGlobal::_nthreads > 0);
        double pct = 0.0;
        if (SMCGlobal::_nthreads == 1) {
            _thread_schedule.clear();
            _thread_schedule.push_back(make_pair(0,(unsigned)_particle_list.size()));
        }
        else {
            while (_particle_list.size() < SMCGlobal::_nthreads)
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

}

