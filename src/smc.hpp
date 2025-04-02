#pragma once

extern proj::StopWatch stopwatch;
extern proj::Lot::SharedPtr rng;

namespace proj {

    class Particle;

    class SMC {
        public:
                                     SMC()  {clear();}
            virtual                  ~SMC() {}
                        
            enum mode_type_t {
                SPECIES_AND_GENE = 0,
                SPECIES_GIVEN_GENE = 1
            };

            void                     setMode(mode_type_t m)     {_mode = m;}
            bool                     isJointMode() const        {return _mode == SPECIES_AND_GENE;}
            bool                     isConditionalMode() const  {return _mode == SPECIES_GIVEN_GENE;}

            void                     setNParticles(unsigned nparticles, unsigned nsubpops)  {_nparticles = nparticles; _nsubpops = nsubpops; }
            
            void                     setData(Data::SharedPtr d) {_data = d;}
            void                     initFromParticle(Particle & p);
            void                     init();
            void                     run();

            double                   filterParticles(unsigned step, int locus);
            double                   filterParticlesWithinSubpops(unsigned step, int locus);

            vector<Particle> &       getParticles()  {return _particles;}
            const vector<Particle> & getParticlesConst() const {return _particles;}
            
            double                   calcLogSpeciesTreePrior(vector<Forest::coalinfo_t> & coalinfo_vect, bool include_join_probs) const;
            void                     classifyCounts(vector<unsigned> & zeros, vector<unsigned> & nonzeros, const vector<unsigned> & counts) const;
            double                   computeEffectiveSampleSize(const vector<double> & probs) const;
            bool                     compareToReferenceTrees(vector<Particle> particles, map<string, tuple<unsigned, double, double, double, double> > & m);
            void                     outputAnnotatedNexusTreefile(string fn, const vector<tuple<unsigned, double, string, string, string> > & treeinfo) const;
            void                     saveAllSpeciesTrees(string fn, const vector<Particle> & particles, unsigned compression_level = 2);
            void                     saveAllGeneTrees(unsigned gene, string fn, vector<Particle> & particles, unsigned compression_level = 2);
            void                     summarize();
            void                     dumpParticles(SMC & ensemble, vector<unsigned> & kept);
            void                     clear();
            unsigned                 countDistinctGeneTreeTopologies();
            
            typedef shared_ptr<SMC>  SharedPtr;
            
        private:
        
            unsigned                 _mode;
            Data::SharedPtr          _data;
            vector<unsigned>         _counts;
            
            unsigned                 _nparticles;
            unsigned                 _nsubpops;

            // _subpop_id[k] is subpop in which kth particle belongs
            vector<unsigned>         _subpop_id;
            
            unsigned                 _nsteps;
            vector<double>           _starting_log_likelihoods;

            vector<Particle>         _particles;
            double                   _log_marg_like;
            vector<double>           _log_weights;
    };
        
    inline void SMC::clear() {
        _mode = SPECIES_AND_GENE;
        _nsteps = 0;
        _log_marg_like = 0.0;
        _nparticles = 0;
        _nsubpops = 1;
        _counts.clear();
        _subpop_id.clear();
    }
    
    inline void SMC::initFromParticle(Particle & p) {
        assert(_nparticles > 0);
        
#if defined(DEBUGGING_INITFROMPARTICLE)
        // output first level gene trees and species tree
        output(format("\nSpecies forest:\n%s\n") % p.getSpeciesForest().makeNewick(9, true, false));
        auto gene_forests = p.getGeneForestsConst();
        for (auto & gf : gene_forests) {
            output(format("\nGene forest:\n%s\n") % gf.makeNewick(9, true, false));
        }
#endif
        _data = p.getData();
        
        // This function should only be used when estimating species
        // trees from complete gene trees (obtained from p)
        assert(_mode == SPECIES_GIVEN_GENE);
        
        // The _counts vector stores the number of darts that hit
        // each of the _nparticles in multinomial sampling.
        _counts.clear();
        _counts.resize(_nparticles);
        
        // The _subpop_id vector stores the assignments of particles to subpops
        _subpop_id.clear();
        _subpop_id.resize(_nparticles);
        unsigned first = 0;
        unsigned particles_per_subpop = _nparticles/_nsubpops;
        for (unsigned i = 0; i < _nsubpops; i++) {
            unsigned last = first + particles_per_subpop;
            for (unsigned j = first; j < last; j++)
                _subpop_id[j] = i;
        }
    
        // Create log_weights vector that stores log weight of each particle
        _log_weights.resize(_nparticles);
        
        // Replace existing species tree with trivial forest
        p.resetSpeciesForest();
        
        // Compute starting log coalescent likelihood
        vector<Forest::coalinfo_t> coalinfo_vect;
        p.recordAllForests(coalinfo_vect);
#if defined(DEBUG_COALLIKE)
        p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/true);
#else
        p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
#endif
        p.resetPrevLogCoalLike();
            
        // Initialize particle list with _nparticles Particles, each of which
        // is a copy of template_particle
        assert(_particles.size() == 0);
        _particles.resize(_nparticles, p);
        
        // Determine total number of steps required to build the species tree
        assert(G::_nspecies > 1);
        _nsteps = G::_nspecies - 1;

        output(format("\nSMC (level 2) will require %d steps.\n") % _nsteps, G::LogCateg::SECONDLEVEL);

        output(format("\n%12s %12s %24s %12s %12s\n") % "Step" % "ESS" % "logml" % "secs" % "wait", G::LogCateg::SECONDLEVEL);
    }
    
    inline void SMC::init() {
        assert(_nparticles > 0);
        assert(_data);
        
        // Stores number of darts that hit each of the _nparticles in
        // multinomial sampling.
        _counts.clear();
        _counts.resize(_nparticles);
        
        // The _subpop_id vector stores the assignments of particles to subpops
        _subpop_id.clear();
        _subpop_id.resize(_nparticles);
        unsigned first = 0;
        unsigned particles_per_subpop = _nparticles/_nsubpops;
        for (unsigned i = 0; i < _nsubpops; i++) {
            unsigned last = first + particles_per_subpop;
            for (unsigned j = first; j < last; j++) {
                _subpop_id[j] = i;
            }
            first = last;
        }
    
        // Create log_weights vector that stores log weight of each particle
        _log_weights.resize(_nparticles);
        
        // Initialize template particle with trivial forests
        Particle template_particle;
        template_particle.setData(_data);
        template_particle.resetSpeciesForest();
        template_particle.resetGeneForests(/*compute_partials*/true);
        //template_particle.resetDistanceMatrix();

        // Compute initial log likelihood
        _starting_log_likelihoods.resize(G::_nloci, 0.0);
        template_particle.calcLogLikelihood();
        template_particle.resetAllPrevLogLikelihood();
        const vector<GeneForest> & gene_forests = template_particle.getGeneForestsConst();
        for (unsigned g = 0; g < G::_nloci; g++) {
            _starting_log_likelihoods[g] = gene_forests[g].getPrevLogLikelihood();
        }

        // Initialize particles vector
        assert(_particles.size() == 0);
        _particles.resize(_nparticles, template_particle);
                        
        // Determine total number of steps required to build all gene trees
        _nsteps = G::_nloci*(G::_ntaxa - 1);
        output(format("\nSMC (level 1) will require %d steps.\n") % _nsteps, 1);
        
        G::showSettings();
        
        // Display header for progress table
        output(format("\n%12s %12s  %24s %12s %12s\n") % "Step" % "ESS" % "logml" % "secs" % "wait");
    }

    inline void SMC::run() {
        assert(_nsteps > 0);
                
        _log_marg_like = 0.0;
        double cum_secs = 0.0;
        for (unsigned step = 0; step < _nsteps; ++step) {
            stopwatch.start();
            unsigned locus = G::_nloci;
            
            if (isJointMode()) {
                
#if defined(RANDOM_LOCUS_ORDERING)
#               error random locus ordering not yet implemented except for sim
#else
                locus = step % G::_nloci;
                assert(locus < G::_nloci);
#endif
            
                for (unsigned i = 0; i < _nparticles; i++) {
                    // Advance each particle by one coalescent event in one locus
                    _particles[i].proposeCoalescence(step, locus);
                }
            }
            else if (isConditionalMode()) {
                for (unsigned i = 0; i < _nparticles; i++) {
                    // Advance each particle by one speciation event
                    _particles[i].proposeSpeciation(step);
                }
            }
            else {
                throw XProj(format("Unknown SMC mode encountered (%d) in ParallelPolicyNone<T>::particleLoop function") % _mode);
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            if (isJointMode()) {
                if (_nsubpops == 1)
                    ess = filterParticles(step, locus);
                else
                    ess = filterParticlesWithinSubpops(step, locus);
            }
            else {
                ess = filterParticles(step, -1);
            }
            
            double secs = stopwatch.stop();
            
            // Calculate waiting time until done
            cum_secs += secs;
            double avg_per_step = cum_secs/(step + 1);
            unsigned steps_to_go = _nsteps - (step + 1);
            double wait = avg_per_step*steps_to_go;

            if (!isConditionalMode()) {
                VALGRIND_PRINTF("~~> post 1st-level step %d at time %d\n", step, (unsigned)clock());
                VALGRIND_MONITOR_COMMAND(str(format("detailed_snapshot stepsnaps-%d.txt") % step).c_str());
                //ps.debugReport();
            }

            unsigned verbosity = isConditionalMode() ? G::LogCateg::SECONDLEVEL : G::LogCateg::INFO;
            output(format("%12d %12.3f %24.6f %12.3f %12.3f\n") % (step+1) % ess % _log_marg_like % secs % wait, verbosity);
        }

        ps.debugReport();
    }
    
    inline void SMC::summarize() {
        assert(_particles.size() > 0);
        
        string prefix = "1st";
        if (isConditionalMode())
            prefix = "2nd";
                    
        // //temporary!
        // Particle & p = *(_particles.begin());
        // vector<Forest::coalinfo_t> coalinfo_vect;
        // p.recordAllForests(coalinfo_vect);
        // double log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        // output(format("log(coalescent likelihood) = %.9f\n") % log_coallike, 2);
        
        output("\nSpecies trees saved to file \"final-species-trees.tre\"\n", G::LogCateg::VERBOSE);
        string sfn = str(format("%s-final-species-trees") % prefix);
        saveAllSpeciesTrees(sfn, _particles, G::_treefile_compression);
        
        if (isJointMode()) {
            for (unsigned g = 0; g < G::_nloci; g++) {
                output(format("Gene trees for locus %d saved to file \"final-gene%d-trees.tre\"\n") % (g+1) % (g+1), G::LogCateg::VERBOSE);
                
                string fnprefix = str(format("%s-final-gene%d-trees") % prefix % (g+1));
                saveAllGeneTrees(g, fnprefix, _particles, G::_treefile_compression);
            }
        }
        
        if (!isConditionalMode()) {
            unsigned long min_partials_needed = 0;
            unsigned long k = G::_ntaxa - 1;
#if defined(UPGMA_WEIGHTS)
            // ntaxa = 4, k = 3
            //       step 1    step 2    step 3
            // ---------------------------------------------
            // SMC   \/ | |    \/ / |    \/ / /
            //                  \/  |     \/ /
            //                             \/
            //           1   +    1    +    1    = k
            // ---------------------------------------------
            // UPGMA   \/ /       \/
            //          \/
            //           2   +     1   +    0    = k*(k-1)/2
            // ---------------------------------------------
            min_partials_needed = 0.5*k*(k + 1)*G::_nparticles*G::_nloci;
#else
            min_partials_needed = k*G::_nparticles*G::_nloci;
#endif
        
            output(str(format("%20d %s\n") % G::_nspecies % "Species"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % G::_ntaxa % "Taxa"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % G::_nloci % "Loci"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % _nsteps % "Steps"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % G::_npartials_calculated % "Number of partials calculated"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % min_partials_needed % "Minimum partials needed"), G::LogCateg::INFO);
            output(str(format("%20.5f %s\n") % _log_marg_like % "Log marginal likelihood"), G::LogCateg::INFO);
        }
        
        map<string, tuple<unsigned, double, double, double, double> > m;
        bool ok = compareToReferenceTrees(_particles, m);
        
        if (ok) {
            double ref_height = 0.0;
            double mean_height = 0.0;
            double meanKF = 0.0;
            double minKF = G::_infinity;
            double maxKF = 0.0;
            double meanRF = 0.0;
            double minRF = G::_infinity;
            double maxRF = 0.0;
            unsigned total_count = 0;
            for (auto & kv : m) {
                tuple<unsigned, double, double, double, double> & t = kv.second;
                unsigned c = get<0>(t);
                ref_height  = get<1>(t);
                double test_height  = get<2>(t);
                double kf = get<3>(t);
                double rf = get<4>(t);
                total_count += c;
                mean_height += test_height*c;
                meanKF += kf*c;
                if (kf < minKF)
                    minKF = kf;
                if (kf > maxKF)
                    maxKF = kf;
                meanRF += rf*c;
                if (rf < minRF)
                    minRF = rf;
                if (rf > maxRF)
                    maxRF = rf;
            }
            meanKF /= total_count;
            meanRF /= total_count;
            mean_height /= total_count;

            string fn = str(format("%s-report.txt") % prefix);
            ofstream outf(fn);
            outf << str(format("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s")
                % "Loci"
                % "Steps"
                % "Distinct"
                % "RefH"
                % "MeanH"
                % "MeanKF"
                % "MinKF"
                % "MaxKF"
                % "SpreadKF"
                % "MeanRF"
                % "MinRF"
                % "MaxRF"
                % "SpreadRF"
                % "logML") << endl;
            outf << str(format("%12d %12d %12d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f")
                % G::_nloci
                % _nsteps
                % m.size()
                % ref_height
                % mean_height
                % meanKF
                % minKF
                % maxKF
                % (maxKF - minKF)
                % meanRF
                % minRF
                % maxRF
                % (maxRF - minRF)
                % _log_marg_like) << endl;
            outf.close();
            
            string longfn = str(format("%s-report-long-format.txt") % prefix);
            ofstream longoutf(longfn);
            longoutf << str(format("%20d %s\n") % G::_nspecies % "Species");
            longoutf << str(format("%20d %s\n") % G::_ntaxa % "Taxa");
            longoutf << str(format("%20d %s\n") % G::_nloci % "Loci");
            longoutf << str(format("%20d %s\n") % _nsteps % "Steps");
            longoutf << str(format("%20d %s\n") % m.size() % "Distinct topologies in species trees posterior");
            longoutf << str(format("%20d %s\n") % G::_npartials_calculated % "Number of partials calculated");
            longoutf << str(format("%20d %s\n") % ((G::_ntaxa - 1)*G::_nparticles*G::_nloci) % "Maximum partials needed");
            longoutf << str(format("%20.5f %s\n") % ref_height % "Reference species tree height");
            longoutf << str(format("%20.5f %s\n") % mean_height % "Posterior mean species tree height");
            longoutf << str(format("%20.5f %s\n") % meanKF % "Posterior mean Kuhner-Felsenstein species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % minKF % "Minimum Kuhner-Felsenstein species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % maxKF % "Maximum Kuhner-Felsenstein species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % (maxKF - minKF) % "Maximum - minimum Kuhner-Felsenstein species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % meanRF % "Posterior mean Robinson-Foulds species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % minRF % "Minimum Robinson-Foulds species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % maxRF % "Maximum Robinson-Foulds species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % (maxRF - minRF) % "Maximum - minimum Robinson-Foulds species tree distance to reference");
            longoutf << str(format("%20.5f %s\n") % _log_marg_like % "Log marginal likelihood");
            longoutf.close();
            
        }
    }
    
    inline double SMC::filterParticles(unsigned step, int locus) {
        // Sanity checks
        assert(_counts.size() == _nparticles);
        
        vector<double> log_weights(_nparticles, 0.0);
        for (unsigned i = 0; i < _particles.size(); i++) {
#if defined(USE_HEATING)
            log_weights[i] = _particles[i].getPrevLogWeight() + G::_heating_power*_particles[i].getLogWeight();
#else
            log_weights[i] = _particles[i].getLogWeight();
#endif
        }
                        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = G::calcLogSum(log_weights);
        vector<double> probs(_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
                
        // Compute component of the log marginal likelihood due to this step
        _log_marg_like += log_sum_weights - log(_nparticles);
        if (locus > -1 && step < G::_nloci) {
            _log_marg_like += _starting_log_likelihoods[locus];
        }
        
        // Compute effective sample size
        double ess = computeEffectiveSampleSize(probs);
        
        // Zero vector of counts storing number of darts hitting each particle
        _counts.assign(_nparticles, 0);
                
        // Store indexes of particles with non-zero counts in vector nonzeros
        vector<unsigned> zeros;
        zeros.reserve(_nparticles);
        vector<unsigned> nonzeros;
        nonzeros.reserve(_nparticles);

#if defined(SYSTEMATIC_FILTERING)
        assert(probs.size() == _nparticles);
        double cump = probs[0];
        double n = _nparticles;
        double delta = rng->uniform()/n;
        unsigned c = (unsigned)(floor(1.0 + n*(cump - delta)));
        if (c > 0)
            nonzeros.push_back(0);
        else
            zeros.push_back(0);
        _counts[0] = c;
        unsigned prev_cum_count = c;
        for (unsigned i = 1; i < _nparticles; ++i) {
            cump += probs[i];
            double cum_count = floor(1.0 + n*(cump - delta));
            if (cum_count > n)
                cum_count = n;
            unsigned c = (unsigned)cum_count - prev_cum_count;
            if (c > 0)
                nonzeros.push_back(i);
            else
                zeros.push_back(i);
            _counts[i] = c;
            prev_cum_count = cum_count;
        }
#else
        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());
        
        // Throw _nparticles darts
        for (unsigned i = 0; i < _nparticles; ++i) {
            double u = ::rng->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            _counts[which]++;
        }
        
        // //temporary!
        // string countstr = "";
        // for (unsigned i = 0; i < _counts.size(); i++)
        //     countstr += to_string(_counts[i]) + " ";
        // output(format("\ncounts: %s\n") % countstr, G::LogCateg::DEBUGGING);
                        
        classifyCounts(zeros, nonzeros, _counts);
#endif
        
        // Example of following code that replaces dead
        // particles with copies of surviving particles:
        //             0  1  2  3  4  5  6  7  8  9
        // _counts  = {0, 2, 0, 0, 0, 8, 0, 0, 0, 0}  size = 10
        // zeros    = {0, 2, 3, 4, 6, 7, 8, 9}        size =  8
        // nonzeros = {1, 5}                          size =  2
        //
        //  next_zero   next_nonzero   k   copy action taken
        //  --------------------------------------------------------------
        //      0             0        0   _particles[1] --> _particles[0]
        //  --------------------------------------------------------------
        //      1             1        0   _particles[5] --> _particles[2]
        //      2             1        1   _particles[5] --> _particles[3]
        //      3             1        2   _particles[5] --> _particles[4]
        //      4             1        3   _particles[5] --> _particles[6]
        //      5             1        4   _particles[5] --> _particles[7]
        //      6             1        5   _particles[5] --> _particles[8]
        //      7             1        6   _particles[5] --> _particles[9]
        //  --------------------------------------------------------------
        unsigned next_zero = 0;
        unsigned next_nonzero = 0;
        while (next_nonzero < nonzeros.size()) {
            double index_survivor = nonzeros[next_nonzero];
            unsigned ncopies = _counts[index_survivor] - 1;
            for (unsigned k = 0; k < ncopies; k++) {
                double index_nonsurvivor = zeros[next_zero++];
                
                // Replace non-survivor with copy of survivor
                _particles[index_nonsurvivor] = _particles[index_survivor];
            }
            
            ++next_nonzero;
        }
                                    
        return ess;
    }
    
    inline double SMC::filterParticlesWithinSubpops(unsigned step, int locus) {
        // Sanity checks
        assert(_counts.size() == _nparticles);
        assert(_subpop_id.size() == _nparticles);
        
        // Zero vector of counts storing number of darts hitting each particle
        _counts.assign(_nparticles, 0);

        // Determine the size of each subpopulation subset of particles
        unsigned nparticles_per_subpop = _nparticles/_nsubpops;

        // Create a vector next_index in which next_index[s] for subpopulation s
        // equals the next available position in the particles vector
        vector<unsigned> next_index(_nsubpops, 0);
        unsigned first = 0;
        for (unsigned i = 0; i < _nsubpops; i++) {
            next_index[i] = first;
            first += nparticles_per_subpop;
        }
        
        // Create a vector of pointers to particles. This vector stores particles
        // within a given subpopulation together.
        // Example: _nsubpops = 2, _nparticles = 10
        //            +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
        // subpop_id  |   0 |   1 |   1 |   0 |   0 |   1 |   0 |   1 |   0 |   1 |
        //            +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
        // _particles |  p0 |  p1 |  p2 |  p3 |  p4 |  p5 |  p6 |  p7 |  p8 |  p9 |
        //            +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
        // particles  | &p0 | &p3 | &p4 | &p6 | &p8 | &p1 | &p2 | &p5 | &p7 | &p9 |
        //            +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
        vector<Particle *> particles(_nparticles);
        for (unsigned i = 0; i < _particles.size(); i++) {
            unsigned s = _subpop_id[i];
            unsigned j = next_index[s];
            particles[j] = &_particles[i];
            next_index[s]++;
        }

        // Normalize log_weights (within subpopulations)
        // to create discrete probability distribution
        vector<double> ess(_nsubpops, 0.0);
        vector<double> log_weights(nparticles_per_subpop, 0.0);
        vector<double> probs(nparticles_per_subpop, 0.0);
        first = 0;
        for (unsigned subpop = 0; subpop < _nsubpops; subpop++) {
            unsigned last = first + nparticles_per_subpop;

            // Record log-weight for each particle
            unsigned j = 0;
            for (unsigned i = first; i < last; i++) {
                log_weights[j++] = particles[i]->getLogWeight();
            }
                        
            double log_sum_weights = G::calcLogSum(log_weights);
            transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
            
            // Compute component of the log marginal likelihood due to this step
            //TODO: probably not correct for multiple subpopulations
            _log_marg_like += log_sum_weights - log(nparticles_per_subpop);
            if (locus > -1 && step < G::_nloci) {
                _log_marg_like += _starting_log_likelihoods[locus];
            }

            // Compute effective sample size
            //TODO: probably not correct for multiple subpopulations
            ess[subpop] = computeEffectiveSampleSize(probs);
            
            // Compute cumulative probabilities
            partial_sum(probs.begin(), probs.end(), probs.begin());
            
#if defined(SYSTEMATIC_FILTERING)
            // This section needs to be written
            assert(false);
#else
            // Throw _nparticles darts
            for (unsigned i = 0; i < nparticles_per_subpop; ++i) {
                double u = ::rng->uniform();
                auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
                assert(it != probs.end());
                unsigned which = _nsubpops*subpop + (unsigned)distance(probs.begin(), it);
                _counts[which]++;
            }
#endif
                            
            first = last;
        }
                                        
        // Store indices of particles with non-zero counts in vector nonzeros
        vector<unsigned> zeros;
        vector<unsigned> nonzeros;
        classifyCounts(zeros, nonzeros, _counts);
        
        // Example of following code that replaces dead
        // particles with copies of surviving particles:
        //             0  1  2  3  4  5  6  7  8  9
        // _counts  = {0, 2, 0, 0, 0, 8, 0, 0, 0, 0}  size = 10
        // zeros    = {0, 2, 3, 4, 6, 7, 8, 9}        size =  8
        // nonzeros = {1, 5}                          size =  2
        //
        //  next_zero   next_nonzero   copy action taken
        //  ----------------------------------------------------------
        //      0             0        *particles[1] --> *particles[0]
        //  ----------------------------------------------------------
        //      1             1        *particles[5] --> *particles[2]
        //      2             1        *particles[5] --> *particles[3]
        //      3             1        *particles[5] --> *particles[4]
        //      4             1        *particles[5] --> *particles[6]
        //      5             1        *particles[5] --> *particles[7]
        //      6             1        *particles[5] --> *particles[8]
        //      7             1        *particles[5] --> *particles[9]
        //  ----------------------------------------------------------
        unsigned next_zero = 0;
        unsigned next_nonzero = 0;
        while (next_nonzero < nonzeros.size()) {
            double index_survivor = nonzeros[next_nonzero];
            unsigned ncopies = _counts[index_survivor] - 1;
            for (unsigned k = 0; k < ncopies; k++) {
                double index_nonsurvivor = zeros[next_zero++];
                
                // Replace non-survivor with copy of survivor
                *particles[index_nonsurvivor] = *particles[index_survivor];
            }
            
            ++next_nonzero;
        }
        
        // Randomly shuffle assignment of particles to subpopulations
        // First create tmp, a vector of tuples, with each tuple comprising
        // a uniform random variate and the index into _subpop_id
        vector< pair<float, unsigned> > tmp(_nparticles);
        for (unsigned i = 0; i < _nparticles; i++) {
            tmp[i] = make_pair(rng->uniform(), i);
        }
        
        // Sort tmp, which leaves the indexes in random order
        sort(tmp.begin(), tmp.end());
        
        // Make a copy of _subpop_id called tmpid
        vector<unsigned> tmpid(_subpop_id.begin(), _subpop_id.end());
        
        // Finally, recreate _subpop_id by copying subpop assignment from tmpid
        // in the order specified by tmp
        for (unsigned i = 0; i < _nparticles; i++) {
            unsigned index = tmp[i].second;
            unsigned subpop_assigned = tmpid[index];
            _subpop_id[i] = subpop_assigned;
        }
        
        double avg_ess = accumulate(ess.begin(), ess.end(), 0.0)/_nsubpops;
        return avg_ess;
    }
    
    inline unsigned SMC::countDistinctGeneTreeTopologies() {
        set<Forest::treeid_t> unique_topologies;
        for (auto & p : _particles) {
            for (unsigned g = 0; g < G::_nloci; g++) {
                GeneForest & gf = p.getGeneForest(g);
                Forest::treeid_t id;
                gf.storeSplits(id);
                unique_topologies.insert(id);
            }
        }
        return (unsigned)unique_topologies.size();
    }
        
    struct SpeciesTreeDetails {
        unsigned _count;
        map<G::species_t, double> _theta_map;
        double _log_coallike;
        
        SpeciesTreeDetails() : _count(0), _log_coallike(0.0) {}
    };

    inline void SMC::saveAllSpeciesTrees(string fnprefix, const vector<Particle> & particles, unsigned compression_level) {
        assert(compression_level >= 0 && compression_level <= 2);
        typedef tuple<unsigned, double, string, string, string> treeinfo_t;
        treeinfo_t treeinfo;
        vector<treeinfo_t> treeinfo_vect;
        if (compression_level == 2) {
            // Save unique species trees along with their frequency in <fnprefix>.tre
            // Save coalescent log likelihood and thetas in <fnprefix>.txt
            map<string,vector<SpeciesTreeDetails> > tree_info;
            for (const Particle & p : particles) {
                SpeciesTreeDetails info;
                
                // Get count for this particle
                info._count = 1;
                
                // Get newick tree description for this species tree
                string newick = p.getSpeciesForestConst().makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                
                if (isConditionalMode()) {
                    //temporary! while debugging only want to call calcLogCoalescentLikelihood while in 2nd level
                    // Calculate log coalescent likelihood for this species tree
                    vector<Forest::coalinfo_t> coalinfo_vect;
                    p.recordAllForests(coalinfo_vect);
#if defined(DEBUG_COALLIKE)
                    info._log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/false, /*verbose*/true);
#else
                    info._log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/false, /*verbose*/false);
#endif
                }
                
                // Record everything for this particle
                tree_info[newick].push_back(info);
            }
            
            ofstream tmpf(str(format("%s.txt") % fnprefix));
            unsigned i = 0;
            for (auto it = tree_info.begin(); it != tree_info.end(); ++it) {
                const string & newick = it->first;
                vector<SpeciesTreeDetails> & details_vect = it->second;

                // Calculate total count (i.e. freq)
                unsigned total_c = 0;
                for (auto v : details_vect) {
                    total_c += v._count;
                }
                
                tmpf << total_c << " <-- " << newick << endl;
                unsigned k = 0;
                for (auto v : details_vect) {
                    unsigned c = v._count;
                    tmpf << "  particle: " << k << endl;
                    tmpf << "    count: " << c << endl;
                    tmpf << "    log_coal_like: " << str(format("%.9f") % v._log_coallike) << endl;
                    for (auto spp_theta : v._theta_map) {
                        tmpf << "    theta for species " << spp_theta.first << ": " << spp_theta.second << endl;
                    }
                    tmpf << endl;
                    k++;
                }
                
                double pct = 100.0*total_c/_nparticles;
                string note = str(format("freq = %d") % total_c);
                string treename = str(format("'tree%d-freq%d'") % i % total_c);
                treeinfo = make_tuple(total_c, pct, note, treename, newick);
                treeinfo_vect.push_back(treeinfo);
                i++;
            }
            tmpf.close();
        }
        else {
            // Save all species trees (might involve saving the same newick string many times)
            // Some compression still performed (despite compress being false) because
            // a particle with count 100 will only be saved once and labeled
            // as having freq=100
            unsigned i = 0;
            for (const Particle & p : particles) {
                unsigned c = 1;
                string newick = p.getSpeciesForestConst().makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                if (compression_level == 1) {
                    double pct = 100.0*c/_nparticles;
                    string note = str(format("freq = %d") % c);
                    string treename = str(format("'tree%d-freq%d'") % i % c);
                    treeinfo = make_tuple(c, pct, note, treename, newick);
                    treeinfo_vect.push_back(treeinfo);
                    i++;
                }
                else {
                    double pct = 100.0/_nparticles;
                    string note = "freq = 1";
                    for (unsigned j = 0; j < c; j++) {
                        string treename = str(format("'tree%d-freq1'") % i);
                        treeinfo = make_tuple(c, pct, note, treename, newick);
                        treeinfo_vect.push_back(treeinfo);
                        i++;
                    }
                }
            }
        }
        string fn = str(format("%s.tre") % fnprefix);
        outputAnnotatedNexusTreefile(fn, treeinfo_vect);
    }
    
    struct GeneTreeDetails {
        unsigned _count;
        double _log_likelihood;
    };

    inline void SMC::saveAllGeneTrees(unsigned gene_index, string fnprefix, vector<Particle> & particles, unsigned compression_level) {
        assert(compression_level >= 0 && compression_level <= 2);
        typedef tuple<unsigned, double, string, string, string> treeinfo_t;
        treeinfo_t treeinfo;
        vector<treeinfo_t> treeinfo_vect;
        if (compression_level == 2) {
            // Save only unique newick strings
            map<string,vector<GeneTreeDetails> > tree_info;
            for (Particle & p : particles) {
                GeneTreeDetails info;
                
                // Get count for this particle
                info._count = 1;
                
                // Get newick tree description for this gene tree
                assert(gene_index < p.getGeneForestsConst().size());
                GeneForest & gf = p.getGeneForests()[gene_index];
                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                
                // Calculate log-likelihood for this gene tree
#if defined(UPGMA_WEIGHTS)
                info._log_likelihood = gf.calcLogLikelihood(p.getSpeciesForestConst());
#else
                info._log_likelihood = gf.calcLogLikelihood();
#endif
                
                // Record everything for this particle
                tree_info[newick].push_back(info);
            }
            
            ofstream tmpf(str(format("%s.txt") % fnprefix));
            unsigned i = 0;
            for (auto it = tree_info.begin(); it != tree_info.end(); ++it) {
                const string & newick = it->first;
                vector<GeneTreeDetails> & details_vect = it->second;

                // Calculate total count (i.e. freq)
                unsigned total_c = 0;
                for (auto v : details_vect) {
                    total_c += v._count;
                }
                
                tmpf << total_c << " <-- " << newick << endl;
                unsigned k = 0;
                for (auto v : details_vect) {
                    unsigned c = v._count;
                    tmpf << "  particle: " << k << endl;
                    tmpf << "    count: " << c << endl;
                    tmpf << "    log_like: " << str(format("%.9f") % v._log_likelihood) << endl;
                    tmpf << endl;
                    k++;
                }
                
                double pct = 100.0*total_c/_nparticles;
                string note = str(format("freq = %d") % total_c);
                string treename = str(format("'tree%d-freq%d'") % i % total_c);
                treeinfo = make_tuple(total_c, pct, note, treename, newick);
                treeinfo_vect.push_back(treeinfo);
                ++i;
            }
            tmpf.close();
        }
        else {
            unsigned i = 0;
            for (const Particle & p : particles) {
                unsigned c = 1;
                assert(gene_index < p.getGeneForestsConst().size());
                const GeneForest & gf = p.getGeneForestsConst()[gene_index];
                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                if (compression_level == 1) {
                    double pct = 100.0*c/_nparticles;
                    string note = str(format("freq = %d") % c);
                    string treename = str(format("'tree%d-freq%d'") % i % c);
                    treeinfo = make_tuple(c, pct, note, treename, newick);
                    treeinfo_vect.push_back(treeinfo);
                    i++;
                }
                else {
                    double pct = 100.0/_nparticles;
                    string note = "freq = 1";
                    for (unsigned j = 0; j < c; j++) {
                        string treename = str(format("'tree%d-freq1'") % i);
                        treeinfo = make_tuple(c, pct, note, treename, newick);
                        treeinfo_vect.push_back(treeinfo);
                        i++;
                    }
                }
            }
        }
        string fn = str(format("%s.tre") % fnprefix);
        outputAnnotatedNexusTreefile(fn, treeinfo_vect);
    }
            
    inline void SMC::outputAnnotatedNexusTreefile(string fn, const vector<tuple<unsigned, double, string, string, string> > & treeinfo) const {
        ofstream streef(fn);
        streef << "#NEXUS\n\n";
        streef << "begin trees;\n";
        unsigned t = 0;
        for (auto tinfo : treeinfo) {
            //unsigned  count = get<0>(tinfo);
            //double      pct = get<1>(tinfo);
            string     note = get<2>(tinfo);
            string treename = get<3>(tinfo);
            string   newick = get<4>(tinfo);
            streef << str(format("  tree %s = [%s] [&R] %s;\n") % treename % note % newick);
            ++t;
        }
        streef << "end;\n";
        streef.close();
    }

    inline bool SMC::compareToReferenceTrees(vector<Particle> particles, map<string, tuple<unsigned, double, double, double, double> > & m) {
        // Bail out if no reference tree was specified
        if (G::_species_tree_ref_file_name.empty())
            return false;
            
        // Read in the reference species tree
        G::_nexus_taxon_map.clear();
        vector<string> tree_names;
        vector<string> newicks;
        SpeciesForest::readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, G::_nexus_taxon_map, tree_names, newicks);

        // Build only the first species tree in the file
        SpeciesForest ref;
        ref.buildFromNewick(newicks[0]);
        double ref_height = ref.getHeight();
        
        for (Particle & p : particles) {
            unsigned c = 1;
            string newick = p.getSpeciesForest().makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
            if (m.count(newick) == 1) {
                // newick already in m, so just update count
                tuple<unsigned, double, double, double, double> & t = m[newick];
                unsigned cnew = c + get<0>(t);
                double ref_height = get<1>(t);
                double test_height = get<2>(t);
                double kf = get<3>(t);
                double rf = get<4>(t);
                m[newick] = make_tuple(cnew, ref_height, test_height, kf, rf);
            }
            else {
                // newick is distinct, so calculate KF,RF distances and root height
                SpeciesForest & test = p.getSpeciesForest();
                double test_height = test.getHeight();
                pair<double,double> kf_rf = Forest::calcTreeDistances(ref, test);
                m[newick] = make_tuple(c, ref_height, test_height, kf_rf.first, kf_rf.second);
            }
        }
        return true;
    }

    inline double SMC::computeEffectiveSampleSize(const vector<double> & probs) const {
        double ss = 0.0;
        for_each(probs.begin(), probs.end(), [&ss](double w){ss += w*w;});
        double ess = 1.0/ss;
        return ess;
    }

    inline void SMC::classifyCounts(vector<unsigned> & zeros, vector<unsigned> & nonzeros, const vector<unsigned> & counts) const {
        zeros.clear();
        nonzeros.clear();
        for (unsigned k = 0; k < counts.size(); k++) {
            if (counts[k] > 0) {
                nonzeros.push_back(k);
            }
            else {
                zeros.push_back(k);
            }
        }
    }

    double SMC::calcLogSpeciesTreePrior(vector<Forest::coalinfo_t> & coalinfo_vect, bool include_join_probs) const {
        double log_prior = 0.0;
        double tprev = 0.0;
        double t = 0.0;
        unsigned n = G::_nspecies;
        double lambda = G::_lambda;
        for (auto cinfo : coalinfo_vect) {
            unsigned gene_plus_1 = get<1>(cinfo);
            if (gene_plus_1 == 0) {
                // Species tree join
                t = get<0>(cinfo);
                double dt = t - tprev;
                
                // Prior for increment
                log_prior += log(lambda*n) - lambda*n*dt;
                
                // Prior for join
                if (include_join_probs)
                    log_prior += log(2.) - log(n) - log(n-1);
                
                tprev = t;
                n--;
            }
        }
        
        return log_prior;
    }

    void SMC::dumpParticles(SMC & ensemble, vector<unsigned> & kept) {
        // Dump _particles into ensemble, leaving _particles empty
        //move(_particles.begin(), _particles.end(), back_inserter(ensemble._particles));
        for (auto k : kept) {
            ensemble._particles.push_back(_particles[k]);
        }
        _particles.clear();
    }
    
}
