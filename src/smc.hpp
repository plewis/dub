#pragma once

extern proj::StopWatch stopwatch;
extern proj::Lot::SharedPtr rng;

namespace proj {

    class Particle;

    class SMC {
        public:
                                     SMC();
            virtual                  ~SMC();
                        
            enum mode_type_t {
                SPECIES_AND_GENE = 0,
                SPECIES_GIVEN_GENE = 1
            };

            void                     setMode(mode_type_t m)     {_mode = m;}
            bool                     isJointMode() const        {return _mode == SPECIES_AND_GENE;}
            bool                     isConditionalMode() const  {return _mode == SPECIES_GIVEN_GENE;}

            double                   calcPercentSpeciesTreeCompletion() const;

            void                     setNParticles(unsigned nparticles, unsigned nsubpops)  {_nparticles = nparticles; _nsubpops = nsubpops; }

#if defined(USING_MULTITHREADING)
            void                     runFirstLevelMultithread();
            void                     runSecondLevelMultithread();
            void                     advanceParticleRange(unsigned step, unsigned locus, bool rebuild_species_tree, unsigned first_particle, unsigned last_particle);
#else
            void                     runFirstLevelSerial();
            void                     runSecondLevelSerial();
#endif

            Lot::SharedPtr           getLot();
            void                     setRandomNumberSeed(unsigned rnseed);
    
            void                     setData(Data::SharedPtr d) {_data = d;}
            void                     initFromParticle(const Particle & p, Data::SharedPtr data);
            void                     init();
            void                     run();
            void                     keepSecondLevel(vector<unsigned> & kept);
            
            double                   getLogMarginalLikelihood() const;

#if defined(LAZY_COPYING)
            void                    buildNonzeroMap(unsigned locus,
                                        map<const void *, list<unsigned> > & nonzero_map,
                                        const vector<unsigned> & nonzeros);
#endif
            double                   calcLogSum(const vector<double> & log_values) const;
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
            void                     saveReport();
            void                     saveGeneTrees();
            void                     saveSpeciesTrees();
            void                     summarize();
            void                     dumpParticles(SMC & ensemble, vector<unsigned> & kept);
            void                     clear();
            unsigned                 countDistinctGeneTreeTopologies();
            void                     saveParamsForLoRaD(string prefix);
            
            typedef shared_ptr<SMC>  SharedPtr;
            
        private:
        
            unsigned                 _mode;
            Data::SharedPtr          _data;
            Lot::SharedPtr           _lot;
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
        
    inline SMC::SMC() {
        _lot.reset(new Lot);
        clear();
    }
    
    inline SMC::~SMC() {
    }
    
    inline void SMC::clear() {
        _mode = SPECIES_AND_GENE;
        _nsteps = 0;
        _log_marg_like = 0.0;
        _nparticles = 0;
        _nsubpops = 1;
        _counts.clear();
        _subpop_id.clear();
    }
    
    inline Lot::SharedPtr SMC::getLot() {
        return _lot;
    }
    
    inline void SMC::setRandomNumberSeed(unsigned rnseed) {
        _lot->setSeed(rnseed);
    }
    
    inline void SMC::initFromParticle(const Particle & parent_particle, Data::SharedPtr data) {
        assert(_nparticles > 0);
        
#if defined(DEBUGGING_INITFROMPARTICLE)
        // output first level gene trees and species tree
        output(format("\nSpecies forest:\n%s\n") % parent_particle.getSpeciesForest().makeNewick(9, true, false));
        auto gene_forests = p.getGeneForestsConst();
        for (auto & gf : gene_forests) {
            output(format("\nGene forest:\n%s\n") % gf.makeNewick(9, true, false));
        }
#endif
        _data = data;
        
        // This function should only be used when estimating species
        // trees from complete gene trees (obtained from parent_particle)
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
        
        // Make a copy of parent_particle
        Particle p = parent_particle;
        assert(p.isEnsembleCoalInfo());
        
        // Replace existing species tree with trivial forest
        //TODO: this is better done before we start second level
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
            
        // Initialize particle list with _nparticles Particles,
        // each of which is a copy of p
        assert(_particles.size() == 0);
        _particles.resize(_nparticles, p);
        
        // Determine total number of steps required to build the species tree
        assert(G::_nspecies > 1);
        _nsteps = G::_nspecies - 1;
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
#if defined(LAZY_COPYING)
        const vector<GeneForest::SharedPtr> & gene_forest_ptrs = template_particle.getGeneForestPtrsConst();
        for (unsigned g = 0; g < G::_nloci; g++) {
            _starting_log_likelihoods[g] = gene_forest_ptrs[g]->getLogLikelihood();
        }
#else
        const vector<GeneForest> & gene_forests = template_particle.getGeneForestsConst();
        for (unsigned g = 0; g < G::_nloci; g++) {
            _starting_log_likelihoods[g] = gene_forests[g].getLogLikelihood();
        }
#endif

        // Initialize particles vector
        assert(_particles.size() == 0);
        _particles.resize(_nparticles, template_particle);
                        
        // Determine total number of steps required to build all gene trees
        _nsteps = G::_nloci*(G::_ntaxa - 1);
    }
        
    inline double SMC::calcPercentSpeciesTreeCompletion() const {
        unsigned ncomplete = 0;
        for (unsigned i = 0; i < _nparticles; i++) {
            const SpeciesForest & sf = _particles[i].getSpeciesForestConst();
            if (sf.isComplete())
                ncomplete++;
        }
        double percent = 100.0*ncomplete/_nparticles;
        return percent;
    }
        
#if defined(USING_MULTITHREADING)
    inline void SMC::advanceParticleRange(unsigned step, unsigned locus, bool rebuild_species_tree, unsigned first_particle, unsigned last_particle) {
        for (unsigned i = first_particle; i < last_particle; i++) {
            // Advance each particle by one coalescent event in one locus
            _particles[i].proposeCoalescence(
                step,
                locus,
                G::_seed_bank[i],
                rebuild_species_tree
            );
        }
    }
#endif

#if defined(USING_MULTITHREADING)
    //**********************************************
    //********** runFirstLevelMultithread **********
    //**********************************************
    inline void SMC::runFirstLevelMultithread() {
        G::buildThreadSchedule(_nparticles, "particle");
                
        output(format("\nSMC (level 1) will require %d steps.\n") % _nsteps, G::LogCateg::INFO);
        
        G::showSettings();
        
        // Display header for progress table
        output(format("\n%12s %12s  %24s %12s %12s\n") % "Step" % "ESS" % "logml" % "secs" % "wait");
        
        _log_marg_like = 0.0;
        double cum_secs = 0.0;
        for (unsigned step = 0; step < _nsteps; ++step) {
            stopwatch.start();
            unsigned locus = 0;
            
#if defined(RANDOM_LOCUS_ORDERING)
#error random locus ordering not yet implemented except for sim
            //TODO: if implemented, make sure species tree rebuilt for first locus considered, not just locus == 0
#else
            locus = step % G::_nloci;
            assert(locus < G::_nloci);
#endif
            // Replace seeds in seed bank
            G::generateUpdateSeeds(_nparticles);

            if (G::_nthreads > 1) {
                vector<thread> threads;
                for (unsigned i = 0; i < G::_nthreads; i++) {
                    threads.push_back(thread(&SMC::advanceParticleRange,
                        this,
                        step,
                        locus,
                        locus == 0,
                        G::_thread_sched[i].first,
                        G::_thread_sched[i].second)
                    );
                }

                // The join function causes this loop to pause until
                // the ith thread finishes
                for (unsigned i = 0; i < threads.size(); i++) {
                    threads[i].join();
                }
            }
            else {
                advanceParticleRange(step, locus, locus == 0, 0, _nparticles);
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            if (_nsubpops == 1)
                ess = filterParticles(step, locus);
            else
                ess = filterParticlesWithinSubpops(step, locus);

            VALGRIND_PRINTF("~~> post 1st-level step %d at time %d\n", step, (unsigned)clock());
            VALGRIND_MONITOR_COMMAND(str(format("detailed_snapshot stepsnaps-%d.txt") % step).c_str());

            // Calculate waiting time until done
            double secs = stopwatch.stop();
            cum_secs += secs;
            double avg_per_step = cum_secs/(step + 1);
            unsigned steps_to_go = _nsteps - (step + 1);
            double wait = avg_per_step*steps_to_go;
            
            output(format("%12d %12.3f %24.6f %12.3f %12.3f\n") % (step+1) % ess % _log_marg_like % secs % wait, G::LogCateg::INFO);
        }
    }

    //***********************************************
    //********** runSecondLevelMultithread **********
    //***********************************************
    inline void SMC::runSecondLevelMultithread() {
        _log_marg_like = 0.0;
        for (unsigned step = 0; step < _nsteps; ++step) {
            unsigned locus = 0;
            
            for (unsigned i = 0; i < _nparticles; i++) {
                // Advance each particle by one speciation event
                _particles[i].proposeSpeciation(step, _lot);
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            //TODO: should second level use subpops?
            double ess = filterParticles(step, -1);
        }
    }
#else
    //*****************************************
    //********** runFirstLevelSerial **********
    //*****************************************
    inline void SMC::runFirstLevelSerial() {
        output(format("\nSMC (level 1) will require %d steps.\n") % _nsteps, G::LogCateg::INFO);
        
        G::showSettings();
        
        // Display header for progress table
        output(format("\n%12s %12s %12s %24s %12s %12s\n") % "Step" % "ESS" % "spp. tree %" % "logml" % "secs" % "wait");
        
        _log_marg_like = 0.0;
        double cum_secs = 0.0;
        for (unsigned step = 0; step < _nsteps; ++step) {
            stopwatch.start();
            unsigned locus = 0;
            
#if defined(RANDOM_LOCUS_ORDERING)
#error random locus ordering not yet implemented except for sim
            //TODO: if implemented, make sure species tree rebuilt for first locus considered, not just locus == 0
#else
            locus = step % G::_nloci;
            assert(locus < G::_nloci);
#endif
            // Replace seeds in seed bank
            G::generateUpdateSeeds(_nparticles);

#if defined(LAZY_COPYING) && defined(PLOT_INCREMENT_DISTRIBUTIONS)
            // Propose increment from prior for each particle
#error begin again here
            double log_mean = -1.0;
            double log_sd = -1.0;
            for (unsigned i = 0; i < _nparticles; i++) {
                // Advance each particle by one coalescent event in one locus
                _particles[i].proposeCoalescence(
                    step,
                    locus,
                    G::_seed_bank[i],
                    log_mean,
                    log_sd,
                    /*rebuild_species_tree*/ locus == 0
                );
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            if (_nsubpops == 1)
                ess = filterParticles(step, locus, log_mean, log_sd);
            else {
                throw XProj("subpops not yet implemented for PLOT_INCREMENT_DISTRIBUTIONS");
                ess = filterParticlesWithinSubpops(step, locus);
            }

            // Propose increment from empirical distribution for each particle
            
            for (unsigned i = 0; i < _nparticles; i++) {
                // Advance each particle by one coalescent event in one locus
                _particles[i].proposeCoalescence(
                    step,
                    locus,
                    G::_seed_bank[i],
                    log_mean,
                    log_sd,
                    locus == 0 /*rebuild_species_tree*/
                );
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            if (_nsubpops == 1)
                ess = filterParticles(step, locus, log_mean, log_sd);
            else
                ess = filterParticlesWithinSubpops(step, locus);
#else
            for (unsigned i = 0; i < _nparticles; i++) {
                // Advance each particle by one coalescent event in one locus
                _particles[i].proposeCoalescence(
                    step,
                    locus,
                    G::_seed_bank[i],
                    locus == 0 /*rebuild_species_tree*/
                );
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            if (_nsubpops == 1)
                ess = filterParticles(step, locus);
            else
                ess = filterParticlesWithinSubpops(step, locus);
#endif

            VALGRIND_PRINTF("~~> post 1st-level step %d at time %d\n", step, (unsigned)clock());
            VALGRIND_MONITOR_COMMAND(str(format("detailed_snapshot stepsnaps-%d.txt") % step).c_str());

            // Calculate waiting time until done
            double secs = stopwatch.stop();
            cum_secs += secs;
            double avg_per_step = cum_secs/(step + 1);
            unsigned steps_to_go = _nsteps - (step + 1);
            double wait = avg_per_step*steps_to_go;
            
            double spp_tree_pct = calcPercentSpeciesTreeCompletion();
            output(format("%12d %12.3f %12.3f %24.6f %12.3f %12.3f\n") % (step+1) % ess % spp_tree_pct % _log_marg_like % secs % wait, G::LogCateg::INFO);
        }
    }

    //******************************************
    //********** runSecondLevelSerial **********
    //******************************************
    inline void SMC::runSecondLevelSerial() {
        output(format("\nSMC (level 2) will require %d steps.\n") % _nsteps, G::LogCateg::SECONDLEVEL);
        output(format("\n%12s %12s %24s %12s %12s\n") % "Step" % "ESS" % "logml" % "secs" % "wait", G::LogCateg::SECONDLEVEL);
        
        _log_marg_like = 0.0;
        double cum_secs = 0.0;
        for (unsigned step = 0; step < _nsteps; ++step) {
            stopwatch.start();
            
            for (unsigned i = 0; i < _nparticles; i++) {
                // Advance each particle by one speciation event
                _particles[i].proposeSpeciation(step, _lot);
            }

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            //TODO: should second level use subpops?
            ess = filterParticles(step, -1);

            // Calculate waiting time until done
            double secs = stopwatch.stop();
            cum_secs += secs;
            double avg_per_step = cum_secs/(step + 1);
            unsigned steps_to_go = _nsteps - (step + 1);
            double wait = avg_per_step*steps_to_go;
            
            output(format("%12d %12.3f %24.6f %12.3f %12.3f\n") % (step+1) % ess % _log_marg_like % secs % wait, G::LogCateg::SECONDLEVEL);
        }
    }
#endif

    inline void SMC::run() {
        assert(_nsteps > 0);
        bool first_level = isJointMode();
        bool second_level = isConditionalMode();
        assert(first_level || second_level);
        
#if defined(USING_MULTITHREADING)
        if (first_level)
            runFirstLevelMultithread();
        else
            runSecondLevelMultithread();
#else
        if (first_level)
            runFirstLevelSerial();
        else
            runSecondLevelSerial();
#endif
    }
    
    inline void SMC::saveGeneTrees() {
        assert(_particles.size() > 0);
        
        if (G::_save_gene_trees && isJointMode()) {
            string prefix = "1st";
            if (isConditionalMode())
                prefix = "2nd";

            // Save all gene trees for each locus in a separate *.tre file
            // and, if G::_treefile_compression == 2, in a long-format *.txt file
            // that reports the frequency and log likelihood of each
            for (unsigned g = 0; g < G::_nloci; g++) {
                output(format("Gene trees for locus %d saved to file \"%s-final-gene%d-trees.tre\"\n") % prefix % (g+1) % (g+1), G::LogCateg::VERBOSE);
                
                string fnprefix = str(format("%s-final-gene%d-trees") % prefix % (g+1));
                saveAllGeneTrees(g, fnprefix, _particles, G::_treefile_compression);
            }
        }
    }
    
    inline void SMC::saveSpeciesTrees() {
        assert(_particles.size() > 0);
        
        if (G::_save_species_trees) {
            string prefix = "1st";
            if (isConditionalMode())
                prefix = "2nd";

            // Save all species trees as both *.tre file and,
            // if G::_treefile_compression == 2, as a long-format *.txt file
            // that reports the frequency and log coalescent likelihood of each
            output("\nSpecies trees saved to file \"final-species-trees.tre\"\n", G::LogCateg::VERBOSE);
            string sfn = str(format("%s-final-species-trees") % prefix);
            saveAllSpeciesTrees(sfn, _particles, G::_treefile_compression);
        }
    }
    
    inline void SMC::saveReport() {
        string prefix = "1st";
        if (isConditionalMode())
            prefix = "2nd";

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
    
    inline void SMC::summarize() {
        assert(_particles.size() > 0);

        if (!isConditionalMode()) {
            unsigned long min_partials_needed = 0;
            unsigned long k = G::_ntaxa - 1;
            min_partials_needed = k*G::_nparticles*G::_nloci;
        
            output(str(format("%20d %s\n") % G::_nspecies % "Species"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % G::_ntaxa % "Taxa"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % G::_nloci % "Loci"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % _nsteps % "Steps"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % G::_npartials_calculated % "Number of partials calculated"), G::LogCateg::INFO);
            output(str(format("%20d %s\n") % min_partials_needed % "Minimum partials needed"), G::LogCateg::INFO);
            output(str(format("%20.5f %s\n") % _log_marg_like % "Log marginal likelihood"), G::LogCateg::INFO);
        }
    }

#if defined(LAZY_COPYING)
    inline void SMC::buildNonzeroMap(unsigned locus,
        map<const void *, list<unsigned> > & nonzero_map,
        const vector<unsigned> & nonzeros) {
        for (auto i : nonzeros) {
            void * ptr = _particles[i].getGeneForestPtr(locus).get();
            if (nonzero_map.count(ptr) > 0) {
                nonzero_map[ptr].push_back(i);
            }
            else {
                nonzero_map[ptr] = {i};
            }
        }
    }
#endif

    inline double SMC::calcLogSum(const vector<double> & log_values) const {
        double max_logv = *max_element(log_values.begin(), log_values.end());
        
        double factored_sum = 0.0;
        for (auto & logv : log_values) {
            factored_sum += exp(logv - max_logv);
        }
        double log_sum_values = max_logv + log(factored_sum);
        return log_sum_values;
    }
    
    inline double SMC::filterParticles(unsigned step, int locus) {
        bool first_level = (locus > -1);
        Lot::SharedPtr lot = (first_level ? ::rng : _lot);
        
        // Sanity checks
        assert(_counts.size() == _nparticles);
        
        // Build vector log_weights
        vector<double> log_weights(_nparticles, 0.0);
#if defined(LAZY_COPYING) && defined(PLOT_INCREMENT_DISTRIBUTIONS)
        vector<double> incr(_nparticles, 0.0);
#endif
        for (unsigned i = 0; i < _particles.size(); i++) {
#if defined(USE_HEATING)
            log_weights[i] = _particles[i].getPrevLogWeight() + G::_heating_power*_particles[i].getLogWeight();
#else
            log_weights[i] = _particles[i].getLogWeight();
#endif
#if defined(LAZY_COPYING) && defined(PLOT_INCREMENT_DISTRIBUTIONS)
            incr[i] = _particles[i].getProposedIncrement(locus);
#endif
        }
                        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = calcLogSum(log_weights);
        vector<double> probs(_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
                
        // Compute component of the log marginal likelihood due to this step
        _log_marg_like += log_sum_weights - log(_nparticles);
        if (first_level && step < G::_nloci) {
            _log_marg_like += _starting_log_likelihoods[locus];
        }
        
        // Compute effective sample size
        double ess = computeEffectiveSampleSize(probs);
        
#if defined(LAZY_COPYING) && defined(PLOT_INCREMENT_DISTRIBUTIONS)
        // Compute weighted mean increment length
        double mean_incr = 0.0;
        vector<pair<double,double> > tmp;
        for (unsigned i = 0; i < _particles.size(); i++) {
            mean_incr += probs[i]*incr[i];
            tmp.push_back(make_pair(incr[i], probs[i]));
        }

        ofstream tmpf(str(format("incr-%d-%d.R") % step % locus));
        
        // save increments
        tmpf << "incr <- c(";
        tmpf << tmp[0].first;
        for (unsigned i = 1; i < tmp.size(); i++) {
            tmpf << "," << tmp[i].first;
        }
        tmpf << ")\n";
        
        // save probs
        tmpf << "prob <- c(";
        tmpf << tmp[0].second;
        for (unsigned i = 1; i < tmp.size(); i++) {
            tmpf << "," << tmp[i].second;
        }
        tmpf << ")\n";
        
        // plot kernel density
        tmpf << "plot(density(incr, from=0.0))\n";
        
        // plot points with prob > 0.0001 as vertical dotted lines
        double log_sum_incr = 0.0;
        double log_sumsq_incr = 0.0;
        double log_n_incr = 0.0;
        for (unsigned i = 0; i < tmp.size(); i++) {
            if (tmp[i].second > 0.0001) {
                double log_incr = log(tmp[i].first);
                log_n_incr += 1.0;
                log_sum_incr += log_incr;
                log_sumsq_incr += pow(log_incr, 2.0);
                tmpf << "abline(v=" << tmp[i].first << ", lwd=2, lty=\"dotted\", col=\"navy\")\n";
            }
        }
        
        // Calculate parameters of empirical lognormal distribution
        if (log_n_incr > 1) {
            double log_mean_incr = log_sum_incr/log_n_incr;
            double log_var_incr = (log_sumsq_incr - log_n_incr*pow(log_mean_incr,2.0))/(log_n_incr - 1.0);
            double log_sd_incr = sqrt(log_var_incr);
            tmpf << "x <- seq(0.0, max(incr), 0.01)\n";
            tmpf << "curve(dlnorm(x, meanlog=" << log_mean_incr << ", sdlog=" << log_sd_incr << "), lwd=2, lty=\"solid\", col=\"navy\", add=T)\n";
        }
        
        tmpf.close();

        output(format("~~> step %d, locus %d, mean incr = %.5f\n") % step % locus % mean_incr, G::LogCateg::INFO);
#endif

        // Zero vector of counts storing number of darts hitting each particle
        _counts.assign(_nparticles, 0);
                
        // Store indices of particles with zero counts in vector zeros
        vector<unsigned> zeros;
        zeros.reserve(_nparticles);

        // Store indices of particles with non-zero counts in vector nonzeros
        vector<unsigned> nonzeros;
        nonzeros.reserve(_nparticles);

#if defined(SYSTEMATIC_FILTERING)
        assert(probs.size() == _nparticles);
        double cump = probs[0];
        double n = _nparticles;
        double delta = lot->uniform()/n;
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
            double u = lot->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            _counts[which]++;
        }
                                
        classifyCounts(zeros, nonzeros, _counts);
#endif

#if defined(LAZY_COPYING)
        // Create map (nonzero_map) in which the key for an element
        // is the memory address of a gene forest and
        // the value is a vector of indices of non-zero counts.
        // This map is used to determine which of the nonzeros
        // that need to be copied (last nonzero count for any
        // memory address does not need to be copied and can be
        // modified in place).
        map<const void *, list<unsigned> > nonzero_map;
        if (first_level) {
            buildNonzeroMap(locus, nonzero_map, nonzeros);
        }
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
#if defined(LAZY_COPYING)
            if (first_level) {
                _particles[index_survivor].finalizeLatestJoin(locus, index_survivor, nonzero_map);
            }
#endif
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
        // If used in second level, need to replace rng when multithreading
        assert(isJointMode());
        
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
                        
            double log_sum_weights = calcLogSum(log_weights);
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
                Forest::treeid_t id;
#if defined(LAZY_COPYING)
                GeneForest::SharedPtr gfp = p.getGeneForestPtr(g);
                gfp->storeSplits(id);
#else
                GeneForest & gf = p.getGeneForest(g);
                gf.storeSplits(id);
#endif
                unique_topologies.insert(id);
            }
        }
        return (unsigned)unique_topologies.size();
    }
        
    inline void SMC::saveParamsForLoRaD(string prefix) {
        output("Saving parameters for LoRaD...\n", G::LogCateg::INFO);
        output("  Parameter file is \"params-lorad.txt\"\n", G::LogCateg::INFO);
        output("  Run LoRaD using \"rscript lorad.R\"\n", G::LogCateg::INFO);
        unsigned num_distinct = countDistinctGeneTreeTopologies();
        if (num_distinct == 1) {
            output("  Topology identical for all gene trees in all particles\n", G::LogCateg::INFO);
        }
        else {
            output(format("  *** warning ***: %d distinct gene tree topologies (LoRaD assumes 1)\n") % num_distinct, G::LogCateg::INFO);
        }
        
        // Open the file into which parameters will be saved
        // in a format useful to LoRaD
        ofstream loradf(str(format("%s-params-lorad.txt") % prefix));

        // Save column names on first line
        vector<string> names;
        names.push_back("index");
        //names.push_back("particle");
        names.push_back("logFelsL");
        names.push_back("logCoalL");
        names.push_back("logPrior");
        G::getAllParamNames(names);
        loradf << boost::join(names, "\t") << endl;
        
        unsigned particle_object_index = 0;
        unsigned particle_index = 0;
        
        // //temporary!
        //unsigned tmp_nparticle_objects = (unsigned)_particle_list.size();
        
        for (auto & p : _particles) {
            // Stores tuple (height, gene + 1, left species, right species)
            // for each join in either species tree or any gene tree.
            // To get 0-offset gene index, subtract
            // one from second element of tuple (0 means tuple represents a
            // species tree join).
            vector<Forest::coalinfo_t> coalinfo_vect;
            
            // Add species tree joins to coalinfo_vect
            SpeciesForest & sf = p.getSpeciesForest();
            sf.heightsInternalsPreorders();
            sf.buildCoalInfoVect();
            sf.saveCoalInfo(coalinfo_vect);

            // Add gene tree joins to coalinfo_vect
            double log_likelihood = 0.0;
            for (unsigned g = 0; g < G::_nloci; g++) {
                GeneForest::SharedPtr gfp = p.getGeneForestPtr(g);
                gfp->heightsInternalsPreorders();
                gfp->buildCoalInfoVect();
                //gfp->forgetSpeciesTreeAbove(0.0);
                gfp->saveCoalInfo(coalinfo_vect);
                log_likelihood += gfp->calcLogLikelihood();
            }
            log_likelihood *= G::_phi;
            
            // Sort coalinfo_vect from smallest to largest height
            sort(coalinfo_vect.begin(), coalinfo_vect.end());
            
            //Forest::debugShowCoalInfo("SMC::saveParamsForLoRaD", coalinfo_vect);
            
            // Calculate log coalescent likelihood (joint prior for gene trees)
            double log_gene_tree_prior = p.calcLogCoalescentLikelihood(coalinfo_vect,
                /*integrate_out_thetas*/true, /*verbose*/false);
                        
            // Calculate log species tree prior
            double log_species_tree_prior = calcLogSpeciesTreePrior(coalinfo_vect, /*include_join_probs*/false);
                        
            // This vector will hold all increments from the species tree
            // and all gene trees (in order)
            vector<string> params;
            
            // Calculate increments in species tree
            double hprev = 0.0;
            for (auto cinfo : coalinfo_vect) {
                unsigned gene_plus_1 = get<1>(cinfo);
                if (gene_plus_1 == 0) {
                    double h = get<0>(cinfo);
                    params.push_back(str(format("%.9f") % (h - hprev)));
                    hprev = h;
                }
            }
            
            for (unsigned g = 0; g < G::_nloci; g++) {
                // Map to store previous height for each species
                //map<G::species_t, double> height_map;
                hprev = 0.0;
                for (auto cinfo : coalinfo_vect) {
                    unsigned gene_plus_1 = get<1>(cinfo);
                    double               h      = get<0>(cinfo);
                    //G::species_t sleft  = get<2>(cinfo);
                    //G::species_t sright = get<3>(cinfo);
                    if (gene_plus_1 > 0) {
                        unsigned gene_index = gene_plus_1 - 1;
                        if (gene_index == g) {
                            // coalescence event in gene g
                            //assert(sleft == sright);
                            params.push_back(str(format("%.9f") % (h - hprev)));
                            hprev = h;
                        }
                    }
                }
            }

            // Second element of pair is number of copies of particle
            //unsigned n = p.getCount();
            //for (unsigned i = 0; i < n; i++) {
                particle_index++;
                //loradf << str(format("%d\t%d\t%.9f\t%.9f\t%.9f\t%s")
                //    % particle_index
                //    % particle_object_index
                //    % log_likelihood
                //    % log_gene_tree_prior
                //    % log_species_tree_prior
                //    % boost::join(params, "\t")) << endl;
                loradf << str(format("%d\t%.9f\t%.9f\t%.9f\t%s")
                    % particle_index
                    % log_likelihood
                    % log_gene_tree_prior
                    % log_species_tree_prior
                    % boost::join(params, "\t")) << endl;
            //}
            //particle_object_index++;
        }
        loradf.close();
        
        ofstream loradRf("lorad.R");
        loradRf << "library(lorad)\n";
        loradRf << "colspec <- c(";
        loradRf << "\"index\" = \"iteration\",";
        //loradRf << "\"particle\" = \"ignore\",";
        loradRf << " \"logFelsL\" = \"posterior\",";
        loradRf << " \"logCoalL\" = \"posterior\",";
        loradRf << " \"logPrior\" = \"posterior\",";
        for (unsigned i = 4; i < names.size(); i++) {
            loradRf << " \"" << names[i] << "\" = \"positive\"";
            if (i < names.size() - 1)
                loradRf << ",";
        }
        loradRf << ")\n";
        loradRf << "params <- read.table(\"params-lorad.txt\", header=TRUE)\n";
        loradRf << "results <- lorad_estimate(params, colspec, 0.5, \"random\", 0.1)\n";
        loradRf << "lorad_summary(results)\n";
        
        loradRf.close();
        
        // ofstream conf("loradml.conf");
        // conf << "paramfile = params-lorad.txt\n";
        // conf << "plotfprefix = lnlplot\n";
        // conf << "trainingfrac = 0.5\n";
        // conf << "coverage = 0.1\n";
        // conf << "colspec = iteration index\n";
        // //conf << "colspec = ignore particle\n";
        // conf << "colspec = posterior logFelsL\n";
        // conf << "colspec = posterior logCoalL\n";
        // conf << "colspec = posterior logPrior\n";
        // for (unsigned i = 4; i < names.size(); i++) {
        //     conf << "colspec = positive " << names[i] << "\n";
        // }
        // conf.close();
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
#if defined(LAZY_COPYING)
                assert(gene_index < p.getGeneForestPtrsConst().size());
                GeneForest::SharedPtr gfp = p.getGeneForestPtrs()[gene_index];
                string newick = gfp->makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);

                // Calculate log-likelihood for this gene tree
                info._log_likelihood = gfp->calcLogLikelihood();
#else
                assert(gene_index < p.getGeneForestsConst().size());
                GeneForest & gf = p.getGeneForests()[gene_index];
                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                
                // Calculate log-likelihood for this gene tree
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
#if defined(LAZY_COPYING)
                assert(gene_index < p.getGeneForestPtrsConst().size());
                const GeneForest::SharedPtr gfp = p.getGeneForestPtrsConst()[gene_index];
                string newick = gfp->makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
#else
                assert(gene_index < p.getGeneForestsConst().size());
                const GeneForest & gf = p.getGeneForestsConst()[gene_index];
                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
#endif
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

    inline double SMC::calcLogSpeciesTreePrior(vector<Forest::coalinfo_t> & coalinfo_vect, bool include_join_probs) const {
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

    inline void SMC::dumpParticles(SMC & ensemble, vector<unsigned> & kept) {
        // Dump _particles into ensemble, leaving _particles empty
        assert(kept.size() > 0);
        for (auto k : kept) {
            ensemble._particles.push_back(_particles[k]);
        }
        assert(ensemble._particles.size() > 0);
        _particles.clear();
    }

    inline void SMC::keepSecondLevel(vector<unsigned> & kept) {
        assert(isConditionalMode());
        assert(_particles.size() == G::_nparticles2);
        assert(_lot);
        kept.reserve(_particles.size());
        for (unsigned k = 0; k < G::_nkept2; k++) {
            unsigned i = _lot->randint(0, G::_nparticles2 - 1);
            kept.push_back(i);
        }
    }
    
    inline double SMC::getLogMarginalLikelihood() const {
        return _log_marg_like;
    }
}
