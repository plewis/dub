#pragma once

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
extern proj::StopWatch stopwatch;
//POLWAS extern proj::Lot rng;
extern proj::Lot::SharedPtr rng;

namespace proj {

    inline void SMC::clear() {
        _mode = SPECIES_AND_GENE;
        _nparticles = 0;
        _nsteps = 0;
        _log_marg_like = 0.0;
        _counts.clear();
        _update_seeds.clear();
    }
    
    inline void SMC::initFromParticle(Particle & p) {
        assert(_nparticles > 0);
        _data = p.getData();
        
        // This function should only be used when estimating species
        // trees from complete gene trees (obtained from p)
        assert(_mode == SPECIES_GIVEN_GENE);
        
        // The _counts vector stores the number of darts that hit
        // each of the _nparticles in multinomial sampling.
        _counts.clear();
        _counts.resize(_nparticles);
        
        // Create vector of random number seeds to be reused each step.
        _update_seeds.clear();
        _update_seeds.resize(_nparticles);

        // Create log_weights vector that stores log weight of each particle
        _log_weights.resize(_nparticles);
        
        // Initialize first particle with a trivial species forest
        // but actual gene trees
        Particle template_particle;
        template_particle.setSMC(this);
        //template_particle.setCount(1);
        template_particle.setCount(_nparticles);
        template_particle.setData(_data);
        template_particle.resetSpeciesForest();
        template_particle.setThetas();
        template_particle.setGeneTrees(p.getGeneForests());

        // Compute starting log coalescent likelihood
        vector<Forest::coalinfo_t> coalinfo_vect;
        template_particle.recordAllForests(coalinfo_vect);
        double log_coallike = template_particle.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        template_particle.setPrevLogCoalLike(log_coallike);
            
        // Initialize particle list with _nparticles Particles, each of which
        // is a copy of template_particle
        _particle_list.clear();
        //_particle_list.resize(_nparticles, template_particle);
        _particle_list.resize(1, template_particle);
        
        // Determine total number of steps required to build the species tree
        assert(G::_nspecies > 1);
        _nsteps = G::_nspecies - 1;

        output(format("\nSMC (level 2) will require %d steps.\n") % _nsteps, 3);

#if defined(LOG_MEMORY)
        output(format("\n%12s %12s %12s %12s %24s %12s %12s %12s %12s\n") % "Step" % "ESS" % "minlogwt" % "maxlogwt" % "logml" % "in-use" % "stored" % "secs" % "wait", 3);
#else
        output(format("\n%12s %12s %12s %12s %24s %12s %12s\n") % "Step" % "ESS" % "minlogwt" % "maxlogwt" % "logml" % "secs" % "wait", 3);
#endif
    
    }
    
    inline void SMC::init() {
        assert(_nparticles > 0);
        assert(_data);
        
        // Stores number of darts that hit each of the _nparticles in
        // multinomial sampling.
        _counts.clear();
        _counts.resize(_nparticles);
        
        // Create vector of random number seeds to be reused each step.
        _update_seeds.clear();
        _update_seeds.resize(_nparticles);

        // Create log_weights vector that stores log weight of each particle
        _log_weights.resize(_nparticles);
        
        // Initialize template particle with trivial forests
        Particle template_particle;
        template_particle.setSMC(this);
        template_particle.setCount(_nparticles);
        template_particle.setData(_data);
        template_particle.resetSpeciesForest();
        template_particle.resetGeneForests(/*compute_partials*/true);

        // Compute initial log likelihood
        _starting_log_likelihood = template_particle.calcLogLikelihood(/*initializing*/true);

        // Initialize particle map with one particle having count _nparticles
        _particle_list.clear();
        _particle_list.push_back(template_particle);
                        
        // Determine total number of steps required to build all
        // gene trees and the species tree
        _nsteps = G::_ngenes*(G::_ntaxa - 1);
        unsigned verbosity = (_mode == SPECIES_AND_GENE ? 2 : 3);
        if (_mode == SPECIES_AND_GENE) {
            output(format("\nSMC (level 1) will require %d steps.\n") % _nsteps, verbosity);
        }
        
        G::showSettings();
        
#if defined(LOG_MEMORY)
        output(format("\n%12s %12s %12s %12s %24s %12s %12s %12s %12s\n") % "Step" % "ESS" % "minlogwt" % "maxlogwt" % "logml" % "in-use" % "stored" % "secs" % "wait", verbosity);
#else
        output(format("\n%12s %12s %12s %12s %24s %12s %12s\n") % "Step" % "ESS" % "minlogwt" % "maxlogwt" % "logml" % "secs" % "wait", verbosity);
#endif
    
    }

    inline void SMC::run() {
        assert(_nsteps > 0);
        
        vector< pair<double, unsigned> > locus_ordering(G::_ngenes);
        
        _log_marg_like = 0.0;
        double cum_secs = 0.0;
        for (unsigned step = 0; step < _nsteps; ++step) {
            //MARK: main step loop
            stopwatch.start();

            // Behavior of balanceThreads is governed by ParallelPolicy
            // (e.g. ParallelPolicyNone<T> does nothing)
            balanceThreads();

            // Initialize vector to store log weights
            _log_weights.assign(_nparticles, 0.0);
            
            // Initialize vectors to store proposed states
            //_proposed_gene.assign(_nparticles, 0);
            //_proposed_spp.assign(_nparticles, 0);
            //_proposed_first.assign(_nparticles, 0);
            //_proposed_second.assign(_nparticles, 0);
            //_proposed_species_tree_lineages.assign(_nparticles, 0);
            
            // Assign _nparticles random number seeds to use for generating the next step
            G::generateUpdateSeeds(_update_seeds);
            
            // Advance each particle by one coalescent event (if SPECIES_AND_GENE mode)
            // or one speciation event (if SPECIES_GIVEN_GENE mode)
            // Implementation provided by parallelization policy
            unsigned step_modulus = step % G::_ngenes;
            assert(step_modulus < G::_ngenes);
            if (step_modulus == 0) {
                // Time to choose a new random locus ordering
                for (unsigned i = 0; i < G::_ngenes; i++) {
                    double u = rng->uniform();
                    locus_ordering[i] = make_pair(u, i);
                }
                sort(locus_ordering.begin(), locus_ordering.end());
                particleLoop(step, locus_ordering[step_modulus].second, _update_seeds, /*refresh_species_tree*/true);
            }
            else {
                particleLoop(step, locus_ordering[step_modulus].second, _update_seeds, /*refresh_species_tree*/false);
            }

            // Determine minimum and maximum log weight
            double minlogw = *min_element(_log_weights.begin(), _log_weights.end());
            double maxlogw = *max_element(_log_weights.begin(), _log_weights.end());

            // Filter particles using normalized weights and multinomial sampling
            // The _counts vector will be filled with the number of times each particle was chosen.
            double ess = _nparticles;
            if (step_modulus == 0) {
                ess = filterParticles(step, locus_ordering[step_modulus].second, _particle_list, _log_weights, _counts, _update_seeds, /*refresh_species_tree*/true);
            }
            else {
                ess = filterParticles(step, locus_ordering[step_modulus].second, _particle_list, _log_weights, _counts, _update_seeds, /*refresh_species_tree*/false);
            }
            
#if defined(EST_THETA)
            debugSaveThetas(step);
#endif
            
            //debugCheckCountMap(count_map, _nparticles);
            double secs = stopwatch.stop();
            
            // Calculate waiting time until done
            cum_secs += secs;
            double avg_per_step = cum_secs/(step + 1);
            unsigned steps_to_go = _nsteps - (step + 1);
            double wait = avg_per_step*steps_to_go;
                    
        unsigned verbosity = (_mode == SPECIES_AND_GENE ? 2 : 3);
#if defined(LOG_MEMORY)
            //npartials = ps.getNumberConstructed();
            npartials_inuse = ps.getInUse();
            npartials_stored = ps.getStored();
            
            output(format("%12d %12.3f %12.3f %12.3f %24.6f %12d %12d %12.3f %12.3f\n") % (step+1) % ess % minlogw % maxlogw % _log_marg_like % npartials_inuse % npartials_stored % secs % wait, verbosity);
#else
            output(format("%12d %12.3f %12.3f %12.3f %24.6f %12.3f %12.3f\n") % (step+1) % ess % minlogw % maxlogw % _log_marg_like % secs % wait, verbosity);
#endif
        }
        
        // //temporary!
        // Particle & first = *_particle_list.begin();
        // vector<Forest::coalinfo_t> sppinfo_vect;
        // first.getSpeciesForest().saveCoalInfo(sppinfo_vect, /*cap*/true);
        // sort(sppinfo_vect.begin(), sppinfo_vect.end());
        // Forest::debugShowCoalInfo("===== sppinfo for first particle =====", sppinfo_vect);
    }
    
#if defined(EST_THETA)
    inline void SMC::debugSaveThetas(unsigned step) const {
        if (step == 0) {
            //    particle          0        1        2   ...  _nsteps
            //      0         0.01234  0.00213  0.02009   ...  0.02137
            //      1         0.01234  0.00019  0.02009   ...  0.00479
            //      .               .        .        .
            //      .               .        .        .
            //      .               .        .        .
            // _nparticles    0.00157  0.01234  0.00313   ...  0.00112
            
            _debug_theta_distr.clear();
            _debug_theta_distr.resize(_nparticles);
        }
        
        unsigned i = 0;
        for (const auto & p : _particle_list) {
            unsigned c = p.getCount();
            double theta = p.getSpeciesForestConst().getThetaMean();
            for (unsigned j = 0; j < c; ++j) {
                assert(_debug_theta_distr[i].size() == step);
                _debug_theta_distr[i++].push_back(theta);
            }
        }
        assert(i == _nparticles);
    }
    
    inline double SMC::calcPosteriorMeanTheta() const {
        double sum_weights = 0.0;
        double sum_weighted_thetas = 0.0;
        for (const auto & p : _particle_list) {
            const SpeciesForest & sf = p.getSpeciesForestConst();
            double theta_mean = sf.getThetaMean();
            double count = (double)p.getCount();
            sum_weighted_thetas += theta_mean*count;
            sum_weights += count;
        }
        assert(sum_weights > 0.0);
        double posterior_mean = sum_weighted_thetas/sum_weights;
        return posterior_mean;
    }
#endif
    
    inline void SMC::summarize() {
        string prefix = "1st";
        if (isConditionalMode())
            prefix = "2nd";
                    
        output(format("log(marginal likelihood) = %.6f\n") % _log_marg_like, 1);
        //output("\nSpecies trees saved to file \"final-species-trees.tre\"\n", 1);
        
        //temporary!
        Particle & p = *(_particle_list.begin());
        vector<Forest::coalinfo_t> coalinfo_vect;
        p.recordAllForests(coalinfo_vect);
        double log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        output(format("log(coalescent likelihood) = %.9f\n") % log_coallike, 2);
        
        // //temporary!
        // Read in the reference species tree
        //G::_nexus_taxon_map.clear();
        //map<unsigned,unsigned> species_taxon_map;
        //vector<string> species_tree_names;
        //vector<string> species_newicks;
        //// if G::_species_names is empty, readTreeFile will populate it using the taxa block
        //Forest::readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, species_taxon_map, //species_tree_names, species_newicks);
        //p.setSpeciesTree(species_newicks[0]);
        //coalinfo_vect.clear();
        //p.recordAllForests(coalinfo_vect);
        //log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        //output(format("log(coalescent likelihood) for true species tree = %.9f\n") % log_coallike, 2);
        
        string sfn = str(format("%s-final-species-trees") % prefix);
        saveAllSpeciesTrees(sfn, _particle_list, G::_treefile_compression);
        
#if defined(EST_THETA)
        output(format("theta (posterior mean) = %.6f\n") % calcPosteriorMeanTheta(), 1);
        
        // Create a file that can be viewed in Tracer that shows
        // the theta values for each particle in the population for
        // each step in the SMC. The steps (columns) take the place
        // of parameters and the particles (rows) take the place of
        // generations/iterations
        ofstream tracef("tracefile.txt");
        
        // Output row of headers
        tracef << "particle";
        for (unsigned j = 0; j < _nsteps; ++j) {
            tracef << '\t' << "step" << j;
        }
        tracef << endl;
        
        // Output _nparticles rows each comprising _nsteps theta values
        for (unsigned i = 0; i < _debug_theta_distr.size(); ++i) {
            tracef << i;
            for (unsigned j = 0; j < _debug_theta_distr[i].size(); ++j) {
                tracef << '\t' << _debug_theta_distr[i][j];
            }
            tracef << endl;
        }
        tracef.close();
#endif
        
        //saveParamsForLoRaD(prefix);
        
        if (isJointMode()) {
            for (unsigned g = 0; g < G::_ngenes; g++) {
                output(format("Gene trees for locus %d saved to file \"final-gene%d-trees.tre\"\n") % (g+1) % (g+1), 1);
                
                string fnprefix = str(format("%s-final-gene%d-trees") % prefix % (g+1));
                saveAllGeneTrees(g, fnprefix, _particle_list, G::_treefile_compression);
            }
        }
        
        map<string, tuple<unsigned, double, double, double, double> > m;
        bool ok = compareToReferenceTrees(_particle_list, m);
        
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
                % "Genes"
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
                % G::_ngenes
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
        }
    }
    
    inline double SMC::filterParticles(unsigned step, unsigned locus, list<Particle> & particle_list, vector<double> & log_weights, vector<unsigned> & counts, vector<unsigned> & rnseeds, bool refresh_species_tree) {
        //MARK: Proj::filterParticles
        // Sanity checks
        assert(counts.size() == _nparticles);
        assert(log_weights.size() == _nparticles);
                        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = G::calcLogSum(log_weights);
        vector<double> probs(_nparticles, 0.0);
        transform(log_weights.begin(), log_weights.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
                
        // Compute component of the log marginal likelihood due to this step
        _log_marg_like += log_sum_weights - log(_nparticles);
        if (step == 0)
            _log_marg_like += _starting_log_likelihood;
        
        // Compute effective sample size
        double ess = computeEffectiveSampleSize(probs);
        
        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());
        
        // Zero vector of counts storing number of darts hitting each particle
        counts.assign(_nparticles, 0);
                
        // Throw _nparticles darts
        for (unsigned i = 0; i < _nparticles; ++i) {
            //POLWAS double u = rng.uniform();
            double u = rng->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)distance(probs.begin(), it);
            counts[which]++;
        }
                
        // //temporary!
        //unsigned maxn = 0;
        
        // Suppose _nparticles = 20 and counts looks like this:
        //                   <----------- 1st particle -----------><-- 2nd particle --->
        //         counts = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 12}
        //                   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18  19
        //                              nonzeros = [3,6,11]            nonzeros = [19]
        //
        // Here, _particle_list originally contains 2 particles:
        //
        //   Particle  Count
        //          0     12
        //          1      8
        //
        // After filtering, _particle_list will contain the following particle counts:
        //
        //                    Copied   Using
        //   Particle  Count    from    seed
        //          0      1       0       3
        //          1      2       0       6
        //          2      5       0      11
        //          3     12       1      19
        unsigned num_parent_particles = (unsigned)particle_list.size();
        unsigned i = 0;
        for (unsigned j = 0; j < num_parent_particles; j++) {
            auto pit = particle_list.begin();
            advance(pit, j);
            Particle & p = *pit;
            
            // Next n elements of counts belong to particle p
            unsigned n = p.getCount();
            
            // Find all non-zero counts associated with particle p
            vector<unsigned> nonzeros;
            findNonZeroCountsInRange(nonzeros, counts, i, i + n);
                        
            for (auto k : nonzeros) {
                // Make a copy of p and push it onto particle_list
                particle_list.push_back(p);
                
                // Get reference to particle just created
                Particle & plast = *(particle_list.rbegin());
                                                    
                //TODO: these lines have no effect because _proposed_gene, _proposed_spp, _proposed_first, and _proposed_second will all be overwritten
                //plast.setLastProposedGene(_proposed_gene[k]);
                //plast.setLastProposedSpecies(_proposed_spp[k]);
                //plast.setLastProposedFirstIndex(_proposed_first[k]);
                //plast.setLastProposedSecondIndex(_proposed_second[k]);

                if (_mode == SMC::SPECIES_AND_GENE) {
                    // Advance to same coalescence event created when
                    // log weight was calculated but this time keep it
                    // by setting make_permanent to true.
                    double starting_height = -1.0; // negative value means don't rebuild species tree
                    if (refresh_species_tree) {
                        starting_height = p.getMaxGeneTreeHeight();
                    }
                    plast.proposeCoalescence(rnseeds[k], step, locus, k,  /*compute_partial*/true, /*make_permanent*/true, starting_height);
                }
                else if (_mode == SMC::SPECIES_GIVEN_GENE) {
                    // Advance to same speciation event created when
                    // log weight was calculated but this time keep it
                    // by setting make_permanent to true.
                    auto proposed = plast.proposeSpeciation(rnseeds[k], step, k, /*make_permanent*/true);
                    
                    //temporary!
                    assert(G::_speclog.count(step) > 0);
                    assert(G::_speclog.at(step).size() > 0);
                    G::SpecLog & speclog_element = *(G::_speclog.at(step).rbegin());
                    speclog_element._freq = counts[k];
                }
                else {
                    throw XProj(format("Unknown SMC mode encountered (%d) in SMC::filterParticles function") % _mode);
                }
                                
                // Set count for new particle
                plast.setCount(counts[k]);

                // //temporary! Useful for debugging filtering
                //output(format("~~~| Copying %d (%d) -> %d (%d)\n") % j % p.getCount() % (particle_list.size()-1) % plast.getCount(), 2);
            }
                            
            // Flag original particle for deletion
            p.setCount(0);
            
            i += n;
        }
        
        // Eliminate all particles with a count of 0
        pruneParticles(particle_list);

        // //temporary!
        //cerr << "~~> maxn = " << maxn << endl;

        return ess;
    }

    inline unsigned SMC::countDistinctGeneTreeTopologies() {
        set<Forest::treeid_t> unique_topologies;
        for (auto & p : _particle_list) {
            for (unsigned g = 0; g < G::_ngenes; g++) {
                GeneForest & gf = p.getGeneForest(g);
                Forest::treeid_t id;
                gf.storeSplits(id);
                unique_topologies.insert(id);
            }
        }
        return (unsigned)unique_topologies.size();
    }
    
    inline void SMC::saveParamsForLoRaD(string prefix) {
        output("Saving parameters for LoRaD...\n", 1);
        output("  Parameter file is \"params-lorad.txt\"\n", 1);
        output("  Run LoRaD using \"rscript lorad.R\"\n", 1);
        unsigned num_distinct = countDistinctGeneTreeTopologies();
        if (num_distinct == 1) {
            output("  Topology identical for all gene trees in all particles\n",1);
        }
        else {
            output("  *** warning ***: more than one distinct gene tree topology\n",1);
        }
        
        // Open the file into which parameters will be saved
        // in a format useful to LoRaD
        ofstream loradf(str(format("%s-params-lorad.txt") % prefix));

        // Save column names on first line
        vector<string> names;
        names.push_back("index");
        names.push_back("particle");
        names.push_back("logFelsL");
        names.push_back("logCoalL");
        names.push_back("logPrior");
        G::getAllParamNames(names);
        loradf << boost::join(names, "\t") << endl;
        
        unsigned particle_object_index = 0;
        unsigned particle_index = 0;
        
        // //temporary!
        //unsigned tmp_nparticle_objects = (unsigned)_particle_list.size();
        
        for (auto & p : _particle_list) {
            // Stores tuple (height, gene + 1, left species, right species)
            // for each join in either species tree or any gene tree.
            // To get 0-offset gene index, subtract
            // one from second element of tuple (0 means tuple represents a
            // species tree join).
            vector<Forest::coalinfo_t> coalinfo_vect;
            
            // Add species tree joins to coalinfo_vect
            //SpeciesForest & sf = p.getSpeciesForest();
            //sf.heightsInternalsPreorders();
            //sf.saveCoalInfo(coalinfo_vect);

            // Add gene tree joins to coalinfo_vect
            double log_likelihood = 0.0;
            for (unsigned g = 0; g < G::_ngenes; g++) {
                GeneForest & gf = p.getGeneForest(g);
                gf.heightsInternalsPreorders();
                //gf.forgetSpeciesTreeAbove(0.0);
                //gf.saveCoalInfo(coalinfo_vect);
                log_likelihood += gf.calcLogLikelihood(false);
            }
            log_likelihood *= G::_phi;
            
            // Sort coalinfo_vect from smallest to largest height
            //sort(coalinfo_vect.begin(), coalinfo_vect.end());
            
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
            
            for (unsigned g = 0; g < G::_ngenes; g++) {
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
            unsigned n = p.getCount();
            for (unsigned i = 0; i < n; i++) {
                particle_index++;
                loradf << str(format("%d\t%d\t%.9f\t%.9f\t%.9f\t%s")
                    % particle_index
                    % particle_object_index
                    % log_likelihood
                    % log_gene_tree_prior
                    % log_species_tree_prior
                    % boost::join(params, "\t")) << endl;
            }
            particle_object_index++;
        }
        loradf.close();
        
        ofstream loradRf("lorad.R");
        loradRf << "library(lorad)\n";
        loradRf << "colspec <- c(";
        loradRf << "\"index\" = \"iteration\",";
        loradRf << "\"particle\" = \"ignore\",";
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
        
        ofstream conf("loradml.conf");
        conf << "paramfile = params-lorad.txt\n";
        conf << "plotfprefix = lnlplot\n";
        conf << "trainingfrac = 0.5\n";
        conf << "coverage = 0.1\n";
        conf << "colspec = iteration index\n";
        conf << "colspec = ignore particle\n";
        conf << "colspec = posterior logFelsL\n";
        conf << "colspec = posterior logCoalL\n";
        conf << "colspec = posterior logPrior\n";
        for (unsigned i = 4; i < names.size(); i++) {
            conf << "colspec = positive " << names[i] << "\n";
        }
        conf.close();
    }
    
    struct SpeciesTreeDetails {
        unsigned _count;
        map<G::species_t, double> _theta_map;
        double _log_coallike;
        
        SpeciesTreeDetails() : _count(0), _log_coallike(0.0) {}
    };

    inline void SMC::saveAllSpeciesTrees(string fnprefix, const list<Particle> & particle_list, unsigned compression_level) {
        assert(compression_level >= 0 && compression_level <= 2);
        typedef tuple<unsigned, double, string, string, string> treeinfo_t;
        treeinfo_t treeinfo;
        vector<treeinfo_t> treeinfo_vect;
        if (compression_level == 2) {
            // Save unique species trees along with their frequency in <fnprefix>.tre
            // Save coalescent log likelihood and thetas in <fnprefix>.txt
            map<string,vector<SpeciesTreeDetails> > tree_info;
            for (const Particle & p : particle_list) {
                SpeciesTreeDetails info;
                
                // Get count for this particle
                info._count = p.getCount();
                
                // Get newick tree description for this species tree
                string newick = p.getSpeciesForestConst().makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                
                // Calculate log coalescent likelihood for this species tree
                vector<Forest::coalinfo_t> coalinfo_vect;
                p.recordAllForests(coalinfo_vect);
                bool integrate_out_thetas = _mode==SPECIES_GIVEN_GENE;
                                
                info._log_coallike = p.calcLogCoalescentLikelihood(coalinfo_vect, integrate_out_thetas, /*verbose*/false);
                
                // Record map of thetas for each species
                info._theta_map = p.getSpeciesForestConst().getThetaMapConst();
                
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
            for (const Particle & p : particle_list) {
                unsigned c = p.getCount();
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

    inline void SMC::saveAllGeneTrees(unsigned gene_index, string fnprefix, list<Particle> & particle_list, unsigned compression_level) {
        assert(compression_level >= 0 && compression_level <= 2);
        typedef tuple<unsigned, double, string, string, string> treeinfo_t;
        treeinfo_t treeinfo;
        vector<treeinfo_t> treeinfo_vect;
        if (compression_level == 2) {
            // Save only unique newick strings
            map<string,vector<GeneTreeDetails> > tree_info;
            for (Particle & p : particle_list) {
                GeneTreeDetails info;
                
                // Get count for this particle
                info._count = p.getCount();
                
                // Get newick tree description for this gene tree
                assert(gene_index < p.getGeneForestsConst().size());
                GeneForest & gf = p.getGeneForests()[gene_index];
                string newick = gf.makeNewick(/*precision*/9, /*use names*/true, /*coalunits*/false);
                
                // Calculate log-likelihood for this gene tree
                info._log_likelihood = gf.calcLogLikelihood(false);
                
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
            for (const Particle & p : particle_list) {
                unsigned c = p.getCount();
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

    inline bool SMC::compareToReferenceTrees(list<Particle> particle_list, map<string, tuple<unsigned, double, double, double, double> > & m) {
        // Bail out if no reference tree was specified
        if (G::_species_tree_ref_file_name.empty())
            return false;
            
        // Read in the reference species tree
        G::_nexus_taxon_map.clear();
        vector<string> tree_names;
        vector<string> newicks;
        SpeciesForest::readTreefile(G::_species_tree_ref_file_name, /*skip*/0, G::_species_names, G::_nexus_taxon_map, tree_names, newicks);
                
        SpeciesForest ref;
        ref.buildFromNewick(newicks[0]);
        double ref_height = ref.getHeight();
        
        for (Particle & p : particle_list) {
            unsigned c = p.getCount();
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
#if defined(DEBUGGING_SANITY_CHECK)
                test.heightsInternalsPreorders();
                double hcheck = test.getHeight();
                assert(fabs(hcheck - test_height) < 0.001);
#endif
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

    inline void SMC::findNonZeroCountsInRange(vector<unsigned> & nonzeros, const vector<unsigned> & counts, unsigned begin_index, unsigned end_index) const {
        nonzeros.clear();
        for (unsigned k = begin_index; k < end_index; k++) {
            if (counts[k] > 0) {
                nonzeros.push_back(k);
            }
        }
    }

    inline void SMC::pruneParticles(list<Particle> & particle_list) {
        list<Particle>::iterator it = particle_list.begin();
        while (it != particle_list.end()) {
            if (it->getCount() == 0) {
                // Iterator post-increment returns previous iterator,
                // which can now be erased since the it has moved on
                particle_list.erase(it++);
            }
            else {
                // Leave this Particle intact because it has non-zero count
                ++it;
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

    void SMC::dumpParticles(SMC & ensemble) {
        // Dump _particle_list into ensemble, leaving _particle_list empty
        move(_particle_list.begin(), _particle_list.end(), back_inserter(ensemble._particle_list));
        _particle_list.clear();
    }
    
}
