#pragma once

namespace proj {

    inline Particle::Particle() : _gene_trees(nullptr) {
        _lot.reset(new Lot());
        clear();
    }

    inline Particle::Particle(const Particle & other) : _gene_trees(nullptr) {
        _lot.reset(new Lot());
        copyParticleFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _smc = nullptr;
        _prev_log_coallike = 0.0;
        _log_weight = 0.0;
        _gene_forests.clear();
        _species_forest.clear();
        //_nspeciations = 0;
        _last_event = LAST_EVENT_UNDEFINED;
    }
    
    inline void Particle::setSMC(SMC * smc) {
        _smc = smc;
    }
            
    inline void Particle::setData(Data::SharedPtr data) {
        _data = data;
    }
            
    inline void Particle::resetSpeciesForest() {
        _species_forest.createTrivialForest();
        assert(_smc);
        if (_smc->isJointMode()) {
            _species_forest.setJointEstimation(true);
        }
        else {
            _species_forest.setJointEstimation(false);
        }
    }
    
    inline void Particle::setGeneTrees(vector<GeneForest> & gtvect) {
        // //temporary!
        // cerr << "Before reassignment:" << endl;
        // cerr << "  gtvect.size() is        " << gtvect.size() << endl;
        // cerr << "  gtvect is at            " << G::memoryAddressAsString(&gtvect) << endl;
        // cerr << "  &gtvect[0] is at        " << G::memoryAddressAsString(&gtvect[0]) << endl;
        // cerr << "  _gene_forests.size() is " << _gene_forests.size() << endl;
        // cerr << "  _gene_forests is at     " << G::memoryAddressAsString(&_gene_forests[0]) << endl;
        // cerr << "  _gene_trees is at       " << G::memoryAddressAsString(_gene_trees) << endl;
        
        _gene_trees = &gtvect;

        // //temporary!
        // cerr << "After reassignment:" << endl;
        // cerr << "  gtvect.size() is         " << gtvect.size() << endl;
        // cerr << "  gtvect is at             " << G::memoryAddressAsString(&gtvect) << endl;
        // cerr << "  &gtvect[0] is at         " << G::memoryAddressAsString(&gtvect[0]) << endl;
        // cerr << "  _gene_forests.size() is  " << _gene_forests.size() << endl;
        // cerr << "  _gene_forests is at      " << G::memoryAddressAsString(&_gene_forests[0]) << endl;
        // cerr << "  _gene_trees is at        " << G::memoryAddressAsString(_gene_trees) << endl;
        // cerr << "  &(*_gene_trees)[0] is at " << G::memoryAddressAsString(&(*_gene_trees)[0]) << endl;

        debugCheckAllPrevSpeciesStacksEmpty();
    }

    //inline void Particle::forgetSpeciesTreeAbove(double height) {
    //    assert(_gene_trees);
    //    for (auto & gt : *_gene_trees) {
    //        gt.refreshAllHeightsAndPreorders(); //TODO: can do without this call?
    //        gt.forgetSpeciesTreeAbove(height);
    //    }
    //}

    inline void Particle::resetGeneForests(bool compute_partials) {
        assert(G::_ntaxa > 0);
        assert(G::_ngenes > 0);
        _gene_forests.clear();
        _gene_forests.resize(G::_ngenes);
        unsigned g = 0;
        for (auto & gf : _gene_forests) {
            gf.setParticle(this);
            gf.setData(_data);
            gf.setRelRate(G::_relrate_for_gene[g]);
            gf.setGeneIndex(g++);
            gf.createTrivialForest(compute_partials);
        }
    }

    inline void Particle::threadComputePartials(unsigned first, unsigned last) {
        for (unsigned k = first; k < last; k++) {
            _gene_forests[k].computeAllPartials();
        }
    }
    
    inline double Particle::calcLogLikelihood() {
        double log_likelihood = 0.0;
        for (unsigned g = 0; g < G::_ngenes; g++) {
            log_likelihood += _gene_forests[g].calcLogLikelihood();
        }
        return log_likelihood;
    }
    
    struct bitless {
        bool operator()(const G::species_t a, const G::species_t b) const {
            bool returned_value = ((a & b) > 0 ? false : a < b);
            return returned_value;
        }
    };
    
    //struct biteq {
    //    bool operator==(const G::species_t a, const G::species_t b) const {
    //        bool returned_value = ((a & b) > 0 ? true : false);
    //        return returned_value;
    //    }
    //};
    
    inline void Particle::recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) {
        // Record gene trees
        for (unsigned g = 0; g < G::_ngenes; g++) {
            GeneForest & gf = (*_gene_trees)[g];
            gf.saveCoalInfo(coalinfo_vect);
        }
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Record species tree
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Modify coalescent elements according to species tree
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);
    }
    
    inline double Particle::calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) {
        //MARK: calcLogCoalescentLikelihood
        // This function assumes gene forests are complete gene trees (not partial states) and
        // that preorders and heights have been precalculated.
        
        // The symbol b is used for a branch (i.e. segment) of the species tree to match the use of b
        // by Graham Jones (2017). Thus, a branch b is synonymous with a a leaf or ancestral species.
        
        // Uses eq. 4, p. 454, in G. Jones. 2017. Algorithmic improvements
        // to species delimitation and phylogeny estimation under the
        // multispecies coalescent. J. Math. Biol. 74:447-467.
        // doi: 10.1007/s00285-016-1034-0
        
        // eq. 4: log p(G | alpha, beta, sigma) = sum_b { log(r_b)
        //          + alpha log(sigma beta) - (alpha + q_b) log(sigma beta + gamma_b)
        //          + log Gamma(alpha + q_b) - log Gamma(alpha) }
        // where (eq. 2)
        //    p_j      = ploidy (2 = diploid, 1 = haploid) for gene j
        //    k_jb     = no. coal. events for gene j, species tree edge b
        //    c_jbi    = sojourn i to i+1 coalescence for gene j, edge b
        //    q_b      = sum_j k_jb
        //    log(r_b) = sum_j -k_jb log(p_j)
        //    gamma_b  = sum_j (1/p_j) sum_i^{k_jb} binom{n_jb - i}{2} c_jbi
        //
        // If p_j = 2 for all genes:
        //    log(r_b) = -log(2) q_b
        //    gamma_b  = 0.5 sum_j sum_i^{k_jb} binom{n_jb - i}{2} c_jbi
        //
        // For compatibility with starbeast3, we assume:
        //    alpha = 2
        //    beta  = (mean of theta)/4
        //    sigma = 1
        
        // //temporary!
        //Forest::debugShowCoalInfo("Particle::calcLogCoalescentLikelihood", coalinfo_vect);
        
        // Create maps to hold quantities needed by Graham Jones' (2017) formula
        // for coalescent likelihood integrating out theta.
        //typedef map<G::species_t, unsigned, bitless> nmap;
        typedef map<G::species_t, double> dmap;
        
        // n_jb[g][b] holds current number of lineages for gene g, branch b
        vector<dmap> n_jb(G::_ngenes);
        
        // log_r_b[b] holds sum of log(r_b) for branch b
        dmap log_r_b;

        // q_b[b] holds number of coalescent events for branch b
        dmap q_b;

        // gamma_b[b] holds cumulative gamma_b for branch b
        dmap gamma_b;
        
        // Initialize n_jb for gene 0
        for (auto & nm : G::_taxon_names) {
            unsigned i = G::_taxon_to_species[nm];
            G::species_t b = (G::species_t)1 << i;
            if (n_jb[0].count(b) == 0)
                n_jb[0][b] = 1;
            else
                n_jb[0][b]++;
        }
        
        // Initialize n_jb for other genes by copying n_jb for gene 0
        for (unsigned g = 1; g < G::_ngenes; g++) {
            n_jb[g] = n_jb[0];
        }
        
        // prev_height[g][b] holds previous height for gene g and branch b
        vector<map<G::species_t, double> > prev_height(G::_ngenes);
        
        // Ploidy for each gene (currently all genes assumed to be diploid)
        vector<double> p_j(G::_ngenes, 2.0);
        
        // Vector branches stores branches (including ancestral ones)
        vector<G::species_t> branches;
        for (unsigned i = 0; i < G::_nspecies; i++) {
            G::species_t b = (G::species_t)1 << i;
            branches.push_back(b);
        }
        
        // Vector bavail stores branches that are in the species tree at the current height
        set<G::species_t, bitless> bavail(branches.begin(), branches.end());
        
        if (verbose) {
            output("\nParticle::calcLogCoalescentLikelihood:\n", 1);
            output(format("%12s %12s %24s %12s %12s %12s %12s\n") % "height" % "gene" % "spp" % "n_jbk" % "c_jbk" % "log(r_b)" % "q_b", 1);
        }

        // Walk down coalinfo_vect, accumulating Graham Jones r_b, q_b, and gamma_b
        for (auto & cinfo : coalinfo_vect) {
            double               height = get<0>(cinfo);
            unsigned             gene_plus_one   = get<1>(cinfo);
            int                  gene   = gene_plus_one - 1;
            
            ostringstream oss;
            vector<G::species_t> & b_vect = get<2>(cinfo);
            copy(b_vect.begin(), b_vect.end(), ostream_iterator<G::species_t>(oss, " "));
            
            G::species_t banc = 0;
            for (auto b : b_vect)
                banc |= b;
                
            if (gene_plus_one == 0) {
                // Speciation event
                
                // If banc represents a new branch in the species tree,
                // add it to the vector of branches
                if (find(branches.begin(), branches.end(), banc) == branches.end()) {
                    branches.push_back(banc);
                }
                    
                for (unsigned g = 0; g < G::_ngenes; g++) {
                    // Account for non-coalescence since the last coalescent
                    unsigned nsum = 0;
                    for (auto b : b_vect) {
                        if (prev_height[g].count(b) == 0) {
                            prev_height[g][b] = 0.0;
                        }
                        double nb = n_jb[g][b];
                        nsum += nb;
                        double cb = height - prev_height[g].at(b);
                        double gammab = 2.0*nb*(nb-1)*cb/p_j[g];
                        gamma_b[b] += gammab;
                    }
                    
                    // Record the beginning height of the new ancestral branch banc
                    if (prev_height[g].count(banc) == 0) {
                        prev_height[g][banc] = height;
                    }
                    
                    // Pool lineages from descendant species to create
                    // counts for the new ancestral species
                    n_jb[g][banc] = 0;
                    for (auto b : b_vect) {
                        n_jb[g].erase(b);
                    }
                    n_jb[g][banc] = nsum; //nleft + nright;
                    
                }
                      
                for (auto b : b_vect) {
                    bavail.erase(b);
                }
                bavail.insert(banc);

                if (verbose) {
                    output(format("%12.5f %12d %24s %12s %12s %12s %12s\n") % height % 0 % oss.str() % "-" % "-" % "-" % "-", 1);
                }
            }
            else {
                // Coalescent event
                //auto it = bavail.find(banc);
                auto it = find_if(bavail.begin(), bavail.end(), [banc](G::species_t v){return (banc & v) > 0;});
                
                assert(it != bavail.end());
                G::species_t b = *it;
                
                if (prev_height[gene].count(b) == 0) {
                    prev_height[gene][b] = 0.0;
                }
                double height0 = prev_height[gene].at(b);
                
                // Coalescent event: update Jones quantities
                log_r_b[b] += log(4/p_j[gene]);
                q_b[b]++;
                double n_jbk = n_jb[gene][b];
                double c_jbk = height - height0;
                gamma_b[b] += 2.0*n_jbk*(n_jbk-1)*c_jbk/p_j[gene];

                prev_height[gene][b] = height;
                n_jb[gene][b]--;

                if (verbose) {
                    output(format("%12.5f %12d %24s %12d %12.5f %12.5f %12d\n") % height % gene_plus_one % oss.str() % n_jbk % c_jbk % gamma_b[b] % q_b[b], 1);
                }
            }
        }

        if (verbose) {
            output(format("\n%12s %12s %12s %12s %12s\n") % "b" % "q_b" % "log(r_b)" % "gamma_b" % "logL", 1);
        }
    
        double log_likelihood = 0.0;
        
        if (integrate_out_thetas) {
            double alpha = 2.0; // G::_invgamma_shape;
#if defined(EST_THETA)
            double beta = _species_forest.getThetaMean();
#else
            double beta = G::_theta;
#endif
            double B = (double)branches.size();
            log_likelihood  = B*alpha*log(beta);
            log_likelihood -= B*boost::math::lgamma(alpha);
            for (auto b : branches) {
                log_likelihood += log_r_b[b];
                log_likelihood -= (alpha + q_b[b])*log(beta + gamma_b[b]);
                log_likelihood += boost::math::lgamma(alpha + q_b[b]);
                
                if (verbose) {
                    output(format("%12d %12d %12.5f %12.5f %12.5f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % log_likelihood, 1);
                }
            }

            if (verbose) {
                output(format("\nalpha = %.9f\n") % alpha, 1);
                output(format("beta  = %.9f\n") % beta, 1);
                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood, 1);
            }
        }
        else {
#if defined(EST_THETA)
            double theta_b = _species_forest.getThetaMean();
#else
            double theta_b = G::_theta;
#endif
            double sum_log_rb = 0.0;
            double sum_gamma_b = 0.0;
            unsigned sum_qb = 0;
            for (auto b : branches) {
                sum_qb += q_b[b];
                sum_log_rb += log_r_b[b];
                sum_gamma_b += gamma_b[b];
                log_likelihood += log_r_b[b] - gamma_b[b]/theta_b - q_b[b]*log(theta_b);
                
                if (verbose) {
                    output(format("%12d %12d %12.5f %12.5f %12.5f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % log_likelihood, 1);
                }
            }
            
            if (verbose) {
                output(format("\nsum log(r_b) = %.9f\n") % sum_log_rb, 1);
                output(format("sum q_b      = %d\n") % sum_qb, 1);
                output(format("sum gamma_b  = %d\n") % sum_gamma_b, 1);
                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood, 1);
            }
        }
        
        //_prev_log_coallike = log_likelihood;
                    
        return log_likelihood;
    }
    
    inline void Particle::computeAllPartials() {
//#if defined(USING_MULTITHREADING)
//        vector<thread> threads;
//        for (unsigned i = 0; i < G::_nthreads; i++) {
//            threads.push_back(thread(&Particle::threadComputePartials,
//                this,
//                G::_thread_first_gene[i],
//                G::_thread_last_gene[i])
//            );
//        }
//
//        // The join function causes this loop to pause until the ith thread finishes
//        for (unsigned i = 0; i < threads.size(); i++) {
//            threads[i].join();
//        }
//#else
        for (auto & gf : _gene_forests) {
            gf.computeAllPartials();
        }
//#endif
    }
    
    inline void Particle::stowAllPartials() {
        for (auto & gf : _gene_forests) {
            gf.stowAllPartials();
        }
    }
    
    //inline void Particle::setLastProposedGene(unsigned g) {
    //    _last_proposed_gene = g;
    //}

    //inline unsigned Particle::getLastProposedGene() {
    //    return _last_proposed_gene;
    //}
    
    //inline void Particle::setLastProposedSpecies(G::species_t s) {
    //    _last_proposed_spp = s;
    //}

    //inline G::species_t Particle::getLastProposedSpecies() {
    //    return _last_proposed_spp;
    //}
    
    //inline void Particle::setLastProposedFirstIndex(unsigned f) {
    //    _last_proposed_first_index = f;
    //}

    //inline unsigned Particle::getLastProposedFirstIndex() {
    //    return _last_proposed_first_index;
    //}
    
    //inline void Particle::setLastProposedSecondIndex(unsigned s) {
    //  _last_proposed_second_index = s;
    //}
    
    //inline unsigned Particle::getLastProposedSecondIndex() {
    //    return _last_proposed_second_index;
    //}

    inline PartialStore::partial_t Particle::pullPartial(unsigned gene) {
        assert(gene < _gene_forests.size());

        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
#if defined(USING_MULTITHREADING)
        {
            lock_guard<mutex> guard(G::_mutex);
            ptr = ps.getPartial(gene);
        }
#else
        ptr = ps.getPartial(gene);
#endif
        return ptr;
    }

    inline void Particle::stowPartial(unsigned gene, Node * nd) {
        assert(gene < _gene_forests.size());
        assert(nd);
        assert(nd->_partial);
#if defined(USING_MULTITHREADING)
        {
            lock_guard<mutex> guard(G::_mutex);
            ps.putPartial(gene, nd->_partial);

            // Decrement shared pointer reference count
            nd->_partial.reset();
        }
#else
        ps.putPartial(gene, nd->_partial);
        
        // Decrement shared pointer reference count
        nd->_partial.reset();
#endif
    }

    inline void Particle::finalizeProposalSpeciesForest() {
        //MARK: Particle::finalizeProposalSpeciesForest
        _species_forest.finalizeProposal();
        debugCheckAllPrevSpeciesStacksEmpty();
    }
    
    inline void Particle::finalizeProposalAllForests() {
        //MARK: Particle::finalizeProposalAllForests
        _species_forest.finalizeProposal();
        assert(_smc);
        if (_smc->isJointMode()) {
            for (auto & gf : _gene_forests) {
                gf.finalizeProposal();
            }
        }
        
        debugCheckAllPrevSpeciesStacksEmpty();
    }
    
    inline void Particle::debugCheckAllPrevSpeciesStacksEmpty() const {
#if defined(DEBUGGING_SANITY_CHECK)
        unsigned node_index = 0;
        
        // Ensure that all species tree nodes have an empty _prev_species_stack
        for (auto & nd : _species_forest._nodes) {
            if (nd.canRevertSpecies()) {
                cerr << str(format("Node %d in species forest has a non-empty _prev_species_stack") % node_index) << endl;
            }
            assert(!nd.canRevertSpecies());
            node_index++;
        }

        // Ensure that all gene tree nodes have an empty _prev_species_stack
        for (auto & gene_forest : _gene_forests) {
            node_index = 0;
            for (auto & nd : gene_forest._nodes) {
                if (nd.canRevertSpecies()) {
                    cerr << str(format("Node %d in gene forest %d has nonempty species stack") % node_index % gene_forest.getGeneIndex()) << endl;
                }
                assert(!nd.canRevertSpecies());
                node_index++;
            }
        }
#endif
    }
    
    double Particle::findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect) {
        double min_height = G::_infinity;
        for (auto & ci : coalinfo_vect) {
            double h = get<0>(ci);
            if (h >= hstart) {
                vector<G::species_t> & v = get<2>(ci);
                assert(v.size() == 2);
                if (v[0] != v[1]) {
                    if (h < min_height)
                        min_height = h;
                    break;
                }
            }
        }
        
        // for (auto & gf : _gene_forests) {
        //     for (auto ci : gf._coalinfo) {
        //         double h = get<0>(ci);
        //         if (h >= hstart) {
        //             vector<G::species_t> & v = get<2>(ci);
        //             assert(v.size() == 2);
        //             if (v[0] != v[1]) {
        //                 if (h < min_height)
        //                     min_height = h;
        //                 break;
        //             }
        //         }
        //     }
        // }
        
        assert(min_height < G::_infinity);
        return min_height;
    }
    
    inline pair<double, unsigned> Particle::proposeSpeciation(unsigned seed, unsigned step, unsigned pindex, bool make_permanent) {
        //MARK: Particle::proposeSpeciation
        // This function is only used for proposing speciation events when there are
        // complete gene trees available. It thus draws increments from a truncated
        // exponential distribution where the trunction point is the next height at
        // which at least one coalescent event combines lineages from two different
        // species.
        
        // //temporary!
        // if (pindex == 9) {
        //     cerr << "pindex is 9" << endl;
        // }
        
        setSeed(seed);
        
        // Stores tuple (height, 0, vector of species) for each join in the current species forest.
        // Do not cap with ancestral species at this point.
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect);
        
        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Stores tuple (height, gene + 1, vector of species)
        // for each join in any gene tree.
        // To get 0-offset gene index, subtract 1 from 2nd element of tuple
        // (if 2nd element is 0, then tuple represents a species tree join)
        vector<Forest::coalinfo_t> coalinfo_vect;
        
        // Add gene tree joins to coalinfo_vect
        // Just need coalescent events at this point in order to choose
        // limits for species tree increments
        assert(_gene_trees);
        for (unsigned g = 0; g < G::_ngenes; g++) {
            GeneForest & gf = (*_gene_trees)[g];
            gf.saveCoalInfo(coalinfo_vect);
        }
        
        // Sort coalinfo_vect from smallest to largest height
        sort(coalinfo_vect.begin(), coalinfo_vect.end());
        
        double max_height = get<0>((*coalinfo_vect.rbegin()));

        if (step > 0) {
            // Create speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(_lot, left_spp, right_spp, anc_spp, /*mark*/!make_permanent);
            
            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());

            // Adjust elements of coalinfo_vect affected by species tree joins
            _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);
        }
        
        // Draw a speciation increment dt.
        double h = findHeightNextCoalescentEvent(_species_forest.getHeight(), coalinfo_vect);
        assert(h < max_height);
        double log_weight_factor = chooseTruncatedSpeciesForestIncrement(h, /*mark*/!make_permanent);

        unsigned num_species_tree_lineages = _species_forest.getNumLineages();
        if (num_species_tree_lineages == 2) {
            // Create final speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(_lot, left_spp, right_spp, anc_spp, /*mark*/!make_permanent);
            
            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());

            // Adjust elements of coalinfo_vect affected by species tree joins
            _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);
            
            // // Advise all gene trees of the change in the species tree
            // // Nodes that are reassigned save their previous state
            // // to allow reversion
            //double species_forest_height = _species_forest.getHeight();
            //for (auto & gene_forest : _gene_forests) {
            //    gene_forest.mergeSpecies(species_forest_height, left_spp, right_spp, anc_spp);
            //}
        }

        _last_event = Particle::LAST_EVENT_SPECIATION;
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear();
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Adjust elements of coalinfo_vect affected by species tree joins
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

        // //temporary!
        //Forest::debugShowCoalInfo("===== Particle::proposeSpeciation (after) =====", coalinfo_vect);

        // Compute coalescent likelihood and log weight
        //double prev_log_coallike = _prev_log_coallike;
        double log_coallike = calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        double log_weight = log_coallike - _prev_log_coallike + log_weight_factor;

        if (make_permanent) {
            _prev_log_coallike = log_coallike;
            finalizeProposalAllForests();
        }
        else {
            revertToMarkAllForests();
        }

        return make_pair(log_weight, num_species_tree_lineages);
    }
    
    inline pair<double, unsigned> Particle::proposeCoalescence(unsigned seed, unsigned step, unsigned pindex, bool compute_partial, bool make_permanent) {
        //MARK: Particle::proposeCoalescence
        setSeed(seed);
        
        //clearMarkAllForests();
#if defined(EST_THETA)
        double log_weight_modifier = 0.0;
        if (step == 0) {
            if (G::_theta_mean_fixed > 0.0) {
                _species_forest.setThetaMean(G::_theta_mean_fixed);
            }
            else {
                // Draw mean theta from Exponential(r) proposal distribution,
                // where r = 1/G::_theta_proposal_mean
                _species_forest.drawThetaMean(1.0/G::_theta_proposal_mean);
                
                // Calculate weight modifier that takes account of the fact
                // that the theta mean proposal distribution differs from its
                // prior distribution
                double prior_rate    = 1.0/G::_theta_prior_mean;
                double proposal_rate = 1.0/G::_theta_proposal_mean;
                log_weight_modifier = log(prior_rate) - log(proposal_rate) - (prior_rate - proposal_rate)*_species_forest.getThetaMean();
            }
            
            // Draw values of theta for each leaf species from
            // an InverseGamma(2, _theta_mean) distribution
            // (unless G::_theta_mean_fixed was specified and
            // G::_theta_mean_frozen is true, in which case
            // make every species have the same theta value
            // (equal to G::_theta_mean_fixed)
            _species_forest.drawLineageSpecificThetas();
        }
#endif
        
        advance(step, pindex, compute_partial, make_permanent);
        while (lastEventSpeciation()) {
            advance(step, pindex, compute_partial, make_permanent);
        }
        
        unsigned num_species_tree_lineages = _species_forest.getNumLineages();
        
        if (make_permanent)
            finalizeProposalAllForests();
        else
            revertToMarkAllForests();

        double log_weight = log_weight_modifier + getLogWeight();
        return make_pair(log_weight, num_species_tree_lineages);
    }
    
    //inline double Particle::calcTotalCoalRate(double speciation_increment) {
    inline double Particle::calcTotalCoalRate() {
        double total_rate = 0.0;
        for (auto & gene_forest : _gene_forests) {
            //total_rate += gene_forest.calcTotalRate(_species_tuples, speciation_increment);
            total_rate += gene_forest.calcTotalRate(_species_tuples);
        }
        return total_rate;
    }
    
    inline void Particle::clearMarkAllForests() {
        _species_forest.clearMark();
        assert(_smc);
        if (_smc->isJointMode()) {
            for (auto & gf : _gene_forests) {
                gf.clearMark();
            }
        }
    }

    inline void Particle::revertToMarkSpeciesForest() {
        _species_forest.revertToMark();
    }
    
    inline void Particle::revertToMarkAllForests() {
        _species_forest.revertToMark();
        assert(_smc);
        if (_smc->isJointMode()) {
            for (auto & gf : _gene_forests) {
                gf.revertToMark();
            }
        }
    }
    
    inline void Particle::advanceAllLineagesBy(double dt, bool mark) {
        _species_forest.advanceAllLineagesBy(dt, mark);
        for (auto & gene_forest : _gene_forests) {
            gene_forest.advanceAllLineagesBy(dt, mark);
        }
    }
    
    inline double Particle::priorPrior(unsigned step, unsigned pindex, double total_rate, bool compute_partial) {
        double log_weight = 0.0;
        
        // _species_tuples has already been filled
        // Each tuple entry stores:
        //  0. number of lineages
        //  1. gene index
        //  2. species within gene
        //  3. vector<Node *> holding lineage roots for gene/species combination
        //  4. theta (if #define EST_THETA)
        
        // Prior-prior chooses a locus and species in which to have a coalescence
        // from the prior. The probability of coalescence in one particular
        // gene-species pair i is
        //   p(i) = r_i/total_rate
        // where
        //   r_i = n_i*(n_i-1)/theta
        //   total_rate = sum_i r_i
        //   n_i = number of lineages in pair i.
        vector<double> probs(_species_tuples.size());
        transform(_species_tuples.begin(), _species_tuples.end(), probs.begin(), [total_rate](Node::species_tuple_t & spp_tuple){
            double n = (double)get<0>(spp_tuple);
#if defined(EST_THETA)
            double theta = (double)get<4>(spp_tuple);
            return n*(n - 1.0)/(theta*total_rate);
#else
            return n*(n - 1.0)/(G::_theta*total_rate);
#endif
        });
                
#if defined(DEBUGGING_SANITY_CHECK)
        double check = accumulate(probs.begin(), probs.end(), 0.0);
        assert(fabs(check - 1.0) < 0.0001);
#endif
                                
        // Choose which gene-species pair in which the coalescence will happen
        unsigned which = G::multinomialDraw(_lot, probs);
        assert(which < probs.size());
        
        // Get a reference to the relevant gene forest
        unsigned g = get<1>(_species_tuples[which]);
        GeneForest & gf = _gene_forests[g];
        
        // Get the species involved
        G::species_t spp = get<2>(_species_tuples[which]);
            
        // Pull next available node
        Node * anc_node = gf.pullNode();
        if (compute_partial)
            anc_node->_partial = pullPartial(g);
                
        // Create variables in which to store indexes (relative
        // to spp) of nodes joined
        unsigned first_index = 0;
        unsigned second_index = 0;
                
        // Get vector of nodes in the specified species spp
        unsigned i = 0;
        unsigned j = 0;
        Node * first_node = nullptr;
        Node * second_node = nullptr;
        if (gf._lineages_within_species.count(spp) == 0)
            throw XProj(gf.lineagesWithinSpeciesKeyError(spp));
        else {
            auto & node_vect = gf._lineages_within_species.at(spp);
            unsigned n = (unsigned)node_vect.size();
            assert(n > 1);
            
            // Choose a random pair of lineages to join
            pair<unsigned,unsigned> chosen_pair = _lot->nchoose2(n);
            i = chosen_pair.first;
            j = chosen_pair.second;
            
            // Return the chosen pair in the supplied reference variables
            first_index = i;
            second_index = j;
            
            // Join the chosen pair of lineages
            first_node  = node_vect[i];
            second_node = node_vect[j];
            gf.joinLineagePair(anc_node, first_node, second_node);
        }

        anc_node->setSpecies(spp);
        
        assert(first_node->getSpecies() == spp);
        assert(second_node->getSpecies() == spp);
        
        if (compute_partial) {
            // Compute partial likelihood array of ancestral node
            assert(_data);
            assert(anc_node->_left_child);
            assert(anc_node->_left_child->_right_sib);
            assert(anc_node->_partial);

            log_weight = gf.calcPartialArray(anc_node);
            assert(!isnan(log_weight));
            assert(!isinf(log_weight));
        }
        
        gf.addCoalInfoElem(anc_node, gf._mark_coalinfo);
        
        gf._mark_anc_nodes.push(anc_node);
        gf._mark_left_right_pos.push(make_pair(first_index, second_index));
        
        //_last_proposed_gene = g;
        //_last_proposed_spp = spp;
        //_last_proposed_first_index = first_index;
        //_last_proposed_second_index = second_index;
        
        return log_weight;
    }

    inline double Particle::priorPost(unsigned step, unsigned pindex, bool make_permanent) {
    //    // Prior-post tries every possible join within each locus-species pair,
    //    // then chooses one of these possibilities using a multinomial draw from
    //    // the normalized weights. The log weight of the proposal is the log of
    //    // the normalizing constant for the multinomial draw.
        double log_weight = 0.0;
        throw XProj("Particle::priorPost not yet implemented");
    //    
    //    if (make_permanent) {
    //        unsigned                       g = _last_proposed_gene;
    //        G::species_t         spp = _last_proposed_spp;
    //        unsigned             first_index = _last_proposed_first_index;
    //        unsigned            second_index = _last_proposed_second_index;
    //        
    //        // Get reference to correct gene forest
    //        GeneForest & gf = _gene_forests[g];
    //        
    //        // Pull next available node
    //        Node * anc_node = gf.pullNode();
    //        anc_node->_partial = pullPartial(g);
    //        
    //        // Join the chosen pair of lineages
    //        if (gf._lineages_within_species.count(spp) == 0)
    //            throw XProj(gf.lineagesWithinSpeciesKeyError(spp));
    //        else {
    //            auto & node_vect = gf._lineages_within_species.at(spp);
    //            Node * first_node  = node_vect[first_index];
    //            Node * second_node = node_vect[second_index];
    //            gf.joinLineagePair(anc_node, first_node, second_node);
    //        }
    //        
    //        anc_node->setSpecies(spp);
    //
    //        // Recalculate partial for anc_node
    //        double logw = gf.calcPartialArray(anc_node);
    //        assert(!isnan(logw));
    //        assert(!isinf(logw));
    //
    //        gf._mark_anc_nodes.push(anc_node);
    //        gf._mark_left_right_pos.push(make_pair(first_index, second_index));
    //    }
    //    else {
    //        // Create a vector to hold the log weight from every possible join.
    //        vector<double> log_weights;
    //        
    //        // Create a map relating an index into log_weights (key) to a tuple (value)
    //        // containing:
    //        //  <0> the index into _species_tuples
    //        //  <1> the first index of the join, and
    //        //  <2> the second index of the join.
    //        map<unsigned, tuple<unsigned, unsigned, unsigned> > index_map;
    //        
    //        // The element at index 3 of _species_tuples provides
    //        // pointers to the root nodes of each lineage in that
    //        // locus-species combination
    //        unsigned log_weight_index = 0;
    //        unsigned species_tuples_index = 0;
    //        Node * first_node = nullptr;
    //        Node * second_node = nullptr;
    //        for (auto t : _species_tuples) {
    //            // Get number of lineages
    //            unsigned n = get<0>(t);
    //            
    //            // Get locus
    //            unsigned g = get<1>(t);
    //            
    //            // Get species
    //            G::species_t spp = get<2>(t);
    //            
    //            // Get vector of nodes in the specified species spp
    //            Node::ptr_vect_t & node_vect = get<3>(t);
    //            assert(n == node_vect.size());
    //            assert(n > 1);
    //            
    //            // Get reference to correct gene forest
    //            GeneForest & gf = _gene_forests[g];
    //            
    //            // Pull next available node
    //            Node * anc_node = gf.pullNode();
    //            anc_node->_partial = pullPartial(g);
    //
    //            // Visit each possible pair of nodes to join and compute
    //            // log weight of that join
    //            for (unsigned i = 0; i < n - 1; i++) {
    //                for (unsigned j = i + 1; j < n; j++) {
    //                    index_map[log_weight_index] = make_tuple(species_tuples_index, i, j);
    //                    
    //                    // Join the chosen pair of lineages
    //                    first_node  = node_vect[i];
    //                    second_node = node_vect[j];
    //
    //                    // Coalescent events should not cross species boundaries
    //                    assert(first_node->getSpecies() == spp);
    //                    assert(second_node->getSpecies() == spp);
    //                    
    //                    gf.joinLineagePair(anc_node, first_node, second_node);
    //
    //                    // Compute partial likelihood array of ancestral node
    //                    assert(anc_node->_left_child);
    //                    assert(anc_node->_left_child->_right_sib);
    //                    assert(anc_node->_partial);
    //
    //                    double logw = gf.calcPartialArray(anc_node);
    //                    assert(!isnan(logw));
    //                    assert(!isinf(logw));
    //                    log_weights.push_back(logw);
    //                    assert(log_weights.size() == log_weight_index + 1);
    //                    
    //                    // Unjoin the chosen pair of lineages
    //                    gf.unjoinLineagePair(anc_node, first_node, second_node);
    //                    
    //                    log_weight_index++;
    //                }   // j loop
    //            } // i loop
    //            
    //            // Return partial
    //            gf.stowPartial(anc_node);
    //
    //            // Return anc node
    //            gf.stowNode(anc_node);
    //            anc_node = nullptr;
    //            
    //            species_tuples_index++;
    //        } // species_tuples loop
    //        
    //        // Compute log of sum of weights
    //        double sum_log_weights = G::calcLogSum(log_weights);
    //                          
    //        // Create multinomial probability distribution by normalizing
    //        // the log weights
    //        vector<double> probs(log_weights.size(), 0.0);
    //        transform(log_weights.begin(), log_weights.end(), probs.begin(), [sum_log_weights](double logw){return exp(logw - sum_log_weights);});
    //        
    //        // Choose one pair
    //        unsigned which = G::multinomialDraw(_lot, probs);
    //        auto chosen = index_map[which];
    //        species_tuples_index  = get<0>(chosen);
    //        unsigned first_index  = get<1>(chosen);
    //        unsigned second_index = get<2>(chosen);
    //
    //        unsigned             g         = get<1>(_species_tuples[species_tuples_index]);
    //        G::species_t spp       = get<2>(_species_tuples[species_tuples_index]);
    //        Node::ptr_vect_t &   node_vect = get<3>(_species_tuples[species_tuples_index]);
    //
    //        log_weight = sum_log_weights;
    //
    //        // Get reference to correct gene forest
    //        GeneForest & gf = _gene_forests[g];
    //        
    //        // Pull next available node
    //        Node * anc_node = gf.pullNode();
    //        anc_node->_partial = pullPartial(g);
    //        
    //        // Join the chosen pair of lineages
    //        first_node  = node_vect[first_index];
    //        second_node = node_vect[second_index];
    //        gf.joinLineagePair(anc_node, first_node, second_node);
    //        
    //        anc_node->setSpecies(spp);
    //
    //        //TODO: necessary?
    //        double logw = gf.calcPartialArray(anc_node);
    //        assert(!isnan(logw));
    //        assert(!isinf(logw));
    //
    //        gf._mark_anc_nodes.push(anc_node);
    //        gf._mark_left_right_pos.push(make_pair(first_index, second_index));
    //
    //        _last_proposed_gene = g;
    //        _last_proposed_spp = spp;
    //        _last_proposed_first_index = first_index;
    //        _last_proposed_second_index = second_index;
    //    }
        return log_weight;
    }

    inline void Particle::advance(unsigned step, unsigned pindex, bool compute_partial, bool make_permanent) {
        //MARK: advance
        // Clear species_tuples vector. Each tuple entry stores:
        //  0. number of lineages
        //  1. gene index
        //  2. species within gene
        //  3. vector<Node *> holding lineage roots for gene/species combination
        //  4. theta for species (if #define EST_THETA)
        _species_tuples.clear();
        
        // Draw a speciation increment Delta. Note: speciation_increment will
        // equal "infinity" if species tree is complete.
        auto incr_rate = _species_forest.drawIncrement(_lot);
        double Delta = incr_rate.first;
        //double speciation_rate = incr_rate.second;
        
        // Visit each species within each locus, storing the number of lineages in
        // _species_tuples and computing the total coalescence rate (total_rate).
        //double total_rate = calcTotalCoalRate(Delta);
        double total_rate = calcTotalCoalRate();
        
        // Draw coalescence increment delta ~ Exponential(total_rate)
        double delta = (total_rate > 0.0 ? -log(1.0 - _lot->uniform())/total_rate : G::_infinity);
        
        // If delta < Delta, then a coalescent event happened; otherwise a speciation event happened.
        bool is_speciation = (delta > Delta);
                
        // Advance all forests by increment dt
        double dt = is_speciation ? Delta : delta;

        // Forests save dt to allow reversion
        advanceAllLineagesBy(dt, /*mark*/!make_permanent);
        
        if (is_speciation) {
            _last_event = Particle::LAST_EVENT_SPECIATION;
            
            // Create speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(_lot, left_spp, right_spp, anc_spp, /*mark*/!make_permanent);
            double species_forest_height = _species_forest.getHeight();
                                    
#if defined(EST_THETA)
            _species_forest.drawThetaForSpecies(anc_spp);
#endif

            // Advise all gene trees of the change in the species tree
            // Nodes that are reassigned save their previous state
            // to allow reversion
            for (auto & gene_forest : _gene_forests) {
                gene_forest.mergeSpecies(species_forest_height, left_spp, right_spp, anc_spp);
            }
        }
        else {
            _last_event = Particle::LAST_EVENT_COALESCENCE;

            if (G::_prior_post) {
                _log_weight = priorPost(step, pindex, make_permanent);
            }
            else {
                _log_weight = priorPrior(step, pindex, total_rate, compute_partial);
            }
        }
    }

#if defined(DEBUGGING_SANITY_CHECK)
    inline void Particle::debugStoreForests() {
        _debug_sfbefore = _species_forest.makeNewick(6,true,false);
        _debug_gfbefore.clear();
        for (auto gf : _gene_forests) {
            _debug_gfbefore.push_back(gf.makeNewick(6,true,false));
        }
    }
    
    inline void Particle::debugCheckForests() {
        string sfafter = _species_forest.makeNewick(6,true,false);
        assert(_debug_sfbefore == sfafter);
        unsigned g = 0;
        for (auto gf : _gene_forests) {
            string gfafter = gf.makeNewick(6,true,false);
            assert(gfafter == _debug_gfbefore[g++]);
        }
    }
#endif

    inline void Particle::refreshHeightsInternalsPreorders() {
        _species_forest.heightsInternalsPreorders();
        for (auto & gf : _gene_forests) {
            gf.heightsInternalsPreorders();
        }
    }
    
    inline SpeciesForest & Particle::getSpeciesForest() {
        return _species_forest;
    }
    
    inline const SpeciesForest & Particle::getSpeciesForestConst() const {
        return _species_forest;
    }
    
    inline GeneForest & Particle::getGeneForest(unsigned gene) {
        assert(gene < G::_ngenes);
        return _gene_forests[gene];
    }

    inline vector<GeneForest> & Particle::getGeneForests() {
        return _gene_forests;
    }
    
    inline const vector<GeneForest> & Particle::getGeneForestsConst() const {
        return _gene_forests;
    }
    
    inline const GeneForest & Particle::getGeneForest(unsigned gene) const {
        assert(gene < G::_ngenes);
        return _gene_forests[gene];
    }

    inline void Particle::debugShowMarkVariables(string title) const {
        output(format("\n*** %s ***\n") % title, 3);
        // Show species forest mark variables
        output("  Species forest:\n",3);
        if (!_species_forest._mark_coalinfo.empty()) {
            for (auto & ci : _species_forest._mark_coalinfo) {
                output(format("    coalinfo: h=%.5f (l=%d, r=%d)") % get<0>(ci) % get<2>(ci)[0] % get<2>(ci)[1], 3);
            }
        }
        if (!_species_forest._mark_increments.empty()) {
            output(format("    increments: %d") % _species_forest._mark_increments.size(), 3);
        }
        if (!_species_forest._mark_left_right_pos.empty()) {
            output(format("    left_right_pos: %d") % _species_forest._mark_left_right_pos.size(), 3);
        }
        if (!_species_forest._mark_anc_nodes.empty()) {
            output(format("    anc_nodes: %d") % _species_forest._mark_anc_nodes.size(), 3);
        }
    }
    
    inline unsigned Particle::debugCountNumCoalEvents() const {
        // Returns number of coalescent events over all gene trees
        unsigned total = 0;
        for (auto & gf : _gene_forests) {
            unsigned n = gf.getNumLineages();
            total += G::_ntaxa - n;
        }
        return total;
    }
            
    inline void Particle::debugCheckPartials() const {
        for (auto & gf : _gene_forests) {
            gf.debugCheckPartials();
        }
    }
    
    inline void Particle::debugShowAllGeneForests() const {
        for (auto & gf : _gene_forests) {
            cerr << gf.makeNewick(0, true, false) << endl;
        }
    }
    
    inline void Particle::copyParticleFrom(const Particle & other) {
        clear();

        // Performs a deep copy of other to this particle
        _data = other._data;
        
        // Performs a deep copy of other to this particle
        _smc = other._smc;
        
        if (other._gene_trees) {
            // Copy pointer to complete gene trees
            _gene_trees = other._gene_trees;
        }
        else {
            // Copy gene forests
            _gene_trees = nullptr;
            assert(G::_ngenes == other._gene_forests.size());
            _gene_forests.resize(G::_ngenes);
            for (unsigned i = 0; i < G::_ngenes; i++) {
                _gene_forests[i] = other._gene_forests[i];
                _gene_forests[i].setParticle(this);
            }
        }
        
        // Copy species forest
        _species_forest = other._species_forest;
        
        _prev_log_coallike = other._prev_log_coallike;
        //_nspeciations = other._nspeciations;
        _last_event = other._last_event;
        _count = other._count;
#if defined(USING_MULTITHREADING)
        _xtra = other._xtra;
        _begin_index = other._begin_index;
#endif
        // no need to copy these as they are just temporary workspaces
        //_last_proposed_gene = 0;
        //_last_proposed_spp = 0;
        //_last_proposed_first_index   = 0;
        //_last_proposed_second_index  = 0;
        
        // No need to copy _log_weight
        // No need to copy _lot
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
    inline double Particle::chooseTruncatedSpeciesForestIncrement(double truncate_at, bool mark) {
        double upper_bound = truncate_at - _species_forest.getHeight();
        auto incr_rate = _species_forest.drawTruncatedIncrement(_lot, truncate_at);
        double incr = incr_rate.first;
        double rate = incr_rate.second;
        assert(rate > 0.0);
        _species_forest.advanceAllLineagesBy(incr, mark);
        double log_factor = log(1.0 - exp(-rate*upper_bound));
        return log_factor;
    }
}
