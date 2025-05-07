#pragma once

namespace proj {

    class Particle {
            
        public:
        
            Particle();
            Particle(const Particle & other);
            ~Particle();
            
            void clear();
            
            void setData(Data::SharedPtr data);
            Data::SharedPtr getData() {return _data;}
            Data::SharedPtrConst getDataConst() const {return _data;}
                        
            double setThetas();
            
            void resetSpeciesForest();
            void resetGeneForests(bool compute_partials);

            //void threadComputePartials(unsigned first, unsigned last);
            void computeAllPartials();
            
            bool isEnsembleCoalInfo();
            void buildEnsembleCoalInfo();
            void rebuildCoalInfo();
            void recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) const;
            
            double calcLogLikelihood();
            double calcTotalPrevLogLikelihood();
            void   resetAllPrevLogLikelihood();

            double calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) const;
            double calcPrevLogCoalLike();
            void   resetPrevLogCoalLike();
            
            double getMaxGeneTreeHeight() const;
                                    
#if defined(LAZY_COPYING)
#else
            double calcTotalCoalRate(unsigned locus) const;
            pair<Node *, Node *> chooseNodesToJoin(const Node::ptr_vect_t & node_ptr_vect) const;
            void calcProbCoalescenceWithinSpecies(vector<double> & probs, double total_rate);
#endif

            void   advanceAllLineagesBy(double dt);
            double findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
            double calcMaxGeneForestHeight() const;

            unsigned proposeSpeciation(unsigned step, Lot::SharedPtr lot);
            void proposeCoalescence(unsigned step, unsigned locus, unsigned rnseed, bool rebuild_species_tree);

            bool advanceByOneCoalescence(unsigned step, unsigned locus, bool first_attempt);

#if defined(LAZY_COPYING)
#else
            double priorPrior(unsigned step, unsigned locus, double total_rate);
#endif
            
            double getLogWeight() const {return _log_weight;}

#if defined(USE_HEATING)
            double getPrevLogWeight() const {return _prev_log_weight;}
#endif
            
            //Lot::SharedPtr getLot() const {return _lot;}
            //void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            
            void refreshHeightsInternalsPreorders();

            void                  setSpeciesTree(string species_tree_newick);
            void                  copySpeciesForestFrom(SpeciesForest & sf);
            SpeciesForest       & getSpeciesForest();
            const SpeciesForest & getSpeciesForestConst() const;
            
#if defined(LAZY_COPYING)
            vector<GeneForest::SharedPtr> & getGeneForestPtrs();
            const vector<GeneForest::SharedPtr> & getGeneForestPtrsConst() const;
            GeneForest::SharedPtr getGeneForestPtr(unsigned locus);
            const GeneForest::SharedPtr getGeneForestPtrConst(unsigned locus) const;
            void finalizeLatestJoin(int locus, unsigned index, map<const void *, list<unsigned> > & nonzero_map);
#else
            vector<GeneForest>       & getGeneForests();
            const vector<GeneForest> & getGeneForestsConst() const;
            GeneForest &          getGeneForest(unsigned gene);
            const GeneForest &    getGeneForestConst(unsigned gene) const;
#endif
            
            void setGeneTrees(vector<GeneForest> & gtvect);
            pair<double,double> chooseTruncatedSpeciesForestIncrement(double truncate_at, Lot::SharedPtr lot);
            
            unsigned debugCountNumCoalEvents() const;
            void debugCheckPartials() const;
            void debugShowAllGeneForests() const;
            void debugCheckAllPrevSpeciesStacksEmpty() const;
            
            void copyParticleFrom(const Particle & other);
            void operator=(const Particle & other);
            
            void setLastProposedGene(unsigned g);
            unsigned getLastProposedGene();
            
            void setLastProposedSpecies(G::species_t s);
            G::species_t getLastProposedSpecies();
            
            void setLastProposedFirstIndex(unsigned f);
            unsigned getLastProposedFirstIndex();
            
            void setLastProposedSecondIndex(unsigned s);
            unsigned getLastProposedSecondIndex();
            
            void stowAllPartials(unsigned locus);
            
            void setRandomNumberSeed(unsigned rnseed);
            
            typedef shared_ptr<Particle> SharedPtr;
                                
        protected:
                
            PartialStore::partial_t pullPartial(unsigned gene);
            void stowPartial(unsigned gene, Node * nd);

            Data::SharedPtr            _data;

            // Even though it is a shared pointer, _lot is a private random number
            // generator not shared with any other particle and has nothing to
            // to do with the global Lot shared_ptr ::rng defined in main.cpp
            mutable Lot::SharedPtr     _lot;
            
            vector<double>             _prev_log_likelihoods;

#if defined(LAZY_COPYING)
            mutable vector<GeneForestExtension> _gene_forest_extensions;
            vector<GeneForest::SharedPtr> _gene_forest_ptrs;
#else
            vector<GeneForest>         _gene_forests;
            mutable vector<Node::species_tuple_t> _species_tuples;
#endif
            SpeciesForest              _species_forest;
            
            vector<Forest::coalinfo_t> _ensemble_coalinfo;  // second level only
            
            mutable double              _log_coal_like;
            mutable double              _prev_log_coal_like;
            mutable double              _log_weight;

#if defined(USE_HEATING)
            mutable double              _prev_log_weight;
#endif
    };

    void Particle::setSpeciesTree(string species_tree_newick) {
        _species_forest.clear();
        _species_forest.buildFromNewick(species_tree_newick);
    }
    
    void Particle::copySpeciesForestFrom(SpeciesForest & sf) {
        _species_forest = sf;
    }

    inline Particle::Particle() {
        _lot.reset(new Lot());
        clear();
    }

    inline Particle::Particle(const Particle & other) {
        _lot.reset(new Lot());
        copyParticleFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _log_coal_like = 0.0;
        _prev_log_coal_like = 0.0;
        _log_weight = 0.0;
#if defined(USE_HEATING)
        _prev_log_weight = 0.0;
#endif
        _prev_log_likelihoods.clear();
#if defined(LAZY_COPYING)
        _gene_forest_ptrs.clear();
#else
        _gene_forests.clear();
#endif
        _species_forest.clear();
    }
    
    inline void Particle::setData(Data::SharedPtr data) {
        _data = data;
    }
            
    inline void Particle::resetSpeciesForest() {
        _species_forest.createTrivialForest();
    }
    
    inline void Particle::resetGeneForests(bool compute_partials) {
        assert(G::_ntaxa > 0);
        assert(G::_nloci > 0);
        _prev_log_likelihoods.clear();
        _prev_log_likelihoods.resize(G::_nloci);
#if defined(LAZY_COPYING)
        _gene_forest_ptrs.clear();
        _gene_forest_ptrs.resize(G::_nloci);
        _gene_forest_extensions.resize(G::_nloci);
        for (unsigned g = 0; g < G::_nloci; g++) {
            _prev_log_likelihoods[g] = 0.0;
            _gene_forest_ptrs[g] = GeneForest::SharedPtr(new GeneForest());
            GeneForest::SharedPtr gfp = _gene_forest_ptrs[g];
            gfp->setData(_data);
            gfp->setRelRate(G::_relrate_for_gene[g]);
            gfp->setGeneIndex(g);
            gfp->createTrivialForest(compute_partials);
        }
#else
        _gene_forests.clear();
        _gene_forests.resize(G::_nloci);
        unsigned g = 0;
        for (auto & gf : _gene_forests) {
            _prev_log_likelihoods[g] = 0.0;
            gf.setData(_data);
            gf.setRelRate(G::_relrate_for_gene[g]);
            gf.setGeneIndex(g++);
            gf.createTrivialForest(compute_partials);
        }
#endif
    }

//    inline void Particle::threadComputePartials(unsigned first, unsigned last) {
//        for (unsigned k = first; k < last; k++) {
//#if defined(LAZY_COPYING)
//            _gene_forest_ptrs[k]->computeAllPartials();
//#else
//            _gene_forests[k].computeAllPartials();
//#endif
//        }
//    }
    
    inline double Particle::calcLogLikelihood() {
        double log_likelihood = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
#if defined(LAZY_COPYING)
            log_likelihood += _gene_forest_ptrs[g]->calcLogLikelihood();
#else
            log_likelihood += _gene_forests[g].calcLogLikelihood();
#endif
        }
        return log_likelihood;
    }
    
    inline void Particle::resetAllPrevLogLikelihood() {
#if defined(LAZY_COPYING)
        for (unsigned g = 0; g < G::_nloci; g++) {
            _prev_log_likelihoods[g] = _gene_forest_ptrs[g]->getLogLikelihood();
        }
#else
        for (unsigned g = 0; g < G::_nloci; g++) {
            _prev_log_likelihoods[g] = _gene_forests[g].getLogLikelihood();
        }
#endif
    }
    
    inline void Particle::resetPrevLogCoalLike() {
        _prev_log_coal_like = _log_coal_like;
    }
    
    inline double Particle::calcTotalPrevLogLikelihood() {
        double prev_log_likelihood = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            prev_log_likelihood += _prev_log_likelihoods[g];
        }
        return prev_log_likelihood;
    }
    
    struct bitless {
        bool operator()(const G::species_t a, const G::species_t b) const {
            bool returned_value = ((a & b) > 0 ? false : a < b);
            return returned_value;
        }
    };
    
    inline void Particle::rebuildCoalInfo() {
        // Rebuild coal info vectors, stripping effect of previous species tree
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            // Each gene forest's _coalinfo vector stores a tuple for each internal node:
            // <1> height
            // <2> gene_index + 1
            // <3> left child species
            // <4> right child species
            gfp->buildCoalInfoVect();
        }
#else
        for (auto & gf : _gene_forests) {
            // Each gene forest's _coalinfo vector stores a tuple for each internal node:
            // <1> height
            // <2> gene_index + 1
            // <3> left child species
            // <4> right child species
            gf.buildCoalInfoVect();
        }
#endif
        _species_forest.buildCoalInfoVect();
    }
    
    inline bool Particle::isEnsembleCoalInfo() {
        return (bool)(_ensemble_coalinfo.size() > 0);
    }
    
    inline void Particle::buildEnsembleCoalInfo() {
        // Store a tuple (height, gene + 1, vector of species)
        // for each join in any gene tree in _ensemble_coalinfo.
        // To get 0-offset gene index,
        // subtract 1 from 2nd element of tuple (if 2nd element
        // is 0, then tuple represents a species tree join)
        _ensemble_coalinfo.clear();

        // Add gene tree joins to _ensemble_coalinfo_vect
#if defined(LAZY_COPYING)
        for (unsigned g = 0; g < G::_nloci; g++) {
            GeneForest::SharedPtr gfp = _gene_forest_ptrs[g];
            gfp->saveCoalInfo(_ensemble_coalinfo);
        }
#else
        for (unsigned g = 0; g < G::_nloci; g++) {
            GeneForest & gf = _gene_forests[g];
            gf.saveCoalInfo(_ensemble_coalinfo);
        }
#endif
        
        // Sort coalinfo_vect from smallest to largest height
        sort(_ensemble_coalinfo.begin(), _ensemble_coalinfo.end());
    }
    
    inline void Particle::recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) const {
    
        // Record gene trees if not already done
        if (_ensemble_coalinfo.empty()) {
            for (unsigned g = 0; g < G::_nloci; g++) {
#if defined(LAZY_COPYING)
                GeneForest::SharedPtr gfp = _gene_forest_ptrs[g];
                assert(gfp->getNumLineages() == 1);
                gfp->saveCoalInfo(coalinfo_vect);
#else
                const GeneForest & gf = _gene_forests[g];
                assert(gf.getNumLineages() == 1);
                gf.saveCoalInfo(coalinfo_vect);
#endif
            }
            sort(coalinfo_vect.begin(), coalinfo_vect.end());
        }
        else {
            // Copy precalculated ensemble coalinfo vector
            coalinfo_vect = _ensemble_coalinfo;
        }

        // Record species tree
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Modify coalescent elements according to species tree
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect, /*capstone*/false);
        
        //BUG lines below fix 2nd level bug 2024-06-19
        // see also Particle::proposeSpeciation
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());
    }
    
    inline double Particle::calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) const {
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
        
#if defined(DEBUG_COALLIKE)
        if (!coalinfo_vect.empty()) {
            output("\nContents of coalinfo_vect:\n");
            output(format("%12s %6s %s\n") % "height" % "locus" % "species", G::LogCateg::DEBUGGING);
            for (auto & cinfo : coalinfo_vect) {
                double      height = get<0>(cinfo);
                unsigned     locus = get<1>(cinfo);
                auto & sppvect = get<2>(cinfo);
                vector<string> s;
                for (auto x : sppvect) {
                    s.push_back(to_string(x));
                }
                output(format("%12.9f %6d %s\n") % height % locus % boost::algorithm::join(s, ","), G::LogCateg::DEBUGGING);
            }
        }
#endif
                
        // Create maps to hold quantities needed by Graham Jones' (2017) formula
        // for coalescent likelihood integrating out theta.
        typedef map<G::species_t, double> dmap;
        
        // n_jb[g][b] holds starting (leaf-end) number of lineages for locus g, branch b
        vector<dmap> n_jb(G::_nloci);
        
        // log_r_b[b] holds sum of log(r_b) for branch b
        dmap log_r_b;

        // q_b[b] holds number of coalescent events for branch b
        dmap q_b;

        // gamma_b[b] holds cumulative gamma_b for branch b
        dmap gamma_b;
        
        // Initialize n_jb for locus 0
        for (auto & nm : G::_taxon_names) {
            // Find index of species corresponding to taxon name nm
            unsigned i = G::_taxon_to_species.at(nm);
            
            // Get species from index
            G::species_t b = (G::species_t)1 << i;
            
            // Increment count of species b in locus 0
            if (n_jb[0].count(b) == 0)
                n_jb[0][b] = 1;
            else
                n_jb[0][b]++;
        }
        
        // Initialize n_jb for other loci by copying n_jb for locus 0
        // Assumes the same number of individuals have been sampled from each species for all loci
        for (unsigned g = 1; g < G::_nloci; g++) {
            n_jb[g] = n_jb[0];
        }
        
        // prev_height[g][b] holds previous height for locus g and branch b
        vector<map<G::species_t, double> > prev_height(G::_nloci);
        
        // Ploidy for each locus (currently all loci assumed to be diploid)
        vector<double> p_j(G::_nloci, 2.0);

#if defined(HACK_FOR_SNAKE_ATP)
        if (G::_hack_atp_index < 0) {
            // Find index of ATP gene
            for (unsigned i = 0; i < G::_nloci; i++) {
                if (G::_gene_names[i] == "ATP") {
                    G::_hack_atp_index = i;
                    break;
                }
            }
        }
        assert(G::_hack_atp_index > -1);
        p_j[G::_hack_atp_index] = 1.0;
#endif
        
        // Vector branches stores branches (including ancestral ones)
        vector<G::species_t> branches;
        for (unsigned i = 0; i < G::_nspecies; i++) {
            G::species_t b = (G::species_t)1 << i;
            branches.push_back(b);
        }
        
        // Vector bavail stores branches that are in the species tree at the current height
        set<G::species_t, bitless> bavail(branches.begin(), branches.end());
        
        if (verbose) {
            output("\nParticle::calcLogCoalescentLikelihood:\n");
            output("Key:\n");
            output("  height:   height above leaf-level of this coalescence or speciation event\n");
            output("  locus:    0 if speciation event, locus index for coalescence events\n");
            output("  spp:      species of the two lineages involved in a coalescence/speciation event\n");
            output("  b:        the species tree branch (i.e. species) to which this event contributes\n");
            output("  n_jbk:    the number of coalescences for branch b\n");
            output("  c_jbk:    the sojourn time prior to coalescence k\n");
            output("  log(r_b): the log of 4/pj, where pj is the ploidy of locus j\n");
            output("  gamma_b:  cumulative 2*n*(n-1)*c_jbk/pj over all loci for branch b\n");
            output("  q_b:      the cumulative number of coalescences over all loci for branch b\n");
            output(format("%12s %12s %25s %12s %12s %12s %12s %12s\n") % "height" % "locus" % "spp" % "b" % "n_jbk" % "c_jbk" % "gamma_b" % "q_b");
        }

        // Walk down coalinfo_vect, accumulating Graham Jones r_b, q_b, and gamma_b
        for (auto & cinfo : coalinfo_vect) {
            double               height = get<0>(cinfo);
            unsigned             locus_plus_one   = get<1>(cinfo);
            int                  locus   = locus_plus_one - 1;
            
            vector<G::species_t> & b_vect = get<2>(cinfo);

            string spp_joined;
            if (verbose) {
                ostringstream oss;
                copy(b_vect.begin(), b_vect.end(), ostream_iterator<G::species_t>(oss, " "));
                string spp_joined = oss.str();
                boost::algorithm::trim_right(spp_joined);
            }
            
            G::species_t banc = 0;
            for (auto b : b_vect)
                banc |= b;
                
            if (locus_plus_one == 0) {
                // Speciation event
                
                // If banc represents a new branch in the species tree,
                // add it to the vector of branches
                if (find(branches.begin(), branches.end(), banc) == branches.end()) {
                    branches.push_back(banc);
                }
                    
                for (unsigned g = 0; g < G::_nloci; g++) {
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
                    output(format("%12.9f %12d %25s %12s %12s %12s %12s %12s\n") % height % 0 % spp_joined % banc % "-" % "-" % "-" % "-");
                }
            }
            else {
                // Coalescent event

                // Determine b, the edge of the species tree we're dealing with
                // Note that b may be an edge that is not yet in the species tree (i.e. the root edge)
                auto it = find_if(bavail.begin(), bavail.end(), [banc](G::species_t v){return (banc & v) > 0;});
                assert(it != bavail.end());
                G::species_t b = *it;
                
                if (prev_height[locus].count(b) == 0) {
                    prev_height[locus][b] = 0.0;
                }
                double height0 = prev_height[locus].at(b);
                
                // Coalescent event: update Jones quantities
                log_r_b[b] += log(4/p_j[locus]);
                q_b[b]++;
                double n_jbk = n_jb[locus][b];
                double c_jbk = height - height0;
                gamma_b[b] += 2.0*n_jbk*(n_jbk-1)*c_jbk/p_j[locus];

                prev_height[locus][b] = height;
                n_jb[locus][b]--;

                if (verbose) {
                    output(format("%12.9f %12d %25s %12d %12d %12.9f %12.9f %12d\n") % height % locus_plus_one % spp_joined % b % n_jbk % c_jbk % gamma_b[b] % q_b[b]);
                }
            }
        }

        if (verbose) {
            output(format("\n%12s %12s %12s %12s %12s %12s\n") % "b" % "q_b" % "log(r_b)" % "gamma_b" % "theta_b" % "logL");
        }
    
        double log_likelihood = 0.0;
        
        if (integrate_out_thetas) {
            double alpha = 2.0; // G::_invgamma_shape;
            double beta = G::_theta;
            assert(beta > 0.0);
            double B = (double)branches.size();
            log_likelihood  = B*alpha*log(beta);
            log_likelihood -= B*boost::math::lgamma(alpha);
            for (auto b : branches) {
                log_likelihood += log_r_b[b];
                log_likelihood -= (alpha + q_b[b])*log(beta + gamma_b[b]);
                log_likelihood += boost::math::lgamma(alpha + q_b[b]);
                
                if (verbose) {
                    output(format("%12d %12d %12.9f %12.9f %12.9f %12.9f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % beta % log_likelihood);
                }
            }

            if (verbose) {
                output(format("\nalpha = %.9f\n") % alpha);
                output(format("beta  = %.9f\n") % beta);
                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood);
            }
        }
        else {
            double theta_b = G::_theta;
            assert(theta_b > 0.0);
            double sum_log_rb = 0.0;
            double sum_gamma_b = 0.0;
            unsigned sum_qb = 0;
            for (auto b : branches) {
                sum_qb += q_b[b];
                sum_log_rb += log_r_b[b];
                sum_gamma_b += gamma_b[b];
                log_likelihood += log_r_b[b] - gamma_b[b]/theta_b - q_b[b]*log(theta_b);
                
                if (verbose) {
                    output(format("%12d %12d %12.9f %12.9f %12.9f %12.9f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % theta_b % log_likelihood, 1);
                }
            }
            
            if (verbose) {
                output(format("\nsum log(r_b) = %.9f\n") % sum_log_rb);
                output(format("sum q_b      = %d\n") % sum_qb);
                output(format("sum gamma_b  = %d\n") % sum_gamma_b);
                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood);
            }
        }

        _log_coal_like = log_likelihood;
        return log_likelihood;
    }
    
    inline void Particle::computeAllPartials() {
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            gfp->computeAllPartials();
        }
#else
        for (auto & gf : _gene_forests) {
            gf.computeAllPartials();
        }
#endif
    }
        
    inline PartialStore::partial_t Particle::pullPartial(unsigned gene) {
#if defined(USING_MULTITHREADING)
        lock_guard<mutex> guard(G::_mutex);
#endif
#if defined(LAZY_COPYING)
        assert(gene < _gene_forest_ptrs.size());
#else
        assert(gene < _gene_forests.size());
#endif

        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
        ptr = ps.getPartial(gene);
        return ptr;
    }

    inline void Particle::stowPartial(unsigned gene, Node * nd) {
#if defined(USING_MULTITHREADING)
        lock_guard<mutex> guard(G::_mutex);
#endif
#if defined(LAZY_COPYING)
        assert(gene < _gene_forest_ptrs.size());
#else
        assert(gene < _gene_forests.size());
#endif
        assert(nd);
        assert(nd->_partial);
        ps.putPartial(gene, nd->_partial);
        
        // Decrement shared pointer reference count
        nd->_partial.reset();
    }
    
    inline void Particle::setRandomNumberSeed(unsigned rnseed) {
        _lot->setSeed(rnseed);
    }

    inline void Particle::stowAllPartials(unsigned locus) {
#if defined(LAZY_COPYING)
        assert(locus < _gene_forest_ptrs.size());
        _gene_forest_ptrs[locus]->stowAllPartials();
#else
        assert(locus < _gene_forests.size());
        _gene_forests[locus].stowAllPartials();
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
        
        assert(min_height < G::_infinity);
        return min_height;
    }
    
#if defined(DEBUG_COALLIKE)
    inline unsigned Particle::proposeSpeciation(unsigned step, Lot::SharedPtr lot) {
        // This function is only used for proposing speciation events when there are
        // complete gene trees available. It thus draws increments from a truncated
        // exponential distribution where the trunction point is the next height at
        // which at least one coalescent event combines lineages from two different
        // species.
        unsigned num_species_tree_lineages = 0;
        
        // Store tuple (height, 0, vector of species) for
        // each join in the current species forest.
        // Do not cap with ancestral species at this point.
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Store a tuple (height, gene + 1, vector of species)
        // for each join in any gene tree. To get 0-offset gene index,
        // subtract 1 from 2nd element of tuple (if 2nd element
        // is 0, then tuple represents a species tree join)
        vector<Forest::coalinfo_t> coalinfo_vect;
        
        // Add gene tree joins to coalinfo_vect
        // Just need coalescent events at this point in order to choose
        // limits for species tree increments
#if defined(LAZY_COPYING)
        for (unsigned g = 0; g < G::_nloci; g++) {
            GeneForest::SharedPtr gfp = _gene_forest_ptrs[g];
            gfp->saveCoalInfo(coalinfo_vect);
        }
#else
        for (unsigned g = 0; g < G::_nloci; g++) {
            GeneForest & gf = _gene_forests[g];
            gf.saveCoalInfo(coalinfo_vect);
        }
#endif
        
        // Sort coalinfo_vect from smallest to largest height
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Get maximum height of any gene tree
        double max_height = get<0>((*coalinfo_vect.rbegin()));
        
        if (step > 0) {
            output("\nSpecies tree before creating speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);

            // Create speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(lot, left_spp, right_spp, anc_spp);
            
            output("\nSpecies tree after creating speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);

            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _species_forest.buildCoalInfoVect();
            _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());

            // Show coalinfo_vect before fixing up
            _species_forest.debugShowCoalInfo("sppinfo_vect after speciation event", sppinfo_vect, /*fn*/"");
            _species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
                
            // Adjust elements of coalinfo_vect affected by species tree joins
            _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

            // Show coalinfo_vect before fixing up
            _species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
        }
                
        // Draw a speciation increment dt.
        double forest_height = _species_forest.getHeight();
        //speclog_element._height = forest_height;
        double h = findHeightNextCoalescentEvent(forest_height, coalinfo_vect);
        assert(h <= max_height);
        
        pair<double,double> tmp = chooseTruncatedSpeciesForestIncrement(h);
        double log_weight_factor = tmp.first;
        //double increment = tmp.second;
        
        num_species_tree_lineages = _species_forest.getNumLineages();
        if (num_species_tree_lineages == 2) {
            output("\nSpecies tree before creating FINAL speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);
            // Create final speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(lot, left_spp, right_spp, anc_spp);
            
            output("\nSpecies tree after creating FINAL speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);
        }
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear();
        _species_forest.buildCoalInfoVect();
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Show coalinfo_vect before fixing up
        _species_forest.debugShowCoalInfo("sppinfo_vect", sppinfo_vect, /*fn*/"");
        _species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
                
        // Adjust elements of coalinfo_vect affected by species tree joins
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

        // Show coalinfo_vect after fixing up
        _species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");

        // Add speciations into coalinfo_vect
        //BUG lines below fix 2nd level bug 2024-06-19
        // see also Particle::recordAllForests
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Compute coalescent likelihood and log weight
        //double prev_log_coallike = _prev_log_coallike;
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/true);

#if defined(USE_HEATING)
        _prev_log_weight = _log_weight;
#endif
        _log_weight = _log_coal_like - _prev_log_coal_like + log_weight_factor;

        resetPrevLogCoalLike();
        
        return num_species_tree_lineages;
    }
#else
    // DEBUG_COALLIKE not defined
    inline unsigned Particle::proposeSpeciation(unsigned step, Lot::SharedPtr lot) {
        // This function is only used for proposing speciation events when there are
        // complete gene trees available. It thus draws increments from a truncated
        // exponential distribution where the trunction point is the next height at
        // which at least one coalescent event combines lineages from two different
        // species.
        unsigned num_species_tree_lineages = 0;
        
        // Store tuple (height, 0, vector of species) for
        // each join in the current species forest.
        // Do not cap with ancestral species at this point.
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Make a copy of _ensemble_coalinfo;
        vector<Forest::coalinfo_t> coalinfo_vect = _ensemble_coalinfo;

        // Get maximum height of any gene tree
        double max_height = get<0>((*coalinfo_vect.rbegin()));
        
        if (step > 0) {
            // Create speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(lot, left_spp, right_spp, anc_spp);
            
            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _species_forest.buildCoalInfoVect();
            _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());

            // Adjust elements of coalinfo_vect affected by species tree joins
            _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect, /*capstone*/false);
        }
                
        // Draw a speciation increment dt.
        double forest_height = _species_forest.getHeight();
        double h = findHeightNextCoalescentEvent(forest_height, coalinfo_vect);
        assert(h <= max_height);
        
        pair<double,double> tmp = chooseTruncatedSpeciesForestIncrement(h, lot);
        double log_weight_factor = tmp.first;
        
        num_species_tree_lineages = _species_forest.getNumLineages();
        if (num_species_tree_lineages == 2) {
            // Create final speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(lot, left_spp, right_spp, anc_spp);
        }
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear();
        _species_forest.buildCoalInfoVect();
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
                
        // Adjust elements of coalinfo_vect affected by species tree joins
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect, /*capstone*/true);

        // Add speciations into coalinfo_vect
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Compute coalescent likelihood and log weight
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
#if defined(USE_HEATING)
        _prev_log_weight = _log_weight;
#endif
        _log_weight = _log_coal_like - _prev_log_coal_like + log_weight_factor;

        resetPrevLogCoalLike();
        
        return num_species_tree_lineages;
    }
#endif // DEBUG_COALLIKE not defined

    inline double Particle::calcMaxGeneForestHeight() const {
        double max_gene_forest_height = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
#if defined(LAZY_COPYING)
            double h = _gene_forest_ptrs[g]->getHeight();
#else
            double h = _gene_forests[g].getHeight();
#endif
            if (h > max_gene_forest_height)
                max_gene_forest_height = h;
        }
        return max_gene_forest_height;
    }
    
    inline void Particle::proposeCoalescence(unsigned step, unsigned locus, unsigned rnseed, bool rebuild_species_tree) {
        _lot->setSeed(rnseed);
        
        // Rebuild species tree starting from the height of the
        // tallest gene forest over all loci
        if (rebuild_species_tree) {
            double max_height = calcMaxGeneForestHeight();
            _species_forest.rebuildStartingFromHeight(_lot, max_height);
        }

        // Record previous log likelihood
#if defined(LAZY_COPYING)
        _prev_log_likelihoods[locus] = _gene_forest_ptrs[locus]->getLogLikelihood();
#else
        _prev_log_likelihoods[locus] = _gene_forests[locus].getLogLikelihood();
#endif

        // Advance locus gene forest by one coalescent event.
        // Function advanceByOneCoalescence returns
        // - true if a coalescence event was created
        // - false if the next speciation event occurred first
        bool done = advanceByOneCoalescence(step, locus, /*first_attempt*/true);
        while (!done) {
            done = advanceByOneCoalescence(step, locus, /*first_attempt*/false);
        }
    }
    
#if defined(LAZY_COPYING)
#else
    inline double Particle::calcTotalCoalRate(unsigned locus) const {
        double total_rate = 0.0;
        total_rate += _gene_forests[locus].calcTotalRate(_species_tuples);
        return total_rate;
    }
#endif
    
    inline void Particle::advanceAllLineagesBy(double dt) {
        _species_forest.advanceAllLineagesBy(dt);
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            gfp->advanceAllLineagesBy(dt);
        }
#else
        for (auto & gene_forest : _gene_forests) {
            gene_forest.advanceAllLineagesBy(dt);
        }
#endif
    }
    
#if defined(LAZY_COPYING)
#else
    inline void Particle::calcProbCoalescenceWithinSpecies(vector<double> & probs, double total_rate) {
        unsigned nspecies = (unsigned)_species_tuples.size();
        assert(nspecies > 0);
        probs.clear();
        probs.resize(nspecies);
        
        // The probability of coalescence in one particular species i is
        //   probs[i] =                 r_i / total_rate
        //            = (n_i*(n_i-1)/theta) / (sum_i r_i)
        // where n_i = number of lineages in species i.
        transform(
            _species_tuples.begin(),
            _species_tuples.end(),
            probs.begin(),
            [total_rate](Node::species_tuple_t & spp_tuple) {
                double n = (double)get<0>(spp_tuple);
                return n*(n - 1.0)/(G::_theta*total_rate);
            }
        );
    }
#endif
    
#if defined(LAZY_COPYING)
#else
    inline pair<Node *, Node *> Particle::chooseNodesToJoin(const Node::ptr_vect_t & node_ptr_vect) const {
        unsigned n = (unsigned)node_ptr_vect.size();
        assert(n > 1);
        
        // Choose a random pair of lineages to join
        pair<unsigned,unsigned> chosen_pair = lot->nchoose2(n);
        unsigned i = chosen_pair.first;
        unsigned j = chosen_pair.second;
        
        // Record chosen pair of lineages
        Node * first  = node_ptr_vect[i];
        Node * second = node_ptr_vect[j];
        return make_pair(first, second);
    }
#endif
    
#if defined(LAZY_COPYING)
#else
    inline double Particle::priorPrior(unsigned step, unsigned locus, double total_rate) {
        double log_weight = 0.0;
        
        // _species_tuples has already been filled with entries specific to locus
        // Each tuple entry stores:
        //  0. number of lineages (unsigned)
        //  1. gene index (unsigned)
        //  2. species within gene (G::species_t)
        //  3. lineages (vector<Node *>)
        
        // Get the gene forest extension involved in this proposed coalescence
        GeneForest & gf = _gene_forests[locus];
        
        // Calculate probabilities of coalescence within each species
        vector<double> probs;
        calcProbCoalescenceWithinSpecies(probs, total_rate);
                
        // Choose the species in which the coalescence will happen
        unsigned which = G::multinomialDraw(lot, probs);
        assert(which < probs.size());
        G::species_t spp = get<2>(_species_tuples[which]);
                    
        // Choose the two nodes to join
        pair<Node *, Node *> chosen = chooseNodesToJoin(gf._lineages_within_species.at(spp));
        assert(chosen.first->getSpecies() == spp);
        assert(chosen.second->getSpecies() == spp);

        // Pull an ancestral node
        Node * anc_node = gf.pullNode();
        anc_node->_edge_length = 0.0;
        anc_node->setSpecies(spp);
        if (!G::_simulating) {
            anc_node->_partial = pullPartial(locus);
        }
        
        // Set anc_node split to union of the two child splits
        anc_node->_split.resize(G::_ntaxa);
        anc_node->_split += chosen.first->_split;
        anc_node->_split += chosen.second->_split;

        // Set ancestral node height
        double h1 = chosen.first->_height + chosen.first->_edge_length;
        double h2 = chosen.second->_height + chosen.second->_edge_length;
        anc_node->_height = 0.5*(h1 + h2);

        // Perform the join
        gf.joinLineagePair(anc_node, chosen.first, chosen.second);
        gf.removeTwoAddOne(gf._lineages, chosen.first, chosen.second, anc_node);
        gf.refreshAllPreorders();

        if (!G::_simulating) {
            // Compute partial likelihood array of ancestral node
            assert(_data);
            assert(anc_node->_left_child);
            assert(anc_node->_left_child->_right_sib);
            assert(anc_node->_partial);
            
            log_weight = gf.calcPartialArray(anc_node);

            assert(!isnan(log_weight));
            assert(!isinf(log_weight));

#if defined(DEBUG_CHECK_WEIGHTS)
            // Check log weight
            double logL = gf.calcLogLikelihood();
            double prevLogL = _prev_log_likelihoods[locus];
            double check_log_weight = logL - prevLogL;
            assert(fabs(check_log_weight - log_weight) < G::_small_enough);
#endif
        }
                
        return log_weight;
    }
#endif
    
    inline bool Particle::advanceByOneCoalescence(unsigned step, unsigned locus, bool first_attempt) {
#if defined(LAZY_COPYING)
        if (first_attempt) {
            // Create temporary gene forest extending existing forest
            // without touching existing forest (which may be used
            // by many particles)
            assert(_gene_forest_extensions.size() == G::_nloci);
            _gene_forest_extensions[locus].dock(_gene_forest_ptrs[locus], pullPartial(locus), _lot);
        }
        GeneForestExtension & gfx = _gene_forest_extensions[locus];
#else
        GeneForest & gf = _gene_forests[locus];

        // Clear species_tuples vector. Each tuple entry stores:
        //  0. number of lineages
        //  1. gene index
        //  2. species within gene
        //  3. vector<Node *> holding lineage roots for gene/species combination
        _species_tuples.clear();
#endif
        
#if defined(LAZY_COPYING)
        // Get current height of this locus' gene forest extension
        double gene_forest_height = gfx.getHeight();
#else
        // Get current height of this locus' gene forest
        double gene_forest_height = gf.getHeight();
#endif

        // Find height of next speciation event
        auto speciation_tuple = _species_forest.findNextSpeciationEvent(gene_forest_height);
        double speciation_height = get<0>(speciation_tuple);
        double speciation_delta  = (speciation_height == G::_infinity)
                                    ? G::_infinity
                                    : speciation_height - gene_forest_height;
        G::species_t left_spp    = get<1>(speciation_tuple);
        G::species_t right_spp   = get<2>(speciation_tuple);

#if defined(LAZY_COPYING)
#else
        G::species_t anc_spp     = get<3>(speciation_tuple);
#endif
                        
        // Visit each species within the specified locus, computing
        // the total coalescence rate (total_rate).
#if defined(LAZY_COPYING)
        double total_rate = gfx.calcTotalRate(G::_theta);
#else
        double total_rate = gf.calcTotalRate(_species_tuples);
#endif
        
        // Draw coalescence increment delta ~ Exponential(total_rate)
        double delta = total_rate > 0.0
                       ? -log(1.0 - _lot->uniform())/total_rate
                       : G::_infinity;
        
        // Sanity check: delta and speciation_delta cannot both
        // be equal to infinity
        assert(!(isinf(delta) && isinf(speciation_delta)));

        bool done = false;
        if (delta > speciation_delta) {
            // Speciation event comes before coalescence event
            
#if defined(LAZY_COPYING)
            // Bring gene forest extension up to height of speciation event
            gfx.advanceProposedDeltaBy(speciation_delta);
                                                
            // Advise gene forest extension of the change in the species tree
            gfx.mergeSpecies(left_spp, right_spp);
#else
            // Bring gene forest up to height of speciation event
            gf.advanceAllLineagesBy(speciation_delta);
                                                
            // Advise gene tree of the change in the species tree
            // Nodes that are reassigned save their previous state
            // to allow reversion
            gf.mergeSpecies(speciation_height, left_spp, right_spp, anc_spp);
#endif
        }
        else {
            // Coalescence event comes before speciation event
            done = true;

#if defined(LAZY_COPYING)
            // Advance gene forest height by an amount delta
            gfx.advanceProposedDeltaBy(delta);
            gfx.coalesce(total_rate);
            _log_weight = gfx.getLogWeight();
#else
            // Advance gene forest height by an amount delta
            gf.advanceAllLineagesBy(delta);

            // Perform a random join and compute the log weight
            _log_weight = priorPrior(step, locus, total_rate);
#endif

#if defined(USE_HEATING)
            _prev_log_weight = _log_weight;
#endif
        }
        
        return done;
    }

    inline void Particle::refreshHeightsInternalsPreorders() {
        _species_forest.heightsInternalsPreorders();
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            gfp->heightsInternalsPreorders();
        }
#else
        for (auto & gf : _gene_forests) {
            gf.heightsInternalsPreorders();
        }
#endif
    }
    
    inline SpeciesForest & Particle::getSpeciesForest() {
        return _species_forest;
    }
        
    inline const SpeciesForest & Particle::getSpeciesForestConst() const {
        return _species_forest;
    }

#if defined(LAZY_COPYING)
    inline void Particle::finalizeLatestJoin(int locus, unsigned index, map<const void *, list<unsigned> > & nonzero_map) {
        // Makes join closest to leaf-level in _gene_forest_extensions[locus]
        // permanent, then undocks _gene_forest_extensions[locus]
        
        // Get reference to gene forest extension for this locus
        GeneForestExtension & gfx = _gene_forest_extensions[locus];
        
        // Get pointer to gene forest for this locus
        GeneForest::SharedPtr gfp = _gene_forest_ptrs[locus];
        
        // If we are not finalizing the last particle for this
        // gene forest object, make a copy that can be modified
        // without affecting other surviving particles
        unsigned nz = (unsigned)nonzero_map[gfp.get()].size();
        if (nz > 1) {
            // Remove the element corresponding to index
            list<unsigned> & v = nonzero_map[gfp.get()];
            auto it = find(v.begin(), v.end(), index);
            assert(it != v.end());
            v.erase(it);
            
            // Make a copy of the object pointed to by gfp
            GeneForest::SharedPtr gfcpy = GeneForest::SharedPtr(new GeneForest());
            *gfcpy = *gfp;
            _gene_forest_ptrs[locus] = gfcpy;
            
            // Let gpf point to the copy
            gfp = gfcpy;
        }
        
        // Copy log likelihood
        gfp->setLogLikelihood(_prev_log_likelihoods[locus] + gfx.getLogWeight());
                        
        // Get splits for children of _proposed_anc
        const Node * anc = gfx.getProposedAnc();
        assert(anc);
        const Node * lchild = gfx.getProposedLChild();
        assert(lchild);
        const Node * rchild = gfx.getProposedRChild();
        assert(rchild);
        Split lsplit = lchild->_split;
        Split rsplit = rchild->_split;
        
        assert(anc->_split.isEquivalent(lsplit + rsplit));
        
        // Recreate extension's join in the actual gene forest
        double incr = gfx.getProposedDelta();
        assert(incr > 0.0);
        
        gfp->addIncrAndJoin(incr, lsplit, rsplit, gfx);
        
        //gfp->debugCheckBleedingEdge("after finalizeLatestJoin", _gene_forest_extensions[locus].getHeight());
        
        // Can now get rid of extension
        _gene_forest_extensions[locus].undock();
    }
#endif

#if defined(LAZY_COPYING)
    inline GeneForest::SharedPtr Particle::getGeneForestPtr(unsigned locus) {
        assert(locus < G::_nloci);
        assert(_gene_forest_ptrs[locus]);
        return _gene_forest_ptrs[locus];
    }
#else
    inline GeneForest & Particle::getGeneForest(unsigned locus) {
        assert(locus < G::_nloci);
        return _gene_forests[locus];
    }
#endif

#if defined(LAZY_COPYING)
    inline vector<GeneForest::SharedPtr> & Particle::getGeneForestPtrs() {
        return _gene_forest_ptrs;
    }
#else
    inline vector<GeneForest> & Particle::getGeneForests() {
        return _gene_forests;
    }
#endif
    
#if defined(LAZY_COPYING)
    inline const vector<GeneForest::SharedPtr> & Particle::getGeneForestPtrsConst() const {
        return _gene_forest_ptrs;
    }
#else
    inline const vector<GeneForest> & Particle::getGeneForestsConst() const {
        return _gene_forests;
    }
#endif
    
#if defined(LAZY_COPYING)
    inline const GeneForest::SharedPtr Particle::getGeneForestPtrConst(unsigned locus) const {
        assert(locus < G::_nloci);
        return _gene_forest_ptrs[locus];
    }
#else
    inline const GeneForest & Particle::getGeneForestConst(unsigned locus) const {
        assert(locus < G::_nloci);
        return _gene_forests[locus];
    }
#endif

    inline unsigned Particle::debugCountNumCoalEvents() const {
        // Returns number of coalescent events over all gene trees
        unsigned total = 0;
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            unsigned n = gfp->getNumLineages();
            total += G::_ntaxa - n;
        }
#else
        for (auto & gf : _gene_forests) {
            unsigned n = gf.getNumLineages();
            total += G::_ntaxa - n;
        }
#endif
        return total;
    }
            
    inline void Particle::debugCheckPartials() const {
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            gfp->debugCheckPartials();
        }
#else
        for (auto & gf : _gene_forests) {
            gf.debugCheckPartials();
        }
#endif
    }
    
    inline void Particle::debugShowAllGeneForests() const {
#if defined(LAZY_COPYING)
        for (auto gfp : _gene_forest_ptrs) {
            cerr << gfp->makeNewick(0, true, false) << endl;
        }
#else
        for (auto & gf : _gene_forests) {
            cerr << gf.makeNewick(0, true, false) << endl;
        }
#endif
    }
    
    inline void Particle::copyParticleFrom(const Particle & other) {
        clear();

        // Performs a deep copy of other to this particle
        _data = other._data;
        _prev_log_likelihoods = other._prev_log_likelihoods;
        _ensemble_coalinfo = other._ensemble_coalinfo;
        
#if defined(LAZY_COPYING)
        // Ensure that _gene_forest_extensions is allocated
        if (_gene_forest_extensions.size() == 0) {
            _gene_forest_extensions.resize(G::_nloci);
        }
        else {
            assert(_gene_forest_extensions.size() == G::_nloci);
            
            // Undock all gene forest extensions
            for_each(_gene_forest_extensions.begin(), _gene_forest_extensions.end(),
                [](GeneForestExtension & p){p.undock();});
        }

        // Copy gene forest pointers only
        assert(G::_nloci == other._gene_forest_ptrs.size());
        _gene_forest_ptrs.resize(G::_nloci);
        for (unsigned i = 0; i < G::_nloci; i++) {
            _gene_forest_ptrs[i] = other._gene_forest_ptrs[i];
        }
#else
        assert(G::_nloci == other._gene_forests.size());
        _gene_forests.resize(G::_nloci);
        for (unsigned i = 0; i < G::_nloci; i++) {
            _gene_forests[i] = other._gene_forests[i];
        }
#endif
        
        // Copy species forest
        _species_forest     = other._species_forest;
        
        _log_coal_like      = other._log_coal_like;
        _prev_log_coal_like = other._prev_log_coal_like;
        _log_weight         = other._log_weight;
#if defined(USE_HEATING)
        _prev_log_weight    = other._prev_log_weight;
#endif
        
#if defined(LAZY_COPYING)
#else
        _species_tuples     = other._species_tuples;
#endif
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
    inline pair<double,double> Particle::chooseTruncatedSpeciesForestIncrement(double truncate_at, Lot::SharedPtr lot) {
        double upper_bound = truncate_at - _species_forest.getHeight();
        auto incr_rate_cum = _species_forest.drawTruncatedIncrement(lot, upper_bound);
        double incr = get<0>(incr_rate_cum);
        double rate = get<1>(incr_rate_cum);
        assert(rate > 0.0);
        _species_forest.advanceAllLineagesBy(incr);
        double cum = get<2>(incr_rate_cum);
        assert(cum > 0.0);
        double log_factor = log(cum);
        return make_pair(log_factor,incr);
    }
    
    inline double Particle::getMaxGeneTreeHeight() const {
        double maxh = 0.0;
        for (unsigned i = 0; i < G::_nloci; i++) {
#if defined(LAZY_COPYING)
            double h = _gene_forest_ptrs[i]->getHeight();
#else
            double h = _gene_forests[i].getHeight();
#endif
            if (h > maxh)
                maxh = h;
        }
        return maxh;
    }
}
