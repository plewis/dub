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
                        
            double setThetas();
            
            void resetSpeciesForest();
            void resetGeneForests(bool compute_partials);

//#if defined(UPGMA_WEIGHTS)
//            void resetDistanceMatrix();
//#endif
            
            void threadComputePartials(unsigned first, unsigned last);
            void computeAllPartials();
            
            void rebuildCoalInfo();
            void recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) const;
            
            double calcLogLikelihood();
            double calcPrevLogLikelihood();
            void   resetAllPrevLogLikelihood();

            double calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) const;
            double calcPrevLogCoalLike();
            void   resetPrevLogCoalLike();
            
            double getMaxGeneTreeHeight() const;
                                    
            void   advanceAllLineagesBy(double dt);
            double findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
            double calcTotalCoalRate(unsigned locus);

            unsigned proposeSpeciation(unsigned step);
            void proposeCoalescence(unsigned step, unsigned locus);

            bool advanceByOneCoalescence(unsigned step, unsigned locus);

            double priorPrior(unsigned step, unsigned locus, double total_rate);
            
            double getLogWeight() const {return _log_weight;}
            double getPrevLogWeight() const {return _prev_log_weight;}
            
            //Lot::SharedPtr getLot() const {return _lot;}
            //void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            
            void refreshHeightsInternalsPreorders();

            void                  setSpeciesTree(string species_tree_newick);
            void                  copySpeciesForestFrom(SpeciesForest & sf);
            SpeciesForest       & getSpeciesForest();
            const SpeciesForest & getSpeciesForestConst() const;

            vector<GeneForest>       & getGeneForests();
            const vector<GeneForest> & getGeneForestsConst() const;
            
            GeneForest       & getGeneForest(unsigned gene);
            const GeneForest & getGeneForestConst(unsigned gene) const;
            
            void setGeneTrees(vector<GeneForest> & gtvect);
            pair<double,double> chooseTruncatedSpeciesForestIncrement(double truncate_at);
            
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
            
            typedef shared_ptr<Particle> SharedPtr;
                                
        protected:
                
            PartialStore::partial_t pullPartial(unsigned gene);
            void stowPartial(unsigned gene, Node * nd);

            Data::SharedPtr            _data;

            vector<GeneForest>         _gene_forests;
            SpeciesForest              _species_forest;
            
//#if defined(UPGMA_WEIGHTS)
//            vector<vector<double> >     _dmatrix;
//            vector<Split>               _dmatrix_rows;
//#endif
            
            vector<Node::species_tuple_t> _species_tuples;
            
            mutable double              _log_coal_like;
            mutable double              _prev_log_coal_like;

            mutable double              _log_weight;
            mutable double              _prev_log_weight;
            
            // Even though a shared pointer, _lot is a private random number
            // generator not shared with any other particle and has nothing to
            // to do with the global Lot shared_ptr ::rng defined in main.cpp
            //mutable Lot::SharedPtr  _lot;
    };

    void Particle::setSpeciesTree(string species_tree_newick) {
        _species_forest.clear();
        _species_forest.buildFromNewick(species_tree_newick);
    }
    
    void Particle::copySpeciesForestFrom(SpeciesForest & sf) {
        _species_forest = sf;
    }

    inline Particle::Particle() {
        //_lot.reset(new Lot());
        clear();
    }

    inline Particle::Particle(const Particle & other) {
        //_lot.reset(new Lot());
        copyParticleFrom(other);
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _log_coal_like = 0.0;
        _prev_log_coal_like = 0.0;
        _log_weight = 0.0;
        _prev_log_weight = 0.0;
        _gene_forests.clear();
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
        _gene_forests.clear();
        _gene_forests.resize(G::_nloci);
        unsigned g = 0;
        for (auto & gf : _gene_forests) {
            gf.setData(_data);
            gf.setRelRate(G::_relrate_for_gene[g]);
            gf.setGeneIndex(g++);
            gf.createTrivialForest(compute_partials);
            gf.setParticle(this);
        }
    }

    inline void Particle::threadComputePartials(unsigned first, unsigned last) {
        for (unsigned k = first; k < last; k++) {
            _gene_forests[k].computeAllPartials();
        }
    }
    
    inline double Particle::calcLogLikelihood() {
        double log_likelihood = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
#if defined(UPGMA_WEIGHTS)
            log_likelihood += _gene_forests[g].calcLogLikelihood(_species_forest);
#else
            log_likelihood += _gene_forests[g].calcLogLikelihood();
#endif
        }
        return log_likelihood;
    }
    
    inline void Particle::resetAllPrevLogLikelihood() {
        for (unsigned g = 0; g < G::_nloci; g++) {
            _gene_forests[g].resetPrevLogLikelihood();
        }
    }
    
    inline void Particle::resetPrevLogCoalLike() {
        _prev_log_coal_like = _log_coal_like;
    }
    
    inline double Particle::calcPrevLogLikelihood() {
        double prev_log_likelihood = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            prev_log_likelihood += _gene_forests[g].getPrevLogLikelihood();
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
        for (auto & gf : _gene_forests) {
            // Each gene forest's _coalinfo vector stores a tuple for each internal node:
            // <1> height
            // <2> gene_index + 1
            // <3> left child species
            // <4> right child species
            gf.buildCoalInfoVect();
        }
        _species_forest.buildCoalInfoVect();
    }
    
    inline void Particle::recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) const {
        // Record gene trees
        for (unsigned g = 0; g < G::_nloci; g++) {
            const GeneForest & gf = _gene_forests[g];
            assert(gf.getNumLineages() == 1);
            gf.saveCoalInfo(coalinfo_vect);
        }
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Record species tree
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // Modify coalescent elements according to species tree
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);
        
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
        for (auto & gf : _gene_forests) {
            gf.computeAllPartials();
        }
    }
        
    inline PartialStore::partial_t Particle::pullPartial(unsigned gene) {
        assert(gene < _gene_forests.size());

        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
        ptr = ps.getPartial(gene);
        return ptr;
    }

    inline void Particle::stowPartial(unsigned gene, Node * nd) {
        assert(gene < _gene_forests.size());
        assert(nd);
        assert(nd->_partial);
        ps.putPartial(gene, nd->_partial);
        
        // Decrement shared pointer reference count
        nd->_partial.reset();
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
    
    inline unsigned Particle::proposeSpeciation(unsigned step) {
        // This function is only used for proposing speciation events when there are
        // complete gene trees available. It thus draws increments from a truncated
        // exponential distribution where the trunction point is the next height at
        // which at least one coalescent event combines lineages from two different
        // species.
        unsigned num_species_tree_lineages = 0;
        
        // Stores tuple (height, 0, vector of species) for each join in the current species forest.
        // Do not cap with ancestral species at this point.
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfo(sppinfo_vect);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // coalinfo_vect stores a tuple (height, gene + 1, vector of species)
        // for each join in any gene tree. To get 0-offset gene index,
        // subtract 1 from 2nd element of tuple (if 2nd element is 0, then
        // tuple represents a species tree join)
        vector<Forest::coalinfo_t> coalinfo_vect;
        
        // Add gene tree joins to coalinfo_vect
        // Just need coalescent events at this point in order to choose
        // limits for species tree increments
        for (unsigned g = 0; g < G::_nloci; g++) {
            GeneForest & gf = _gene_forests[g];
            gf.saveCoalInfo(coalinfo_vect);
        }
        
        // Sort coalinfo_vect from smallest to largest height
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // //temporary!
        // G::SpecLog speclog_element;
        // speclog_element._seed = seed;
        
        // Get maximum height of any gene tree
        double max_height = get<0>((*coalinfo_vect.rbegin()));
        
        if (step > 0) {
#if defined(DEBUG_COALLIKE)
            output("\nSpecies tree before creating speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);
#endif

            // Create speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(rng, left_spp, right_spp, anc_spp);
            
#if defined(DEBUG_COALLIKE)
            output("\nSpecies tree after creating speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);
#endif

            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _species_forest.buildCoalInfoVect();
            _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            _species_forest.debugShowCoalInfo("sppinfo_vect after speciation event", sppinfo_vect, /*fn*/"");
            _species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
#endif
                
            // Adjust elements of coalinfo_vect affected by species tree joins
            _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            _species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
#endif
                
            // //temporary!
            // speclog_element._left = left_spp;
            // speclog_element._right = right_spp;
            // speclog_element._anc = anc_spp;
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
#if defined(DEBUG_COALLIKE)
            output("\nSpecies tree before creating FINAL speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);
#endif
            // Create final speciation event
            G::species_t left_spp, right_spp, anc_spp;
            _species_forest.speciationEvent(rng, left_spp, right_spp, anc_spp);
            
#if defined(DEBUG_COALLIKE)
            output("\nSpecies tree after creating FINAL speciation:\n", G::LogCateg::DEBUGGING);
            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), G::LogCateg::DEBUGGING);
#endif
            // Let sppinfo_vect reflect current state of species forest
            //sppinfo_vect.clear();
            //_species_forest.buildCoalInfoVect();
            //_species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            //sort(sppinfo_vect.begin(), sppinfo_vect.end());

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            //_species_forest.debugShowCoalInfo("sppinfo_vect after FINAL speciation event", sppinfo_vect, /*fn*/"");
            //_species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
#endif
            // Adjust elements of coalinfo_vect affected by species tree joins
            //_species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            //_species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
#endif
        }
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear();
        _species_forest.buildCoalInfoVect();
        _species_forest.saveCoalInfo(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
#if defined(DEBUG_COALLIKE)
        // Show coalinfo_vect before fixing up
        _species_forest.debugShowCoalInfo("sppinfo_vect", sppinfo_vect, /*fn*/"");
        _species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
#endif
                
        // Adjust elements of coalinfo_vect affected by species tree joins
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

#if defined(DEBUG_COALLIKE)
        // Show coalinfo_vect after fixing up
        _species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
#endif

        // Add speciations into coalinfo_vect
        //BUG lines below fix 2nd level bug 2024-06-19
        // see also Particle::recordAllForests
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Compute coalescent likelihood and log weight
        //double prev_log_coallike = _prev_log_coallike;
#if defined(DEBUG_COALLIKE)
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/true);
#else
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
#endif
        _prev_log_weight = _log_weight;
        _log_weight = _log_coal_like - _prev_log_coal_like + log_weight_factor;

        resetPrevLogCoalLike();
        
        // //temporary!
        // speclog_element._logw = log_weight;
        // speclog_element._filtered = true;
        // G::_speclog[step].push_back(speclog_element);

        return num_species_tree_lineages;
    }
    
    inline void Particle::proposeCoalescence(unsigned step, unsigned locus) {
        // If first locus, rebuild species tree starting from
        // the height of the tallest gene forest over all loci
        if (locus == 0) {
            double max_gene_forest_height = 0.0;
            for (unsigned g = 0; g < G::_nloci; g++) {
                double h = _gene_forests[g].getHeight();
                if (h > max_gene_forest_height)
                    max_gene_forest_height = h;
            }
            _species_forest.rebuildStartingFromHeight(max_gene_forest_height);
            
            // //temporary!
            // cerr << _species_forest.makeNewick(9, /*use_names*/true, /*coalunitstrue*/false) << endl;
        }
        
        // Copy _log_likelihood to _prev_log_likelihood
        _gene_forests[locus].resetPrevLogLikelihood();
        
        // Advance locus gene forest by one coalescent event.
        // Returns true if a coalescence event was created, and
        // returns false if the next speciation event occurred first
        bool done = advanceByOneCoalescence(step, locus);
        while (!done) {
            done = advanceByOneCoalescence(step, locus);
        }
    }
    
    inline double Particle::calcTotalCoalRate(unsigned locus) {
        double total_rate = 0.0;
        total_rate += _gene_forests[locus].calcTotalRate(_species_tuples);
        return total_rate;
    }
    
    inline void Particle::advanceAllLineagesBy(double dt) {
        _species_forest.advanceAllLineagesBy(dt);
        for (auto & gene_forest : _gene_forests) {
            gene_forest.advanceAllLineagesBy(dt);
        }
    }
    
    inline double Particle::priorPrior(unsigned step, unsigned locus, double total_rate) {
        double log_weight = 0.0;
        
        // _species_tuples has already been filled
        // Each tuple entry stores:
        //  0. number of lineages
        //  1. gene index
        //  2. species within gene
        //  3. vector<Node *> holding lineage roots for gene/species combination

        // Prior-prior chooses species in which to have a coalescence event drawn
        // from the prior. The probability of coalescence in one particular
        // species i is
        //   p(i) = r_i/total_rate
        // where
        //   r_i = n_i*(n_i-1)/theta
        //   total_rate = sum_i r_i
        //   n_i = number of lineages in i.
        vector<double> probs(_species_tuples.size());
        transform(_species_tuples.begin(), _species_tuples.end(), probs.begin(), [total_rate,locus](Node::species_tuple_t & spp_tuple){
            double n = (double)get<0>(spp_tuple);
            assert(locus == (unsigned)get<1>(spp_tuple));
            return n*(n - 1.0)/(G::_theta*total_rate);
        });
                
        // Choose the species in which the coalescence will happen
        assert(probs.size() > 0);
        unsigned which = G::multinomialDraw(::rng, probs);
        assert(which < probs.size());
        
        GeneForest & gf = _gene_forests[locus];
        
        // Get the species involved
        G::species_t spp = get<2>(_species_tuples[which]);
            
        // Pull next available node
        Node * anc_node = gf.pullNode();
        if (!G::_simulating) {
            anc_node->_partial = pullPartial(locus);
        }
                
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
            pair<unsigned,unsigned> chosen_pair = ::rng->nchoose2(n);
            i = chosen_pair.first;
            j = chosen_pair.second;
            
            // Return the chosen pair in the supplied pointer variables
            first_index = i;
            second_index = j;
            
            // Join the chosen pair of lineages
            first_node  = node_vect[i];
            second_node = node_vect[j];
            gf.joinLineagePair(anc_node, first_node, second_node);
        }

        anc_node->setSpecies(spp);
        
        // Set anc_node split to union of the two child splits
        anc_node->_split.resize(G::_ntaxa);
        anc_node->_split += first_node->_split;
        anc_node->_split += second_node->_split;
        
        assert(first_node->getSpecies() == spp);
        assert(second_node->getSpecies() == spp);

        // Fix up _lineages
        gf.removeTwoAddOne(gf._lineages, first_node, second_node, anc_node);
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

#if defined(UPGMA_WEIGHTS)
            // Update _dmatrix and _dmatrix_rows
            Split s1 = first_node->getSplit();
            Split s2 = second_node->getSplit();
            
            //output(format("s1 = %s\n") % s1.createPatternRepresentation());
            //output(format("s2 = %s\n") % s2.createPatternRepresentation());
            //output("_dmatrix_rows:\n");
            //unsigned z = 0;
            //for (auto s : _dmatrix_rows) {
            //    output(format("%6d = %s\n") % z++ % s.createPatternRepresentation());
            //}
            
            G::mergeDMatrixPair(gf._dmatrix_rows, gf._dmatrix, s1, s2);

            double logL = gf.calcLogLikelihood(_species_forest);
            double prevLogL = gf.getPrevLogLikelihood();
            log_weight = logL - prevLogL;

            //output(format("\nLocus %d gene forest (logL = %g, prevLogL = %g, logw = %g):\n") % locus % logL % prevLogL % log_weight);
            //output(format("%s\n") % gf.makeNewick(9, true, false));
#else
            double logL = gf.calcLogLikelihood();
            double prevLogL = gf.getPrevLogLikelihood();
            double check_log_weight = logL - prevLogL;
            assert(fabs(check_log_weight - log_weight) < G::_small_enough);
#endif
        }
                
        return log_weight;
    }
    
    inline bool Particle::advanceByOneCoalescence(unsigned step, unsigned locus) {
        // Clear species_tuples vector. Each tuple entry stores:
        //  0. number of lineages
        //  1. gene index
        //  2. species within gene
        //  3. vector<Node *> holding lineage roots for gene/species combination
        _species_tuples.clear();
        
        // Get current height of this locus' gene forest
        double gene_forest_height = _gene_forests[locus].getHeight();
        
        // Find height of next speciation event
        auto speciation_tuple = _species_forest.findNextSpeciationEvent(gene_forest_height);
        double speciation_height = get<0>(speciation_tuple);
        double speciation_delta  = (speciation_height == G::_infinity) ? G::_infinity : speciation_height - gene_forest_height;
        G::species_t left_spp    = get<1>(speciation_tuple);
        G::species_t right_spp   = get<2>(speciation_tuple);
        G::species_t anc_spp     = get<3>(speciation_tuple);
                        
        // Visit each species within the specified locus, storing
        // the number of lineages in _species_tuples and computing
        // the total coalescence rate (total_rate).
        double total_rate = calcTotalCoalRate(locus);
        
        // Draw coalescence increment delta ~ Exponential(total_rate)
        double delta = (total_rate > 0.0 ? -log(1.0 - ::rng->uniform())/total_rate : G::_infinity);
        
        bool done = false;
        if (delta > speciation_delta) {
            // Speciation event comes before coalescence event
            
            // Forests save delta to allow reversion
            _gene_forests[locus].advanceAllLineagesBy(speciation_delta);
                                                
            // Advise gene tree of the change in the species tree
            // Nodes that are reassigned save their previous state
            // to allow reversion
            _gene_forests[locus].mergeSpecies(speciation_height, left_spp, right_spp, anc_spp);
        }
        else {
            // Coalescence event comes before speciation event
            done = true;

            // Forests save delta to allow reversion
            _gene_forests[locus].advanceAllLineagesBy(delta);

            _prev_log_weight = _log_weight;
            _log_weight = priorPrior(step, locus, total_rate);
        }
        
        return done;
    }

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
        assert(gene < G::_nloci);
        return _gene_forests[gene];
    }

    inline vector<GeneForest> & Particle::getGeneForests() {
        return _gene_forests;
    }
    
    inline const vector<GeneForest> & Particle::getGeneForestsConst() const {
        return _gene_forests;
    }
    
    inline const GeneForest & Particle::getGeneForestConst(unsigned gene) const {
        assert(gene < G::_nloci);
        return _gene_forests[gene];
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
        
        // Copy gene forests
        assert(G::_nloci == other._gene_forests.size());
        _gene_forests.resize(G::_nloci);
        for (unsigned i = 0; i < G::_nloci; i++) {
            _gene_forests[i] = other._gene_forests[i];
            _gene_forests[i].setParticle(this);
        }
        
//#if defined(UPGMA_WEIGHTS)
//        _dmatrix = other._dmatrix;
//        _dmatrix_rows = other._dmatrix_rows;
//#endif
        
        // Copy species forest
        _species_forest     = other._species_forest;
        
        _log_coal_like      = other._log_coal_like;
        _prev_log_coal_like = other._prev_log_coal_like;
        _log_weight         = other._log_weight;
        _prev_log_weight    = other._prev_log_weight;
        
        _species_tuples     = other._species_tuples;        
    }
    
    inline void Particle::operator=(const Particle & other) {
        copyParticleFrom(other);
    }
        
    inline pair<double,double> Particle::chooseTruncatedSpeciesForestIncrement(double truncate_at) {
        double upper_bound = truncate_at - _species_forest.getHeight();
        auto incr_rate_cum = _species_forest.drawTruncatedIncrement(::rng, upper_bound);
        double incr = get<0>(incr_rate_cum);
        double rate = get<1>(incr_rate_cum);
        assert(rate > 0.0);
//#if defined(ADHOC_REMOVE_BEFORE_RELEASE)
//        //temporary!
//        incr = 0.15;
//#endif
        _species_forest.advanceAllLineagesBy(incr);
        double cum = get<2>(incr_rate_cum);
        assert(cum > 0.0);
        double log_factor = log(cum);
        return make_pair(log_factor,incr);
    }
    
    inline double Particle::getMaxGeneTreeHeight() const {
        double maxh = 0.0;
        for (unsigned i = 0; i < G::_nloci; i++) {
            double h = _gene_forests[i].getHeight();
            if (h > maxh)
                maxh = h;
        }
        return maxh;
    }
    
//#if defined(UPGMA_WEIGHTS)
//    inline void Particle::resetDistanceMatrix() {
//        _dmatrix = G::_dmatrix;
//        _dmatrix_rows = G::_dmatrix_rows;
//    }
//#endif
}
