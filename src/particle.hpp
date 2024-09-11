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
            
            void threadComputePartials(unsigned first, unsigned last);
            void computeAllPartials();
            
            void recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) const;
            
            double calcLogLikelihood();
            double calcPrevLogLikelihood();
            void   resetAllPrevLogLikelihood();
            double calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) const;
            
            double getMaxGeneTreeHeight() const;
                                    
            void   advanceAllLineagesBy(double dt);
            double findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
            double calcTotalCoalRate(unsigned locus);

            void proposeCoalescence(unsigned step, unsigned locus);

            bool advanceByOneCoalescence(unsigned step, unsigned locus);

            double priorPrior(unsigned step, unsigned locus, double total_rate);
            
            double getLogWeight() const {return _log_weight;}
            
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
            const GeneForest & getGeneForest(unsigned gene) const;
            
            void setGeneTrees(vector<GeneForest> & gtvect);
            pair<double,double> chooseTruncatedSpeciesForestIncrement(double truncate_at);
            
            unsigned debugCountNumCoalEvents() const;
            void debugCheckPartials() const;
            void debugShowAllGeneForests() const;
            void debugCheckAllPrevSpeciesStacksEmpty() const;
            
            void copyParticleFrom(const Particle & other);
            void operator=(const Particle & other);
            
            void stowAllPartials();

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
            
            vector<Node::species_tuple_t> _species_tuples;
            
            mutable double          _log_weight;
            
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
        _log_weight = 0.0;
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
#if defined(UPGMA_CONSTRAINED)
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
        //cerr << "*** inside Particle::calcLogCoalescentLikelihood ***" << endl;
        
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
        vector<dmap> n_jb(G::_nloci);
        
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
        for (unsigned g = 1; g < G::_nloci; g++) {
            n_jb[g] = n_jb[0];
        }
        
        // prev_height[g][b] holds previous height for gene g and branch b
        vector<map<G::species_t, double> > prev_height(G::_nloci);
        
        // Ploidy for each gene (currently all genes assumed to be diploid)
        vector<double> p_j(G::_nloci, 2.0);
        
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
            output("Key:\n", 1);
            output("  height:   height above leaf-level of this coalescence or speciation event\n", 1);
            output("  gene:     0 if speciation event, gene index for coalescence events\n", 1);
            output("  spp:      species of the two lineages involved in a coalescence/speciation event\n", 1);
            output("  b:        the species tree branch (i.e. species) to which this event contributes\n", 1);
            output("  n_jbk:    the number of coalescences for branch b\n", 1);
            output("  c_jbk:    the sojourn time prior to coalescence k\n", 1);
            output("  log(r_b): the log of 4/pj, where pj is the ploidy of gene j\n", 1);
            output("  q_b:      the cumulative number of coalescences over all genes for branch b\n", 1);
            output(format("%12s %12s %24s %12s %12s %12s %12s %12s\n") % "height" % "gene" % "spp" % "b" % "n_jbk" % "c_jbk" % "log(r_b)" % "q_b", 1);
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
                    output(format("%12.5f %12d %24s %12s %12s %12s %12s %12s\n") % height % 0 % banc % oss.str() % "-" % "-" % "-" % "-", 1);
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
                    output(format("%12.5f %12d %24s %12d %12d %12.5f %12.5f %12d\n") % height % gene_plus_one % oss.str() % b % n_jbk % c_jbk % gamma_b[b] % q_b[b], 1);
                }
            }
        }

        if (verbose) {
            output(format("\n%12s %12s %12s %12s %12s %12s\n") % "b" % "q_b" % "log(r_b)" % "gamma_b" % "theta_b" % "logL", 1);
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
                    output(format("%12d %12d %12.5f %12.5f %12.5f %12.5f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % beta % log_likelihood, 1);
                }
            }

            if (verbose) {
                output(format("\nalpha = %.9f\n") % alpha, 1);
                output(format("beta  = %.9f\n") % beta, 1);
                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood, 1);
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
                    output(format("%12d %12d %12.5f %12.5f %12.5f %12.5f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % theta_b % log_likelihood, 1);
                }
            }
            
            if (verbose) {
                output(format("\nsum log(r_b) = %.9f\n") % sum_log_rb, 1);
                output(format("sum q_b      = %d\n") % sum_qb, 1);
                output(format("sum gamma_b  = %d\n") % sum_gamma_b, 1);
                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood, 1);
            }
        }
                            
        return log_likelihood;
    }
    
    inline void Particle::computeAllPartials() {
        for (auto & gf : _gene_forests) {
            gf.computeAllPartials();
        }
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
    
    inline void Particle::proposeCoalescence(unsigned step, unsigned locus) {
        if (locus == 0) {
            double max_gene_forest_height = 0.0;
            for (unsigned g = 0; g < G::_nloci; g++) {
                double h = _gene_forests[g].getHeight();
                if (h > max_gene_forest_height)
                    max_gene_forest_height = h;
            }
            _species_forest.rebuildStartingFromHeight(max_gene_forest_height);
        }
        
        _gene_forests[locus].resetPrevLogLikelihood();
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

        // Prior-prior chooses a locus and species in which to have a coalescence
        // from the prior. The probability of coalescence in one particular
        // gene-species pair i is
        //   p(i) = r_i/total_rate
        // where
        //   r_i = n_i*(n_i-1)/theta
        //   total_rate = sum_i r_i
        //   n_i = number of lineages in pair i.
        vector<double> probs(_species_tuples.size());
        transform(_species_tuples.begin(), _species_tuples.end(), probs.begin(), [total_rate,locus](Node::species_tuple_t & spp_tuple){
            double n = (double)get<0>(spp_tuple);
            assert(locus == (unsigned)get<1>(spp_tuple));
            return n*(n - 1.0)/(G::_theta*total_rate);
        });
                
        // Choose which gene-species pair in which the coalescence will happen
        if (probs.size() == 0) {
            cerr << "debug stop" << endl;
        }
        unsigned which = G::multinomialDraw(::rng, probs);
        assert(which < probs.size());
        
        GeneForest & gf = _gene_forests[locus];
        
        // Get the species involved
        G::species_t spp = get<2>(_species_tuples[which]);
            
        // Pull next available node
        Node * anc_node = gf.pullNode();
        anc_node->_partial = pullPartial(locus);
                
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

        // Fix up _lineages
        gf.removeTwoAddOne(gf._lineages, first_node, second_node, anc_node);
        gf.refreshAllPreorders();

        // Compute partial likelihood array of ancestral node
        assert(_data);
        assert(anc_node->_left_child);
        assert(anc_node->_left_child->_right_sib);
        assert(anc_node->_partial);
        
        log_weight = gf.calcPartialArray(anc_node);
        assert(!isnan(log_weight));
        assert(!isinf(log_weight));

#if defined(UPGMA_WEIGHTS)
#   if defined(UPGMA_CONSTRAINED)
        double logL = gf.calcLogLikelihood(_species_forest);
#   else
        double logL = gf.calcLogLikelihood();
#   endif
        double prevLogL = gf.getPrevLogLikelihood();
        log_weight = logL - prevLogL;

        //output(format("\nLocus %d gene forest (logL = %g, prevLogL = %g, logw = %g):\n") % locus % logL % prevLogL % log_weight, 0);
        //output(format("%s\n") % gf.makeNewick(9, true, false), 0);
#else
        double logL = gf.calcLogLikelihood();
        double prevLogL = gf.getPrevLogLikelihood();
        double check_log_weight = logL - prevLogL;
        assert(fabs(check_log_weight - log_weight) < G::_small_enough);
#endif
                
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
            // Forests save delta to allow reversion
            _gene_forests[locus].advanceAllLineagesBy(speciation_delta);
                                                
            // Advise gene tree of the change in the species tree
            // Nodes that are reassigned save their previous state
            // to allow reversion
            _gene_forests[locus].mergeSpecies(speciation_height, left_spp, right_spp, anc_spp);
        }
        else {
            done = true;

            // Forests save delta to allow reversion
            _gene_forests[locus].advanceAllLineagesBy(delta);

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
    
    inline const GeneForest & Particle::getGeneForest(unsigned gene) const {
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
        
        // Copy species forest
        _species_forest = other._species_forest;
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
        _species_forest.advanceAllLineagesBy(incr);
        double log_factor = get<2>(incr_rate_cum);
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
    
}
