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
            
            void setSMC(SMC * smc);
            
            double setThetas();
            
            void resetSpeciesForest();
            void resetGeneForests(bool compute_partials);
            
            //void incrementSpeciations() {_nspeciations++;}
            //unsigned getSpeciations() const {return _nspeciations;}
            //void clearSpeciations() {_nspeciations = 0;}
            
            void threadComputePartials(unsigned first, unsigned last);
            void computeAllPartials();
            
            void recordAllForests(vector<Forest::coalinfo_t> & coalinfo_vect) const;
            
            double calcLogLikelihood();
            double calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) const;
            
            double getPrevLogCoalLike() {return _prev_log_coallike;}
            void setPrevLogCoalLike(double lnL) {_prev_log_coallike = lnL;}
                        
            //double calcTotalCoalRate(double speciation_increment);
            double calcTotalCoalRate();
            
            void clearMarkAllForests();
            void revertToMarkSpeciesForest();
            void revertToMarkAllForests();
            void advanceAllLineagesBy(double dt, bool mark);
            double findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
            pair<double, unsigned> proposeSpeciation(unsigned seed, unsigned step, unsigned pindex, bool make_permanent);
            pair<double, unsigned> proposeCoalescence(unsigned seed, unsigned step, unsigned pindex, bool compute_partial, bool make_permanent);
            void finalizeProposalSpeciesForest();
            void finalizeProposalAllForests();
            double priorPrior(unsigned step, unsigned pindex, double total_rate, bool compute_partial);
            double priorPost(unsigned step, unsigned pindex, bool make_permanent);
            void advance(unsigned step, unsigned pindex, bool compute_partial, bool make_permanent);
            
            unsigned getCount() const {return _count;}
            void setCount(unsigned cnt) {_count = cnt;}
            
#if defined(USING_MULTITHREADING)
            unsigned getXtra() const {return _xtra;}
            void setXtra(unsigned x) {_xtra = x;}
            
            unsigned getBeginIndex() const {return _begin_index;}
            void setBeginIndex(unsigned x) {_begin_index = x;}
#endif

#if defined(DEBUGGING_SANITY_CHECK)
            void debugStoreForests();
            void debugCheckForests();
#endif

            double getLogWeight() const {return _log_weight;}
            
            Lot::SharedPtr getLot() const {return _lot;}
            void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            
            enum last_event_t {
                LAST_EVENT_UNDEFINED = 0,
                LAST_EVENT_SPECIATION = 1,
                LAST_EVENT_COALESCENCE = 2
            };
            
            void setLastEvent(last_event_t last) {_last_event = last;}
            bool lastEventUndefined() const {return _last_event == LAST_EVENT_UNDEFINED;}
            bool lastEventSpeciation() const {return _last_event == LAST_EVENT_SPECIATION;}
            bool lastEventCoalescence() const {return _last_event == LAST_EVENT_COALESCENCE;}
            
            void refreshHeightsInternalsPreorders();

            void                  setSpeciesTree(string species_tree_newick);
            void                  copySpeciesForestFrom(SpeciesForest & sf);
            SpeciesForest       & getSpeciesForest();
            const SpeciesForest & getSpeciesForestConst() const;

            vector<GeneForest>       & getGeneForests();
            const vector<GeneForest> & getGeneForestsConst() const;
            
            vector<GeneForest> * getGeneTrees() {return _gene_trees;}
            const vector<GeneForest> * getGeneTreesConst() const {return _gene_trees;}
            
            GeneForest       & getGeneForest(unsigned gene);
            const GeneForest & getGeneForest(unsigned gene) const;
            
            void setGeneTrees(vector<GeneForest> & gtvect);
            //void forgetSpeciesTreeAbove(double height);
            pair<double,double> chooseTruncatedSpeciesForestIncrement(double truncate_at, bool mark);
            
            void debugShowMarkVariables(string title) const;
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

            SMC *                      _smc;
            Data::SharedPtr            _data;

            vector<GeneForest>         _gene_forests;
            
            // points to nullptr by default
            vector<GeneForest> *       _gene_trees;

            SpeciesForest              _species_forest;
            //unsigned                 _nspeciations;
            last_event_t               _last_event;
            unsigned                   _count;
            
#if defined(USING_MULTITHREADING)
            unsigned               _xtra;
            unsigned               _begin_index;
#endif

            vector<Node::species_tuple_t> _species_tuples;
            
            //unsigned                _last_proposed_gene;
            //G::species_t            _last_proposed_spp;
            //unsigned                _last_proposed_first_index;
            //unsigned                _last_proposed_second_index;

#if defined(DEBUGGING_SANITY_CHECK)
            string                  _debug_sfbefore;
            vector<string>          _debug_gfbefore;
#endif

            // Only used if SMC mode is SPECIES_GIVEN_GENE
            double                  _prev_log_coallike;
            
            mutable double          _log_weight;
            
            // Even though a shared pointer, _lot is a private random number
            // generator not shared with any other particle and has nothing to
            // to do with the global Lot shared_ptr rng defined in main.cpp
            mutable Lot::SharedPtr  _lot;
    };

    void Particle::setSpeciesTree(string species_tree_newick) {
        _species_forest.clear();
        _species_forest.buildFromNewick(species_tree_newick);
    }
    
    void Particle::copySpeciesForestFrom(SpeciesForest & sf) {
        _species_forest = sf;
    }

}
