#pragma once

namespace proj {

    class Particle;

    class SMC : public ParallelPolicyNone<SMC> {
        friend class ParallelPolicyNone<SMC>;
        public:
                                SMC()  {clear();}
            virtual             ~SMC() {}
            
            enum mode_type_t {
                SPECIES_AND_GENE = 0,
                SPECIES_GIVEN_GENE = 1
            };
            bool isJointMode() const {return _mode == SPECIES_AND_GENE;}
            bool isConditionalMode() const {return _mode == SPECIES_GIVEN_GENE;}
            
            unsigned            getNParticles() {return _nparticles;}
            list<Particle> &    getParticles()  {return _particle_list;}
            
            void setMode(mode_type_t m)     {_mode = m;}
            void setNParticles(unsigned n)  {_nparticles = n;}
            void setData(Data::SharedPtr d) {_data = d;}
            void init();
            void initFromParticle(Particle & p);
            void run();
            double filterParticles(unsigned step, list<Particle> & particle_list, vector<double> & log_weights, vector<unsigned> & counts, vector<unsigned> & rnseeds);
            double calcLogSpeciesTreePrior(vector<Forest::coalinfo_t> & coalinfo_vect, bool include_join_probs) const;
            void pruneParticles(list<Particle> & particle_list);
            void findNonZeroCountsInRange(vector<unsigned> & nonzeros, const vector<unsigned> & counts, unsigned begin_index, unsigned end_index) const;
            double computeEffectiveSampleSize(const vector<double> & probs) const;
            bool compareToReferenceTrees(list<Particle> particle_list, map<string, tuple<unsigned, double, double, double, double> > & m);
            void outputAnnotatedNexusTreefile(string fn, const vector<tuple<unsigned, double, string, string, string> > & treeinfo) const;
            void saveAllSpeciesTrees(string fn, const list<Particle> & particle_list, unsigned compression_level = 2);
            void saveAllGeneTrees(unsigned gene, string fn, const list<Particle> & particle_list, unsigned compression_level = 2);
            //void                         saveUniqueSpeciesTrees(string fn, const vector<Particle> particles, const vector<unsigned> & counts);
            void summarize();
            void dumpParticles(SMC & ensemble);
            void clear();
            unsigned countDistinctGeneTreeTopologies();
            void saveParamsForLoRaD(string prefix);
            
            typedef shared_ptr<SMC> SharedPtr;
            
        private:
        
            Data::SharedPtr              _data;
            unsigned                     _mode;
            unsigned                     _nparticles;
            vector<unsigned>             _counts;
            vector<unsigned>             _update_seeds;
            
            unsigned                     _nsteps;
            //double                       _starting_log_likelihood;

            list<Particle>               _particle_list;
            double                       _log_marg_like;
            vector<double>               _log_weights;
            //vector<unsigned>             _proposed_gene;
            //vector<G::species_t>         _proposed_spp;
            //vector<unsigned>             _proposed_first;
            //vector<unsigned>             _proposed_second;
            //vector<unsigned>             _proposed_species_tree_lineages;
            
    };
        
}
