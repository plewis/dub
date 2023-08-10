#pragma once

extern void output(string msg);
extern void output(string msg, unsigned level);

namespace proj {

    struct Epoch  {
        typedef map<Node::species_t, unsigned> lineage_counts_t;

        enum epoch_t {
            init_epoch       = 0x01,
            coalescent_epoch = 0x02,
            speciation_epoch = 0x04
        };
            
        Epoch(epoch_t t, double h);
            
        virtual bool isInitEpoch()       const {return _type == init_epoch;}
        virtual bool isCoalescentEpoch() const {return _type == coalescent_epoch;}
        virtual bool isSpeciationEpoch() const {return _type == speciation_epoch;}
        
        void reset() {_valid = true;}
        
        bool operator==(const Epoch & other) const;
        
        string leftSpeciesAsStr()  const {return speciesSetToStr(_left_species);}
        string rightSpeciesAsStr() const {return speciesSetToStr(_right_species);}
        string ancSpeciesAsStr()   const {return speciesSetToStr(_anc_species);}
        string speciesSetAsStr()   const {return speciesSetToStr(_species);}
            
        int                 _type;
        double              _height;
        bool                _valid;
        
        // Used only for init_epoch
        lineage_counts_t    _lineage_counts;    // counts of the number of lineages in each species just prior to coalescence
        
        // Used only for coalescent_epoch
        Node *              _coalescence_node;   // the coalescent node created by the coalescence event
        Node::species_t _species;               // species of coalescent node
        
        // Used only for init_epoch and coalescent_epoch
        int                 _gene;              // index of gene in which coalescence occurred
        
        // Used only for speciation_epoch
        //unsigned _left_species;
        //unsigned _right_species;
        //unsigned _anc_species;
        Node::species_t _left_species;
        Node::species_t _right_species;
        Node::species_t _anc_species;
        
        static string speciesSetToStr(const Node::species_t & s);
    };
    
    typedef list<Epoch> epoch_list_t;
    
    inline Epoch::Epoch(epoch_t t, double h) {
        _type                       = t;
        _height                     = h;
        _valid                      = false;
        _gene                       = -1;
        _coalescence_node           = nullptr;
        //_left_species               = -1;
        //_right_species              = -1;
        //_anc_species                = -1;
        // _lineage_counts is empty by default and does not need to be initialized
        // _species_set is empty by default and does not need to be initialized
    }
    
    inline bool Epoch::operator==(const Epoch & other) const {
        bool eq = true;
        eq = eq && _type                    == other._type;
        eq = eq && _height                  == other._height;
        eq = eq && _valid                   == other._valid;
        eq = eq && _gene                    == other._gene;
        eq = eq && _coalescence_node        == other._coalescence_node;
        eq = eq && _left_species            == other._left_species;
        eq = eq && _right_species           == other._right_species;
        eq = eq && _anc_species             == other._anc_species;
        eq = eq && _species             == other._species;
        eq = eq && _lineage_counts          == other._lineage_counts;
        return eq;
    }
    
    inline string Epoch::speciesSetToStr(const Node::species_t & s) {
        ostringstream oss;
        copy(s.begin(), s.end(), ostream_iterator<unsigned>(oss, "+"));

        // Remove the trailing '+'
        string returned_str(oss.str());
        returned_str.pop_back();
        
        return returned_str;
    }

    inline bool epochLess (const Epoch & i, const Epoch & j) {
        return (i._height < j._height);
    }
    
    inline epoch_list_t::iterator pushFrontEpoch(epoch_list_t & epochs, const Epoch & e) {
        return epochs.insert(epochs.begin(), e);
    }
    
    inline epoch_list_t::iterator insertEpochBefore(epoch_list_t & epochs, const Epoch & e, epoch_list_t::iterator & it) {
        return epochs.insert(it, e);
    }
    
    inline epoch_list_t::iterator pushBackEpoch(epoch_list_t & epochs, const Epoch & e) {
        epochs.push_back(e);
        auto rit = epochs.rbegin();
        rit++;
        return rit.base();
    }
    
    inline void removeEpoch(epoch_list_t & epochs, const Epoch & e) {
        auto it = find(epochs.begin(), epochs.end(), e);
        assert(it != epochs.end());
        epochs.erase(it);
    }
    
    inline void removeEpochAt(epoch_list_t & epochs, epoch_list_t::iterator & it) {
        epochs.erase(it);
    }
    
    inline void resetAllEpochs(epoch_list_t & epochs) {
        for_each(epochs.begin(), epochs.end(), [](Epoch & e){e.reset();});
    }

    inline bool checkEpochs(const epoch_list_t & epochs, int valid_types) {
        double prev_height = -1.0;
        for (auto it = epochs.begin(); it != epochs.end(); it++) {
            if (it->_height < prev_height) {
                output(str(format("epoch height (%.9f) not less than previous epoch's height (%.9f)\n") % it->_height % prev_height));
                return false;
            }
            bool ok = (bool)(it->_type & valid_types);
            if (!ok) {
                output(str(format("epoch type (%d) not among valid types (%d)\n") % it->_type % valid_types));
                return false;
            }
            prev_height = it->_height;
        }
        return true;
    }
}
