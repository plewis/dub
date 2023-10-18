#pragma once    

namespace proj {

    class TreeManip;
    class Likelihood;
    class Updater;
    class Epoch;
    class CoalescentEpoch;
    class SpeciationEpoch;
    class Forest;
    class SpeciesForest;
    class GeneForest;
    class Particle;

    class Node {

        friend class TreeManip;
        friend class Likelihood;
        friend class Updater;
        friend class Epoch;
        friend class CoalescentEpoch;
        friend class SpeciationEpoch;
        friend class Forest;
        friend class SpeciesForest;
        friend class GeneForest;
        friend class Particle;

        public:
                                        Node();
                                        ~Node();

                    typedef vector<Node *>  ptr_vect_t;
#if defined(SPECIES_IS_BITSET)
                    typedef unsigned long   species_t;
#else
                    typedef set<unsigned>   species_t;
#endif
        
                    Node *              getParent()                 {return _parent;}
                    const Node *        getParent() const           {return _parent;}

                    Node *              getLeftChild()              {return _left_child;}
                    const Node *        getLeftChild() const        {return _left_child;}

                    Node *              getRightSib()               {return _right_sib;}
                    const Node *        getRightSib() const         {return _right_sib;}

                    string              getName()                   {return _name;}
                    const string        getName() const             {return _name;}

                    int                 getNumber() const           {return _number;}
                    Split               getSplit()                  {return _split;}
        
                    bool                isSelected()                {return _flags & Flag::Selected;}
                    void                select()                    {_flags |= Flag::Selected;}
                    void                deselect()                  {_flags &= ~Flag::Selected;}

                    bool                isSelPartial()              {return _flags & Flag::SelPartial;}
                    void                selectPartial()             {_flags |= Flag::SelPartial;}
                    void                deselectPartial()           {_flags &= ~Flag::SelPartial;}

                    bool                isSelTMatrix()              {return _flags & Flag::SelTMatrix;}
                    void                selectTMatrix()             {_flags |= Flag::SelTMatrix;}
                    void                deselectTMatrix()           {_flags &= ~Flag::SelTMatrix;}

                    bool                isAltPartial()              {return _flags & Flag::AltPartial;}
                    void                setAltPartial()             {_flags |= Flag::AltPartial;}
                    void                clearAltPartial()           {_flags &= ~Flag::AltPartial;}

                    bool                isAltTMatrix()              {return _flags & Flag::AltTMatrix;}
                    void                setAltTMatrix()             {_flags |= Flag::AltTMatrix;}
                    void                clearAltTMatrix()           {_flags &= ~Flag::AltTMatrix;}
                    
                    void                flipTMatrix()               {isAltTMatrix() ? clearAltTMatrix() : setAltTMatrix();}
                    void                flipPartial()               {isAltPartial() ? clearAltPartial() : setAltPartial();}

                    double              getEdgeLength() const       {return _edge_length;}
                    void                setEdgeLength(double v);
                    
                    const species_t &   getSpecies() const {return _species;}
                    species_t &         getSpecies() {return _species;}
                    void                setSpecies(const species_t & other);
                    void                setSpeciesToUnion(const species_t & left, const species_t & right);
#if !defined(SPECIES_IS_BITSET)
                    void                setSpeciesFromUnsigned(unsigned spp);
#endif
        
                    double              getHeight() const       {return _height;}
                    void                setHeight(double v);
        
                    unsigned            countChildren() const;

                    void                clearPointers()             {_left_child = _right_sib = _parent = 0;}
                    
                    static string       taxonNameToSpeciesName(string taxon_name);

#if defined(SPECIES_IS_BITSET)
                    static string       speciesStringRepresentation(species_t species);
                    static void         setSpeciesMask(species_t & mask, unsigned nspecies);
                    static void         setSpeciesBits(species_t & to_species, const species_t & from_species, bool init_to_zero_first);
                    static void         unsetSpeciesBits(species_t & to_species, const species_t & from_species);
                    static void         setSpeciesBit(species_t & to_species, unsigned i, bool init_to_zero_first);
                    static void         unsetSpeciesBit(species_t & to_species, unsigned i);
#endif
                                        
            static const double _smallest_edge_length;

        private:
        
            enum Flag {
                Selected   = (1 << 0),
                SelPartial = (1 << 1),
                SelTMatrix = (1 << 2),
                AltPartial = (1 << 3),
                AltTMatrix = (1 << 4)
            };

            void                clear();
            
            // Adding data members? Be sure to get them copied in Forest::operator=(const Forest & other)

            Node *          _left_child;
            Node *          _right_sib;
            Node *          _parent;
            int             _number;
            string          _name;
            double          _edge_length;
            double          _height;         // distance from node to any leaf
            species_t       _species;    // set of species (indices) compatible with this node
            Split           _split;
            int             _flags;
            
            PartialStore::partial_t _partial;
    };
    
#if defined(SPECIES_IS_BITSET)
    inline string Node::speciesStringRepresentation(species_t species) {
        species_t species_copy = species;
        unsigned bits_avail = (unsigned)sizeof(species_t);
        string s;
        for (unsigned i = 0; i < bits_avail; ++i) {
            species_t bitmask = ((species_t)1 << i);
            bool bit_is_set = ((species_copy & bitmask) > (species_t)0);
            if (bit_is_set) {
                // Add species i to the string
                s += to_string(i);
                
                // Zero that bit so we know when we are done
                species_copy &= ~bitmask;
            }
            if (!species_copy) {
                // If species_copy is zero, there are no more bits set
                break;
            }
        }
        return s;
    }
    
    inline void Node::setSpeciesMask(species_t & mask, unsigned nspecies) {
        mask = (species_t)0;
        for (unsigned i = 0; i < nspecies; ++i) {
            mask |= ((species_t)1 << i);
        }
    }
    
    inline void Node::setSpeciesBits(species_t & to_species, const species_t & from_species, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (species_t)0;

        // Copy bits in from_species to to_species
        to_species |= from_species;
    }
    
    inline void Node::unsetSpeciesBits(species_t & to_species, const species_t & from_species) {
        // Zero from_species' bits in to_species
        to_species &= ~from_species;
    }
    
    inline void Node::setSpeciesBit(species_t & to_species, unsigned i, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (species_t)0;
            
        // Set ith bit in to_species
        to_species |= ((species_t)1 << i);
    }
    
    inline void Node::unsetSpeciesBit(species_t & to_species, unsigned i) {
        // Unset ith bit in to_species
        to_species &= ~((species_t)1 << i);
    }
#endif
    
    inline Node::Node() {
        clear();
    }

    inline Node::~Node() {
    }

    inline void Node::clear() { 
        _flags = 0;
        clearPointers();    
        //_left_child = 0;
        //_right_sib = 0;
        //_parent = 0;      
        _number = -1;
        _name = "";
        _edge_length = _smallest_edge_length;
        _height = 0.0;
#if defined(SPECIES_IS_BITSET)
        _species = (species_t)0;
#else
        _species.clear();
#endif
        _partial = nullptr;
    }   

    inline void Node::setEdgeLength(double v) {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

    inline void Node::setHeight(double h) {
        _height = h;
    }

    inline unsigned Node::countChildren() const {
        unsigned n_children = 0;
        for (Node * child = _left_child; child; child=child->_right_sib) {
            n_children ++;
        }
        return n_children;
    }
    
    inline void Node::setSpecies(const species_t & other) {
        _species = other;
    }
    
    inline void Node::setSpeciesToUnion(const species_t & left, const species_t & right) {
        _species = left;
#if defined(SPECIES_IS_BITSET)
        Node::setSpeciesBits(_species, right, /*init_to_zero_first*/false);
#else
        _species.insert(right.begin(), right.end());
#endif
        
        // Internal nodes can have species that are unions and hence have size > 1, but
        // leaf nodes should be assigned to just a single species. This function should only
        // be called for internal nodes.
        assert(_left_child);
    }
    
#if defined(SPECIES_IS_BITSET)
    // setSpeciesFromUnsigned not defined
#else
    inline void Node::setSpeciesFromUnsigned(unsigned spp) {
        _species = {spp};
        
        // Internal nodes can have species that are unions and hence have size > 1, but
        // leaf nodes should be assigned to just a single species. This function should only
        // be called for leaf nodes.
        assert(!_left_child);
    }
#endif
    
    inline string Node::taxonNameToSpeciesName(string tname) {
        vector<string> before_after;
        split(before_after, tname, boost::is_any_of("^"));
        if (before_after.size() != 2)
            throw XProj(format("Expecting taxon name to conform to taxon^species pattern: %s") % tname);
        return before_after[1];
    }

}
