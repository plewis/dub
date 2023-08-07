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
                    typedef set<unsigned>   species_t;
        
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
                    void                setSpeciesFromUnsigned(unsigned spp);
        
                    double              getHeight() const       {return _height;}
                    void                setHeight(double v);
        
                    unsigned            countChildren() const;

                    void                clearPointers()             {_left_child = _right_sib = _parent = 0;}
                    
                    static string       taxonNameToSpeciesName(string taxon_name);
                                        
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
        _species.clear();
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
        _species.insert(right.begin(), right.end());
        
        // Internal nodes can have species that are unions and hence have size > 1, but
        // leaf nodes should be assigned to just a single species. This function should only
        // be called for internal nodes.
        assert(_left_child);
    }
    
    inline void Node::setSpeciesFromUnsigned(unsigned spp) {
        _species = {spp};
        
        // Internal nodes can have species that are unions and hence have size > 1, but
        // leaf nodes should be assigned to just a single species. This function should only
        // be called for leaf nodes.
        assert(!_left_child);
    }
    
    inline string Node::taxonNameToSpeciesName(string tname) {
        vector<string> before_after;
        split(before_after, tname, boost::is_any_of("^"));
        if (before_after.size() != 2)
            throw XProj(format("Expecting taxon name to conform to taxon^species pattern: %s") % tname);
        return before_after[1];
    }

}
