#pragma once

//extern void output(string msg, proj::G::verbosity_t verb);
//extern void output(format & fmt, proj::G::verbosity_t level);
//extern proj::PartialStore         ps;
//extern proj::StopWatch            stopwatch;
//extern proj::Lot::SharedPtr       rng;
//extern proj::Partition::SharedPtr partition;
//extern proj::Data::SharedPtr      data;

namespace proj {

    // Particle encapsulates generic particle data and methods
    class Particle {
        public:
            Particle() : _is_species_tree(false), _nleaves(0), _index(0), _height(0.0), _next_node_number(0) {}
            ~Particle() {}
            
            void clear();
            
            // Setters
            void setIndex(unsigned i) {_index = i;}

            // Getters
            unsigned getIndex() const {return _index;}
            unsigned getNLeaves() const {return _nleaves;}
            double   getHeight() const {return _height;}
            
            void lineagesRemoveTwoAddOne(Node * del1, Node * del2, Node * add);
            void extendAllLineagesBy(double incr);
            void makeAnc(Node * anc, Node * lchild, Node * rchild);
            
            void buildFromNewick(const string newick, map<unsigned,unsigned> & taxon_map, vector<string> & taxon_names);
            string makeNewick(unsigned precision, bool use_names) const;
            void   buildPreordersVector(vector<Node::ptr_vect_t> & preorders) const;
            void storeEdgelensBySplit(map<Split, double> & edgelenmap) const;
            void storeSplits(Split::treeid_t & splitset) const;
            void refreshSplits();
            Node * findNextPreorder(Node * nd) const;
            
            // Statics
            static pair<double,double> calcTreeDistances(const Particle & ref, const Particle & test);
            
            // Virtuals
            virtual void recordJoinInfo(vector<G::join_info_t> & join_info) const;
            virtual void operator=(const Particle & other);

            // Pure virtuals
            virtual pair<double, double> drawIncrement(Lot::SharedPtr lot) = 0;
            virtual void joinRandomPair(Lot::SharedPtr lot) = 0;
            virtual double calcLogLikelihood() const = 0;
            virtual void createTrivialForest() = 0;
            virtual string info() const = 0;
                        
        protected:
            bool canHaveSibling(Node * nd) const;

            vector<Node>            _nodes;
            vector<Node *>          _lineages;
            unsigned                _nleaves;
            unsigned                _index;
            double                  _height;
            bool                    _is_species_tree;
            unsigned                _next_node_number;
            vector<unsigned>        _returned_node_numbers;
    };
    
    inline void Particle::clear() {
        _lineages.clear();
        _nodes.clear();
        _nleaves = 0;
        _index = 0;
        _height = 0.0;
        _next_node_number = 0;
        _is_species_tree = false;
    }
     
    inline void Particle::extendAllLineagesBy(double incr) {
        _height += incr;
        for (auto nd : _lineages) {
            double new_edge_length = nd->getEdgeLength() + incr;
            nd->setEdgeLength(new_edge_length);
        }
    }
    
    inline void Particle::makeAnc(Node * anc, Node * lchild, Node * rchild) {
        // Sanity checks
        assert(anc);
        assert(lchild);
        assert(rchild);
        assert(lchild->_species > 0);
        assert(rchild->_species > 0);
        
        // Connections
        anc->_left_child = lchild;
        lchild->_right_sib = rchild;
        lchild->_parent = anc;
        rchild->_parent = anc;
        
        // anc height is the current height of the particle
        anc->setHeight(_height);
        
        // Species of anc is union of child species
        anc->_species = lchild->_species | rchild->_species;
        
        // Set split of new internal node
        anc->_split.resize(_nleaves);
        anc->_split.addSplit(lchild->_split);
        anc->_split.addSplit(rchild->_split);
        
        // Adjust _lineages
        lineagesRemoveTwoAddOne(lchild, rchild, anc);
    }
    
    inline void Particle::lineagesRemoveTwoAddOne(Node * del1, Node * del2, Node * add) {
        // Get iterator to first node to be deleted and remove from v
        auto it = find(_lineages.begin(), _lineages.end(), del1);
        assert(it != _lineages.end());
        _lineages.erase(it);
        
        // Get iterator to second node to be deleted and remove from v
        it = find(_lineages.begin(), _lineages.end(), del2);
        assert(it != _lineages.end());
        _lineages.erase(it);
        
        // Add node
        _lineages.push_back(add);
    }

    inline bool Particle::canHaveSibling(Node * nd) const {
        assert(nd);
        if (!nd->_parent) {
            // Trying to give root node a sibling
            return false;
        }

        bool nd_can_have_sibling = true;
        if (nd != nd->_parent->_left_child) {
            // nd is NOT the left child of its parent
            if (nd->_parent->_parent) {
                // Trying to give a sibling to a sibling of nd, and nd's parent is not the root
                nd_can_have_sibling = false;
            }
            else {
                // Trying to add a sibling to the right child of the root node
                nd_can_have_sibling = false;
            }
        }

        return nd_can_have_sibling;
    }

    inline void Particle::buildFromNewick(const string newick, map<unsigned,unsigned> & taxon_map, vector<string> & taxon_names) {
        // Builds strictly-bifurcating ultrametric, rooted, complete (i.e.
        // _lineages contains just one element) forest from the supplied newick
        // tree description. Assumes newick either uses names for leaves or,
        // if it specifies numbers, the numbers correspond to keys in
        // taxon_map, which translates taxon numbers in newick
        // strings to the index of the taxon in taxon_names.
        
        // Remove all nodes and all pointers to nodes
        clear();
        
        unsigned curr_leaf = 0;
        unsigned num_edge_lengths = 0;

        // Set that ensures that no two leaf nodes have the same number
        set<unsigned> used;
        
        // Remove comments from the supplied newick string
        string commentless_newick = newick;
        regex commentexpr("\\[.*?\\]");
        commentless_newick = regex_replace(commentless_newick, commentexpr, string(""));

        // Count leaves in newick description
        regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
        sregex_iterator m1(commentless_newick.begin(), commentless_newick.end(), taxonexpr);
        sregex_iterator m2;
        unsigned nleaves = (unsigned)distance(m1, m2);
        if (nleaves < 2)
            throw XProj("Expecting newick tree description to have at least 2 leaves");

        // Assumes rooted, fully-bifurcating tree
        //
        //  A    B    C    nleaves = 3
        //   \   /   /     nnodes  = 5 = 2*3 - 1
        //    \ /   /
        //     X   /
        //      \ /
        //       Y
        unsigned max_nodes = 2*nleaves - 1;

        // Resize the _nodes vector
        _nodes.resize(max_nodes);
        
        try {
            // Root node is the first one in the _nodes array
            assert(_next_node_number == 0);
            Node * nd   = &_nodes[_next_node_number];
            Node * root = &_nodes[_next_node_number++];
            
            // The _lineages vector will have just one entry because
            // this will be a complete tree
            _lineages.push_back(root);
                        
            // Some flags to keep track of what we did last
            enum {
                Prev_Tok_LParen		= 0x01,	// previous token was a left parenthesis ('(')
                Prev_Tok_RParen		= 0x02,	// previous token was a right parenthesis (')')
                Prev_Tok_Colon		= 0x04,	// previous token was a colon (':')
                Prev_Tok_Comma		= 0x08,	// previous token was a comma (',')
                Prev_Tok_Name		= 0x10,	// previous token was a node name (e.g. '2', 'P._articulata')
                Prev_Tok_EdgeLen	= 0x20	// previous token was an edge length (e.g. '0.1', '1.7e-3')
            };
            unsigned previous = Prev_Tok_LParen;

            // Some useful flag combinations
            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

            // Set to true while reading an edge length
            bool inside_edge_length = false;
            string edge_length_str;
            unsigned edge_length_position = 0;

            // Set to true while reading a node name surrounded by (single) quotes
            bool inside_quoted_name = false;

            // Set to true while reading a node name not surrounded by (single) quotes
            bool inside_unquoted_name = false;

            // Set to start of each node name and used in case of error
            unsigned node_name_position = 0;
            
            // Useful if there is an error to show where in the newick string the problems began
            unsigned position_in_string = 0;

            // Loop through the characters in newick
            for (auto ch : commentless_newick) {
                position_in_string++;

                if (inside_quoted_name) {
                    if (ch == '\'') {
                        // Done reading quoted name
                        inside_quoted_name = false;
                        node_name_position = 0;
                        if (!nd->_left_child) {
                            // Set node number
                            int num = G::extractNodeNumberFromName(nd->_name, used);
                            assert(num > 0);
                            nd->_number = taxon_map[(unsigned)num];
                            
                            // Set node name
                            assert(nd->_number < 2*taxon_names.size() - 1);
                            if (nd->_number < taxon_names.size())
                                nd->_name = taxon_names[nd->_number];
                            
                            // Set node species
                            G::setSpeciesBit(nd->_species, G::_taxon_to_species.at(nd->_name), /*init_to_zero_first*/true);

                            curr_leaf++;
                        }
                        previous = Prev_Tok_Name;
                    }
                    else if (iswspace(ch))
                        nd->_name += ' ';
                    else
                        nd->_name += ch;

                    continue;
                }
                else if (inside_unquoted_name) {
                    if (ch == '(')
                        throw XProj(str(format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));

                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                        inside_unquoted_name = false;

                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(str(format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));

                        if (!nd->_left_child) {
                            // Set node number
                            int num = G::extractNodeNumberFromName(nd->_name, used);
                            assert(num > 0);
                            nd->_number = G::_nexus_taxon_map[(unsigned)num];

                            // Set node name
                            assert(nd->_number < 2*taxon_names.size() - 1);
                            if (nd->_number < taxon_names.size())
                                nd->_name = taxon_names[nd->_number];
                                
                            // Set node species
                            G::setSpeciesBit(nd->_species, G::_taxon_to_species.at(nd->_name), /*init_to_zero_first*/true);

                            curr_leaf++;
                        }

                        previous = Prev_Tok_Name;
                    }
                    else {
                        nd->_name += ch;
                        continue;
                    }
                }
                else if (inside_edge_length) {
                    if (ch == ',' || ch == ')' || iswspace(ch)) {
                        inside_edge_length = false;
                        edge_length_position = 0;
                        double edge_len = G::extractEdgeLen(edge_length_str);
                        nd->setEdgeLength(edge_len);
                        ++num_edge_lengths;
                        previous = Prev_Tok_EdgeLen;
                    }
                    else {
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid)
                            throw XProj(str(format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
                        edge_length_str += ch;
                        continue;
                    }
                }

                if (iswspace(ch))
                    continue;

                switch(ch) {
                    case ';':
                        break;

                    case ')':
                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
                        if (!nd->_parent)
                            throw XProj(str(format("Too many right parentheses at position %d in tree description") % position_in_string));

                        // Expect right paren only after an edge length, a node name, or another right paren
                        if (!(previous & RParen_Valid))
                            throw XProj(str(format("Unexpected right parenthesis at position %d in tree description") % position_in_string));

                        // Go down a level
                        nd = nd->_parent;
                        if (!nd->_left_child->_right_sib)
                            throw XProj(str(format("Internal node has only one child at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_RParen;
                        break;

                    case ':':
                        // Expect colon only after a node name or another right paren
                        if (!(previous & Colon_Valid))
                            throw XProj(str(format("Unexpected colon at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_Colon;
                        break;

                    case ',':
                        // Expect comma only after an edge length, a node name, or a right paren
                        if (!nd->_parent || !(previous & Comma_Valid))
                            throw XProj(str(format("Unexpected comma at position %d in tree description") % position_in_string));

                        // Check for polytomies
                        if (!canHaveSibling(nd)) {
                            throw XProj(str(format("Polytomy found in the following tree description but polytomies are not allowed:\n%s") % newick));
                        }

                        // Create the sibling
                        nd->_right_sib = &_nodes[_next_node_number++];
                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;

                        break;

                    case '(':
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw XProj(str(format("Not expecting left parenthesis at position %d in tree description") % position_in_string));

                        // Create new node above and to the left of the current node
                        assert(!nd->_left_child);
                        nd->_left_child = &_nodes[_next_node_number++];
                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;

                        break;

                    case '\'':
                        // Encountered an apostrophe, which always indicates the start of a
                        // node name (but note that node names do not have to be quoted)

                        // Expect node name only after a:
                        // - left paren (child's name)
                        // - comma (sib's name)
                        // - right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(str(format("Not expecting node name at position %d in tree description") % position_in_string));

                        // Prepare for receiving the node name in upcoming characters
                        nd->_name.clear();
                        inside_quoted_name = true;
                        node_name_position = position_in_string;

                        break;

                    default:
                        // Get here if ch is not one of ();:,'

                        // Expecting either an edge length or an unquoted node name
                        if (previous == Prev_Tok_Colon) {
                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                            inside_edge_length = true;
                            edge_length_position = position_in_string;
                            edge_length_str = ch;
                        }
                        else {
                            // Get the next node name character
                            nd->_name = ch;
                            inside_unquoted_name = true;
                            node_name_position = position_in_string;
                        }
                }   // end of switch statement
            }   // loop over characters in newick string

            if (inside_unquoted_name)
                throw XProj(str(format("Tree description ended before end of node name starting at position %d was found: \"%s\"") % node_name_position % newick));
            if (inside_edge_length && !(nd == root))
                throw XProj(str(format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
            if (inside_quoted_name)
                throw XProj(str(format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));
            _nleaves = curr_leaf;
            refreshSplits();
        }
        catch(XProj x) {
            clear();
            throw x;
        }
    }
    
    inline string Particle::makeNewick(unsigned precision, bool use_names) const  {
        // Place basal polytomy (if there is one) at a height 10% greater than
        // _height
        double basal_polytomy_height = 0.0; //_height*(1.1);
        bool is_complete_tree = (bool)(_lineages.size() == 1);
        
        //const format basal_subtree_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_name_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_number_format( str(format("%%d:%%.%df") % precision) );
        const format internal_node_format( str(format("):%%.%df") % precision) );
        
        // Build preorders vector
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        vector<string> subtree_newicks;
        for (auto & preorder : preorders) {
            string subtree_newick;
            stack<Node *> node_stack;
            for (auto nd : preorder) {
                if (nd->_left_child) {
                    subtree_newick += "(";
                    node_stack.push(nd);
                }
                else {
                    double edge_length = nd->_edge_length;
                    if (!nd->_parent) {
                        // Subtree consists of just this one leaf node
                        edge_length += basal_polytomy_height;
                    }
                    if (use_names) {
                        if (precision > 0)
                            subtree_newick += str(format(tip_node_name_format) % nd->_name % edge_length);
                        else
                            subtree_newick += str(format("%s") % nd->_name);
                    } else {
                        if (precision > 0)
                            subtree_newick += str(format(tip_node_number_format) % (nd->_number + 1) % edge_length);
                        else
                            subtree_newick += str(format("%d") % (nd->_number + 1));
                    }
                    if (nd->_right_sib)
                        subtree_newick += ",";
                    else if (nd->_parent) {
                        Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                        double popped_edge_length = popped->_edge_length;
                        while (popped && !popped->_right_sib) {
                            node_stack.pop();
                            if (node_stack.empty()) {
                                if (is_complete_tree) {
                                    subtree_newick += ")";
                                }
                                else {
                                    // This is the root of one of several subtrees, so
                                    // it is important to preserve its edge length
                                    popped_edge_length += basal_polytomy_height;
                                    if (precision > 0)
                                        subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                    else
                                        subtree_newick += ")";
                                }
                                popped = 0;
                            }
                            else {
                                if (precision > 0)
                                    subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                else
                                    subtree_newick += ")";
                                popped = node_stack.top();
                                popped_edge_length = popped->_edge_length;
                            }
                        }
                        if (popped && popped->_right_sib) {
                            node_stack.pop();
                            if (precision > 0) {
                                subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                subtree_newick += ",";
                            }
                            else
                                subtree_newick += "),";
                        }
                    }
                }
            }
            subtree_newicks.push_back(subtree_newick);
        }
        
        string newick;
        if (is_complete_tree)
            newick = subtree_newicks[0];
        else {
            string separator = str(format(":%.5f,") % basal_polytomy_height);
            string insides = boost::join(subtree_newicks, ",");
            newick = str(format("(%s)") % insides);
        }

        return newick;
    }

    inline void Particle::buildPreordersVector(vector<Node::ptr_vect_t> & preorders) const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        preorders.clear();
        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
            // lineage is a Node::ptr_vect_t (i.e. vector<Node *>)
            // lineage[0] is the first node pointer in the preorder sequence for this lineage
            // Add a new vector to preorders containing, for now, just the root of the subtree
            preorders.push_back({nd});
            
            // Now add the nodes above the root in preorder sequence
            Node::ptr_vect_t & preorder = *preorders.rbegin();
            Node * nd2 = preorder[0];
            while (true) {
                nd2 = findNextPreorder(nd2);
                if (nd2)
                    preorder.push_back(nd2);
                else
                    break;
            }
        }
    }
    
    inline Node * Particle::findNextPreorder(Node * nd) const {
        assert(nd);
        Node * next = 0;
        if (!nd->_left_child && !nd->_right_sib) {
            // nd has no children and no siblings, so next preorder is the right sibling of
            // the first ancestral node that has a right sibling.

            Node * anc = nd->_parent;
            while (anc && !anc->_right_sib)
                anc = anc->_parent;
            if (anc) {
                // We found an ancestor with a right sibling
                next = anc->_right_sib;
            }
            else {
                // nd is last preorder node in the tree
                next = 0;
            }
        }
        else if (nd->_right_sib && !nd->_left_child) {
            // nd has no children (it is a tip), but does have a sibling on its right
            next = nd->_right_sib;
        }
        else if (nd->_left_child && !nd->_right_sib) {
            // nd has children (it is an internal node) but no siblings on its right
            next = nd->_left_child;
        }
        else {
            // nd has both children and siblings on its right
            next = nd->_left_child;
        }
        return next;
    }

    inline void Particle::recordJoinInfo(vector<G::join_info_t> & join_info) const {
        // Visit all internal nodes, for each recording the tuple:
        // height, is_species_tree, left child's species, right child's species.
        join_info.clear();

        // Build preorders vector
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        for (auto & preorder : preorders) {
            for (auto nd : preorder) {
                if (nd->getLeftChild()) {
                    // internal node
                    join_info.push_back({
                        nd->getHeight(),
                        _is_species_tree,
                        nd->getLeftChild()->getSpecies(),
                        nd->getLeftChild()->getRightSib()->getSpecies()
                    });
                }
            }
        }
                
        // Sort by increasing height above leaf level
        sort(join_info.begin(), join_info.end());
    }
    
    inline void Particle::operator=(const Particle & other) {
        _height                     = other._height;
        _next_node_number           = other._next_node_number;
        _index                      = other._index;
        _nleaves                    = other._nleaves;
        _is_species_tree            = other._is_species_tree;
        _returned_node_numbers      = other._returned_node_numbers;
        
        // Create node map: if _nodes[3]._number = 2, then node_map[2] = 3 (i.e. node number 2 is at index 3 in _nodes vector)
        map<unsigned, unsigned> node_map;
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            int other_node_number = other._nodes[i]._number;
            if (other_node_number > -1) {
                node_map[other_node_number] = i;
            }
        }
        
        _nodes.resize(other._nodes.size());
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            unsigned j = 0;
            int n = 0;
            
            // _left_child
            _nodes[i]._left_child = nullptr;
            if (other._nodes[i]._left_child) {
                n = other._nodes[i]._left_child->_number;
                assert(n > -1);
                j = node_map[n];
                _nodes[i]._left_child = &_nodes[j];
            }
            
            // _right_sib
            _nodes[i]._right_sib = nullptr;
            if (other._nodes[i]._right_sib) {
                n = other._nodes[i]._right_sib->_number;
                assert(n > -1);
                j = node_map[n];
                _nodes[i]._right_sib = &_nodes[j];
            }

            // _parent
            _nodes[i]._parent = nullptr;
            if (other._nodes[i]._parent) {
                n = other._nodes[i]._parent->_number;
                assert(n > -1);
                j = node_map[n];
                _nodes[i]._parent = &_nodes[j];
            }
            
            _nodes[i]._number      = other._nodes[i]._number;
            _nodes[i]._name        = other._nodes[i]._name;
            _nodes[i]._edge_length = other._nodes[i]._edge_length;
            _nodes[i]._height      = other._nodes[i]._height;
            _nodes[i]._species     = other._nodes[i]._species;
            _nodes[i]._split       = other._nodes[i]._split;
            _nodes[i]._flags       = other._nodes[i]._flags;
                        
            // _nodes[i]._partial copied by GParticle::operator=
        }

        // Build _lineages
        _lineages.clear();
        for (auto & nd : _nodes) {
            // Assume that nodes that are being used have _number > -1
            // and nodes serving as the root of a lineage have no parent
            if (nd._number > -1 && !nd._parent)
                _lineages.push_back(&nd);
        }
        assert(_lineages.size() == other._lineages.size());
    }

    inline void Particle::storeEdgelensBySplit(map<Split, double> & edgelenmap) const {
        edgelenmap.clear();

        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
            Node * nd2 = nd;
            while (nd2) {
                edgelenmap[nd2->_split] = nd2->_edge_length;
                nd2 = findNextPreorder(nd2);
            }
        }
    }

    inline void Particle::storeSplits(Split::treeid_t & splitset) const {
        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
            Node * nd2 = nd;
            while (nd2) {
                if (nd2->_left_child) {
                    // add this internal node's split to splitset
                    splitset.insert(nd2->_split);
                }
                nd2 = findNextPreorder(nd2);
            }
        }
    }
    
    inline void Particle::refreshSplits() {
        // Resize all splits (also clears all splits)
        for (auto & nd : _nodes) {
            nd._split.resize(_nleaves);
        }

        // Now do a postorder traversal and add the bit corresponding
        // to the current node in its parent node's split
        if (_lineages.size() == 0)
            return;
            
        vector<Node::ptr_vect_t> preorders;
        buildPreordersVector(preorders);
        
        for (auto & preorder : preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child == nullptr) {
                    // set bit corresponding to this leaf node's number
                    nd->_split.setBitAt(nd->_number);
                }

                if (nd->_parent) {
                    // parent's bits are the union of the bits set in all its children
                    nd->_parent->_split.addSplit(nd->_split);
                }
            }
        }
    }

    inline pair<double,double> Particle::calcTreeDistances(const Particle & ref, const Particle & test) {
        unsigned nleaves = ref.getNLeaves();
        assert(nleaves == test.getNLeaves());
        
        // Store splits from reference tree
        Split::treeid_t ref_splits;
        ref.storeSplits(ref_splits);
        
        // Store edge lengths from reference tree
        map<Split, double> ref_edgelen_map;
        ref.storeEdgelensBySplit(ref_edgelen_map);
                
        // Store splits from test tree
        Split::treeid_t test_splits;
        test.storeSplits(test_splits);
        
        // Store edge lengths from test tree
        map<Split, double> test_edgelen_map;
        test.storeEdgelensBySplit(test_edgelen_map);
                
        // Now calculate squares for leaf nodes, storing in KLleaves
        std::vector<double> KLleaves(nleaves);
        Split s;
        s.resize(nleaves);
        Split sroot;
        sroot.resize(nleaves);
        for (unsigned i = 0; i < nleaves; i++) {
            s.clear();
            s.setBitAt(i);
            sroot.setBitAt(i);
            assert(ref_edgelen_map.count(s) == 1);
            assert(test_edgelen_map.count(s) == 1);
            double ref_edge_length  = ref_edgelen_map[s];
            double test_edge_length = test_edgelen_map[s];
            double square = pow(test_edge_length - ref_edge_length, 2.0);
            KLleaves[i] = square;
        }
        
        // Store union of refsplits and testsplits in allsplits
        Split::treeid_t all_splits;
        set_union(
            ref_splits.begin(), ref_splits.end(),
            test_splits.begin(), test_splits.end(),
            std::inserter(all_splits, all_splits.begin()));
        
        // Traverse allsplits, storing squared branch length differences in KLinternals
        std::vector<double> KLinternals(all_splits.size());
        double RFdist = 0.0;
        unsigned i = 0;
        for (auto s : all_splits) {
            if (s == sroot)
                continue;
            bool s_in_ref  = ref_edgelen_map.count(s) == 1;
            bool s_in_test = test_edgelen_map.count(s) == 1;
            assert(s_in_ref || s_in_test);
            if (!s_in_ref) {
                double test_edge_length = test_edgelen_map[s];
                double test_square = pow(test_edge_length, 2.0);
                KLinternals[i++] = test_square;
                RFdist += 1.0;
            }
            else if (!s_in_test) {
                double ref_edge_length = ref_edgelen_map[s];
                double ref_square = pow(ref_edge_length, 2.0);
                KLinternals[i++] = ref_square;
                RFdist += 1.0;
            }
            else {
                double test_edge_length = test_edgelen_map[s];
                double ref_edge_length = ref_edgelen_map[s];
                double square = pow(test_edge_length - ref_edge_length, 2.0);
                KLinternals[i++] = square;
            }
        }
            
        // Calculate KL distance
        double KFSS = 0.0;
        for (auto square : KLinternals) {
            KFSS += square;
        }
        for (auto square : KLleaves) {
            KFSS += square;
        }
        assert(KFSS >= 0.0);
        double KFdist = sqrt(KFSS);
        return make_pair(KFdist, RFdist);
    }
    
}

