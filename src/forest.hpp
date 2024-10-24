#pragma once

using boost::algorithm::join;
extern proj::PartialStore ps;

class Particle;

extern proj::Lot::SharedPtr rng;

namespace proj {

    class Forest {
    
        friend class Particle;
        
        public:
        
            typedef tuple<double, unsigned, vector<G::species_t> >  coalinfo_t;
            typedef set<Split> treeid_t;

                                Forest();
                                ~Forest();
                
            virtual void        clear();
            double              getPrevHeight() const;
            double              getHeight() const;
            unsigned            getNumLineages() const;
            
            Node *              findNodeNumbered(unsigned num);
                        
            static void         readTreefile(const string filename,
                                        unsigned skip,
                                        vector<string> & leaf_names,
                                        map<unsigned,unsigned> & taxon_map,
                                        vector<string> & names,
                                        vector<string> & newicks);
                                        
            string              makeNewick(unsigned precision = 6,
                                        bool use_names = true,
                                        bool coalunits = false) const;
            void                buildFromNewick(const string newick);
            void                storeSplits(Forest::treeid_t & splitset);
            //void                alignTaxaUsing(const map<unsigned, unsigned> & taxon_map);
            
            unsigned            advanceAllLineagesBy(double dt);
            void                heightsInternalsPreorders();

            double              mrcaHeight(G::species_t lspp, G::species_t rspp) const;

            virtual void        createTrivialForest(bool compute_partials) = 0;
            virtual bool        isSpeciesForest() const = 0;
            virtual void        recordHeights(vector<double> & height_vect) const = 0;
            virtual void        setSpeciesFromNodeName(Node * nd) = 0;

            virtual void        addCoalInfoElem(const Node *, vector<coalinfo_t> & recipient) = 0;
            virtual void        saveCoalInfo(vector<coalinfo_t> & coalinfo_vect, bool cap = false) const = 0;
            
            static void         debugShowCoalInfo(string title, vector<coalinfo_t> & coalinfo_vect, string fn = "");
            static pair<double,double> calcTreeDistances(Forest & ref, Forest & test);
            static bool         subsumed(G::species_t test_species, G::species_t subtending_species);

            void        debugShowPreorder() const;
            void        debugShowLineages() const;
            
            Node::ptr_vect_t & getLineages();

            void        removeOne(Node::ptr_vect_t & node_vect, Node * del);
            void        addOne(Node::ptr_vect_t & node_vect, Node * add);
            void        removeTwoAddOne(Node::ptr_vect_t & node_vect, Node * del1, Node * del2, Node * add);
            void        addTwoRemoveOne(Node::ptr_vect_t & node_vect, Node * del1, Node * del2, Node * add);
            void        addTwoRemoveOneAt(Node::ptr_vect_t & node_vect, unsigned pos1, Node * del1, unsigned pos2, Node * del2, Node * add);

            void        refreshAllPreorders() const;

        protected:
        
            int         extractNodeNumberFromName(string node_name, set<unsigned> & used);
            void        extractEdgeLen(Node * nd, string edge_length_string);
            void        setNodeNameFromNumber(Node * nd);
            Node *      findNextPreorder(Node * nd) const;
            bool        canHaveSibling(Node * nd) const;
            void        refreshPreorder(Node::ptr_vect_t & preorder) const;
            void        storeEdgelensBySplit(map<Split, double> & edgelenmap);
            void        resizeAllSplits(unsigned nleaves);
            void        refreshAllHeightsAndPreorders();
            Node *      pullNode();
            void        stowNode(Node * nd);
            void        joinLineagePair(Node * anc, Node * first, Node * second);
            void        unjoinLineagePair(Node * anc, Node * first, Node * second);
            
            void        debugCheckCoalInfoSorted(const vector<coalinfo_t> & coalinfo_vect) const;
            void        debugCheckAllPreorders() const;
            void        debugCheckPreorder(const Node::ptr_vect_t & preorder) const;
            void        debugShowNodeInfo(string title);
            void        debugShowNodePtrVector(Node::ptr_vect_t& v, string title);

            virtual void operator=(const Forest & other);

            // NOTE: any variables added must be copied in operator=
            
            vector<unsigned>        _unused_nodes;

            double                  _forest_height;
            //unsigned                _next_node_index;
            double                  _log_likelihood;
            double                  _prev_log_likelihood;
            double                  _log_prior;
            vector<Node>            _nodes;
            Node::ptr_vect_t        _lineages;
            vector<coalinfo_t>      _coalinfo;
                                                
            // Because these can be recalculated at any time, they should not
            // affect the const status of the Forest object
            mutable unsigned                 _next_node_number;
            mutable vector<Node::ptr_vect_t> _preorders;
};
    
    inline Forest::Forest() {
    }

    inline Forest::~Forest() {
        clear();
    }

    inline void Forest::clear() {
        _unused_nodes.clear();
        _forest_height = 0.0;
        //_next_node_index = 0;
        _next_node_number = 0;
        _log_likelihood = 0.0;
        _prev_log_likelihood = 0.0;
        _log_prior = 0.0;
        _nodes.clear();
        _preorders.clear();
        _lineages.clear();
        _coalinfo.clear();
    }
    
    inline Node::ptr_vect_t & Forest::getLineages() {
        return _lineages;
    }
    
    inline double Forest::getHeight() const {
        return _forest_height;
    }
    
    inline unsigned Forest::getNumLineages() const {
        return (unsigned)_lineages.size();
    }
    
    inline int Forest::extractNodeNumberFromName(string node_name, set<unsigned> & used) {
        // Attempts to convert node_name to a number and throws exception if node's
        // name cannot be converted to an integer or if node number has already
        // been encountered previously
        bool success = true;
        int x;
        try {
            x = stoi(node_name);
        }
        catch(invalid_argument &) {
            // node_name could not be converted to an integer value
            success = false;
        }

        if (success) {
            // conversion succeeded
            // attempt to insert x into the set of node numbers already used
            pair<set<unsigned>::iterator, bool> insert_result = used.insert(x);
            if (insert_result.second) {
                // insertion was made, so x has NOT already been used
                return x;
            }
            else {
                // insertion was not made, so set already contained x
                throw XProj(str(format("leaf number %d used more than once") % x));
            }
        }
        else
            throw XProj(str(format("node name (%s) not interpretable as an integer") % node_name));
    }
    
    inline void Forest::extractEdgeLen(Node * nd, string edge_length_string) {
        assert(nd);
        bool success = true;
        double d = 0.0;
        try {
            d = stod(edge_length_string);
        }
        catch(invalid_argument &) {
            // edge_length_string could not be converted to a double value
            success = false;
        }

        if (success) {
            // conversion succeeded
            nd->setEdgeLength(d);
        }
        else
            throw XProj(str(format("%s is not interpretable as an edge length") % edge_length_string));
    }

    inline void Forest::setNodeNameFromNumber(Node * nd) {
        assert(nd->_number > -1);
        unsigned n = nd->_number;
        if (isSpeciesForest()) {
            assert(n < 2*G::_species_names.size() - 1);
            if (n < G::_species_names.size())
                nd->_name = G::_species_names[n];
        }
        else {
            assert(n < 2*G::_taxon_names.size() - 1);
            if (n < G::_taxon_names.size())
                nd->_name = G::_taxon_names[n];
        }
    }
    
    inline Node * Forest::findNextPreorder(Node * nd) const {
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
    
    inline bool Forest::canHaveSibling(Node * nd) const {
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
    
    inline string Forest::makeNewick(unsigned precision, bool use_names, bool coalunits) const  {
        refreshAllPreorders();
        //refreshAllHeightsAndPreorders();
        // Place basal polytomy (if there is one) at a height 10% greater than
        // the _forest_height
        double basal_polytomy_height = 0.0; //_forest_height*(1.1);
        bool is_complete_tree = (bool)(_preorders.size() == 1);
        double coalfactor = (coalunits ? (2.0/G::_theta) : 1.0);
        
        //const format basal_subtree_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_name_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_number_format( str(format("%%d:%%.%df") % precision) );
        const format internal_node_format( str(format("):%%.%df") % precision) );
        
        vector<string> subtree_newicks;
        for (auto preorder : _preorders) {
            string subtree_newick;
            stack<Node *> node_stack;
            for (auto nd : preorder) {
                assert(nd->_number > -1);
                if (nd->_left_child) {
                    subtree_newick += "(";
                    node_stack.push(nd);
                }
                else {
                    double edge_length = nd->_edge_length*coalfactor;
                    if (!nd->_parent) {
                        // Subtree consists of just this one leaf node
                        edge_length += basal_polytomy_height*coalfactor;
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
                        double popped_edge_length = popped->_edge_length*coalfactor;
                        while (popped && !popped->_right_sib) {
                            node_stack.pop();
                            if (node_stack.empty()) {
                                if (is_complete_tree) {
                                    subtree_newick += ")";
                                }
                                else {
                                    // This is the root of one of several subtrees, so
                                    // it is important to preserve its edge length
                                    popped_edge_length += basal_polytomy_height*coalfactor;
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
                                popped_edge_length = popped->_edge_length*coalfactor;
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
            //if (subtree_newick[0] == '(')
            //    subtree_newicks.push_back(str(format(basal_subtree_format) % subtree_newick % basal_polytomy_height));
            //else
                subtree_newicks.push_back(subtree_newick);
        }
        
        string newick;
        if (is_complete_tree)
            newick = subtree_newicks[0];
        else {
            string separator = str(format(":%.5f,") % basal_polytomy_height);
            string insides = join(subtree_newicks, ",");
            newick = str(format("(%s)") % insides);
        }

        return newick;
    }

    //inline void Forest::alignTaxaUsing(const map<unsigned, unsigned> & taxon_map) {
    //    for (auto & preorder : _preorders) {
    //        for (auto nd : boost::adaptors::reverse(preorder)) {
    //            if (!nd->_left_child) {
    //                // leaf node
    //                if (taxon_map.count(nd->_number + 1) == 0) {
    //                    output("\ntaxon_map:\n",1);
    //                    for (auto t : taxon_map) {
    //                        output(format("  key: %d --> value: %d\n") % t.first % t.second,1);
    //                    }
    //                    throw XProj(format("%d is not a key in taxon_map") % (nd->_number + 1));
    //                }
    //                else {
    //                    unsigned leaf_index = taxon_map.at(nd->_number + 1);
    //                    nd->_number = leaf_index;
    //                    nd->_name = G::_species_names[leaf_index];
    //                }
    //            }
    //        }
    //    }
    //}
    
    inline void Forest::buildFromNewick(const string newick) {
        // Builds strictly-bifurcating ultrametric, rooted, complete (i.e.
        // _lineages contains just one element) forest from the supplied newick
        // tree description. Assumes newick either uses names for leaves or,
        // if it specifies numbers, the numbers correspond to keys in
        // G::_nexus_taxon_map, which translates taxon numbers in newick
        // strings to the index of the taxon in G::_species_names (if
        // building a species tree) or G::_taxon_names (if building a
        // gene tree).
        
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
        
        // Set _number for each node to its index in _nodes vector
        _unused_nodes.clear();
        for (int i = 0; i < _nodes.size(); i++) {
            _nodes[i]._number = -1;
            _nodes[i]._my_index = i;
            _unused_nodes.push_back(i);
        }
        
        try {
            // Root node
            Node * nd   = pullNode();
            Node * root = nd;
            
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
                            // leaf node
                            int num = extractNodeNumberFromName(nd->_name, used);
                            assert(num > 0);
                            nd->_number = G::_nexus_taxon_map[(unsigned)num];
                            setNodeNameFromNumber(nd);
                            setSpeciesFromNodeName(nd);
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
                            // leaf node
                            int num = extractNodeNumberFromName(nd->_name, used);
                            assert(num > 0);
                            nd->_number = G::_nexus_taxon_map[(unsigned)num];
                            setNodeNameFromNumber(nd);
                            setSpeciesFromNodeName(nd);
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
                        extractEdgeLen(nd, edge_length_str);
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
                        nd->_right_sib = pullNode();
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
                        nd->_left_child = pullNode();
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
                throw XProj(str(format("Expecting single quote to end node name at position %d in tree description") % node_name_position));

            heightsInternalsPreorders();
        }
        catch(XProj x) {
            clear();
            throw x;
        }
    }
        
    inline void Forest::debugCheckPreorder(const Node::ptr_vect_t & preorder) const {
        unsigned which = 0;
        Node * nd = preorder[which++];
        while (true) {
            nd = findNextPreorder(nd);
            if (nd) {
                assert(preorder.size() > which);
                assert(nd == preorder[which++]);
            }
            else
                break;
        }
    }
    
    inline void Forest::refreshPreorder(Node::ptr_vect_t & preorder) const {
        // Assumes preorder just contains the root node when this function is called
        // Also assumes that _next_node_number was initialized prior to calling this function
        assert(preorder.size() == 1);
        
        Node * nd = preorder[0];
        while (true) {
            nd = findNextPreorder(nd);
            if (nd) {
                preorder.push_back(nd);
                if (nd->_left_child)
                    nd->_number = _next_node_number++;
            }
            else
                break;
        }
    }
    
    inline void Forest::storeEdgelensBySplit(map<Split, double> & edgelenmap) {
        edgelenmap.clear();
        for (auto & preorder : _preorders) {
            for (auto nd : preorder) {
                edgelenmap[nd->_split] = nd->_edge_length;
            }
        }
    }

    Node * Forest::findNodeNumbered(unsigned num) {
        Node * found = nullptr;
        for (auto & preorder : _preorders) {
            for (auto nd : preorder) {
                if (nd->_number == num) {
                    found = nd;
                }
            }
        }
        if (!found)
            throw XProj(format("Could not find node numbered %d in forest") % num);
        return found;
    }

    inline void Forest::storeSplits(Forest::treeid_t & splitset) {
        // Resize all splits (also clears all splits)
        resizeAllSplits(isSpeciesForest() ? G::_nspecies : G::_ntaxa);

        // Now do a postorder traversal and add the bit corresponding
        // to the current node in its parent node's split
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // add this internal node's split to splitset
                    splitset.insert(nd->_split);
                }
                else {
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

    inline void Forest::resizeAllSplits(unsigned nleaves) {
        for (auto & nd : _nodes) {
            nd._split.resize(nleaves);
        }
    }
    
    inline void Forest::debugShowPreorder() const {
        assert(_lineages.size() > 0);
        unsigned sub = 0;
        for (auto preorder : _preorders) {
            output(format("\nSubtree %d\n") % sub, 1);
            for (auto nd : preorder) {
                output(format("  node %d:\n") % nd->_number, 1);
                output(format("    name    = \"%s\"\n") % nd->_name, 1);
                output(format("    height  = %.5f\n") % nd->_height, 1);
                output(format("    edgelen = %.5f\n") % nd->_edge_length, 1);
                output(format("    par     = %s\n") % (nd->_parent ? to_string(nd->_parent->_number) : "none"), 1);
                output(format("    lchild  = %s\n") % (nd->_left_child ? to_string(nd->_left_child->_number) : "none"), 1);
                output(format("    rchild  = %s\n") % (nd->_right_sib ? to_string(nd->_right_sib->_number) : "none"), 1);
            }
            sub++;
        }
    }

    inline void Forest::debugShowLineages() const {
        assert(_lineages.size() > 0);
        output("\n_lineages:\n", 0);
        for (auto nd : _lineages) {
            output(format("  %d \"%s\" (%s)\n") % nd->_number % nd->_name % G::memoryAddressAsString(nd), 0);
            if (nd->_left_child) {
                output(format("    left child: %d \"%s\" (%s)\n") % nd->_left_child->_number % nd->_left_child->_name % G::memoryAddressAsString(nd->_left_child), 0);
            }
                
            if (nd->_right_sib) {
                output(format("    right sib: %d \"%s\" (%s)\n") % nd->_right_sib->_number % nd->_right_sib->_name %  G::memoryAddressAsString(nd->_right_sib), 0);
            }
                
            if (nd->_parent) {
                output(format("    parent: %d \"%s\" (%s)\n") % nd->_parent->_number % nd->_parent->_name %  G::memoryAddressAsString(nd->_parent), 0);
            }
        }
    }

    inline void Forest::debugCheckAllPreorders() const {
        // For each subtree stored in _lineages, check to ensure that the
        // node pointers stored in _preorders are indeed in preorder sequence
        if (_lineages.size() == 0)
            return;
        
        unsigned which_lineage = 0;
        for (auto nd : _lineages) {
            if (_preorders[which_lineage][0] != nd) {
                output("Forest::debugCheckAllPreorders found a problem\n" , 1);
                output(format("  Node %d (name %s) not equal to _preorders[%d]\n") % nd->_number % nd->_name % which_lineage, 1);
                
                // Print out every lineage root node in _preorders
                unsigned i = 0;
                output(format("%12s %12s %s\n") % "index" % "number" % "name", 1);
                for (const Node::ptr_vect_t & v : _preorders) {
                    output(format("%12d %12d %s\n") % (i++) % v[0]->_number % v[0]->_name, 1);
                }
                throw XProj("Aborting from within Forest::debugCheckAllPreorders");
            }
            
            // Now check that the nodes above nd are in preorder sequence
            debugCheckPreorder(_preorders[which_lineage]);
            
            which_lineage++;
        }
    }
    
    inline void Forest::refreshAllPreorders() const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        _next_node_number = isSpeciesForest() ? G::_nspecies : G::_ntaxa;
        _preorders.clear();
        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
            if (nd->_left_child) {
                nd->_number = _next_node_number++;
            }
            
            // lineage is a Node::ptr_vect_t (i.e. vector<Node *>)
            // lineage[0] is the first node pointer in the preorder sequence for this lineage
            // Add a new vector to _preorders containing, for now, just the root of the subtree
            _preorders.push_back({nd});
            
            // Now add the nodes above the root in preorder sequence
            Node::ptr_vect_t & preorder_vector = *_preorders.rbegin();
            refreshPreorder(preorder_vector);
        }
    }
    
    inline void Forest::refreshAllHeightsAndPreorders() {
        // Refreshes _preorders and then recalculates heights of all nodes
        // and sets _forest_height
        
        resizeAllSplits(isSpeciesForest() ? G::_nspecies : G::_ntaxa);
        refreshAllPreorders();
        _forest_height = 0.0;

        // Set heights for each lineage in turn
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                assert(nd->_number > -1);
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != G::_infinity);
                    if (nd->_height + nd->_edge_length > _forest_height)
                        _forest_height = nd->_height + nd->_edge_length;
                }
                else {
                    // nd is a leaf node
                    nd->_height = 0.0;
                    if (nd->_edge_length > _forest_height)
                        _forest_height = nd->_edge_length;

                    // Set bit corresponding to this leaf node's number
                    nd->_split.setBitAt(nd->_number);
                }
                
                if (nd->_parent) {
                    // Set parent's height if nd is right-most child of its parent
                    bool is_rightmost_child = !nd->_right_sib;
                    double parent_height = nd->_height + nd->_edge_length;
                    if (is_rightmost_child) {
                        nd->_parent->_height = parent_height;
                    }
                    
                    // If nd is not its parent's rightmost child, check ultrametric assumption
                    assert(!is_rightmost_child || fabs(nd->_parent->_height - parent_height) < G::_small_enough);

                    // Parent's bits are the union of the bits set in all its children
                    nd->_parent->_split.addSplit(nd->_split);
                }
            }
        }
    }
    
    inline void Forest::heightsInternalsPreorders() {
        // Recalculates heights of all nodes, refreshes preorder sequences, and renumbers internal nodes
        assert(_lineages.size() == 1);
        refreshAllPreorders();
        _forest_height = 0.0;

        // First internal node number is the number of leaves
        int next_node_number = (unsigned)(isSpeciesForest() ? G::_species_names.size() : G::_taxon_names.size());

        // Renumber internal nodes in postorder sequence for each lineage in turn
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != G::_infinity);
                    if (nd->_height > _forest_height)
                        _forest_height = nd->_height;
                    nd->_number = next_node_number++;
                    nd->_species = nd->_left_child->_species;
                    assert(nd->_left_child->_right_sib);
                    nd->_species |= nd->_left_child->_right_sib->_species;
                    assert(nd->_left_child->_right_sib->_right_sib == nullptr);
                    //addCoalInfoElem(nd);
                }
                else {
                    // nd is a leaf node
                    assert(nd->_number > -1);
                    nd->_height = 0.0;
                }
                                
                if (nd->_parent) {
                    // Set parent's height if nd is right-most child of its parent
                    bool is_rightmost_child = !nd->_right_sib;
                    double parent_height = nd->_height + nd->_edge_length;
                    if (is_rightmost_child) {
                        nd->_parent->_height = parent_height;
                    }
                    
                    // If nd is not its parent's rightmost child, check ultrametric assumption
                    assert(!is_rightmost_child || fabs(nd->_parent->_height - parent_height) < G::_small_enough);
                }
            }
        }
    }
        
    inline void Forest::readTreefile(const string filename,
                                     unsigned skip,
                                     vector<string> & leaf_names,
                                     map<unsigned,unsigned> & taxon_map,
                                     vector<string> & tree_names,
                                     vector<string> & newicks) {
        // Reads file, skipping the first skip trees in each trees block and
        // storing newick tree descriptions (in newicks) and tree names (in tree_names).
        // User-supplied taxon_map will be filled with entries such that taxon_map[k]
        // provides the index (into leaf_names) of the taxon whose number in the taxa block
        // (and thus also in the newick tree description) equals k. (If the tree file does
        // not supply a taxa block, one will be created in which the ordering of taxa is
        // arbitrary (perhaps the order in which taxa are encountered in the first tree?).

        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

        // Read and process the file
        //MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
        MultiFormatReader nexus_reader(-1, NxsReader::IGNORE_WARNINGS);

        // Both of these needed to suppress "storing read block" messages
        // see NxsReader::statusMessage in nxsreader.cpp
        nexus_reader.SetAlwaysReportStatusMessages(false);
        nexus_reader.SetWarningOutputLevel(NxsReader::SUPPRESS_WARNINGS_LEVEL);
        
        try {
            nexus_reader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...) {
            nexus_reader.DeleteBlocksFromFactories();
            throw;
        }

        // Get number of taxa blocks
        int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
        
        // Process the trees blocks associated with each taxa block
        for (int i = 0; i < num_taxa_blocks; ++i) {
            NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(i);
            string taxa_block_title = taxa_block->GetTitle();
            
            // Ensure that the number of leaves accords with expectations
            unsigned num_leaves = (unsigned)taxa_block->GetNumTaxonLabels();
            if (num_leaves != leaf_names.size()) {
                if (leaf_names.empty()) {
                    // Populate leaf_names using taxa block
                    leaf_names.resize(num_leaves);
                    for (unsigned t = 0; t < num_leaves; t++) {
                        leaf_names[t] = taxa_block->GetTaxonLabel(t);
                    }
                    output(format("Assuming %d taxon names are those found in taxa block") % num_leaves, 1);
                }
                else {
                    throw XProj(format("Found %d taxa in taxa block %d but expecting %d") % num_leaves % (i+1) % leaf_names.size());
                }
            }
            
            // It is possible (probable even) that the order of taxa in taxa_block will differ from
            // the order in leaf_names. Therefore, taxon_map is constructed the keys of which are the
            // taxon numbers in the newick tree description (corresponding to taxa_block); the values
            // are the index of the taxon in leaf_names. For example (branch lengths omitted):
            //
            // #NEXUS
            // begin trees;
            //   tree tree1 = [&R] (A,(((D,B),C),E));
            // end;
            //
            // Supplied to readTreeFile function:
            //   leaf_names = {"A","B","C","D","E"}
            //                  0   1   2   3   4
            //
            // Populated by readTreeFile function:
            //   tree_names = {"tree1"}
            //   newicks = {"(1,(((2,3),4),5))"}
            //   taxon_map: (newick)    (leaf_names)
            //                key          value
            //                 1    -->      0 (i.e. "A")
            //                 2    -->      3 (i.e. "D")
            //                 3    -->      1 (i.e. "B")
            //                 4    -->      2 (i.e. "C")
            //                 5    -->      4 (i.e. "E")
            
            taxon_map.clear();
            for (unsigned t = 0; t < num_leaves; t++) {
                string taxon_label = taxa_block->GetTaxonLabel(t);
                auto it = find(leaf_names.begin(), leaf_names.end(), taxon_label);
                if (it == leaf_names.end()) {
                    throw XProj(format("Taxon label \"%s\" was not expected ") % taxon_label);
                }
                unsigned d = (unsigned)(distance(leaf_names.begin(), it));
                taxon_map[t+1] = d;
            }
            
            // Get number of trees blocks for this taxa block
            const unsigned num_trees_blocks = nexus_reader.GetNumTreesBlocks(taxa_block);
            
            // Get all trees for each trees block, skipping the first skip trees
            for (unsigned j = 0; j < num_trees_blocks; ++j) {
                // GetTreeName is inexplicably non-const, necessitating
                // making trees_block non-const
                NxsTreesBlock * trees_block = nexus_reader.GetTreesBlock(taxa_block, j);
                unsigned num_trees = trees_block->GetNumTrees();
                if (skip < num_trees) {
                    for (unsigned t = skip; t < num_trees; ++t) {
                        const NxsFullTreeDescription & d = trees_block->GetFullTreeDescription(t);

                        if (!d.IsRooted()) {
                            throw XProj(format("Tree %d in trees block %d and taxa block %d is unrooted; expecting only rooted trees") % (t+1) % (j+1) % (i+1));
                        }
                        
                        // Store the tree name
                        string tree_name = trees_block->GetTreeName(t);
                        tree_names.push_back(tree_name);
                        
                        // Store the newick tree description
                        string newick = d.GetNewick();
                        newicks.push_back(newick);
                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

        // No longer any need to store raw data from nexus file
        nexus_reader.DeleteBlocksFromFactories();
    }
            
    inline unsigned Forest::advanceAllLineagesBy(double dt) {
        // Add t to the edge length of all lineage root nodes, unless there
        // is just one lineage, in which case do nothing
        unsigned n = (unsigned)_lineages.size();
        if (n > 1) {
            for (auto nd : _lineages) {
                double elen = nd->getEdgeLength() + dt;
                assert(elen >= 0.0 || fabs(elen) < Node::_smallest_edge_length);
                nd->setEdgeLength(elen);
                ++n;
            }
        
            // Add to to the current forest height
            _forest_height += dt;
        }
        
        return n;
    }
    
    inline Node * Forest::pullNode() {
        //if (_next_node_index == _nodes.size()) {
        if (_unused_nodes.empty()) {
            unsigned nleaves = 2*(isSpeciesForest() ? G::_nspecies : G::_ntaxa) - 1;
            //throw XProj(str(format("Forest::pullNode tried to return a node beyond the end of the _nodes vector (_next_node_index = %d equals %d nodes allocated for %d leaves)") % _next_node_index % _nodes.size() % nleaves));
            throw XProj(str(format("Forest::pullNode tried to return a node beyond the end of the _nodes vector (%d nodes allocated for %d leaves)") % _nodes.size() % nleaves));
        }
        double node_index = _unused_nodes.back();
        _unused_nodes.pop_back();
        
        //Node * new_nd = &_nodes[_next_node_index++];
        Node * new_nd = &_nodes[node_index];
        assert(new_nd->getMyIndex() == node_index);
        
        assert(!new_nd->_partial);
        new_nd->clear();
        new_nd->_number = -2;
        new_nd->_split.resize(G::_ntaxa);
        return new_nd;
    }
    
    inline void Forest::stowNode(Node * nd) {
        // Get index of nd in _nodes vector
        int offset = nd->_my_index;
        assert(offset > -1);
        _unused_nodes.push_back(offset);
        //_next_node_index--;
        //_next_node_number--;
        //assert(nd == &_nodes[_next_node_index]);
        nd->clear();
    }
    
    inline void Forest::joinLineagePair(Node * new_nd, Node * subtree1, Node * subtree2) {
        // Note: must call pullNode to obtain new_nd before calling this function
        assert(new_nd);
        assert(subtree1);
        assert(subtree2);
        assert(new_nd->_my_index > -1);
        new_nd->_name        = "anc-" + to_string(new_nd->_my_index);
        new_nd->_left_child  = subtree1;
        new_nd->_edge_length = 0.0;
        new_nd->_species = 0;
                
        // Calculate height of the new node
        // (should be the same whether computed via the left or right child)
        double h1 = subtree1->_height + subtree1->_edge_length;
        double h2 = subtree2->_height + subtree2->_edge_length;

        // //temporary! Probably need to reinstate this assert
        //assert(fabs(h1 - h2) < G::_small_enough);

        new_nd->_height = (h1 + h2)/2.0;

        // Finish connecting new trio of nodes
        subtree1->_right_sib = subtree2;
        subtree1->_parent    = new_nd;
        subtree2->_parent    = new_nd;
    }
    
    inline void Forest::unjoinLineagePair(Node * anc, Node * subtree1, Node * subtree2) {
        // Note: be sure to call stowNode for anc after calling this function
        // Reset members set by joinLineagePair function
        anc->_number = -1;
        anc->_name = "";
        anc->_left_child = nullptr;
        anc->_edge_length = 0.0;
        anc->_height = 0.0;
        anc->_species = 0;
        //anc->_partial = nullptr;
        subtree1->_right_sib = nullptr;
        subtree1->_parent = nullptr;
        subtree2->_parent = nullptr;
    }
    
    inline void Forest::removeOne(Node::ptr_vect_t & v, Node * del) {
        // Get iterator to node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), del);
        assert(it != v.end());
        v.erase(it);
    }
    
    inline void Forest::addOne(Node::ptr_vect_t & v, Node * add) {
        // Add node
        v.push_back(add);
    }
    
    inline void Forest::removeTwoAddOne(Node::ptr_vect_t & v, Node * del1, Node * del2, Node * add) {
        // Get iterator to first node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), del1);
        assert(it != v.end());
        v.erase(it);
        
        // Get iterator to second node to be deleted and remove from v
        it = find(v.begin(), v.end(), del2);
        assert(it != v.end());
        v.erase(it);
        
        // Add node
        v.push_back(add);
    }
    
    inline void Forest::addTwoRemoveOne(Node::ptr_vect_t & v, Node * add1, Node * add2, Node * rem) {
        // Get iterator to node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), rem);
        assert(it != v.end());
        v.erase(it);
        
        // Add node2
        v.push_back(add1);
        v.push_back(add2);
    }
    
    inline void Forest::addTwoRemoveOneAt(Node::ptr_vect_t & v, unsigned pos1, Node * add1, unsigned pos2, Node * add2, Node * rem) {
        // Get iterator to node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), rem);
        assert(it != v.end());
        v.erase(it);
        
        // Insert add1 and add2 into v so that they end up in
        // positions pos1 and pos2, respectively
        if (pos1 < pos2) {
            it = v.begin();
            advance(it, pos1);
            v.insert(it, add1);
            
            it = v.begin();
            advance(it, pos2);
            v.insert(it, add2);
        }
        else {
            it = v.begin();
            advance(it, pos2);
            v.insert(it, add2);
            
            it = v.begin();
            advance(it, pos1);
            v.insert(it, add1);
        }
    }
    
    inline void Forest::debugShowNodePtrVector(Node::ptr_vect_t & v, string title) {
        output(format("\nNode::ptr_vect_t: %s\n") % title, 2);
        output(format("%6s %12s %6s %10s\n") % "index" % "address" % "number" % "species", 2);
        for (unsigned i = 0; i < v.size(); ++i) {
            string memory_address = G::memoryAddressAsString((void *)&_nodes[i]);
            string species_set = G::speciesStringRepresentation(v[i]->getSpecies());
            output(format("%6d %12s %6d %10s\n") % i % memory_address % v[i]->_number % species_set, 2);
        }
        output("\n", 2);
    }

    inline void Forest::debugShowCoalInfo(string title, vector<Forest::coalinfo_t> & coalinfo_vect, string fn) {
        if (fn.size() > 0) {
            ofstream tmpf(fn);
            tmpf << str(format("\n%s:\n") % title);
            tmpf << str(format("%12s %12s %s\n") % "height" % "gene" % "child spp");
            for (auto cinfo : coalinfo_vect) {
                double                                   h = get<0>(cinfo);
                unsigned                       gene_plus_1 = get<1>(cinfo);
                vector<G::species_t> &         spp = get<2>(cinfo);
                //G::species_t sleft       = get<2>(cinfo);
                //G::species_t sright      = get<3>(cinfo);
                
                ostringstream oss;
                copy(spp.begin(), spp.end(), ostream_iterator<G::species_t>(oss, " "));
                if (gene_plus_1 == 0)
                    tmpf << str(format("%12.9f %12s %12s\n") % h % "(0)" % oss.str());
                else
                    tmpf << str(format("%12.9f %12d %12s\n") % h % gene_plus_1 % oss.str());
            }
            tmpf << endl;
            tmpf.close();
        }
        else {
            output(format("\n%s:\n") % title, 1);
            output(format("%12s %12s %s\n") % "height" % "gene" % "child spp", 1);
            for (auto cinfo : coalinfo_vect) {
                double                                   h = get<0>(cinfo);
                unsigned                       gene_plus_1 = get<1>(cinfo);
                vector<G::species_t> &         spp = get<2>(cinfo);
                //G::species_t sleft       = get<2>(cinfo);
                //G::species_t sright      = get<3>(cinfo);
                
                ostringstream oss;
                copy(spp.begin(), spp.end(), ostream_iterator<G::species_t>(oss, " "));
                if (gene_plus_1 == 0)
                    output(format("%12.9f %12s %12s\n") % h % "(0)" % oss.str(), 1);
                else
                    output(format("%12.9f %12d %12s\n") % h % gene_plus_1 % oss.str(), 1);
            }
            output("\n", 1);
        }
    }
    
    inline void Forest::operator=(const Forest & other) {
        _forest_height                  = other._forest_height;
        //_next_node_index                = other._next_node_index;
        _unused_nodes                   = other._unused_nodes;
        _next_node_number               = other._next_node_number;
        _log_likelihood                 = other._log_likelihood;
        _prev_log_likelihood            = other._prev_log_likelihood;
        
        // Create node map: if other._nodes[3]._number = 2, then node_map[2] = 3
        // (i.e. node number 2 is at index 3 in _nodes vector)
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
            _nodes[i]._my_index    = other._nodes[i]._my_index;
            _nodes[i]._name        = other._nodes[i]._name;
            _nodes[i]._edge_length = other._nodes[i]._edge_length;
            _nodes[i]._height      = other._nodes[i]._height;
            _nodes[i]._species     = other._nodes[i]._species;
            _nodes[i]._split       = other._nodes[i]._split;
            _nodes[i]._flags       = other._nodes[i]._flags;
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
        
        _coalinfo = other._coalinfo;
        
        refreshAllPreorders();
    }

    inline bool Forest::subsumed(G::species_t test_species, G::species_t subtending_species) {
        bool not_equal = (test_species != subtending_species);
        bool is_subset = ((test_species & subtending_species) == test_species);
        if (not_equal && is_subset)
            return true;
        else
            return false;
    }

    inline pair<double,double> Forest::calcTreeDistances(Forest & ref, Forest & test) {
        // Determine whether ref and test are species forests or gene forests
        bool ref_is_species_forest  = ref.isSpeciesForest();
        bool test_is_species_forest = test.isSpeciesForest();
        bool species_forest_comparison = ref_is_species_forest && test_is_species_forest;
        bool gene_forest_comparison = !ref_is_species_forest && !test_is_species_forest;
        assert(species_forest_comparison || gene_forest_comparison);
        unsigned nlvs = 0;
        if (species_forest_comparison)
            nlvs = G::_nspecies;
        else
            nlvs = G::_ntaxa;
        
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
        std::vector<double> KLleaves(nlvs);
        Split s;
        s.resize(nlvs);
        Split sroot;
        sroot.resize(nlvs);
        for (unsigned i = 0; i < nlvs; i++) {
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
    
    inline void Forest::debugCheckCoalInfoSorted(const vector<coalinfo_t> & coalinfo_vect) const {
        double h0 = get<0>(coalinfo_vect[0]);
        for (auto & ci : coalinfo_vect) {
            assert(get<0>(ci) >= h0);
        }
    }
            
    inline double Forest::mrcaHeight(G::species_t lspp, G::species_t rspp) const {
        double h = G::_infinity;
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    G::species_t aspp = nd->getSpecies();
                    if ((aspp & lspp) && (aspp & rspp)) {
                        if (nd->_height < h)
                            h = nd->_height;
                    }
                }
            }
        }
        assert(h != G::_infinity);
        return h;
    }

}
