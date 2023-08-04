#pragma once

using boost::algorithm::join;
extern proj::PartialStore ps;

class Particle;

extern proj::Lot rng;

namespace proj {

    class Forest {
    
        friend class Particle;
        
        public:
        
                Forest();
                ~Forest();
                
            double getHeight() const {return _forest_height;}

            virtual void clear();

            virtual void copyEpochsFrom(const epoch_list_t & other);

            double calcLogCoalescentLikelihood(const epoch_list_t & epochs, unsigned gene) const;
            
            static void debugShow(const format & f);
            
            static string unsignedVectToString(const vector<unsigned> & v);
            static string lineageCountsMapToString(const Epoch::lineage_counts_t & lcm);
            static void createDefaultGeneTreeNexusTaxonMap();
            static void createDefaultSpeciesTreeNexusTaxonMap();
            
            void debugCheckSpeciesSets() const;

            static void readTreefile(const string filename,
                                        unsigned skip,
                                        const vector<string> & leaf_names,
                                        map<unsigned,unsigned> & taxon_map,
                                        vector<string> & names,
                                        vector<string> & newicks);
                                        
            double  calcMaxT(epoch_list_t & epochs, double h0, set<Node::species_t> current_species);
            string  makeNewick(unsigned precision = 6, bool use_names = true, bool coalunits = false) const;
            void    buildFromNewick(const string newick);
            void    debugShowEpochs(const epoch_list_t & epochs);
            
            virtual void createTrivialForest(bool compute_partials) = 0;
            virtual bool isSpeciesForest() const = 0;
            
            static double calcLogSum(const vector<double> & log_values);
            static string memoryAddressAsString(const void * ptr);
            static string speciesSetAsString(const Node::species_t & s);
            
            static unsigned multinomialDraw(const vector<double> & probs);

            static map<string, unsigned>    _taxon_to_species;

            static unsigned                 _nstates;
            
            static unsigned                 _ntaxa;
            static vector<string>           _taxon_names;

            static unsigned                 _nspecies;
            static vector<string>           _species_names;
            static map<unsigned,unsigned>   _nexus_taxon_map;

            static unsigned                 _ngenes;
            static vector<string>           _gene_names;

            static double                   _theta;
            static double                   _lambda;
            
            static double                   _theta_prior_mean;
            static double                   _lambda_prior_mean;
            
            static bool                      _update_theta;
            static bool                      _update_lambda;
            
            static double                   _small_enough;
            static double                   _infinity;
            static double                   _negative_infinity;
            
            static vector<double>           _cumprobs;  // workspace used by multinomialDraw
                        
        protected:
        
            int         extractNodeNumberFromName(Node * nd, set<unsigned> & used);
            void        extractEdgeLen(Node * nd, string edge_length_string);
            void        setNodeNameFromNumber(Node * nd);
            void        setSpeciesFromNodeName(Node * nd);
            void        setSpeciesFrommNodeName(Node * nd);
            Node *      findNextPreorder(Node * nd) const;
            bool        canHaveSibling(Node * nd) const;
            void        refreshPreorder(Node::ptr_vect_t & preorder) const;
            void        refreshAllPreorders() const;
            void        refreshAllHeightsAndPreorders();
            void        heightsInternalsPreorders();
            void        advanceAllLineagesBy(double t);
            Node *      joinLineagePair(Node * first, Node * second);
            void        unjoinLineagePair(Node * anc, Node * first, Node * second);
            void        removeTwoAddOne(Node::ptr_vect_t & node_vect, Node * del1, Node * del2, Node * add);
            void        addTwoRemoveOne(Node::ptr_vect_t & node_vect, Node * del1, Node * del2, Node * add);
            void        debugShowNodeInfo(string title);
            void        debugShowNodePtrVector(Node::ptr_vect_t& v, string title);

            virtual epoch_list_t & digest() = 0;
            virtual void operator=(const Forest & other);

            // NOTE: any variables added must be copied in operator=

            double                  _forest_height;
            unsigned                _next_node_index;
            unsigned                _next_node_number;
            
            double                  _prev_log_likelihood;
            double                  _prev_log_coalescent_likelihood;

            vector<Node>            _nodes;
            
            Node::ptr_vect_t        _lineages;
            
            epoch_list_t            _epochs;
                        
            // Because these can be recalulated at any time, they should not
            // affect the const status of the Forest object
            mutable vector<Node::ptr_vect_t> _preorders;
            mutable Epoch::lineage_counts_t _counts_workspace;
            
            // Debugging tool, should not affect const status
            mutable bool _debug_coal_like;
};
    
    inline Forest::Forest() {
    }

    inline Forest::~Forest() {
        clear();
    }

    inline void Forest::clear() {
        _forest_height = 0.0;
        _next_node_index = 0;
        _next_node_number = 0;
        _prev_log_likelihood = Forest::_negative_infinity;
        _prev_log_coalescent_likelihood = Forest::_negative_infinity;
        _nodes.clear();
        _preorders.clear();
        _lineages.clear();
        _epochs.clear();
        _counts_workspace.clear();
        _debug_coal_like = false;
    }
    
    inline void Forest::copyEpochsFrom(const epoch_list_t & other) {
        _epochs = other;
    }
        
    inline int Forest::extractNodeNumberFromName(Node * nd, set<unsigned> & used) {
        // Attempts to convert node name to a number and throws exception if node's
        // name cannot be converted to an integer or if node number has already
        // been encountered previously
        assert(nd);
        bool success = true;
        int x;
        try {
            x = stoi(nd->_name);
        }
        catch(invalid_argument &) {
            // node name could not be converted to an integer value
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
            throw XProj(str(format("node name (%s) not interpretable as an integer") % nd->_name));
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

    inline void Forest::setSpeciesFromNodeName(Node * nd) {
        nd->_species = {Forest::_taxon_to_species[nd->_name]};
    }
    
    inline void Forest::setNodeNameFromNumber(Node * nd) {
        unsigned n = nd->_number;
        if (isSpeciesForest()) {
            assert(n < 2*_species_names.size() - 1);
            if (n < _species_names.size())
                nd->_name = _species_names[n];
        }
        else {
            assert(n < 2*_taxon_names.size() - 1);
            if (n < _taxon_names.size())
                nd->_name = _taxon_names[n];
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
        double coalfactor = (coalunits ? (2.0/Forest::_theta) : 1.0);
        
        //const format basal_subtree_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_name_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_number_format( str(format("%%d:%%.%df") % precision) );
        const format internal_node_format( str(format("):%%.%df") % precision) );
        
        vector<string> subtree_newicks;
        for (auto preorder : _preorders) {
            string subtree_newick;
            stack<Node *> node_stack;
            for (auto nd : preorder) {
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
                        subtree_newick += str(format(tip_node_name_format)
                            % nd->_name
                            % edge_length);
                    } else {
                        subtree_newick += str(format(tip_node_number_format)
                            % (nd->_number + 1)
                            % edge_length);
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
                                    subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                }
                                popped = 0;
                            }
                            else {
                                subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                popped = node_stack.top();
                                popped_edge_length = popped->_edge_length*coalfactor;
                            }
                        }
                        if (popped && popped->_right_sib) {
                            node_stack.pop();
                            subtree_newick += str(format(internal_node_format) % popped_edge_length);
                            subtree_newick += ",";
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

    inline void Forest::buildFromNewick(const string newick) {
        // Builds strictly-bifurcating ultrametric, rooted, complete (i.e.
        // _lineages contains just one element) forest from the supplied newick
        // tree description. Assumes newick either uses names for leaves or,
        // if it specifies numbers, the numbers correspond to keys in
        // Forest::_nexus_taxon_map, which translates taxon numbers in newick
        // strings to the index of the taxon in Forest::_species_names (if
        // building a species tree) or Forest::_taxon_names (if building a
        // gene tree).

        //cerr << "Forest::buildFromNewick:" << endl;
        
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
            Node * nd   = &_nodes[_next_node_index];
            Node * root = &_nodes[_next_node_index++];
            
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
                            int num = extractNodeNumberFromName(nd, used);
                            assert(num > 0);
                            nd->_number = Forest::_nexus_taxon_map[(unsigned)num];
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
                            int num = extractNodeNumberFromName(nd, used);
                            assert(num > 0);
                            nd->_number = Forest::_nexus_taxon_map[(unsigned)num];
                            //nd->_number = (unsigned)(num - 1);
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

                        //cerr << "    edgelen " << edge_length_str << endl;
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
                        if (_next_node_index == _nodes.size())
                            throw XProj(str(format("Too many nodes specified by tree description (exceeds the %d nodes allocated for %d leaves)") % _nodes.size() % nleaves));
                        nd->_right_sib = &_nodes[_next_node_index++];
                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;

                        //cerr << "  node " << (_next_node_index-1) << " (rsib)" << endl;

                        break;

                    case '(':
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw XProj(str(format("Not expecting left parenthesis at position %d in tree description") % position_in_string));

                        // Create new node above and to the left of the current node
                        assert(!nd->_left_child);
                        if (_next_node_index == _nodes.size())
                            throw XProj(str(format("malformed tree description (more than %d nodes specified)") % _nodes.size()));
                        nd->_left_child = &_nodes[_next_node_index++];
                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;

                        //cerr << "  node " << (_next_node_index - 1) << " (lchild)" << endl;

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
                throw XProj(str(format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
            if (inside_edge_length && !(nd == root))
                throw XProj(str(format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
            if (inside_quoted_name)
                throw XProj(str(format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));

            heightsInternalsPreorders();
        }
        catch(XProj x) {
            clear();
            throw x;
        }
    }
    
    inline void Forest::refreshPreorder(Node::ptr_vect_t & preorder) const {
        // Assumes preorder just contains the root node when this function is called
        assert(preorder.size() == 1);
        
        Node * nd = preorder[0];
        while (true) {
            nd = findNextPreorder(nd);
            if (nd)
                preorder.push_back(nd);
            else
                break;
        }
    }
    
    inline void Forest::refreshAllPreorders() const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        _preorders.clear();
        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
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
        
        refreshAllPreorders();
        _forest_height = 0.0;

        // Set heights for each lineage in turn
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != Forest::_infinity);
                    if (nd->_height + nd->_edge_length > _forest_height)
                        _forest_height = nd->_height + nd->_edge_length;
                }
                else {
                    // nd is a leaf node
                    nd->_height = 0.0;
                    if (nd->_edge_length > _forest_height)
                        _forest_height = nd->_edge_length;
                }
                
                if (nd->_parent) {
                    // Set parent's height if nd is right-most child of its parent
                    bool is_rightmost_child = !nd->_right_sib;
                    double parent_height = nd->_height + nd->_edge_length;
                    if (is_rightmost_child) {
                        nd->_parent->_height = parent_height;
                    }
                    
                    // If nd is not its parent's rightmost child, check ultrametric assumption
                    assert(!is_rightmost_child || fabs(nd->_parent->_height - parent_height) < Forest::_small_enough);
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
        _next_node_number = (unsigned)(isSpeciesForest() ? _species_names.size() : _taxon_names.size());

        // Renumber internal nodes in postorder sequence for each lineage in turn
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != Forest::_infinity);
                    if (nd->_height > _forest_height)
                        _forest_height = nd->_height;
                    nd->_number = _next_node_number++;
                }
                else {
                    // nd is a leaf node
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
                    assert(!is_rightmost_child || fabs(nd->_parent->_height - parent_height) < Forest::_small_enough);
                }
            }
        }
    }
    
    inline void Forest::debugShow(const format & f) {
#if defined(DEBUGGING)
        cout << str(f) << endl;
#endif
    }
    
    inline string Forest::unsignedVectToString(const vector<unsigned> & v) {
        ostringstream oss;
        
        // Create string with values separated by commas
        copy(v.begin(), v.end(), ostream_iterator<unsigned>(oss, ","));

        // Remove the trailing comma
        oss.seekp(-1, oss.end);
        oss.put('\0');
        
        // Return string
        return oss.str();
    }
    
    inline string Forest::lineageCountsMapToString(const Epoch::lineage_counts_t & lcm) {
        ostringstream oss;
        
        // Output species as headers
        oss << "***   ";
        for (auto iter = lcm.begin(); iter != lcm.end(); ++iter) {
            if (iter->second > 0)
                oss << str(format("%10s") % speciesSetAsString(iter->first));
        }
        oss << endl;

        // Output lineage counts
        oss << "***   ";
        for (auto iter = lcm.begin(); iter != lcm.end(); ++iter) {
            if (iter->second > 0)
                oss << str(format("%10d") % iter->second);
        }
        oss << endl;
        
        // Return string
        return oss.str();
    }
    
    inline void Forest::createDefaultSpeciesTreeNexusTaxonMap() {
        // Build default _nexus_taxon_map used by buildFromNewick
        // that assumes taxon indices in the newick tree description are in
        // the same order as _species_names
        _nexus_taxon_map.clear();
        for (unsigned i = 0; i < Forest::_nspecies; ++i) {
            _nexus_taxon_map[i+1] = i;
        }
    }
    
    inline void Forest::createDefaultGeneTreeNexusTaxonMap() {
        // Build default _nexus_taxon_map used by buildFromNewick
        // that assumes taxon indices in the newick tree description are in
        // the same order as _taxon_names
        _nexus_taxon_map.clear();
        for (unsigned i = 0; i < Forest::_ntaxa; ++i) {
            _nexus_taxon_map[i+1] = i;
        }
    }
    
    inline void Forest::debugCheckSpeciesSets() const {
#if defined(DEBUGGING)
        for (auto & preorder : _preorders) {
            for (auto nd : preorder) {
                for (auto s : nd->_species) {
                    assert(s < Forest::_nspecies);
                }
            }
        }
#endif
    }

    inline void Forest::readTreefile(const string filename,
                                     unsigned skip,
                                     const vector<string> & leaf_names,
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
                throw XProj(format("Found %d taxa in taxa block %d but was expecting %d") % num_leaves % (i+1) % leaf_names.size());
            }
            
            // It is possible (probable even) that the order of taxa in taxa_block will differ from
            // the order in leaf_names. Therefore, taxon_map is constructed the keys of which are the
            // taxon numbers in the newick tree description (corresponding to taxa_block); the values
            // are the index of the taxon in leaf_names. For example:
            //
            //   #nexus              newick = (1,2,(3,4)) --> (A,D,(C,B))
            //   begin trees;        leaf_names: [A,B,C,D]
            //     translate
            //       1 A,            t           = 0  1  2  3
            //       2 D,            taxon_label = A  D  C  B  <-- order in taxa block
            //       3 C,            d           = 0  3  2  1  <-- index into leaf_names
            //       4 B             t+1         = 1  2  3  4
            //     ;                 taxon_map: {1:0, 2:3, 3:2, 4:1}
            //   end;                                 | |
            //                                        | index into leaf_names vector
            //                                        number in newick string
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
                    //cout << "Trees block contains " << ntrees << " tree descriptions.\n";
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
    
    inline double Forest::calcLogSum(const vector<double> & log_values) {
        double max_logv = *max_element(log_values.begin(), log_values.end());
        
        double factored_sum = 0.0;
        for (auto & logv : log_values) {
            factored_sum += exp(logv - max_logv);
        }
        double log_sum_values = max_logv + log(factored_sum);
        return log_sum_values;
    }
    
    inline double Forest::calcMaxT(epoch_list_t & epochs, double h0, set<Node::species_t> current_species) {
        // This function returns the height of the next coalescent epoch in which
        // lineages from different species coalesce. That sets the upper limit
        // (looking backwards in time) for the next speciation event.
        
        debugShowEpochs(epochs);
        
        // Find first epoch with height greater than h0
        auto it0 = find_if(epochs.begin(), epochs.end(), [h0](Epoch & ep){return ep._height > h0;});
        
        // Starting with it0, look for coalescent epoch in which lineages from different current_species coalesce
        auto it = find_if(it0, epochs.end(), [&current_species](Epoch & ep){
            if (ep.isCoalescentEpoch()) {
                bool different_species = current_species.find(ep._species) == current_species.end();
                return different_species;
            }
            else {
                return false;
            }
        });
        
        // We should only be calling this function if another speciation event is possible
        // and, if another speciation event is possible, there must be a coalescent event
        // that combines lineages from different species.
        assert(it != epochs.end());

        return it->_height;
    }
        
    inline void Forest::advanceAllLineagesBy(double t) {
        assert(t >= 0.0);
        assert(t != Forest::_infinity);
        
        // Add t to the edge length of all lineage root nodes
        for (auto nd : _lineages) {
            double elen = nd->getEdgeLength() + t;
            //cout << nd->_name << ": elen = " << elen << endl;
            nd->setEdgeLength(elen);
        }
        
        // Add to to the current forest height
        _forest_height += t;
    }
    
    inline unsigned Forest::multinomialDraw(const vector<double> & probs) {
        // Compute cumulative probababilities
        _cumprobs.resize(probs.size());
        partial_sum(probs.begin(), probs.end(), _cumprobs.begin());

        // Draw a Uniform(0,1) random deviate
        double u = rng.uniform();

        // Find first element in _cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6                <-- u
        //
        // _cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(_cumprobs.begin(), _cumprobs.end(), [u](double cumpr){return cumpr > u;});
                
        assert(it != _cumprobs.end());
        return (unsigned)distance(_cumprobs.begin(), it);
    }
    
    inline Node * Forest::joinLineagePair(Node * subtree1, Node * subtree2) {
        // Get new node to serve as the ancestral node
        Node * new_nd = &_nodes[_next_node_index++];
        new_nd->_number = _next_node_number++;
        new_nd->_name        = "anc-" + to_string(new_nd->_number);
        new_nd->_left_child  = subtree1;
        new_nd->_edge_length = 0.0;
        assert(new_nd->_partial == nullptr);
        
        // Calculate height of the new node
        // (should be the same whether computed via the left or right child)
        double h1 = subtree1->_height + subtree1->_edge_length;
        double h2 = subtree2->_height + subtree2->_edge_length;
        assert(fabs(h1 - h2) < _small_enough);
        new_nd->_height = (h1 + h2)/2.0;

        // Finish connecting new trio of nodes
        subtree1->_right_sib = subtree2;
        subtree1->_parent    = new_nd;
        subtree2->_parent    = new_nd;

        return new_nd;
    }
    
    inline void Forest::unjoinLineagePair(Node * anc, Node * subtree1, Node * subtree2) {
        // Reset members set by joinLineagePair function
        anc->_number = -1;
        anc->_name = "";
        anc->_left_child = nullptr;
        anc->_edge_length = 0.0;
        anc->_height = 0.0;
        //anc->_partial = nullptr;
        subtree1->_right_sib = nullptr;
        subtree1->_parent = nullptr;
        subtree2->_parent = nullptr;
        _next_node_index--;
        _next_node_number--;
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
    
    inline void Forest::debugShowEpochs(const epoch_list_t & epochs) {
#if defined(DEBUG_COAL_LIKE)
        if (!_debug_coal_like)
            return;
            
        cout << "\nEpochs:" << endl;
        for (auto & epoch : epochs) {
            double h = epoch._height;
            if (epoch.isInitEpoch()) {
                ostringstream oss;
                for (auto it = epoch._lineage_counts.begin(); it != epoch._lineage_counts.end(); it++)
                    oss << it->second << " ";
                oss.seekp(-1, oss.cur); // Remove the trailing ' '
                cout << "  " << (epoch._valid ? " " : "x") << h << str(format(": init gene tree %d (lineage counts: %s)") % epoch._gene % oss.str()) << endl;
            }
            else if (epoch.isCoalescentEpoch()) {
                assert(epoch._species.size() > 0);
                string species_set = Forest::speciesSetAsString(epoch._species);
                cout << "  " << (epoch._valid ? " " : "x") << h << str(format(": coalescence in gene tree %d (species %s)") % epoch._gene % species_set) << endl;
            }
            else {
                assert(epoch._left_species.size() > 0);
                assert(epoch._right_species.size() > 0);
                assert(epoch._anc_species.size() > 0);
                cout << "  " << (epoch._valid ? " " : "x") << h << str(format(": speciation event (%s,%s -> %s)") %  epoch.leftSpeciesAsStr() % epoch.rightSpeciesAsStr() % epoch.ancSpeciesAsStr()) << endl;
            }
         }
#endif
    }
    
    inline string Forest::memoryAddressAsString(const void * ptr) {
        ostringstream memory_address;
        memory_address << ptr;
        return memory_address.str();
    }
    
    inline string Forest::speciesSetAsString(const Node::species_t & s) {
        if (s.size() == 0)
            return string();
        else if (s.size() == 1) {
            unsigned species = *(s.begin());
            return to_string(species);
        }
        else {
            vector<string> species;
            for (auto spp : s) {
                species.push_back(to_string(spp));
            }
            return join(species, "+");
        }
    }
    
    inline void Forest::debugShowNodeInfo(string title) {
#if defined(DEBUGGING)
        cout << "\nNode info: " << title << "\n";
        cout << str(format("%6s %12s %6s %10s %10s %6s %6s %6s %12s %12s\n") % "index" % "address" % "number" % "species" % "type" % "lchild" % "rsib" % "parent" % "height" % "brlen");
        for (unsigned i = 0; i < _nodes.size(); ++i) {
            string memory_address = memoryAddressAsString((void *)&_nodes[i]);
            string type   = _nodes[i]._number == -1 ? "unused" : (_nodes[i]._left_child ? "internal" : "leaf");
            string lchild = _nodes[i]._left_child ? to_string(_nodes[i]._left_child->_number) : "null";
            string rsib   = _nodes[i]._right_sib ? to_string(_nodes[i]._right_sib->_number) : "null";
            string parent = _nodes[i]._parent ? to_string(_nodes[i]._parent->_number) : "null";
            string species_set = speciesSetAsString(_nodes[i]._species);
            cout << str(format("%6d %12s %6d %10s %10s %6s %6s %6s %12.5f %12.5f\n") % i % memory_address % _nodes[i]._number % species_set % type % lchild % rsib % parent % _nodes[i]._height % _nodes[i]._edge_length);
        }
        cout << endl;
#endif
    }

    inline void Forest::debugShowNodePtrVector(Node::ptr_vect_t & v, string title) {
#if defined(DEBUGGING)
        cout << "\nNode::ptr_vect_t: " << title << "\n";
        cout << str(format("%6s %12s %6s %10s\n") % "index" % "address" % "number" % "species");
        for (unsigned i = 0; i < v.size(); ++i) {
            string memory_address = memoryAddressAsString((void *)&_nodes[i]);
            string species_set = speciesSetAsString(v[i]->getSpecies());
            cout << str(format("%6d %12s %6d %10s\n") % i % memory_address % v[i]->_number % species_set);
        }
        cout << endl;
#endif
    }

    inline double Forest::calcLogCoalescentLikelihood(const epoch_list_t & epochs, unsigned gene) const {
        // Computes the log of the coalescent likelihood, which is the probability of a
        // gene tree given the species tree. If GeneForest, stop when _forest_height is
        // reached; if SpeciesForest, continue until all epochs are considered but combine
        // all lineages into a single ancestral lineage once _forest_height is exceeded.
#if defined(DEBUG_COAL_LIKE)
        if (_debug_coal_like) {
            cout << str(format("\n*** Calculating log(coalescent likelihood) for gene %d\n") % gene);
            cout << str(format("***   theta = %.5f\n") % Forest::_theta);
        }
#endif

        double log_coal_like = 0.0;
        double prev_height = 0.0;
        
        // If building species forest, everything further in the past than _forest_height
        // is considered a single ancestral species (common_pool = true) and the variable
        // ncommonpool keeps track of the number of lineages in this ancestral species
        bool common_pool = false;
        unsigned ncommonpool = 0;
        
        _counts_workspace.clear();
        for (const Epoch & e : epochs) {
            // Ignore any init or coalescent epochs that do not apply to the gene under consideration
            if ((e.isInitEpoch() || e.isCoalescentEpoch()) && e._gene != gene)
                continue;
                
            if (e._height < 0.0) {
                // If height of epoch is negative, it should only be an init epoch
                assert(e.isInitEpoch());
                
                // Initialize lineage counts
                _counts_workspace = e._lineage_counts;
#if defined(DEBUG_COAL_LIKE)
                if (_debug_coal_like) {
                    cout << str(format("***   Init epoch for gene %d:\n") % e._gene);
                    cout << "***     _counts_workspace:\n";
                    cout << Forest::lineageCountsMapToString(_counts_workspace);
                }
#endif
                continue;
            }
            
            if (e._height > _forest_height) {
                if (isSpeciesForest()) {
                    if (!common_pool) {
                        // The coalescent likelihood is being used in the context of advancing
                        // a species tree, in which case we must now dump all lineages into
                        // one ancestral species
                        common_pool = true;
                        ncommonpool = 0;
                        Node::species_t anc;
                        for (auto & map_pair : _counts_workspace) {
                            anc.insert(map_pair.first.begin(), map_pair.first.end());
                            ncommonpool += map_pair.second;
                        }
                        _counts_workspace.clear();
#if defined(DEBUG_COAL_LIKE)
                        if (_debug_coal_like) {
                            cout << str(format("***   Dumping all %d remaining lineages into a single ancestal species because\n") % ncommonpool);
                            cout << str(format("***      we are building a species tree and the epoch height (%.9f)\n") % e._height);
                            cout << str(format("***      exceeds the forest height (%.9f)\n") % _forest_height);
                        }
#endif
                    }
                }
                else {
                    // The coalescent likelihood is being used in the context of advancing
                    // a gene tree, in which case we are done when we get to the current
                    // height of the gene tree.
#if defined(DEBUG_COAL_LIKE)
                    if (_debug_coal_like) {
                        cout << str(format("***   Stopping because we are building a gene tree and the epoch height (%.9f)\n") % e._height);
                        cout << str(format("***      exceeds forest height (%.9f)\n") % _forest_height);
                    }
#endif
                    break;
                }
            }

            assert(common_pool || !_counts_workspace.empty());
                                
            // This is the time since the previous coalescence or speciation
            double t = e._height - prev_height;
                
#if defined(DEBUG_COAL_LIKE)
            if (_debug_coal_like) {
                cout << str(format("***\n***   t = %.20f = %.20f - %.20f\n") % t % e._height % prev_height);
            }
#endif
            if (e.isCoalescentEpoch()) {
                // This epoch ends in a coalescent event, so we must consider:
                // 1. the probability that no lineages coalesce for time t (e^{-rt})
                // 2. the probability density that a pair of lineages in one
                //    species coalesces exactly at time t (r)
                // 3. the probability that the particular pair of lineages that
                //    coalesced would be chosen out of the n lineages still
                //    separate in the species (1/(n choose 2))
                // where r is the coalescent rate (n*(n-1)/theta).
                double lnL = 0.0;
                
#if defined(DEBUG_COAL_LIKE)
                vector<string> lnLparts;
                if (_debug_coal_like) {
                    cout << str(format("***   Coalescent event in gene %d, species %s, at height %.9f\n") % e._gene % e.speciesSetAsStr() % e._height);
                }
#endif
                
                // Calculate probability of no coalescence over time t in all extant lineages
                if (common_pool) {
                    // ncommonpool is the number of lineages in the common pool
                    double log_prob_no_coal = 0.0;
                    if (ncommonpool > 1) {
                        log_prob_no_coal = -1.0*ncommonpool*(ncommonpool-1)*t/Forest::_theta;
                        lnL += log_prob_no_coal;
                    }
#if defined(DEBUG_COAL_LIKE)
                    if (_debug_coal_like) {
                        lnLparts.push_back(str(format("%.9f") % log_prob_no_coal));
                        cout << str(format("***     %.9f log prob no coal in gene %d, common pool (%d lineages)\n") % log_prob_no_coal % gene % ncommonpool);
                    }
#endif
                }
                else {
                    for (auto it = _counts_workspace.begin(); it != _counts_workspace.end(); ++it) {
                        // n is the number of lineages of a particular species
                        unsigned n = it->second;
                        double log_prob_no_coal = 0.0;
                        if (n > 1) {
                            log_prob_no_coal = -1.0*n*(n-1)*t/Forest::_theta;
                            lnL += log_prob_no_coal;
                        }
#if defined(DEBUG_COAL_LIKE)
                        if (_debug_coal_like) {
                            lnLparts.push_back(str(format("%.9f") % log_prob_no_coal));
                            cout << str(format("***     %.9f log prob no coal in gene %d, spp %s (%d lineages)\n") % log_prob_no_coal % gene % speciesSetAsString(it->first) % n);
                        }
#endif
                    }
                }
                
                // Calculate the probability density of a coalescence event at exactly time t
                unsigned n = (common_pool ? ncommonpool : _counts_workspace[e._species]);
                
                assert(n > 1);
                double log_prob_coal = log(float(n)*(n-1)/Forest::_theta);
                lnL += log_prob_coal;
#if defined(DEBUG_COAL_LIKE)
                if (_debug_coal_like) {
                    lnLparts.push_back(str(format("%.9f") % log_prob_coal));
                    if (common_pool)
                        cout << str(format("***     %.9f log prob coal gene %d, common pool (%d lineages)\n") % log_prob_coal % gene % n);
                    else
                        cout << str(format("***     %.9f log prob coal gene %d, spp %s (%d lineages)\n") % log_prob_coal % gene % speciesSetAsString(e._species) % n);
                }
#endif
                
                // Calculate the probability of choosing the two lineages that coalesced
                double log_prob_join = -log(float(n)*(n-1)/2.0);
                lnL += log_prob_join;
#if defined(DEBUG_COAL_LIKE)
                if (_debug_coal_like) {
                    lnLparts.push_back(str(format("%.9f") % log_prob_join));
                    if (common_pool)
                        cout << str(format("***     %.9f log prob join gene %d, common pool (%d lineages)\n") % log_prob_join % gene % n);
                    else
                        cout << str(format("***     %.9f log prob join gene %d, spp %s (%d lineages)\n") % log_prob_join % gene % speciesSetAsString(e._species) % n);
                }
#endif

                // Adjust lineage counts
                if (common_pool) {
                    --ncommonpool;
#if defined(DEBUG_COAL_LIKE)
                    if (_debug_coal_like) {
                        cout << str(format("***     ncommonpool updated to: %d\n") % ncommonpool);
                    }
#endif
                }
                else {
                    _counts_workspace[e._species] -= 1;
#if defined(DEBUG_COAL_LIKE)
                    if (_debug_coal_like) {
                        cout << "***     _counts_workspace updated to:\n";
                        cout << Forest::lineageCountsMapToString(_counts_workspace);
                    }
#endif
                }
                
#if defined(DEBUG_COAL_LIKE)
                if (_debug_coal_like) {
                    cout << str(format("***     %.9f lnL this epoch (=(%s))\n") % lnL % join(lnLparts, ")+("));
                    cout << str(format("***     %.9f log_coal_like\n") % (log_coal_like + lnL));
                }
#endif
                log_coal_like += lnL;
            }
            else {
                // This is a speciation epoch, so need to compute the probability of no coalescence
                // over time t
                double lnL = 0.0;
                
#if defined(DEBUG_COAL_LIKE)
                if (_debug_coal_like) {
                    cout << str(format("***   Speciation event at height %.9f\n") % e._height);
                    cout << str(format("***      %s is the left species\n") % e.leftSpeciesAsStr());
                    cout << str(format("***      %s is the right species\n") % e.rightSpeciesAsStr());
                    cout << str(format("***      %s is the ancestral species\n") % e.ancSpeciesAsStr());
                }
#endif
                // Calculate probability of no coalescence over time t
                for (auto it = _counts_workspace.begin(); it != _counts_workspace.end(); ++it) {
                    // n is the number of lineages of a particular species
                    unsigned n = it->second;
                    double log_prob_no_coal = 0.0;
                    if (n > 1) {
                        log_prob_no_coal = -1.0*n*(n-1)*t/Forest::_theta;
                        lnL += log_prob_no_coal;
                    }
#if defined(DEBUG_COAL_LIKE)
                    if (_debug_coal_like) {
                        cout << str(format("***     %.9f log prob no coal gene %d, spp %s (%d lineages)\n") % log_prob_no_coal % gene % speciesSetAsString(it->first) % n);
                    }
#endif
                }
                
                // Adjust lineage counts
                _counts_workspace[e._anc_species] += _counts_workspace[e._left_species];
                _counts_workspace[e._anc_species] += _counts_workspace[e._right_species];
                _counts_workspace.erase(e._left_species);
                _counts_workspace.erase(e._right_species);
#if defined(DEBUG_COAL_LIKE)
                if (_debug_coal_like) {
                    cout << "***     _counts_workspace updated to:\n";
                    cout << Forest::lineageCountsMapToString(_counts_workspace);
                    cout << str(format("***   %.9f lnL this epoch\n") % lnL);
                    cout << str(format("***   %.9f log_coal_like\n") % (log_coal_like + lnL));
                }
#endif
                log_coal_like += lnL;
            }
                
            prev_height = e._height;
        }
                        
        if (log_coal_like != log_coal_like) {
            cerr << "debug stop" << endl;
        }
        if (log_coal_like == Forest::_negative_infinity) {
            cerr << "debug stop" << endl;
        }
        return log_coal_like;
    }
    
    inline void Forest::operator=(const Forest & other) {
        _debug_coal_like     = other._debug_coal_like;
        _forest_height       = other._forest_height;
        _next_node_index     = other._next_node_index;
        _next_node_number    = other._next_node_number;
        _prev_log_likelihood = other._prev_log_likelihood;
        _prev_log_coalescent_likelihood = other._prev_log_coalescent_likelihood;
        
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
            
            // _nodes[i]._partial copied by GeneForest::operator=
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
        
        // The _preorders vector is recreated as needed
        _preorders.clear();
        
        // Deep-copy the _epochs vector
        _epochs = other._epochs;

        // No need to copy _counts_workspace because it is rebuilt
        // every time it is used
        _counts_workspace.clear();
    }

}
