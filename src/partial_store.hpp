#pragma once

namespace proj {

    struct Partial {
        Partial(unsigned g, unsigned n);
        ~Partial();
        unsigned        _g; // the gene
        vector<double>  _v; // the partial array: length = _nstates*<no. patterns>
        
#if defined(MEMORY_FRUGAL)
        bool            _in_storage;
#endif
        
#if defined(LOG_MEMORY)
        static vector<unsigned> _nconstructed;  // no. partials allocated for each gene
        static vector<unsigned> _ndestroyed;    // no. partials destroyed for each gene
        static vector<unsigned> _max_in_use;    // maximum no. partials in use at any one time for each gene
        static vector<unsigned long> _bytes_per_partial;    // no. bytes required to store partial array for each gene
        static unsigned long _total_max_in_use;      // total max. in use over all genes
        static unsigned long _total_max_bytes;       // total max. bytes over all genes
        static unsigned _nstates;               // no. of states (4 for DNA)
#endif
    };

    inline Partial::Partial(unsigned g, unsigned n) {
        _g = g;
        _v.resize(n);

#if defined(MEMORY_FRUGAL)
        _in_storage = true;
#endif

#if defined(LOG_MEMORY)
        assert(g < _nconstructed.size());
        _nconstructed[g]++;
        
        if (_nconstructed[g] -_ndestroyed[g] > _max_in_use[g]) {
            _max_in_use[g] = _nconstructed[g] -_ndestroyed[g];
        }
        
        _total_max_in_use = 0;
        _total_max_bytes = 0;
        for (unsigned g = 0; g < _nconstructed.size(); ++g) {
            unsigned in_use = _nconstructed[g] -_ndestroyed[g];
            _total_max_in_use += in_use;
            _total_max_bytes = in_use*_bytes_per_partial[g];
        }
#endif
    }
    
    inline Partial::~Partial() {
#if defined(LOG_MEMORY)
        _ndestroyed[_g]++;
#endif
    }

    class PartialStore {
    
        public:
                                    PartialStore();
                                    ~PartialStore();
                                    
            typedef shared_ptr<Partial>          partial_t;
            
            typedef vector<partial_t>            vect_partial_t;
            typedef vector<vect_partial_t>       leaf_partials_t;
            typedef vector<vect_partial_t>       storage_t;
            
            void            clear();
            void            setNGenes(unsigned ngenes);
            partial_t       getPartial(unsigned gene);
            void            putPartial(unsigned gene, partial_t partial);
            void            setNElements(unsigned nelements, unsigned gene);
            unsigned        getNElements(unsigned gene) const {return _nelements[gene];}
            
#if defined(LOG_MEMORY)
            unsigned        getInUse();
            unsigned        getStored();
            unsigned        getNumberConstructed();
            unsigned        getNumberDestroyed();
            void            memoryReport(ofstream & memf) const;
#endif
            
        private:
        
            vector<unsigned> _nelements;
            storage_t        _storage;
    };

    inline PartialStore::PartialStore() {
    }

    inline PartialStore::~PartialStore() {
        clear();
    }
    
    inline void PartialStore::clear() {
        _nelements.clear();
        _storage.clear();
    }
    
#if defined(LOG_MEMORY)
    inline unsigned PartialStore::getNumberConstructed() {
        return (unsigned)accumulate(Partial::_nconstructed.begin(), Partial::_nconstructed.end(), 0);
    }

    inline unsigned PartialStore::getNumberDestroyed() {
        return (unsigned)accumulate(Partial::_ndestroyed.begin(), Partial::_ndestroyed.end(), 0);
    }

    inline unsigned PartialStore::getStored() {
        unsigned total_stored = 0;
        unsigned ngenes = (unsigned)_storage.size();
        for (unsigned g = 0; g < ngenes; ++g) {
            total_stored += (unsigned)_storage[g].size();
        }
        return total_stored;
    }
    
    inline unsigned PartialStore::getInUse() {
        unsigned total_in_use = 0;
        unsigned ngenes = (unsigned)Partial::_nconstructed.size();
        for (unsigned g = 0; g < ngenes; ++g) {
            unsigned in_use = Partial::_nconstructed[g] - Partial::_ndestroyed[g];
            total_in_use += in_use;
        }
        return total_in_use;
    }
#endif

    inline void PartialStore::setNGenes(unsigned ngenes) {
        // Should be called before any partials are stored
        assert(_nelements.empty());
        assert(_storage.empty());
        
        // Resize both containers
        _nelements.resize(ngenes);
        _storage.resize(ngenes);
        
#if defined(LOG_MEMORY)
        Partial::_nconstructed.clear();
        Partial::_nconstructed.resize(ngenes, 0);
        Partial::_ndestroyed.clear();
        Partial::_ndestroyed.resize(ngenes, 0);
        Partial::_max_in_use.clear();
        Partial::_max_in_use.resize(ngenes, 0);
        Partial::_bytes_per_partial.clear();
        Partial::_bytes_per_partial.resize(ngenes, 0);
#endif
    }

    inline void PartialStore::setNElements(unsigned nelements, unsigned gene) {
        assert(_nelements.size() > gene);
        _nelements[gene] = nelements;
#if defined(LOG_MEMORY)
        Partial::_bytes_per_partial[gene] = (unsigned)nelements*sizeof(double);
#endif
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned gene) {
        // Check to make sure supplied value of gene is valid
        assert(_nelements.size() > gene);
        assert(_nelements[gene] > 0);
        
        partial_t partial;
        if (_storage[gene].empty()) {
            // No stored partials for this gene, so allocate one
            partial = partial_t(new Partial(gene, _nelements[gene]));
        }
        else {
            size_t n = _storage[gene].size();
            partial = _storage[gene].at(n-1);
            _storage[gene].pop_back();
        }

#if defined(MEMORY_FRUGAL)
        partial->_in_storage = false;
#endif

        return partial;
    }
    
    inline void PartialStore::putPartial(unsigned gene, partial_t partial) {
        // Check to make sure supplied value of gene is valid
        assert(_nelements.size() > gene);
        assert(_nelements[gene] > 0);

        // Store the partial for later
        assert(partial->_v.size() == _nelements[gene]);
        partial->_v.assign(_nelements[gene], 0.0);

#if defined(MEMORY_FRUGAL)
        partial->_in_storage = true;
#endif
        
        _storage[gene].push_back(partial);
    }
    
#if defined(LOG_MEMORY)
    inline void PartialStore::memoryReport(ofstream & memf) const {
        memf << "\nPartialStore memory report:\n\n";
        memf << "Note: one partial is computed for each leaf and simply copied when\n";
        memf << "      trivial gene forests are created. Thus, the upper bound for the\n";
        memf << "      number of partials needed is (nleaves - 1)*nparticles + nleaves\n\n";
        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") %         "gene" %        "bytes" %    "allocated" %       "in use" %   "max in use" %    "max bytes");
        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
        unsigned total_allocated = 0;
        unsigned total_in_use = 0;
        unsigned total_max_in_use = 0;
        unsigned long total_max_nbytes = 0L;
        unsigned ngenes = (unsigned)Partial::_nconstructed.size();
        for (unsigned g = 0; g < ngenes; ++g) {
            unsigned nallocated = Partial::_nconstructed[g];
            total_allocated += nallocated;
            
            unsigned in_use = Partial::_nconstructed[g] - Partial::_ndestroyed[g];
            total_in_use += in_use;

            unsigned max_in_use = Partial::_max_in_use[g];
            total_max_in_use += max_in_use;

            unsigned long max_nbytes = _nelements[g]*max_in_use*sizeof(double);
            total_max_nbytes += max_nbytes;
            
            unsigned bytes_each = _nelements[g]*sizeof(double);
            
            memf << str(format("  %12d %12d %12d %12d %12d %12d\n") % g % bytes_each % nallocated % in_use % max_in_use % max_nbytes);
        }
        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
        memf << str(format("  %12s %12s %12s %12d %12d %12d\n") % " " % " " % total_allocated % total_in_use % total_max_in_use % total_max_nbytes);
        memf << str(format("  Total megabytes: %.5f\n") % (1.0*total_max_nbytes/1048576));
        memf << str(format("\n  %12d Maximum partials in use at any one time over all genes\n") % Partial::_total_max_in_use);
        memf << str(format("  %12d Maximum bytes in use by partials at any one time over all genes\n") % Partial::_total_max_bytes);
    }
#endif
    
 }
