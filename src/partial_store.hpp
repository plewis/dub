#pragma once

namespace proj {

    class PartialStore {
    
        public:
                                    PartialStore();
                                    ~PartialStore();
                                    
            typedef shared_ptr<vector<double> >  partial_t;
            typedef vector<partial_t>            vect_partial_t;
            typedef vector<vect_partial_t>       leaf_partials_t;
            typedef vector<vect_partial_t>       storage_t;
            
            void       setNGenes(unsigned ngenes);
            partial_t  getPartial(unsigned gene);
            void       stowPartial(partial_t & p, unsigned gene);
            void       setNElements(unsigned nelements, unsigned gene);
            unsigned   getNElements(unsigned gene) const {return _nelements[gene];}
            unsigned   getNAllocated(unsigned gene) const {return _nallocated[gene];}
            
            void           memoryReport(ofstream & memf) const;
            unsigned long  getTotalBytesAllocated() const;
            
#if defined(DEBUG_PARTIAL_STORE)
            void       showStorage(unsigned gene);
            void       history(string message);
#endif
            
        private:
        
            vector<unsigned> _nelements;
            storage_t        _storage;
            vector<unsigned> _nallocated;
    };

    inline PartialStore::PartialStore() {
#if defined(DEBUG_PARTIAL_STORE)
        history("Created PartialStore object.");
#endif
    }

    inline PartialStore::~PartialStore() {
        _nelements.clear();
        _nallocated.clear();
        _storage.clear();
        
#if defined(DEBUG_PARTIAL_STORE)
        history("Destroyed PartialStore object.");
#endif
    }
    
#if defined(DEBUG_PARTIAL_STORE)
    inline void PartialStore::history(string message) {
        ofstream histf("pshistory.txt", ios::out | ios::app);
        histf << message << endl;
        histf.close();
    }
#endif

    inline void PartialStore::setNGenes(unsigned ngenes) {
        // Should be called before any partials are stored
        assert(_nelements.empty());
        assert(_nallocated.empty());
        assert(_storage.empty());
        
        // Resize both containers
        _nelements.resize(ngenes);
        _nallocated.resize(ngenes, 0);
        _storage.resize(ngenes);
        
#if defined(DEBUG_PARTIAL_STORE)
        history(str(format("Setting number of genes to %d") % ngenes));
#endif
    }

    inline void PartialStore::setNElements(unsigned nelements, unsigned gene) {
        assert(_nelements.size() > gene);
        _nelements[gene] = nelements;
        
#if defined(DEBUG_PARTIAL_STORE)
        history(str(format("Setting number of elements for gene %d to %d") % gene % nelements));
#endif
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned gene) {
        assert(_nelements.size() > gene);
        assert(_nallocated.size() > gene);
        assert(_storage.size() > gene);
        assert(_nelements[gene] > 0);
        if (_storage[gene].size() > 0) {
            partial_t last = _storage[gene].back();
            _storage[gene].pop_back();
#if defined(DEBUG_PARTIAL_STORE)
            ostringstream memory_address;
            memory_address << last.get();
            history(str(format("Popping partial %s (%d)") % memory_address.str() % last.use_count()));
#endif
            // //temporary!
            //cerr << str(format("~~~> using stored partial for gene %d\n") % gene);

            return last;
        }
        else {
            partial_t ptr = partial_t(new vector<double>(_nelements[gene]));
            _nallocated[gene]++;
#if defined(DEBUG_PARTIAL_STORE)
            ostringstream memory_address;
            memory_address << ptr.get();
            history(str(format("Allocating partial %s (%d)") % memory_address.str() % ptr.use_count()));
#endif

            // //temporary!
            //cerr << str(format("~~~> allocating new partial for gene %d\n") % gene);

            return ptr;
        }
    }
    
    inline void PartialStore::stowPartial(partial_t & p, unsigned gene) {
        assert(_storage.size() > gene);
        fill(p->begin(), p->end(), 0.0);
        
#if defined(DEBUG_PARTIAL_STORE)
        ostringstream memory_address;
        memory_address << p.get();
        history(str(format("Pushing partial %s (%d)") % memory_address.str() % p.use_count()));
#endif

        // //temporary!
        //cerr << str(format("~~~> stowing partial for gene %d\n") % gene);

        _storage[gene].push_back(p);
    }
    
#if defined(DEBUG_PARTIAL_STORE)
    inline void PartialStore::showStorage(unsigned gene) {
        history(str(format("Storage for gene %d") % gene));
        for (unsigned i = 0; i < _storage[gene].size(); ++i) {
            ostringstream memory_address;
            memory_address << _storage[gene][i].get();
            history(str(format("  stored partial %s (%d)") % memory_address.str() % _storage[gene][i].use_count()));
        }
    }
#endif

    inline unsigned long PartialStore::getTotalBytesAllocated() const {
        unsigned long total_bytes = 0L;
        for (unsigned g = 0; g < _nallocated.size(); ++g) {
            total_bytes += _nelements[g]*_nallocated[g];
        }
        return total_bytes;
    }
    
    inline void PartialStore::memoryReport(ofstream & memf) const {
        memf << "\nPartialStore memory report:\n\n";
        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") %         "gene" %         "size" %       "in use" %       "stored" %        "total" %        "bytes");
        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
        unsigned long total_bytes = 0L;
        unsigned total_in_use = 0;
        unsigned total_stored = 0;
        unsigned total_allocated = 0;
        for (unsigned g = 0; g < _nallocated.size(); ++g) {
            unsigned long nbytes = _nelements[g]*_nallocated[g];
            unsigned nallocated = _nallocated[g];
            total_allocated += nallocated;
            unsigned nstored = (unsigned)_storage[g].size();
            total_stored += nstored;
            unsigned in_use = nallocated - nstored;
            total_in_use += in_use;
            memf << str(format("  %12d %12d %12d %12d %12d %12d\n") % g % _nelements[g] % in_use % nstored % nallocated % nbytes);
            total_bytes += nbytes;
        }
        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
        memf << str(format("  %12s %12s %12s %12d %12d %12d\n") % " " % " " % total_in_use % total_stored % total_allocated % total_bytes);
        memf << str(format("  Total megabytes: %.5f\n") % (1.0*total_bytes/1048576));
    }
    
 }
