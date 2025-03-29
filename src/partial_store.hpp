#pragma once

namespace proj {

    struct Partial {
        Partial(unsigned g, unsigned n);
        ~Partial();
        unsigned        _g; // the locus
        vector<double>  _v; // the partial array: length = _nstates*<no. patterns>
    };

    inline Partial::Partial(unsigned g, unsigned n) {
        _g = g;
        _v.resize(n);
    }
    
    inline Partial::~Partial() {
    }

    class PartialStore {
    
        public:
                                    PartialStore();
                                    ~PartialStore();
                                    
            typedef shared_ptr<Partial>          partial_t;
            
            typedef vector<partial_t>            vect_partial_t;
            typedef vector<vect_partial_t>       leaf_partials_t;
            typedef vector<vect_partial_t>       storage_t;
            
            void            setNLoci(unsigned nloci);
            partial_t       getPartial(unsigned locus);
            void            putPartial(unsigned locus, partial_t partial);
            void            setNElements(unsigned nelements, unsigned locus);
            unsigned        getNElements(unsigned locus) const {return _nelements[locus];}
            
            void            debugReport() const;
            
        private:
        
            vector<unsigned> _nelements;
            storage_t        _storage;
            
            unsigned         _total_partials_created;
            unsigned         _total_elements_created;
            unsigned         _total_partials_reused;
            unsigned         _total_elements_reused;
    };

    inline PartialStore::PartialStore() {
        _total_partials_created = 0;
        _total_elements_created = 0;
        _total_partials_reused = 0;
        _total_elements_reused = 0;
    }

    inline PartialStore::~PartialStore() {
        _nelements.clear();
        _storage.clear();
    }
    
    inline void PartialStore::setNLoci(unsigned nloci) {
        // Should be called before any partials are stored
        assert(_nelements.empty());
        assert(_storage.empty());
        
        // Resize both containers
        _nelements.resize(nloci);
        _storage.resize(nloci);
    }

    inline void PartialStore::setNElements(unsigned nelements, unsigned locus) {
        assert(_nelements.size() > locus);
        _nelements[locus] = nelements;
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned locus) {
        // Check to make sure supplied value of locus is valid
        assert(_nelements.size() > locus);
        assert(_nelements[locus] > 0);
        
        partial_t partial;
#if defined(REUSE_PARTIALS)
        if (_storage[locus].empty()) {
            // No stored partials for this locus, so allocate one
            partial = partial_t(new Partial(locus, _nelements[locus]));
            _total_partials_created++;
            _total_elements_created += _nelements[locus];
        }
        else {
            size_t n = _storage[locus].size();
            partial = _storage[locus].at(n-1);
            _storage[locus].pop_back();
            _total_partials_reused++;
            _total_elements_reused += _nelements[locus];
        }
#else
        // No partials are being stored in this version of the program
        // so allocate a new one
        partial = partial_t(new Partial(locus, _nelements[locus]));
        _total_partials_created++;
        _total_elements_created += _nelements[locus];
#endif
        return partial;
    }
    
    inline void PartialStore::putPartial(unsigned locus, partial_t partial) {
#if defined(REUSE_PARTIALS)
        // Check to make sure supplied value of locus is valid
        assert(_nelements.size() > locus);
        assert(_nelements[locus] > 0);

        // Store the partial for later
        assert(partial->_v.size() == _nelements[locus]);
        partial->_v.assign(_nelements[locus], 0.0);
        _storage[locus].push_back(partial);
#endif
    }
        
    inline void PartialStore::debugReport() const {
        output("\nPartialStore report:\n", G::LogCateg::DEBUGGING);
        output(format("  %12d partials created\n") % _total_partials_created, G::LogCateg::DEBUGGING);
        output(format("  %12d partials reused\n") % _total_partials_reused, G::LogCateg::DEBUGGING);
        output(format("  %12d elements created\n") % _total_elements_created, G::LogCateg::DEBUGGING);
        output(format("  %12d elements reused\n") % _total_elements_reused, G::LogCateg::DEBUGGING);
        //if (_total_partials_reused > 0) {
            output(format("  %12s %12s %12s\n") % "Locus" % "Stored" % "Size", G::LogCateg::DEBUGGING);
            unsigned partials_in_storage = 0;
            unsigned elements_in_storage = 0;
            for (unsigned i = 0; i < _storage.size(); i++) {
                unsigned n = (unsigned)_storage[i].size();
                unsigned nelem = n*_nelements[i];
                partials_in_storage += n;
                elements_in_storage += nelem;
                output(format("  %12d %12d %12d\n") % (i+1) % n % nelem, G::LogCateg::DEBUGGING);
            }
            output(format("  %12s %12d %12d\n\n") % "" % partials_in_storage % elements_in_storage, G::LogCateg::DEBUGGING);
        //}
    }
 }
