#pragma once

namespace proj {

    class PartialStore {
    
        public:
                                    PartialStore();
                                    ~PartialStore();
            
            typedef shared_ptr<vector<double> > partial_t;
            typedef vector<partial_t>           vect_partial_t;
            
            partial_t  getPartial();
            void       getVectorOfPartials(unsigned n, PartialStore::vect_partial_t & v);
            void       stowPartial(partial_t p);
            void       setNElements(unsigned nelements);
            unsigned   getNElements() const {return _nelements;}
            
        private:
        
            unsigned _nelements;
            vector<partial_t> storage;
    };

    inline PartialStore::PartialStore() {
        _nelements = 0;
    }

    inline PartialStore::~PartialStore() {
        storage.clear();
    }

    inline void PartialStore::setNElements(unsigned nelements) {
        assert(nelements > 0);
        _nelements = nelements;
    }

    inline PartialStore::partial_t PartialStore::getPartial() {
        assert(_nelements > 0);
        if (storage.size() > 0) {
            partial_t last = storage.back();
            storage.pop_back();
            return last;
        }
        else
            return partial_t(new vector<double>(_nelements));
    }
    
    inline void PartialStore::getVectorOfPartials(unsigned n, PartialStore::vect_partial_t & v) {
        assert(_nelements > 0);
        v.resize(n);
        unsigned i = 0;
        
        // First use up any stored partials
        while (storage.size() > 0) {
            partial_t last = storage.back();
            storage.pop_back();
            v[i++] = last;
        }
        
        // Allocate remaining partials if necessary
        while (i < n) {
            v[i++] = partial_t(new vector<double>(_nelements));
        }
    }
    
    inline void PartialStore::stowPartial(partial_t p) {
        fill(p->begin(), p->end(), 0.0);
        storage.push_back(p);
    }
    
 }
