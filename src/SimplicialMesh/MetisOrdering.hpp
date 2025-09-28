public:

    virtual cref<Permutation<Int>> MetisOrdering() const override
    {
        std::string tag ("MetisOrdering");
        
        if( !this->InPersistentCacheQ(tag))
        {
            TOOLS_PTIMER(timer,ClassName()+"::"+tag);

            Permutation<Int> perm;
            
#ifdef TENSORS_HAS_METIS
            
            auto & A = H1Metric();
            
            perm = Sparse::Metis<Int>()(
                A.Outer().data(), A.Inner().data(), A.RowCount(), ThreadCount()
            );
#else
            wprint(ClassName()+"::MetisOrdering: Sparse::Metis not available. Using default NestedDissectionOrdering.");
                    
            perm = NestedDissectionOrdering();
#endif
            
            this->SetPersistentCache( tag, std::move(perm) );
        }
        
        return this->template GetPersistentCache<Permutation<Int>>(tag);
    }
