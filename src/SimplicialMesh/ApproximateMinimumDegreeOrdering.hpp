public:

    virtual cref<Permutation<Int>> ApproximateMinimumDegreeOrdering() const override
    {
        std::string tag ("ApproximateMinimumDegreeOrdering");
        
        if( !this->InPersistentCacheQ(tag))
        {
            ptic(ClassName()+"::"+tag);
            
            Permutation<Int> perm;
            
#ifdef TENSORS_HAS_AMD
            
            auto & A = H1Metric();
            
            perm = Sparse::ApproximateMinimumDegree<Int>()(
                A.Outer().data(), A.Inner().data(), A.RowCount(), ThreadCount()
            );
#else
            wprint(ClassName()+"::ApproximateMinimumDegreeOrdering: Sparse::ApproximateMinimumDegree not available. Using default NestedDissectionOrdering.");
                    
            perm = NestedDissectionOrdering();
#endif
            
            this->SetPersistentCache( tag, std::move(perm) );
            
            ptoc(ClassName()+"::"+tag);
        }
        
        return this->template GetPersistentCache<Permutation<Int>>(tag);
    }
