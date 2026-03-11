public:

    virtual cref<Permutation_T> ApproximateMinimumDegreeOrdering() const override
    {
        std::string tag ("ApproximateMinimumDegreeOrdering");
        
        if( !this->InPersistentCacheQ(tag))
        {
            TOOLS_PTIMER(timer,ClassName()+"::"+tag);
            
            Permutation_T perm;
            
#ifdef TENSORS_HAS_AMD
            
            auto & A = H1Metric();
            
            perm = Sparse::ApproximateMinimumDegree<Int,Parallel>()(
                A.Outer().data(), A.Inner().data(), A.RowCount(), ThreadCount()
            );
#else
            wprint(ClassName()+"::ApproximateMinimumDegreeOrdering: Sparse::ApproximateMinimumDegree not available. Using default NestedDissectionOrdering.");
                    
            perm = NestedDissectionOrdering();
#endif
            
            this->SetPersistentCache( tag, std::move(perm) );
        }
        
        return this->template GetPersistentCache<Permutation_T>(tag);
    }
