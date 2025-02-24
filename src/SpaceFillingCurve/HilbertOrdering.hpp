public:

    template<typename Real, typename ExtInt>
    void HilbertOrdering( cptr<Real> X, mptr<ExtInt> perm, const Int point_count )
    {
        TOOLS_PTIC(ClassName()+"::HilbertOrdering");

        Tensor1<HilbertCode_T,Int> M = CoordsToHilbert( X, point_count );
        
        TOOLS_PTIC("Ordering");
        Ordering( M.data(), perm, point_count );
        TOOLS_PTOC("Ordering");
        
        TOOLS_PTOC(ClassName()+"::HilbertOrdering");
    }


    template<typename Real>
    Tensor1<Int,Int> HilbertOrdering( cptr<Real> X, const Int point_count )
    {
        Tensor1<Int,Int> perm ( point_count );
        
        HilbertOrdering( X, perm.data(), point_count );
        
        return perm;
    }
