public:

    template<typename Real, typename ExtInt>
    void HilbertOrdering( cptr<Real> X, mptr<ExtInt> perm, const Int point_count )
    {
        TOOLS_PTIMER(timer,ClassName()+"::HilbertOrdering");
        Tensor1<HilbertCode_T,Int> M = CoordsToHilbert( X, point_count );
        Ordering( M.data(), perm, point_count );
    }


    template<typename Real>
    Tensor1<Int,Int> HilbertOrdering( cptr<Real> X, const Int point_count )
    {
        Tensor1<Int,Int> perm ( point_count );
        HilbertOrdering( X, perm.data(), point_count );
        return perm;
    }
