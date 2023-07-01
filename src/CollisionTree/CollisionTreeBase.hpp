#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_, bool symmetricQetric>
    class CollisionTreeBase
    {
        ASSERT_FLOAT(Real_   );
        ASSERT_INT  (Int_    );
        ASSERT_INT  (LInt_   );
        ASSERT_FLOAT(SReal_  );
        ASSERT_FLOAT(ExtReal_);
        
    public:
        
        using Real              = Real_;
        using Int               = Int_;
        using LInt              = LInt_;
        using SReal             = SReal_;
        using ExtReal           = ExtReal_;
        
        using ClusterTreeBase_T = ClusterTreeBase<Real,Int,LInt,SReal,ExtReal>;
        
        using CollisionMatrix_T = Sparse::MatrixCSR<SReal,Int,LInt>;
        
        
        CollisionTreeBase () = default;
        
        virtual ~CollisionTreeBase() = default;
        
        CollisionTreeBase ( const ClusterTreeBase_T & restrict S, const ClusterTreeBase_T & restrict T )
        {}

        
    protected:
        
        
    public:
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        static constexpr bool SymmetricQ()
        {
            return symmetricQetric;
        }

        virtual const CollisionMatrix_T & PrimitiveCollisionMatrix() const = 0;
        
        virtual const ClusterTreeBase_T & GetS() const = 0;
        
        virtual const ClusterTreeBase_T & GetT() const = 0;

        virtual ExtReal MaximumSafeStepSize(
            const SReal t_,
            const SReal TOL
        ) const = 0;
        
        virtual std::string ClassName() const
        {
            return  std::string("CollisionTreeBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }

    }; // class CollisionTreeBase
    
} // namespace Repulsor
