#pragma once

namespace Repulsor
{
    template<
        typename Real_, typename Int_, typename LInt_, typename SReal_,
        typename ExtReal_, bool symmetricQetric
    >
    class CollisionTreeBase
    {
        static_assert(IntQ<   Int_>,"");
        static_assert(IntQ<  LInt_>,"");
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
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
        
        CollisionTreeBase ( cref<ClusterTreeBase_T> S, cref<ClusterTreeBase_T> T )
        {}

        
    protected:
        
        
    public:
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        static constexpr bool SymmetricQ()
        {
            return symmetricQetric;
        }

        virtual cref<CollisionMatrix_T> PrimitiveCollisionMatrix() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetS() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetT() const = 0;

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
