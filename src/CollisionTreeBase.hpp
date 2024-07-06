#pragma once

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal, typename ExtInt>
    class CollisionTreeBase
    {
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
        static_assert(IntQ< Int>,"");
        
    public:
        
        using ClusterTreeBase_T = ClusterTreeBase<Real,Int,SReal,ExtReal,ExtInt>;
        
        using CollisionMatrix_T = SparseMatrixCSR<SReal,Int,Int>;
        
        CollisionTreeBase () {}
        
        virtual ~CollisionTreeBase() = default;
        
    protected:
        
        
    public:
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        virtual bool SymmetricQ() const = 0;

        virtual const CollisionMatrix_T & ClusterCollisionMatrix() const = 0;
        
        virtual cref<CollisionMatrix_T> PrimitiveCollisionMatrix() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetS() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetT() const = 0;

        virtual SReal MaximumSafeStepSize() const = 0;
        
//        virtual Int CollisionCount() const = 0;
        
//        virtual std::string Stats() const = 0;
        

        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CollisionTreeBase) + "<"+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }

    }; // CollisionTreeBase
    
} // namespace Repulsor
