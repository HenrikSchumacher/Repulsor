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
        
        // Default constructor
        CollisionTreeBase() = default;
        // Destructor
        virtual ~CollisionTreeBase() = default;
        // Copy constructor
        CollisionTreeBase( const CollisionTreeBase & other ) = default;
        // Copy assignment operator
        CollisionTreeBase & operator=( const CollisionTreeBase & other ) = default;
        // Move constructor
        CollisionTreeBase( CollisionTreeBase && other ) = default;
        // Move assignment operator
        CollisionTreeBase & operator=( CollisionTreeBase && other ) = default;
        
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
            return TOOLS_TO_STD_STRING(CollisionTreeBase) + "<"+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }

    }; // CollisionTreeBase
    
} // namespace Repulsor
