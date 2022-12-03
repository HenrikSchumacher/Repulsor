#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_, bool is_symmetric>
    class CollisionTreeBase
    {
        ASSERT_FLOAT(Real_   );
        ASSERT_INT  (Int_    );
        ASSERT_FLOAT(SReal_  );
        ASSERT_FLOAT(ExtReal_);
        
    public:
        
        using Real              = Real_;
        using Int               = Int_;
        using SReal             = SReal_;
        using ExtReal           = ExtReal_;
        
        using ClusterTreeBase_T = ClusterTreeBase<Real,Int,SReal,ExtReal>;
        
        using CollisionMatrix_T = SparseMatrixCSR<SReal,Int,size_t>;
        
        CollisionTreeBase () = default;
        
        virtual ~CollisionTreeBase() = default;
        
        CollisionTreeBase ( const ClusterTreeBase_T & S, const ClusterTreeBase_T & T )
        {}

        
    protected:
        
        
    public:
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        static constexpr bool IsSymmetric()
        {
            return is_symmetric;
        }

        virtual const CollisionMatrix_T & PrimitiveCollisionMatrix() const = 0;
        
        virtual const ClusterTreeBase_T & GetS() const = 0;
        
        virtual const ClusterTreeBase_T & GetT() const = 0;

        virtual ExtReal MaximumSafeStepSize( const SReal t_) const = 0;
        
        virtual std::string ClassName() const
        {
            return  "CollisionTreeBase<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }

    }; // class CollisionTreeBase
    
} // namespace Repulsor
