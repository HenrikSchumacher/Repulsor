#pragma once

#define CLASS CollisionTreeBase

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal, bool is_symmetric>
    class CLASS
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);        
        
    public:
        
        using ClusterTreeBase_T = ClusterTreeBase<Real,Int,SReal,ExtReal>;
        
        using CollisionMatrix_T = SparseMatrixCSR<SReal,Int,size_t>;
        
        CLASS () = default;
        
        virtual ~CLASS() = default;
        
        CLASS ( const ClusterTreeBase_T & S, const ClusterTreeBase_T & T )
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
        
//        virtual Int CollisionCount() const = 0;
        
//        virtual std::string Stats() const = 0;
        

        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS) + "<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }

    }; // CLASS
    
} // namespace Repulsor

#undef CLASS
