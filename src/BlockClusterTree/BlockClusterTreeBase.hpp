#pragma once

#define CLASS BlockClusterTreeBase

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
        
        using Setting_T         = BlockClusterTreeSettings;

        using Inter_T           = SparseBinaryMatrixCSR<Int>;
        using VeryNear_T        = SparseBinaryMatrixCSR<Int>;
        using Near_T            = SparseBinaryMatrixCSR<Int>;
        using Far_T             = SparseBinaryMatrixCSR<Int>;
    
        using ClusterTreeBase_T = ClusterTreeBase<Real,Int,SReal,ExtReal>;

    protected:

    public:
        
        virtual ~CLASS() = default;
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        virtual Real FarFieldSeparationParameter() const = 0;
        
        virtual Real NearFieldSeparationParameter() const = 0;
        
        static constexpr bool IsSymmetric()
        {
            return is_symmetric;
        }
        
        virtual Int PrimitiveIntersectionCount() const = 0;
        
        virtual Int VeryNearFieldInteractionCount() const = 0;
        
        virtual Int NearFieldInteractionCount() const = 0;

        virtual Int FarFieldInteractionCount() const = 0;
        
        virtual const Inter_T & PrimitiveIntersectionMatrix() const = 0;

        virtual const VeryNear_T & VeryNear() const = 0;
        
        virtual const Near_T & Near() const = 0;
        
        virtual const Far_T & Far() const = 0;
        
        virtual const ClusterTreeBase_T & GetS() const = 0;
        
        virtual const ClusterTreeBase_T & GetT() const = 0;

        virtual const Setting_T & Settings() const = 0;
        
        virtual std::string Stats() const = 0;
        
    public:
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS) + "<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+","+ToString(is_symmetric)+">";
        }
        
    };
    
} //namespace Repulsor

#undef CLASS
