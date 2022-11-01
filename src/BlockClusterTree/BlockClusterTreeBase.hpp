#pragma once

#define CLASS BlockClusterTreeBase

namespace Repulsor
{
    
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_, bool is_symmetric>
    class CLASS
    {
        ASSERT_FLOAT(Real_   );
        ASSERT_INT  (Int_    );
        ASSERT_FLOAT(SReal_  );
        ASSERT_FLOAT(ExtReal_);
        
    public:
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        
        using LInt    = size_t;
//        using LInt    = uint32_t;
        
        using Setting_T         = BlockClusterTreeSettings;

        using Pattern_T         = SparseBinaryMatrixCSR<Int,LInt>;
        
        using Values_T          = Tensor2<Real,LInt>;
        using ValueContainer_T  = std::unordered_map<std::string,Values_T>;
        
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
        
        virtual const Pattern_T & PrimitiveIntersectionMatrix() const = 0;

        virtual const Pattern_T & VeryNear() const = 0;
        
        virtual const Pattern_T & Near() const = 0;
        
        virtual const Pattern_T & Far() const = 0;
        
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
