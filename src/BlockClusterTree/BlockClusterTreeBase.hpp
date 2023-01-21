#pragma once

namespace Repulsor
{
    
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_, bool is_symmetric>
    class BlockClusterTreeBase
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
        
        using Setting_T          = BlockClusterTreeSettings;

        using Inter_Pattern_T    = Sparse::BinaryMatrixCSR<Int,LInt>;
        using VeryNear_Pattern_T = Sparse::BinaryMatrixCSR<Int,LInt>;
        using Near_Pattern_T     = Sparse::BinaryMatrixCSR<Int,LInt>;
        using Far_Pattern_T      = Sparse::BinaryMatrixCSR<Int,LInt>;
        
        using ClusterTreeBase_T  = ClusterTreeBase<Real,Int,SReal,ExtReal>;

    protected:

    public:
        
        virtual ~BlockClusterTreeBase() = default;
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        virtual Real FarFieldSeparationParameter() const = 0;
        
        virtual Real NearFieldSeparationParameter() const = 0;
        
        static constexpr bool IsSymmetric()
        {
            return is_symmetric;
        }
        
        virtual LInt PrimitiveIntersectionCount() const = 0;
        
        virtual LInt VeryNearFieldInteractionCount() const = 0;
        
        virtual LInt NearFieldInteractionCount() const = 0;

        virtual LInt FarFieldInteractionCount() const = 0;
        
        virtual const Inter_Pattern_T & PrimitiveIntersectionMatrix() const = 0;

        virtual const VeryNear_Pattern_T & VeryNear() const = 0;
        
        virtual const Near_Pattern_T & Near() const = 0;
        
        virtual const Far_Pattern_T & Far() const = 0;
        
        virtual const ClusterTreeBase_T & GetS() const = 0;
        
        virtual const ClusterTreeBase_T & GetT() const = 0;

        virtual const Setting_T & Settings() const = 0;
        
        virtual std::string Stats() const = 0;
        
    public:
        
        virtual std::string ClassName() const
        {
            return  "BlockClusterTreeBase<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+","+ToString(is_symmetric)+">";
        }
        
    }; // class BlockClusterTreeBase
    
} //namespace Repulsor
