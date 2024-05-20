#pragma once

namespace Repulsor
{
    
    template<typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_, bool symmetricQ>
    class BlockClusterTreeBase
    {
        static_assert(IntQ< Int_>,"");
        static_assert(IntQ<LInt_>,"");
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
    public:
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = LInt_;
        
        using Setting_T          = BlockClusterTreeSettings;

        using Inter_Pattern_T    = Sparse::BinaryMatrixCSR<Int,LInt>;
        using VeryNear_Pattern_T = Sparse::BinaryMatrixCSR<Int,LInt>;
        using Near_Pattern_T     = Sparse::BinaryMatrixCSR<Int,LInt>;
        using Far_Pattern_T      = Sparse::BinaryMatrixCSR<Int,LInt>;
        
        using ClusterTreeBase_T  = ClusterTreeBase<Real,Int,LInt,SReal,ExtReal>;

    protected:

    public:
        
        virtual ~BlockClusterTreeBase() = default;
        
        virtual Int AmbDim() const = 0;

        virtual Int ThreadCount() const = 0;
        
        virtual Real FarFieldSeparationParameter() const = 0;
        
        virtual Real NearFieldSeparationParameter() const = 0;
        
        static constexpr bool SymmetricQ()
        {
            return symmetricQ;
        }
        
        virtual LInt PrimitiveIntersectionCount() const = 0;
        
        virtual LInt VeryNearFieldInteractionCount() const = 0;
        
        virtual LInt NearFieldInteractionCount() const = 0;

        virtual LInt FarFieldInteractionCount() const = 0;
        
        virtual cref<Inter_Pattern_T> PrimitiveIntersectionMatrix() const = 0;

        virtual cref<VeryNear_Pattern_T> VeryNear() const = 0;
        
        virtual cref<Near_Pattern_T> Near() const = 0;
        
        virtual cref<Far_Pattern_T> Far() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetS() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetT() const = 0;

        virtual cref<Setting_T> Settings() const = 0;
        
        virtual std::string Stats() const = 0;
        
    public:
        
        virtual std::string ClassName() const
        {
            return  std::string("BlockClusterTreeBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+ToString(symmetricQ)+">";
        }
        
    }; // class BlockClusterTreeBase
    
} //namespace Repulsor
