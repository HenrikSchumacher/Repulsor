#pragma once

#define CLASS Traversor_Kernel

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class alignas(OBJECT_ALIGNMENT) CLASS
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
    public:
        
        virtual ~CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_
        )
        :   S(S_)
        ,   T(T_)
        ,   tree_string( S_.ClassName() )
        {
            ptic(className());
            
            ptoc(className());
        } // Constructor
        
    protected:
        
        const ClusterTree_T & S;
        const ClusterTree_T & T;

        const std::string tree_string;
        
    public:

        virtual void LoadClusterS( const Int i ) = 0;
        virtual void LoadClusterT( const Int j ) = 0;
        
        virtual void LoadPrimitiveS( const Int i ) = 0;
        virtual void LoadPrimitiveT( const Int j ) = 0;
        
        virtual bool IsAdmissable() = 0;
        
        virtual bool ClusterScoreS() = 0;
        virtual bool ClusterScoreT() = 0;
        
        virtual void ComputeLeafDiagonal() = 0;
        virtual void ComputeLeaf() = 0;
        virtual void ComputeLeafSwapped() = 0;
        
        virtual void ComputeAdmissable() = 0;
        virtual void ComputeAdmissableSwapped() = 0;

        
    protected:
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    private:
        
        std::string className()
        {
            return TO_STD_STRING(CLASS) + "<"+tree_string+">";
        }
      
        
    };
    
} //namespace Repulsor

#undef CLASS

