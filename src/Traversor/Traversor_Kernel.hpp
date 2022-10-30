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
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        using Int     = typename ClusterTree_T::Int;
        
    public:
        
        CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_
        )
        :   tree_string( S_.ClassName() )
        {}
        
        CLASS(
            const CLASS & other
        )
        :   tree_string ( other.tree_string                )
        {}
        
        virtual ~CLASS() = default;
        
        
        
        friend void swap(CLASS &A, CLASS &B)
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;

            swap( A.tree_string, B.tree_string );
        }
        
//        // Copy assignment
//        CLASS & operator=(CLASS other) // Intentionally no pass-by-reference here!
//        {
//            swap(*this, other);
//
//            return *this;
//        }
        
        // Move constructor
        CLASS( CLASS && other ) noexcept
        :   CLASS()
        {
            swap(*this, other);
        }

        /* Move assignment operator */
        CLASS & operator=( CLASS && other ) noexcept
        {
            if( this != &other )
            {
                swap( *this, other );
            }
            return *this;
        }
        
    protected:
        
        std::string tree_string;
        
    public:

        virtual void LoadClusterS( const Int i ) = 0;
        virtual void LoadClusterT( const Int j ) = 0;
        
        virtual void LoadPrimitiveS( const Int i ) = 0;
        virtual void LoadPrimitiveT( const Int j ) = 0;
        
        virtual bool IsAdmissable() = 0;
        
        virtual SReal ClusterScoreS() = 0;
        virtual SReal ClusterScoreT() = 0;
        
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
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS) + "<"+tree_string+">";
        }
      
        
    };
    
} //namespace Repulsor

#undef CLASS

