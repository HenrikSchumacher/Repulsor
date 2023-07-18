#pragma once

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class alignas(ObjectAlignment) ClusterTreePairTraversor_Kernel
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        using Int     = typename ClusterTree_T::Int;
        using LInt    = typename ClusterTree_T::LInt;
        
    public:
        
        ClusterTreePairTraversor_Kernel() = default;
        
        ClusterTreePairTraversor_Kernel(
            const ClusterTree_T & restrict S_,
            const ClusterTree_T & restrict T_
        )
//        :   tree_string( S_.ClassName() )
        :   tree_string( "..." )
        {}
        
        ClusterTreePairTraversor_Kernel(
            const ClusterTreePairTraversor_Kernel & restrict other
        )
        :   tree_string ( other.tree_string                )
        {}
        
        ~ClusterTreePairTraversor_Kernel() = default;
        
        
        friend void swap(ClusterTreePairTraversor_Kernel &A, ClusterTreePairTraversor_Kernel &B)
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;

            swap( A.tree_string, B.tree_string );
        }
        
//        // Copy assignment
//        ClusterTreePairTraversor_Kernel & operator=(ClusterTreePairTraversor_Kernel other) // Intentionally no pass-by-reference here!
//        {
//            swap(*this, other);
//
//            return *this;
//        }
        
        // Move constructor
        ClusterTreePairTraversor_Kernel( ClusterTreePairTraversor_Kernel && other ) noexcept
        :   ClusterTreePairTraversor_Kernel()
        {
            swap(*this, other);
        }

        /* Move assignment operator */
        ClusterTreePairTraversor_Kernel & operator=( ClusterTreePairTraversor_Kernel && other ) noexcept
        {
            if( this != &other )
            {
                swap( *this, other );
            }
            return *this;
        }
        
    protected:
        
        std::string tree_string;
        
    private:
        
        std::string ClassName() const
        {
            return std::string("ClusterTreePairTraversor_Kernel") + "<"+tree_string+">";
        }
      
        
    }; // class ClusterTreePairTraversor_Kernel
    
} //namespace Repulsor

