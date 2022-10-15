#pragma once

#include "Traversor/Traversor_Kernel.hpp"

#define CLASS Traversor

namespace Repulsor
{
    template<typename Kernel_T, bool is_symmetric, bool leafs_are_singletons >
    class CLASS
    {
        
    public:
        
        using ClusterTree_T = typename Kernel_T::ClusterTree_T;
        using Real          = typename Kernel_T::Real;
        using Int           = typename Kernel_T::Int;
        using SReal         = typename Kernel_T::SReal;
        using ExtReal       = typename Kernel_T::ExtReal;
        
    public:

        // In order to prevent GetS() and GetT() shooting a segfault, we have to initialize S and T here. This is the only case in which CLASS owns these raw pointers.
        
        virtual ~CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_,
            std::vector<Kernel_T> & kernels_
        )
        :   kernels (kernels_)
        ,   thread_count( static_cast<Int>(kernels.size()) )
        ,   S_C_left  ( S_.ClusterLeft().data() )
        ,   S_C_right ( S_.ClusterRight().data() )
        ,   S_C_begin ( S_.Begin().data() )
        ,   S_C_end   ( S_.End().data() )
        ,   T_C_left  ( T_.ClusterLeft().data() )
        ,   T_C_right ( T_.ClusterRight().data() )
        ,   T_C_begin ( T_.Begin().data() )
        ,   T_C_end   ( T_.End().data() )
        {
            ptic(className());
            
            if constexpr ( is_symmetric )
            {
                assert( std::addressof(S_) == std::addressof(T_) );
            }
            
            ptoc(className());
        } // Constructor
        
    protected:

        static constexpr Int zero = static_cast<Int>(0);
        static constexpr Int max_depth = 128;
        static constexpr Int null = static_cast<Int>(0);
    
        std::vector<Kernel_T> kernels;
        
        const Int thread_count = static_cast<Int>(1);
        
        BlockSplitMethod block_split_method = BlockSplitMethod::Parallel;

        std::deque<Int> i_queue;
        std::deque<Int> j_queue;
        
        const Int * restrict const S_C_left;
        const Int * restrict const S_C_right;
        const Int * restrict const S_C_begin;
        const Int * restrict const S_C_end;
        
        const Int * restrict const T_C_left;
        const Int * restrict const T_C_right;
        const Int * restrict const T_C_begin;
        const Int * restrict const T_C_end;
        
    public:
        
        Int ThreadCount() const
        {
            return thread_count;
        }
        
        constexpr bool IsSymmetric() const
        {
            return is_symmetric;
        }
        
    protected:
        
        std::string ClassName() const
        {
            return className();
        }
        
    private:
        
        std::string className()
        {
            return TO_STD_STRING(CLASS) + "<"+kernels[0].ClassName()+">";
        }
      
        
    };
    
} //namespace Repulsor

#undef BASE
#undef CLASS
