#pragma once

#include "Traversor/Traversor_Kernel.hpp"

#define CLASS Traversor

namespace Repulsor
{
    template<typename Kernel_T, bool is_symmetric, bool leaves_are_singletons >
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
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_,
            std::vector<Kernel_T> & kernels_
        )
        :   kernels (kernels_)
        ,   thread_count( static_cast<Int>(kernels.size()) )
        ,   S_C_left  ( S_.ClusterLeft().data()          )
        ,   S_C_right ( S_.ClusterRight().data()         )
        ,   S_C_begin ( S_.ClusterBegin().data()         )
        ,   S_C_end   ( S_.ClusterEnd().data()           )
        ,   T_C_left  ( T_.ClusterLeft().data()          )
        ,   T_C_right ( T_.ClusterRight().data()         )
        ,   T_C_begin ( T_.ClusterBegin().data()         )
        ,   T_C_end   ( T_.ClusterEnd().data()           )
        {
            if constexpr ( is_symmetric )
            {
                assert( std::addressof(S_) == std::addressof(T_) );
            }
        } // Constructor
        
        virtual ~CLASS() = default;
        
    protected:

        static constexpr Int   null      = static_cast<Int>(0);
        static constexpr Int   max_depth = 128;
        static constexpr SReal zero      = static_cast<SReal>(0);
    
        std::vector<Kernel_T> & kernels;
        
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
        
    public:
        
        void Traverse()
        {
            ptic(className()+"::Traverse");
            
            if( thread_count > 1 )
            {
                Traverse_Parallel();
            }
            else
            {
                Traverse_Sequential();
            }
            
            ptoc(className()+"::Traverse");
        }
        
    protected:
        
#include "Traversor/Traverse_DepthFirst.hpp"
        
        void Traverse_Sequential()
        {
            ptic(className()+"::Traverse_Sequential");
            
            Traverse_DepthFirst(0,0);
            
            ptoc(className()+"::Traverse_Sequential");
        }
        
    protected:
        
#include "Traversor/Traverse_BreadthFirst.hpp"
        
        void Traverse_Parallel()
        {
            ptic(className()+"::Traverse_Parallel");

            Traverse_BreadthFirst( 0, 0, static_cast<Int>(4) * ThreadCount() * ThreadCount() );

            #pragma omp parallel for num_threads( thread_count ) schedule( dynamic )
            for( Int k = 0; k < static_cast<Int>(i_queue.size()); ++k )
            {
                Traverse_DepthFirst(i_queue[k], j_queue[k]);
            }
            
            ptoc(className()+"::Traverse_Parallel");
        }

    private:
        
        std::string ClassName() const
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS) + "<"+kernels[0].ClassName()+">";
        }
      
        
    };
    
} //namespace Repulsor

#undef BASE
#undef CLASS
