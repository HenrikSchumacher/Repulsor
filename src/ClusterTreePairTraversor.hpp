#pragma once

#include "ClusterTreePairTraversor/ClusterTreePairTraversor_Kernel.hpp"

namespace Repulsor
{
    template<typename Kernel_T, bool symmetricQ, bool leaves_are_singletons>
    class ClusterTreePairTraversor
    {
    public:
        
        using ClusterTree_T = typename Kernel_T::ClusterTree_T;
        using Real          = typename Kernel_T::Real;
        using Int           = typename Kernel_T::Int;
        using SReal         = typename Kernel_T::SReal;
        using ExtReal       = typename Kernel_T::ExtReal;
        
    public:

        // In order to prevent GetS() and GetT() shooting a segfault, we have to initialize S and T here. This is the only case in which Traversor owns these raw pointers.
        
        ClusterTreePairTraversor(
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
            if constexpr ( symmetricQ )
            {
                if( std::addressof(S_) != std::addressof(T_) )
                {
                    eprint(className()+": symmetricQ == true, bu S != T.");
                }
            }
        } // Constructor
        
        ~ClusterTreePairTraversor() = default;
        
    protected:

        static constexpr Int   null      = static_cast<Int>(0);
        static constexpr Int   max_depth = 128;
        static constexpr SReal zero      = static_cast<SReal>(0);
    
        std::vector<Kernel_T> & kernels;
        
        const Int thread_count = static_cast<Int>(1);
        
        BlockSplitMethod block_split_method = BlockSplitMethod::Parallel;

        std::deque<Int> i_queue;
        std::deque<Int> j_queue;
        
        
        ptr<Int> S_C_left;
        ptr<Int> S_C_right;
        ptr<Int> S_C_begin;
        ptr<Int> S_C_end;
        
        ptr<Int> T_C_left;
        ptr<Int> T_C_right;
        ptr<Int> T_C_begin;
        ptr<Int> T_C_end;
        
    public:
        
        Int ThreadCount() const
        {
            return thread_count;
        }
        
        constexpr bool SymmetricQ() const
        {
            return symmetricQ;
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
        
#include "ClusterTreePairTraversor/Traverse_DepthFirst.hpp"
        
        void Traverse_Sequential()
        {
            ptic(className()+"::Traverse_Sequential");
            
            Traverse_DepthFirst( 0, 0, 0 );
            
            ptoc(className()+"::Traverse_Sequential");
        }
        
    protected:
        
#include "ClusterTreePairTraversor/Traverse_BreadthFirst.hpp"
        
        void Traverse_Parallel()
        {
            ptic(className()+"::Traverse_Parallel");
            
            Traverse_BreadthFirst( 0, 0, 0, static_cast<Int>(16) * ThreadCount() );
            
            
            ParallelDo_Dynamic(
                [this]( const Int thread, const Int k )
                {
                    debug_print(className()+"::Traverse_Parallel: Requesting kernel on thread " + ToString(thread) + " for task " + ToString(k) + "." );
                    
                    Traverse_DepthFirst( thread, i_queue[k], j_queue[k] );
                },
                Scalar::Zero<Int>, static_cast<Int>( i_queue.size() ), Scalar::One<Int>,
                ThreadCount()
            );
            
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
            return "ClusterTreePairTraversor<"+kernels[0].ClassName()+">";
        }
      
        
    }; // class ClusterTreePairTraversor
    
} //namespace Repulsor
