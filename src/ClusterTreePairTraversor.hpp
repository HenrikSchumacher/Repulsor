#pragma once

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

        ClusterTreePairTraversor(
            cref<ClusterTree_T> S_,
            cref<ClusterTree_T> T_,
            mref<std::vector<Kernel_T>> kernels_
        )
        :   kernels (kernels_)
        ,   thread_count( static_cast<Int>(kernels.size()) )
        ,   S_C_L     ( S_.ClusterLeft().data()          )
        ,   S_C_R     ( S_.ClusterRight().data()         )
        ,   S_C_begin ( S_.ClusterBegin().data()         )
        ,   S_C_end   ( S_.ClusterEnd().data()           )
        ,   T_C_L     ( T_.ClusterLeft().data()          )
        ,   T_C_R     ( T_.ClusterRight().data()         )
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
        static constexpr SReal zero      = Scalar::Zero<SReal>;
    
        mref<std::vector<Kernel_T>> kernels;
        
        const Int thread_count = static_cast<Int>(1);
        
        BlockSplitMethod block_split_method = BlockSplitMethod::Parallel;

//        std::deque<Int> i_queue;
//        std::deque<Int> j_queue;
        
        std::deque<std::pair<Int,Int>> queue;
        
        
        cptr<Int> S_C_L;
        cptr<Int> S_C_R;
        cptr<Int> S_C_begin;
        cptr<Int> S_C_end;
        
        cptr<Int> T_C_L;
        cptr<Int> T_C_R;
        cptr<Int> T_C_begin;
        cptr<Int> T_C_end;
        
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
            TOOLS_PTIMER(timer,className()+"::Traverse");
            
            if( thread_count > 1 )
            {
                Traverse_Parallel();
            }
            else
            {
                Traverse_Sequential();
            }
        }
        
    protected:
        
#include "ClusterTreePairTraversor/Traverse_DepthFirst.hpp"
        
        void Traverse_Sequential()
        {
            TOOLS_PTIMER(timer,className()+"::Traverse_Sequential");
            
            Traverse_DepthFirst( Int(0), Int(0), Int(0) );
        }
        
    protected:
        
#include "ClusterTreePairTraversor/Traverse_BreadthFirst.hpp"
        
        void Traverse_Parallel()
        {
            TOOLS_PTIMER(timer,className()+"::Traverse_Parallel");
            
            Traverse_BreadthFirst( Int(0), Int(0), Int(0), static_cast<Int>(16) * ThreadCount() );
            
            
            ParallelDo_Dynamic(
                [this]( const Int thread, const Int k )
                {
                    TOOLS_DEBUG_PRINT(className()+"::Traverse_Parallel: Requesting kernel on thread " + ToString(thread) + " for task " + ToString(k) + "." );
                    
                    auto [i,j] = queue[static_cast<Size_T>(k)];
                    
                    Traverse_DepthFirst( thread, i, j );
                },
                Int(0), static_cast<Int>( queue.size() ), Int(1),
                ThreadCount()
            );
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
