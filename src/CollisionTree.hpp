#pragma once

#include "CollisionTree/CollisionTreeBase.hpp"

#define CLASS CollisionTree
#define BASE CollisionTreeBase<Real,Int,SReal,ExtReal,is_symmetric>

namespace Repulsor
{
    
    template <int AMB_DIM, typename SReal, typename Int>
        inline void AABB_CollisionTimeInterval(
           const SReal * restrict const p_0,
           const SReal * restrict const p_1,
           const SReal * restrict const q_0,
           const SReal * restrict const q_1,
           SReal & t_first,    // returning per reference to avoid std::pair
           SReal & t_last      // returning per reference to avoid std::pair
        )
        {
            // Assuming p_0, p_1, q_0, q_1 are serizalizations in the format of some AABB<AMB_DIM,GJK_Real,Int,SReal>

            
            // Returns interval {t_first,t_last} of those t in [0,1] such that the AABBs defined by
            //     P(t) = (1-t) * p_0 + t * p_1
            // and
            //     Q(t) = (1-t) * q_0 + t * q_1
            // intersect.
            
            static constexpr SReal zero = static_cast<SReal>(0);
            static constexpr SReal one  = static_cast<SReal>(1);
            
            SReal t_0 = zero;
            SReal t_1 = one;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                // Compute the intervals of the four AABB in the k-th coordinate direction.
                
                const Int i = 1+k;
                const Int j = 1+k+AMB_DIM;
                
                const SReal P_min_0 = p_0[i] - p_0[j];
                const SReal P_max_0 = p_0[i] + p_0[j];
                
                const SReal P_min_1 = p_1[i] - p_1[j];
                const SReal P_max_1 = p_1[i] + p_1[j];
                
                const SReal Q_min_0 = q_0[i] - q_0[j];
                const SReal Q_max_0 = q_0[i] + q_0[j];
                
                const SReal Q_min_1 = q_1[i] - q_1[j];
                const SReal Q_max_1 = q_1[i] + q_1[j];
                
                bool initially_intersecting = (P_min_0 <= Q_max_0) && (Q_min_0 <= P_max_0);
    //            bool finally_intersecting   = (P_min_1 <= Q_max_1) && (Q_min_1 <= P_max_1);
                
                const SReal m_0 = Q_max_0 - Q_max_1 - P_min_0 + P_min_1;
                const SReal m_1 = P_max_0 - P_max_1 - Q_min_0 + Q_min_1;
                
                SReal A = ( m_0 != zero ) ? std::min(one,(Q_max_0 - P_min_0)/m_0) : one;
                SReal B = ( m_1 != zero ) ? std::min(one,(P_max_0 - Q_min_0)/m_1) : one;
                
                if( A > B )
                {
                    std::swap(A,B);
                }
                
                SReal a;
                SReal b;
                
    //            GJK_DUMP(initially_intersecting);
    //            GJK_DUMP(A);
    //            GJK_DUMP(B);
                
                if( initially_intersecting )
                {
                    a = zero;
                    b = (A>zero) ? A : ( (B<=zero) ? one : B );
                }
                else
                {
                    // initially_intersecting == false
                    
                    if( /*A < zero &&*/ B < zero )
                    {
                        a = one;
                        b = one;
                    }
                    else
                    {
                        if( A >= zero/* && B >= zero*/ )
                        {
                            a = A;
                            b = B;
                        }
                        else
                        {
                            a = B;
                            b = one;
                        }
                    }
                }
                
                t_0 = std::max(t_0,a);
                t_1 = std::min(t_1,b);
            }
            
            if( t_1 < t_0 )
            {
                // empty intersection interval
                t_first = one;
                t_last  = one;
            }
            else
            {
                t_first = t_0;
                t_last  = t_1;
            }
        }
    
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal, bool is_symmetric>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);        
        
    public:
        
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        using CollisionMatrix_T = typename BASE::CollisionMatrix_T;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using MovingPrimitive_T = typename ClusterTree_T::MovingPrimitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
        
        using PrimitiveCollisionFinder_T  = CollisionFinder<AMB_DIM,Real,Int,SReal>;
        
        
        CLASS( const ClusterTree_T & S_, const ClusterTree_T & T_, const SReal t_init_ )
        :   S( S_ )
        ,   T( T_ )
        ,   t_init(t_init_)
        ,   thread_count( std::min(S_.ThreadCount(), T_.ThreadCount()) )
//        ,   is_symmetric( std::addressof(S_) == std::addressof(T_) )
        ,   collision_i( thread_count )
        ,   collision_j( thread_count )
        ,   collision_t( thread_count )
        {
            ptic(className()+"()");
            
            if constexpr ( is_symmetric)
            {
                assert( std::addressof(S_) == std::addressof(T_) );
            }
            
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                collision_i[thread] = std::vector<Int>();
                collision_j[thread] = std::vector<Int>();
                collision_t[thread] = std::vector<SReal>();
            }
            ptoc(className()+"()");
        }
        
//        // Copy constructor
//        CLASS( const CLASS & other )
//        :   S( other.S_ )
//        ,   T( other.T_ )
//        ,   t_init(other.t_init_)
//        ,   thread_count( other.thread_count )
//        ,   is_symmetric( other.is_symmetric )
//        ,   collision_i( thread_count )
//        ,   collision_j( thread_count )
//        ,   collision_t( thread_count )
//        {}
        
        
        virtual ~CLASS() = default;
        
    protected:
        
        // Not very elegant to use raw pointers here, but maybe acceptable due to constness.
//        const ClusterTree_T * S = nullptr; // "left"  BVH (output side of matrix-vector multiplication)
//        const ClusterTree_T * T = nullptr; // "right" BVH (input  side of matrix-vector multiplication)

        const ClusterTree_T & S; // "left"  BVH (output side of matrix-vector multiplication)
        const ClusterTree_T & T; // "right" BVH (input  side of matrix-vector multiplication)

        
        const SReal t_init = 1;
        
        const Int thread_count = 1;
        
//        const bool is_symmetric = false;
        
        mutable std::deque<Int> i_queue;
        mutable std::deque<Int> j_queue;
        
        mutable std::vector<std::vector<Int>>   collision_i;
        mutable std::vector<std::vector<Int>>   collision_j;
        mutable std::vector<std::vector<SReal>> collision_t;
        
        mutable bool P_collision_matrix_initialized = false;
        mutable CollisionMatrix_T P_collision_matrix;
        mutable bool C_collision_matrix_initialized = false;
        mutable CollisionMatrix_T C_collision_matrix;
        
        static constexpr Int max_depth = 128;
        static constexpr Int null = static_cast<Int>(0);
        
    protected:
        
//####################################################################################################
//      FindCollidingPrimitives
//####################################################################################################
        
        void FindCollidingPrimitives() const
        {
            ptic(className()+"::FindCollidingPrimitives");
            
            auto & C_matrix = ClusterCollisionMatrix();
            
            if( C_matrix.NonzeroCount() <= 0 )
            {
                P_collision_matrix = CollisionMatrix_T(0,0,0,1);
                
                DUMP(P_collision_matrix.Stats());
                
                ptoc(className()+"::FindCollidingPrimitives");
                return;
            }
            
            
            DUMP(is_symmetric);
            
            const Int expected = is_symmetric ?
                    (S.PrimitiveCount() + S.PrimitiveAdjacencyMatrix().NonzeroCount()) : 0;
            
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                collision_i[thread] = std::vector<Int>();
                collision_j[thread] = std::vector<Int>();
                collision_t[thread] = std::vector<SReal>();
                
                collision_i[thread].reserve(expected);
                collision_j[thread].reserve(expected);
                collision_t[thread].reserve(expected);
            }
            
            const auto & job_ptr = C_matrix.JobPtr();
            
            const Int * restrict const S_begin = S.ClusterBegin().data();
            const Int * restrict const S_end   = S.ClusterEnd().data();
            const Int * restrict const T_begin = T.ClusterBegin().data();
            const Int * restrict const T_end   = T.ClusterEnd().data();
            
            const Tensor2<SReal,Int> & p = S.PrimitiveSerializedData();
            const Tensor2<SReal,Int> & u = S.PrimitiveVelocitiesSerializedData();
            const Tensor2<SReal,Int> & q = T.PrimitiveSerializedData();
            const Tensor2<SReal,Int> & v = T.PrimitiveVelocitiesSerializedData();
                
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                const Int b_i_begin = job_ptr[thread  ];
                const Int b_i_end   = job_ptr[thread+1];
                
                const Int * restrict const outer = C_matrix.Outer().data();
                const Int * restrict const inner = C_matrix.Inner().data();
                
                PrimitiveCollisionFinder_T C (
                    S.MovingPrimitivePrototype(),
                    T.MovingPrimitivePrototype()
                );
                
                auto & collision_idx  = collision_i[thread];
                auto & collision_jdx  = collision_j[thread];
                auto & collision_time = collision_t[thread];
                
                if( is_symmetric )
                {
                    const SparseBinaryMatrixCSR<Int> & A = S.PrimitiveAdjacencyMatrix();
                    
                    for( Int b_i = b_i_begin; b_i < b_i_end; ++b_i )
                    {
                        const Int b_k_begin = outer[b_i  ];
                        const Int b_k_end   = outer[b_i+1];
                        
                        for( Int b_k = b_k_begin; b_k < b_k_end; ++b_k )
                        {
                            const Int b_j = inner[b_k];
                            
                            const Int i_begin = S_begin[b_i];
                            const Int i_end   = S_end  [b_i];
                            
                            const Int j_end   = T_end  [b_j];
                            
                            for( Int i = i_begin; i < i_end; ++i )
                            {
                                const Int j_begin = (b_i==b_j) ? i+1 : T_begin[b_j];
                                
                                for( Int j = j_begin; j < j_end; ++j )
                                {
                                    if( A.FindNonzeroPosition(i,j) < 0 )
                                    {
                                        const SReal t = C.FindMaximumSafeStepSize(
                                            p.data(i), u.data(i),
                                            q.data(j), v.data(j),
                                            t_init, false
                                        );
                                        
                                        if( t < t_init )
                                        {
                                            collision_idx.push_back(i);
                                            collision_jdx.push_back(j);
                                            collision_time.push_back(t);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    //TODO: Code for asymmetric matrices has yet to be debugged.
                    wprint(className()+"FindCollidingPrimitives: Code for asymmetric matrices has yet to be debugged.");
                    for( Int b_i = b_i_begin; b_i < b_i_end; ++b_i )
                    {
                        const Int b_k_begin = outer[b_i  ];
                        const Int b_k_end   = outer[b_i+1];
                        
                        for( Int b_k = b_k_begin; b_k < b_k_end; ++b_k )
                        {
                            const Int b_j = inner[b_k];
                            
                            const Int i_begin = S_begin[b_i];
                            const Int i_end   = S_end  [b_i];
                            
                            const Int j_begin = T_begin[b_j];
                            const Int j_end   = T_end  [b_j];
                            
                            for( Int i = i_begin; i < i_end; ++i )
                            {
                                for( Int j = j_begin; j < j_end; ++j )
                                {
                                    const SReal t = C.FindMaximumSafeStepSize(
                                        p.data(i), u.data(i),
                                        q.data(j), v.data(j),
                                        t_init, false
                                    );
                                    
                                    if( t < t_init )
                                    {
                                        collision_idx.push_back(i);
                                        collision_jdx.push_back(j);
                                        collision_time.push_back(t);
                                    }
                                }
                            }
                        }
                    }
                } // if( is_symmetric )
                
            } // #pragma omp parallel for num_threads( job_ptr.Size()-1 )

            
            P_collision_matrix = CollisionMatrix_T(
                collision_i, collision_j, collision_t,
                S.PrimitiveCount(), T.PrimitiveCount(),
                std::min(S.ThreadCount(), T.ThreadCount()),
                false, false
            );
            
            DUMP(P_collision_matrix.Stats());
            // Free memory that is no longer used.
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                collision_i[thread] = std::vector<Int>();
                collision_j[thread] = std::vector<Int>();
                collision_t[thread] = std::vector<SReal>();
            }
            
            ptoc(className()+"::FindCollidingPrimitives");
        }
        
//####################################################################################################
//      FindCollidingClusters
//####################################################################################################
        
        void FindCollidingClusters() const
        {
            ptic(className()+"::FindCollidingClusters");
            
            const Int expected = S.LeafClusterCount() + T.LeafClusterCount();
            
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                collision_i[thread] = std::vector<Int>();
                collision_j[thread] = std::vector<Int>();
                collision_t[thread] = std::vector<SReal>();
                
                collision_i[thread].reserve(expected);
                collision_j[thread].reserve(expected);
                collision_t[thread].reserve(expected);
            }
            
//            FindCollidingClusters_Sequential();
            
            if(
               (S.ClusterSerializedData().Dimension(0) != S.ClusterUpdatedSerializedData().Dimension(0))
               ||
               (T.ClusterSerializedData().Dimension(0) != T.ClusterUpdatedSerializedData().Dimension(0))
               )
            {
                wprint(className()+"::FindCollidingClusters: No update data loaded. Returning zero matrix.");
            }
            else
            {
                FindCollidingClusters_Parallel();
            }
            
            C_collision_matrix = CollisionMatrix_T(
                collision_i, collision_j,collision_t,
                S.ClusterCount(), T.ClusterCount(),
                std::min( S.ThreadCount(), T.ThreadCount() ),
                false, false
            );
            
            DUMP(C_collision_matrix.Stats());
            // Free memory that is no longer used.
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                collision_i[thread] = std::vector<Int>();
                collision_j[thread] = std::vector<Int>();
                collision_t[thread] = std::vector<SReal>();
            }
            
            ptoc(className()+"::FindCollidingClusters");
        }
        
//################################################################################################
//      FindCollidingClusters_Sequential
//################################################################################################

#include "CollisionTree/FindCollidingClusters_DFS.hpp"
        
        void FindCollidingClusters_Sequential() const
        {
            ptic(className()+"::FindCollidingClusters_Sequential");
            
            FindCollidingClusters_DFS(0,0);
            
            ptoc(className()+"::FindCollidingClusters_Sequential");
        }
        
//################################################################################################
//      FindCollidingClusters_Parallel
//################################################################################################
        
#include "CollisionTree/FindCollidingClusters_BFS.hpp"
        
        void FindCollidingClusters_Parallel() const
        {
            ptic(className()+"::FindCollidingClusters_Parallel");
            
            i_queue = std::deque<Int>();
            j_queue = std::deque<Int>();
            

            FindCollidingClusters_BFS( static_cast<Int>(4) * ThreadCount() * ThreadCount() );

            #pragma omp parallel for num_threads( thread_count ) schedule( dynamic )
            for( Int k = 0; k < static_cast<Int>(i_queue.size()); ++k )
            {
                FindCollidingClusters_DFS(i_queue[k], j_queue[k]);
            }

            i_queue = std::deque<Int>();
            j_queue = std::deque<Int>();
            
            ptoc(className()+"::FindCollidingClusters_Parallel");
        }
        

        
        
//####################################################################################################
//      FindMaximumSafeStepSize
//####################################################################################################
        
    public:
        
        SReal MaximumSafeStepSize() const override
        {
            ptic(className()+"::MaximumSafeStepSize");
            
//            SReal t = MaximumSafeStepSize_Sequential();
            
            SReal t = MaximumSafeStepSize_Parallel();
            
            ptoc(className()+"::MaximumSafeStepSize");
            
            return t;
        }
        
//####################################################################################################
//      MaximumSafeStepSize_Sequential
//####################################################################################################
               
#include "CollisionTree/MaximumSafeStepSize_DFS.hpp"
#include "CollisionTree/MaximumSafeStepSize_BFS.hpp"
        
    protected:
        
        SReal MaximumSafeStepSize_Sequential() const
        {
            ptic(className()+"::MaximumSafeStepSize_Sequential");
            
            SReal t;
            
            t = MaximumSafeStepSize_DFS(0,0,t_init);
            
            ptoc(className()+"::MaximumSafeStepSize_Sequential");
            
            return t;
        }        
        
        
//####################################################################################################
//      MaximumSafeStepSize_Parallel
//####################################################################################################
                       
    protected:
        
        SReal MaximumSafeStepSize_Parallel() const
        {
            ptic(className()+"::MaximumSafeStepSize_Parallel");
            
            SReal t_max = t_init;
            
            const Int max_leaves = ( static_cast<Int>(4) * ThreadCount() * ThreadCount() );
            
            i_queue = std::deque<Int>();
            j_queue = std::deque<Int>();
            
            logprint("Peeling the top of tree by BFS search.");
            t_max = MaximumSafeStepSize_BFS( max_leaves, t_max );

            DUMP(t_max);
//                DUMP(ThreadCount());
            
            #pragma omp parallel for num_threads( thread_count ) schedule(dynamic)
            for( Int k = 0; k < static_cast<Int>(i_queue.size()); ++k )
            {
                const SReal t = MaximumSafeStepSize_DFS( i_queue[k], j_queue[k], t_max );
                
                logprint("Chunk "+ToString(k)+" got t_max = "+ToString(t)+".");
                #pragma omp critical(t_max)
                {
                    t_max = std::min( t, t_max );
                }
            }

            i_queue = std::deque<Int>();
            j_queue = std::deque<Int>();
            
            ptoc(className()+"::MaximumSafeStepSize_Parallel");
            
            DUMP(t_max);
            
            return t_max;
        }
        
    public:
        
        virtual Int ThreadCount() const override
        {
            return thread_count;
        }
        
        Int AmbDim() const override
        {
            return AMB_DIM;
        }
        
        const CollisionMatrix_T & ClusterCollisionMatrix() const override
        {
            if( !C_collision_matrix_initialized )
            {
                FindCollidingClusters();
                
                C_collision_matrix_initialized = true;
            }
            
            return C_collision_matrix;
        }
        
        const CollisionMatrix_T & PrimitiveCollisionMatrix() const override
        {
            if( !P_collision_matrix_initialized )
            {
                FindCollidingPrimitives();
                
                P_collision_matrix_initialized = true;
            }
            
            return P_collision_matrix;
        }
        
        const ClusterTree_T & GetS() const override
        {
            return S;
        }
        
        const ClusterTree_T & GetT() const override
        {
            return T;
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS) + "<"+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+","+ToString(is_symmetric)+">";
        }

    }; // CLASS
    
} // namespace Repulsor

#undef BASE
#undef CLASS
