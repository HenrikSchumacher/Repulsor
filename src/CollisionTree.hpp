#pragma once

#include "CollisionTree/CollisionTreeBase.hpp"
#include "CollisionTree/Compute_AABB_CollisionTimeInterval.hpp"
#include "CollisionTree/MaximumSafeStepSize_Kernel.hpp"
#include "CollisionTree/Collision_Kernel.hpp"

namespace Repulsor
{
    template<int AMB_DIM, typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_, bool symmetricQ>
    class CollisionTree : public CollisionTreeBase<Real_,Int_,LInt_,SReal_,ExtReal_,symmetricQ>
    {
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
        static_assert(IntQ< Int_>,"");
        static_assert(IntQ<LInt_>,"");
        
    private:
        
        using Base_T = CollisionTreeBase<Real_,Int_,LInt_,SReal_,ExtReal_,symmetricQ>;
        
    public:
        
        using Real              = Real_;
        using Int               = Int_;
        using SReal             = SReal_;
        using ExtReal           = ExtReal_;
        using LInt              = LInt_;
        
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,LInt,SReal,ExtReal>;
        using CollisionMatrix_T = typename Base_T::CollisionMatrix_T;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using MovingPrimitive_T = typename ClusterTree_T::MovingPrimitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
        
        using PrimitiveCollisionFinder_T = CollisionFinder<AMB_DIM,Real,Int,SReal>;
        
        
        CollisionTree( cref<ClusterTree_T> S_, cref<ClusterTree_T> T_ )
        :   Base_T( S_,T_ )
        ,   S( S_ )
        ,   T( T_ )
        ,   thread_count( Min(S_.ThreadCount(), T_.ThreadCount()) )
        {
            ptic(className()+"()");
            if constexpr ( symmetricQ )
            {
                if( std::addressof(S_) != std::addressof(T_) )
                {
                    eprint(className()+": symmetricQ == true, bu S != T.");
                }
            }
            ptoc(className()+"()");
        }
        
//        // Copy constructor
//        CollisionTree( const CollisionTree & other )
//        :   S( other.S_ )
//        ,   T( other.T_ )
//        ,   thread_count( other.thread_count )
//        {}
        
        
        virtual ~CollisionTree() = default;
        
    protected:
        
        static constexpr Int max_depth = 128;
        static constexpr Int null = static_cast<Int>(0);
        
        cref<ClusterTree_T> S; // "left"  BVH (output side of matrix-vector multiplication)
        cref<ClusterTree_T> T; // "right" BVH (input  side of matrix-vector multiplication)
        
        const Int thread_count = 1;
        
        
        mutable bool P_collision_matrix_initialized = false;
        mutable CollisionMatrix_T P_collision_matrix;
        

//#################################################################################
//      PrimitiveCollisionMatrix
//#################################################################################

    public:
        
        cref<CollisionMatrix_T> PrimitiveCollisionMatrix() const override
        {
            if( !P_collision_matrix_initialized )
            {
                ptic(ClassName()+"PrimitiveCollisionMatrix");
                
                using Kernel_T = Collision_Kernel<ClusterTree_T>;
                
                ptic(ClassName()+"::PrimitiveCollisionMatrix: Prepare kernels");
                std::vector<Kernel_T> kernels;
                kernels.reserve(thread_count);
                
                SReal t = Min( S.UpdateTime(), T.UpdateTime() );
                
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    kernels.emplace_back( S, T, thread, t, static_cast<SReal>(0.0625) );
                }
                ptoc(ClassName()+"::PrimitiveCollisionMatrix: Prepare kernels");
                
                if( S.SplitThreshold()==1 && T.SplitThreshold()==1 )
                {
                    ClusterTreePairTraversor<Kernel_T, symmetricQ, true > traversor (S, T, kernels);
                    
                    traversor.Traverse();
                }
                else
                {
                    ClusterTreePairTraversor<Kernel_T, symmetricQ, false> traversor (S, T, kernels);
                    
                    traversor.Traverse();
                }
                
                std::vector<TripleAggregator<Int,Int,SReal,LInt>> triples (thread_count);
                
                ptic(ClassName()+"::PrimitiveCollisionMatrix: Reduce kernels");

                
                ParallelDo(
                    [&kernels]( const Int thread )
                    {
                        kernels[thread].triples.Finalize();
                    },
                    thread_count
                );
                
                ParallelDo(
                    [&kernels,&triples]( const Int thread )
                    {
                        triples[thread] = std::move( kernels[thread].triples );
                    },
                    thread_count
                );

                ptoc(ClassName()+"::PrimitiveCollisionMatrix: Reduce kernels");
                
                P_collision_matrix = CollisionMatrix_T(
                    triples,
                    S.PrimitiveCount(), T.PrimitiveCount(),
                    thread_count, false, symmetricQ
                );
                
                P_collision_matrix_initialized = true;
                
                ptoc(ClassName()+"PrimitiveCollisionMatrix");
            }
            
            return P_collision_matrix;
        }

//#################################################################################
//      MaximumSafeStepSize
//#################################################################################
        
    public:
        
        ExtReal MaximumSafeStepSize(
            const SReal t_,
            const SReal TOL
        ) const override
        {
            ptic(ClassName()+"::MaximumSafeStepSize");
            
            using Kernel_T = MaximumSafeStepSize_Kernel<ClusterTree_T>;
            
            (void)S.PrimitiveAdjacencyMatrix();
            (void)T.PrimitiveAdjacencyMatrix();
            
            std::vector<Kernel_T> kernels;
            kernels.reserve( thread_count );
            
            
            SReal t = Min( t_, Min( S.UpdateTime(), T.UpdateTime() ) );
            
            ptic(ClassName()+"::MaximumSafeStepSize: Prepare kernels");
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    kernels.emplace_back( S, T, thread,t, TOL );
                }
            ptoc(ClassName()+"::MaximumSafeStepSize: Prepare kernels");
            
            if( S.SplitThreshold()==1 && T.SplitThreshold()==1 )
            {
                ClusterTreePairTraversor<Kernel_T, symmetricQ, true>  traversor (S, T, kernels);
                
                traversor.Traverse();
            }
            else
            {
                ClusterTreePairTraversor<Kernel_T, symmetricQ, false> traversor (S, T, kernels);
                
                traversor.Traverse();
            }
            
            ptic(ClassName()+"::MaximumSafeStepSize: Reduce kernels");
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    SReal s = kernels[thread].MaxTime();
                    t = Min( t, s );
                }
            ptoc(ClassName()+"::MaximumSafeStepSize: Reduce kernels");
            
            
            ptoc(ClassName()+"::MaximumSafeStepSize");
            
            return t;
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
        
        cref<ClusterTree_T> GetS() const override
        {
            return S;
        }
        
        cref<ClusterTree_T> GetT() const override
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
            return  "CollisionTree<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+ToString(symmetricQ)+">";
        }

    }; // class CollisionTree
    
} // namespace Repulsor
