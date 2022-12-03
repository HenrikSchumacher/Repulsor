#pragma once

#include "CollisionTree/CollisionTreeBase.hpp"
#include "CollisionTree/Compute_AABB_CollisionTimeInterval.hpp"
#include "CollisionTree/MaximumSafeStepSize_Kernel.hpp"
#include "CollisionTree/Collision_Kernel.hpp"

namespace Repulsor
{
    template<int AMB_DIM, typename Real_, typename Int_, typename SReal_, typename ExtReal_, bool is_symmetric>
    class CollisionTree : public CollisionTreeBase<Real_,Int_,SReal_,ExtReal_,is_symmetric>
    {
        ASSERT_FLOAT(Real_   );
        ASSERT_INT  (Int_    );
        ASSERT_FLOAT(SReal_  );
        ASSERT_FLOAT(ExtReal_);
        
    public:
        
        using Real              = Real_;
        using Int               = Int_;
        using SReal             = SReal_;
        using ExtReal           = ExtReal_;
        
        using Base_T            = CollisionTreeBase<Real,Int,SReal,ExtReal,is_symmetric>;
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        using CollisionMatrix_T = typename Base_T::CollisionMatrix_T;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using MovingPrimitive_T = typename ClusterTree_T::MovingPrimitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
        
        using PrimitiveCollisionFinder_T = CollisionFinder<AMB_DIM,Real,Int,SReal>;
        
        
        CollisionTree( const ClusterTree_T & S_, const ClusterTree_T & T_ )
        :   Base_T( S_,T_ )
        ,   S( S_ )
        ,   T( T_ )
        ,   thread_count( std::min(S_.ThreadCount(), T_.ThreadCount()) )
        {
            ptic(className()+"()");
            if constexpr ( is_symmetric )
            {
                if( std::addressof(S_) != std::addressof(T_) )
                {
                    eprint(className()+": is_symmetric == true, bu S != T.");
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
        
        const ClusterTree_T & S; // "left"  BVH (output side of matrix-vector multiplication)
        const ClusterTree_T & T; // "right" BVH (input  side of matrix-vector multiplication)
        
        const Int thread_count = 1;
        
        
        mutable bool P_collision_matrix_initialized = false;
        mutable CollisionMatrix_T P_collision_matrix;
        

//####################################################################################################
//      PrimitiveCollisionMatrix
//####################################################################################################

    public:
        
        const CollisionMatrix_T & PrimitiveCollisionMatrix() const override
        {
            if( !P_collision_matrix_initialized )
            {
                ptic(ClassName()+"PrimitiveCollisionMatrix");
                
                using Kernel_T = Collision_Kernel<ClusterTree_T>;
                using LInt     = typename Kernel_T::LInt;
                
                ptic(ClassName()+"::PrimitiveCollisionMatrix: Prepare kernels");
                std::vector<Kernel_T> kernels;
                
                SReal t = std::min( S.UpdateTime(), T.UpdateTime() );
                
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    kernels.emplace_back( S, T, t );
                }
                ptoc(ClassName()+"::PrimitiveCollisionMatrix: Prepare kernels");
                
                if( S.SplitThreshold()==1 && T.SplitThreshold()==1 )
                {
                    Traversor<Kernel_T, is_symmetric, true > traversor (S, T, kernels);
                    
                    traversor.Traverse();
                }
                else
                {
                    Traversor<Kernel_T, is_symmetric, false> traversor (S, T, kernels);
                    
                    traversor.Traverse();
                }
                
                std::vector<TripleAggregator<Int,Int,SReal,LInt>> triples (thread_count);
                
                ptic(ClassName()+"::PrimitiveCollisionMatrix: Reduce kernels");

                #pragma omp parallel for num_threads(thread_count)
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    kernels[thread].triples.Finalize();
                }
                
                #pragma omp parallel for num_threads(thread_count)
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    triples[thread] = std::move( kernels[thread].triples );
                }
                ptoc(ClassName()+"::PrimitiveCollisionMatrix: Reduce kernels");
                
                P_collision_matrix = CollisionMatrix_T(
                    triples,
                    S.PrimitiveCount(), T.PrimitiveCount(),
                    thread_count, false, is_symmetric
                );
                
                P_collision_matrix_initialized = true;
                
                ptoc(ClassName()+"PrimitiveCollisionMatrix");
            }
            
            return P_collision_matrix;
        }

//####################################################################################################
//      MaximumSafeStepSize
//####################################################################################################
        
    public:
        
        ExtReal MaximumSafeStepSize( const SReal t_) const override
        {
            ptic(ClassName()+"::MaximumSafeStepSize");
            
            using Kernel_T = MaximumSafeStepSize_Kernel<ClusterTree_T>;
            
            (void)S.PrimitiveAdjacencyMatrix();
            (void)T.PrimitiveAdjacencyMatrix();
            
            std::vector<Kernel_T> kernels;
            
            SReal t = std::min( t_, std::min( S.UpdateTime(), T.UpdateTime() ) );
            
            ptic(ClassName()+"::MaximumSafeStepSize: Prepare kernels");
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    kernels.emplace_back( S, T, t );
                }
            ptoc(ClassName()+"::MaximumSafeStepSize: Prepare kernels");
            
            if( S.SplitThreshold()==1 && T.SplitThreshold()==1 )
            {
                Traversor<Kernel_T, is_symmetric, true>  traversor (S, T, kernels);
                
                traversor.Traverse();
            }
            else
            {
                Traversor<Kernel_T, is_symmetric, false> traversor (S, T, kernels);
                
                traversor.Traverse();
            }
            
            ptic(ClassName()+"::MaximumSafeStepSize: Reduce kernels");
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    kernels.emplace_back( S, T, t );
                }
            ptoc(ClassName()+"::MaximumSafeStepSize: Reduce kernels");
            
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                SReal s = kernels[thread].MaxTime();
                t = std::min( t, s );
            }
            
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
            return  "CollisionTree<"+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+","+ToString(is_symmetric)+">";
        }

    }; // class CollisionTree
    
} // namespace Repulsor
