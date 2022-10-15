#pragma once

#include "BlockClusterTree/BlockClusterTreeBase.hpp"
#include "BlockClusterTree/BlockSplitter.hpp"

#define CLASS BlockClusterTree
#define BASE  BlockClusterTreeBase<Real,Int,SReal,ExtReal,is_symmetric>

namespace Repulsor
{
    template<int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal, bool is_symmetric>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);
        
    public:
        
        using BlockClusterTreeBase_T = BASE;
        
        using Setting_T = typename BASE::Setting_T;

        using Inter_T           = typename BASE::Inter_T;
        using VeryNear_T        = typename BASE::VeryNear_T;
        using Near_T            = typename BASE::Near_T;
        using Far_T             = typename BASE::Far_T;
        
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
    
        using GJK_T             = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
        using BASE::AmbDim;
        using BASE::ThreadCount;
        using BASE::FarFieldSeparationParameter;
        using BASE::NearFieldSeparationParameter;
        using BASE::Settings;
        
    public:

        // In order to prevent GetS() and GetT() shooting a segfault, we have to initialize S and T here. This is the only case in which CLASS owns these raw pointers.
        
        virtual ~CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_,
            Setting_T settings_ = Setting_T()
        )
        :   S(S_)
        ,   T(T_)
        ,   settings(settings_)
        ,   near_theta2( pow(settings_.near_field_separation_parameter,2) )
        ,    far_theta2( pow(settings_. far_field_separation_parameter,2) )
        ,   thread_count( std::min(S_.ThreadCount(), T_.ThreadCount()) )
//        ,   is_symmetric( std::addressof(S_) == std::addressof(T_) )
        {
            ptic(className());
            
            if constexpr ( is_symmetric )
            {
                assert( std::addressof(S_) == std::addressof(T_) );
            }

            if( std::min( S.PrimitiveCount(), T.PrimitiveCount()) > 0 )
            {
                ComputeBlocks();
            }

            ptoc(className());
        } // Constructor
        
    protected:

        static constexpr Int zero = static_cast<Int>(0);
        
        // Not very elegant to use raw pointers here, but maybe acceptable due to constness.
        const ClusterTree_T & S; // "left"  BVH (output side of matrix-vector multiplication)
        const ClusterTree_T & T; // "right" BVH (input  side of matrix-vector multiplication)
        
        const Setting_T settings;
        
        mutable Inter_T    inter;
        mutable VeryNear_T verynear;
        mutable Near_T     near;
        mutable Far_T      far;
        
        const SReal near_theta2 = static_cast<SReal>(10);
        const SReal  far_theta2 = static_cast<SReal>(0.25);
        
        const Int thread_count = static_cast<Int>(1);
//        const bool is_symmetric;
        
        mutable bool blocks_initialized = false;
        
        mutable std::deque<Int> i_queue;
        mutable std::deque<Int> j_queue;
        
        // Containers for parallel aggregation of {i,j}-index pairs
        mutable std::vector<std::vector<Int>> inter_i;
        mutable std::vector<std::vector<Int>> inter_j;
        mutable std::vector<std::vector<Int>> verynear_i;
        mutable std::vector<std::vector<Int>> verynear_j;
        mutable std::vector<std::vector<Int>> near_i;
        mutable std::vector<std::vector<Int>> near_j;
        mutable std::vector<std::vector<Int>> far_i;
        mutable std::vector<std::vector<Int>> far_j;
        
        static constexpr Int max_depth = 128;
        static constexpr Int null = static_cast<Int>(0);
        
    public:

        Int AmbDim() const override
        {
            return AMB_DIM;
        }
        
        Int ThreadCount() const override
        {
            return thread_count;
        }
        
        Real FarFieldSeparationParameter() const override
        {
            return sqrt(far_theta2);
        }
                        
        Real NearFieldSeparationParameter() const override
        {
            return sqrt(near_theta2);
        }
        
        Int PrimitiveIntersectionCount() const override
        {
            return inter.NonzeroCount();
        }
        
        Int VeryNearFieldInteractionCount() const override
        {
            return verynear.NonzeroCount();
        }
        
        Int NearFieldInteractionCount() const override
        {
            return near.NonzeroCount();
        }
        
        Int FarFieldInteractionCount() const override
        {
            return far.NonzeroCount();
        }

        
        virtual const Inter_T & PrimitiveIntersectionMatrix() const override
        {
            return inter;
        }
        
        const VeryNear_T & VeryNear() const override
        {
            return verynear;
        }
        
        const Near_T & Near() const override
        {
            return near;
        }
        
        const Far_T & Far() const override
        {
            return far;
        }

        
        const ClusterTree_T & GetS() const override
        {
            return S;
        }
        
        const ClusterTree_T & GetT() const override
        {
            return T;
        }

        const Setting_T & Settings() const override
        {
            return settings;
        }
        
        virtual std::string Stats() const override
        {
            std::stringstream s;
            
            s
            << "\n==== "+className()+" Stats ====" << "\n\n"
            << " AmbDim()                        = " <<  AmbDim() << "\n"
            << " FarFieldSeparationParameter()   = " <<  FarFieldSeparationParameter() << "\n"
            << " NearFieldSeparationParameter()  = " <<  NearFieldSeparationParameter() << "\n"
            << " ThreadCount()                   = " <<  ThreadCount() << "\n"
            << "\n"
            << "Tree S \n"
            << S.Stats() << "\n"
            << "Tree T \n"
            << T.Stats() << "\n"
            << " PrimitiveIntersectionCount()    = " <<  PrimitiveIntersectionCount() << "\n"
            << " VeryNearFieldInteractionCount() = " <<  VeryNearFieldInteractionCount() << "\n"
            << " NearFieldInteractionCount()     = " <<  NearFieldInteractionCount() << "\n"
            << " FarFieldInteractionCount()      = " <<  FarFieldInteractionCount() << "\n"

            << "\n---- bool data ----" << "\n"

            << " IsSymmetric()               = " <<  is_symmetric << "\n"

            << "\n==== "+className()+" Stats ====\n" << std::endl;

            return s.str();
        }
        

        
    
//#################################################################################################
//      Initialization
//#################################################################################################

#include "BlockClusterTree/Split_Parallel.hpp"

#include "BlockClusterTree/Split_Sequential_DFS.hpp"
        
#include "BlockClusterTree/Split_Sequential_DFS_1.hpp"
        
    protected:
        
        void ComputeBlocks()
        {
            if( blocks_initialized )
            {
                return;
            }
            
            ptic(className()+"::ComputeBlocks");

            inter_i    = std::vector<std::vector<Int>>(thread_count);
            inter_j    = std::vector<std::vector<Int>>(thread_count);
            verynear_i = std::vector<std::vector<Int>>(thread_count);
            verynear_j = std::vector<std::vector<Int>>(thread_count);
            near_i     = std::vector<std::vector<Int>>(thread_count);
            near_j     = std::vector<std::vector<Int>>(thread_count);
            far_i      = std::vector<std::vector<Int>>(thread_count);
            far_j      = std::vector<std::vector<Int>>(thread_count);
            
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                // TODO: Find a better initial guess for the number of admissable and inadmissable blocks.
                const Int expected = static_cast<Int>(10) * ( S.PrimitiveCount() + T.PrimitiveCount() );

                   inter_i[thread] = std::vector<Int>();
                   inter_j[thread] = std::vector<Int>();
                verynear_i[thread] = std::vector<Int>();
                verynear_j[thread] = std::vector<Int>();
                    near_i[thread] = std::vector<Int>();
                    near_j[thread] = std::vector<Int>();
                     far_i[thread] = std::vector<Int>();
                     far_j[thread] = std::vector<Int>();
                
                near_i[thread].reserve(expected);
                near_j[thread].reserve(expected);
                 far_i[thread].reserve(expected);
                 far_j[thread].reserve(expected);
            }

            switch( settings.block_split_method )
            {
                case BlockSplitMethod::Parallel:
                {
                    Split_Parallel( static_cast<Int>(4) * thread_count * thread_count );
                    break;
                }
                case BlockSplitMethod::Sequential:
                {
                    if( S.SplitThreshold() == 1 && T.SplitThreshold() == 1 )
                    {
                        Split_Sequential_DFS_1( zero, zero );
                    }
                    else
                    {
                        Split_Sequential_DFS( zero, zero );
                    }
                    break;
                }
                default:
                {
                    Split_Parallel( static_cast<Int>(4) * thread_count * thread_count );
                    break;
                }
            }
        
            ptic(className()+"  Far field interaction data");
            
            far = Far_T( far_i, far_j, S.ClusterCount(), T.ClusterCount(),
                    thread_count, false, is_symmetric );
            
            pdump(far.Stats());
            
            ptoc(className()+"  Far field interaction data");

            ptic(className()+"  Near field interaction data");
            
            near = Near_T( near_i, near_j, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            pdump(near.Stats());
            
            ptoc(className()+"  Near field interaction data");
            
            ptic(className()+"  Very near field interaction data");

            verynear = VeryNear_T( verynear_i, verynear_j, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            pdump(verynear.Stats());
            
            ptoc(className()+"  Very near field interaction data");
            
            ptic(className()+"  Primitive intersection data");

            inter = Inter_T( inter_i, inter_j, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            pdump(inter.Stats());
            
            ptoc(className()+"  Primitive intersection data");
            
            blocks_initialized = true;

            // Free memory that is no longer used.
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                   inter_i[thread] = std::vector<Int>();
                   inter_j[thread] = std::vector<Int>();
                verynear_i[thread] = std::vector<Int>();
                verynear_j[thread] = std::vector<Int>();
                    near_i[thread] = std::vector<Int>();
                    near_j[thread] = std::vector<Int>();
                     far_i[thread] = std::vector<Int>();
                     far_j[thread] = std::vector<Int>();
            }
            
            ptoc(className()+"::ComputeBlocks");

        }; // ComputeBlocks

        

        
        std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS) + "<"+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+","+ToString(is_symmetric)+">";
        }
      
        
    };
    
} //namespace Repulsor

#undef BASE
#undef CLASS
