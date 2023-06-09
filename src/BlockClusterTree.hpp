#pragma once

#include "BlockClusterTree/BlockClusterTreeBase.hpp"
#include "BlockClusterTree/BlockSplit_Kernel.hpp"

namespace Repulsor
{
    
    
    template<int AMB_DIM_, typename Real_, typename Int_, typename SReal_, typename ExtReal_, bool is_symmetric_>
    class BlockClusterTree : public BlockClusterTreeBase<Real_,Int_,SReal_,ExtReal_,is_symmetric_>
    {
    public:
        
        
        
        using Real    = Real_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using Int     = Int_;
        
        using BlockClusterTreeBase_T = BlockClusterTreeBase<Real,Int,SReal,ExtReal,is_symmetric_>;
        
        using LInt    = typename BlockClusterTreeBase_T::LInt;
        

        static constexpr Int  AMB_DIM      = AMB_DIM_;
        static constexpr bool is_symmetric = is_symmetric_;
        
        using Setting_T = typename BlockClusterTreeBase_T::Setting_T;

        using Inter_Pattern_T     = typename BlockClusterTreeBase_T::VeryNear_Pattern_T;
        using VeryNear_Pattern_T  = typename BlockClusterTreeBase_T::VeryNear_Pattern_T;
        using Near_Pattern_T      = typename BlockClusterTreeBase_T::Near_Pattern_T;
        using Far_Pattern_T       = typename BlockClusterTreeBase_T::Far_Pattern_T;
        
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        using BlockSplitter_T   = BlockSplit_Kernel<ClusterTree_T, LInt>;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
        
        using GJK_T             = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
    
        using BlockClusterTreeBase_T::AmbDim;
        using BlockClusterTreeBase_T::ThreadCount;
        using BlockClusterTreeBase_T::FarFieldSeparationParameter;
        using BlockClusterTreeBase_T::NearFieldSeparationParameter;
        using BlockClusterTreeBase_T::Settings;
        
    public:

        // In order to prevent GetS() and GetT() shooting a segfault, we have to initialize S and T here. This is the only case in which BlockClusterTree owns these raw pointers.
        
        virtual ~BlockClusterTree() = default;
        
        BlockClusterTree(
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
                if( std::addressof(S_) != std::addressof(T_) )
                {
                    eprint(className()+": is_symmetric == true, bu S != T.");
                }
            }

            if( std::min( S.PrimitiveCount(), T.PrimitiveCount()) > 0 )
            {
                ComputeBlocks();
            }

            ptoc(className());
        } // Constructor
        
    protected:

        static constexpr Int null = static_cast<Int>(0);
        
        // Not very elegant to use raw pointers here, but maybe acceptable due to constness.
        const ClusterTree_T & S; // "left"  BVH (output side of matrix-vector multiplication)
        const ClusterTree_T & T; // "right" BVH (input  side of matrix-vector multiplication)
        
        const Setting_T settings;
        
        mutable Inter_Pattern_T    inter;
        mutable VeryNear_Pattern_T verynear;
        mutable Near_Pattern_T     near;
        mutable Far_Pattern_T      far;
        
        const SReal near_theta2 = static_cast<SReal>(10);
        const SReal  far_theta2 = static_cast<SReal>(0.25);
        
        const Int thread_count = static_cast<Int>(1);
//        const bool is_symmetric;
        
        mutable bool blocks_initialized = false;
        
        static constexpr Int max_depth = 128;
        
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
        
        LInt PrimitiveIntersectionCount() const override
        {
            return inter.NonzeroCount();
        }
        
        LInt VeryNearFieldInteractionCount() const override
        {
            return verynear.NonzeroCount();
        }
        
        LInt NearFieldInteractionCount() const override
        {
            return near.NonzeroCount();
        }
        
        LInt FarFieldInteractionCount() const override
        {
            return far.NonzeroCount();
        }

        
        virtual const Inter_Pattern_T & PrimitiveIntersectionMatrix() const override
        {
            return inter;
        }
        
        virtual const VeryNear_Pattern_T & VeryNear() const override
        {
            return verynear;
        }
        
        virtual const Near_Pattern_T & Near() const override
        {
            return near;
        }
        
        virtual const Far_Pattern_T & Far() const override
        {
            return far;
        }
        
        virtual const ClusterTree_T & GetS() const override
        {
            return S;
        }
        
        virtual const ClusterTree_T & GetT() const override
        {
            return T;
        }

        virtual const Setting_T & Settings() const override
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
        
    protected:
        
        void ComputeBlocks()
        {
            if( blocks_initialized )
            {
                return;
            }

            ptic(className()+"::ComputeBlocks");

            std::vector<BlockSplitter_T> kernels;

            ptic(className()+"::ComputeBlocks: prepare kernels");
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                kernels.emplace_back( S, T, far_theta2, near_theta2 );
            }
            ptoc(className()+"::ComputeBlocks: prepare kernels");
            
            ptic(className()+"::ComputeBlocks: traversal");
            if( (S.SplitThreshold()==1) && (T.SplitThreshold()==1) )
            {
                Traversor<BlockSplitter_T,is_symmetric,true> traversor ( S, T, kernels );

                traversor.Traverse();
            }
            else
            {
                Traversor<BlockSplitter_T,is_symmetric,false> traversor ( S, T, kernels );

                traversor.Traverse();
            }
            ptoc(className()+"::ComputeBlocks: traversal");
            
            
            ptic(className()+"::ComputeBlocks: reduce kernels");
            
            std::vector<PairAggregator<Int,Int,LInt>> inter_idx    (thread_count);
            std::vector<PairAggregator<Int,Int,LInt>> verynear_idx (thread_count);
            std::vector<PairAggregator<Int,Int,LInt>> near_idx     (thread_count);
            std::vector<PairAggregator<Int,Int,LInt>> far_idx      (thread_count);

            
            ParallelDo(
                [&kernels]( const Int thread )
                {
                    kernels[thread].inter_idx.Finalize();
                    kernels[thread].verynear_idx.Finalize();
                    kernels[thread].near_idx.Finalize();
                    kernels[thread].far_idx.Finalize();
                },
                thread_count
            );
            
            size_t intersection_count = 0;
            
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                   inter_idx[thread] = std::move( kernels[thread].inter_idx    );
                verynear_idx[thread] = std::move( kernels[thread].verynear_idx );
                    near_idx[thread] = std::move( kernels[thread].near_idx     );
                     far_idx[thread] = std::move( kernels[thread].far_idx      );

                logprint("Thread "+ToString(thread)+" aggregated this many elements: "
                    + "\n\t intersections:   \t" + ToString(inter_idx[thread].Size())
                    + "\n\t very near field: \t" + ToString(verynear_idx[thread].Size())
                    + "\n\t near field:      \t" + ToString(near_idx[thread].Size())
                    + "\n\t far field:       \t" + ToString(far_idx[thread].Size())
                    + "\n"
                );
                
                intersection_count += static_cast<size_t>(inter_idx[thread].Size());
            }
            
            if( intersection_count > 0 )
            {
                wprint(className()+"::ComputeBlocks detected "+ToString(intersection_count)+" intersections.");
            }
            
            ptoc(className()+"::ComputeBlocks: reduce kernels");
            
            ptoc(className()+"::ComputeBlocks");


            ptic(className()+"  Primitive intersection data");
            
//            inter = Inter_T( inter_idx, inter_jdx, S.PrimitiveCount(), T.PrimitiveCount(),
//                thread_count, false, is_symmetric );
            
            inter = Inter_Pattern_T( inter_idx, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );

            pdump(inter.Stats());

            ptoc(className()+"  Primitive intersection data");
            

            ptic(className()+"  Very near field interaction data");

            verynear = VeryNear_Pattern_T( verynear_idx, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            pdump(verynear.Stats());

            ptoc(className()+"  Very near field interaction data");
            

            ptic(className()+"  Near field interaction data");

            near = Near_Pattern_T( near_idx, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );

            pdump(near.Stats());

            ptoc(className()+"  Near field interaction data");

            
            ptic(className()+"  Far field interaction data");

            far = Far_Pattern_T( far_idx, S.ClusterCount(), T.ClusterCount(),
                    thread_count, false, is_symmetric );

            pdump(far.Stats());

            ptoc(className()+"  Far field interaction data");

            
            blocks_initialized = true;

        }; // ComputeBlocks
        
    private:
        
        static std::string className()
        {
            return  "BlockClusterTree<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+ToString(is_symmetric)+">";
        }
      
    public:
        
        std::string ClassName() const override
        {
            return className();
        }
    
    }; // class BlockClusterTree
    
} //namespace Repulsor
