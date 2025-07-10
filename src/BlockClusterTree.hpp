#pragma once

#include "BlockClusterTree/BlockClusterTreeBase.hpp"
#include "BlockClusterTree/BlockSplit_Kernel.hpp"

namespace Repulsor
{
    
    
    template<int AMB_DIM_, typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_, bool symmetricQ_>
    class BlockClusterTree : public BlockClusterTreeBase<Real_,Int_,LInt_,SReal_,ExtReal_,symmetricQ_>
    {
    public:
        
        using Real    = Real_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using Int     = Int_;
        using LInt    = LInt_;
        
        using BlockClusterTreeBase_T = BlockClusterTreeBase<Real,Int,LInt,SReal,ExtReal,symmetricQ_>;
        

        static constexpr Int  AMB_DIM      = AMB_DIM_;
        static constexpr bool symmetricQ = symmetricQ_;
        
        using Setting_T = typename BlockClusterTreeBase_T::Setting_T;

        using Inter_Pattern_T     = typename BlockClusterTreeBase_T::VeryNear_Pattern_T;
        using VeryNear_Pattern_T  = typename BlockClusterTreeBase_T::VeryNear_Pattern_T;
        using Near_Pattern_T      = typename BlockClusterTreeBase_T::Near_Pattern_T;
        using Far_Pattern_T       = typename BlockClusterTreeBase_T::Far_Pattern_T;
        
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,LInt,SReal,ExtReal>;
        using BlockSplitter_T   = BlockSplit_Kernel<ClusterTree_T>;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
        
        using GJK_T             = GJK<AMB_DIM,GJK_Real,Int>;
    
        using BlockClusterTreeBase_T::AmbDim;
        using BlockClusterTreeBase_T::ThreadCount;
        using BlockClusterTreeBase_T::FarFieldSeparationParameter;
        using BlockClusterTreeBase_T::NearFieldSeparationParameter;
        using BlockClusterTreeBase_T::Settings;
        
    public:
        
        BlockClusterTree(
            cref<ClusterTree_T> S_,
            cref<ClusterTree_T> T_,
            Setting_T settings_ = Setting_T()
        )
        :   S(S_)
        ,   T(T_)
        ,   settings(settings_)
        ,   near_theta2( pow(settings_.near_field_separation_parameter,2) )
        ,    far_theta2( pow(settings_. far_field_separation_parameter,2) )
        ,   thread_count( Min(S_.ThreadCount(), T_.ThreadCount()) )
//        ,   symmetricQ( std::addressof(S_) == std::addressof(T_) )
        {
            TOOLS_PTIC(className());
            
            if constexpr ( symmetricQ )
            {
                if( std::addressof(S_) != std::addressof(T_) )
                {
                    eprint(className()+": symmetricQ == true, bu S != T.");
                }
            }

            if( Min( S.PrimitiveCount(), T.PrimitiveCount()) > 0 )
            {
                ComputeBlocks();
            }

            TOOLS_PTOC(className());
        } // Constructor
        
        
        // No default constructor.
        BlockClusterTree() = delete;
        // Destructor
        virtual ~BlockClusterTree() override = default;
        // Copy constructor
        BlockClusterTree( const BlockClusterTree & other ) = default;
        // Copy assignment operator
        BlockClusterTree & operator=( const BlockClusterTree & other ) = default;
        // Move constructor
        BlockClusterTree( BlockClusterTree && other ) = default;
        // Move assignment operator
        BlockClusterTree & operator=( BlockClusterTree && other ) = default;
        
    protected:

        static constexpr Int null = static_cast<Int>(0);
        
        // Not very elegant to use raw pointers here, but maybe acceptable due to constness.
        cref<ClusterTree_T> S; // "left"  BVH (output side of matrix-vector multiplication)
        cref<ClusterTree_T> T; // "right" BVH (input  side of matrix-vector multiplication)
        
        const Setting_T settings;
        
        mutable Inter_Pattern_T    inter;
        mutable VeryNear_Pattern_T verynear;
        mutable Near_Pattern_T     near;
        mutable Far_Pattern_T      far;
        
        const SReal near_theta2 = static_cast<SReal>(10);
        const SReal  far_theta2 = static_cast<SReal>(0.25);
        
        const Int thread_count = static_cast<Int>(1);
//        const bool symmetricQ;
        
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
            return Sqrt(far_theta2);
        }
                        
        Real NearFieldSeparationParameter() const override
        {
            return Sqrt(near_theta2);
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

        
        virtual cref<Inter_Pattern_T> PrimitiveIntersectionMatrix() const override
        {
            return inter;
        }
        
        virtual cref<VeryNear_Pattern_T> VeryNear() const override
        {
            return verynear;
        }
        
        virtual cref<Near_Pattern_T> Near() const override
        {
            return near;
        }
        
        virtual cref<Far_Pattern_T> Far() const override
        {
            return far;
        }
        
        virtual cref<ClusterTree_T> GetS() const override
        {
            return S;
        }
        
        virtual cref<ClusterTree_T> GetT() const override
        {
            return T;
        }

        virtual cref<Setting_T> Settings() const override
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

            << " SymmetricQ()               = " <<  symmetricQ << "\n"

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

            TOOLS_PTIC(className()+"::ComputeBlocks");

            std::vector<BlockSplitter_T> kernels;

            // Prepare kernels.
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                kernels.emplace_back( S, T, thread, far_theta2, near_theta2 );
            }
            
            // Tree traversal.
            if( (S.SplitThreshold()==1) && (T.SplitThreshold()==1) )
            {
                ClusterTreePairTraversor<BlockSplitter_T,symmetricQ,true > traversor ( S, T, kernels );

                traversor.Traverse();
            }
            else
            {
                ClusterTreePairTraversor<BlockSplitter_T,symmetricQ,false> traversor ( S, T, kernels );

                traversor.Traverse();
            }
            
//            // Reduce kernels.
//            ParallelDo(
//                [&kernels]( const Int thread )
//                {
//                    kernels[thread].inter_idx.Finalize();
//                    kernels[thread].verynear_idx.Finalize();
//                    kernels[thread].near_idx.Finalize();
//                    kernels[thread].far_idx.Finalize();
//                },
//                thread_count
//            );
            
            std::vector<PairAggregator<Int,Int,LInt>> inter_idx    (thread_count);
            std::vector<PairAggregator<Int,Int,LInt>> verynear_idx (thread_count);
            std::vector<PairAggregator<Int,Int,LInt>> near_idx     (thread_count);
            std::vector<PairAggregator<Int,Int,LInt>> far_idx      (thread_count);
            
            Size_T intersection_count = 0;
            
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                intersection_count += static_cast<Size_T>( kernels[thread].inter_idx.Size() );
                
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
                
                
            }
            
            if( intersection_count > 0 )
            {
                wprint(className()+"::ComputeBlocks detected "+ToString(intersection_count)+" intersections.");
            }


            // Primitive intersection data.
            inter = Inter_Pattern_T( inter_idx, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, symmetricQ );

            TOOLS_PDUMP(inter.Stats());
            

            // Very near field interaction data.
            verynear = VeryNear_Pattern_T( verynear_idx, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, symmetricQ );
            
            TOOLS_PDUMP(verynear.Stats());


            // Near field interaction data.
            near = Near_Pattern_T( near_idx, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, symmetricQ );
            
            
            // Far field interaction data.
            far = Far_Pattern_T( far_idx, S.ClusterCount(), T.ClusterCount(),
                    thread_count, false, symmetricQ );
            
            
            blocks_initialized = true;
            
            
            TOOLS_PTOC(className()+"::ComputeBlocks");

        }; // ComputeBlocks
        
    private:
        
        static std::string className()
        {
            return  "BlockClusterTree<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+ToString(symmetricQ)+">";
        }
      
    public:
        
        std::string ClassName() const override
        {
            return className();
        }
    
    }; // class BlockClusterTree
    
} //namespace Repulsor
