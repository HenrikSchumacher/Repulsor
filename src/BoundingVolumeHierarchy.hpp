#pragma once

#include "BoundingVolumeHierarchy/Node.hpp"
#include "BoundingVolumeHierarchy/BoundingVolumeHierarchyBase.hpp"

namespace Repulsor
{
    template<int AMB_DIM_, typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_>
    class BoundingVolumeHierarchy : public BoundingVolumeHierarchyBase<Real_,Int_,LInt_,SReal_,ExtReal_>
    {
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = LInt_;
        
        using Base_T               = BoundingVolumeHierarchyBase<Real,Int,LInt,SReal,ExtReal>;
        
        using SparseMatrix_T       = typename Base_T::SparseMatrix_T;
        using SparseBinaryMatrix_T = typename Base_T::SparseBinaryMatrix_T;
        
        static constexpr Int AMB_DIM = AMB_DIM_;
        
        using Cluster_T         = Cluster<Int>;
        
        using       Primitive_T =       PolytopeBase<AMB_DIM,GJK_Real,Int,SReal>;
        using MovingPrimitive_T = MovingPolytopeBase<AMB_DIM,GJK_Real,Int,SReal>;
        
        using BoundingVolume_T  = AABB<AMB_DIM,GJK_Real,Int,SReal>;
        
        
        using DataContainer_T   = typename Base_T::DataContainer_T;
        using BufferContainer_T = typename Base_T::BufferContainer_T;
        
        // In principle, ThreadTensor3<Real,Int> should have better scaling on multiple socket machines, because I tried to encourages that the thread-local arrays are allocated on local RAM. -- On my tiny Quad Core however, it performs a bit _WORSE_ than Tensor3<Real,Int>.
        using DerivativeContainer_T = typename Base_T::DerivativeContainer_T;
        
        static constexpr Int  FAR_DIM = 1 + AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        static constexpr Int NEAR_DIM = 1 + AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        
    public:
        
        using Base_T::SplitThreshold;
        using Base_T::ThreadCount;
        using Base_T::PrimitiveCount;
        using Base_T::ClusterCount;
        using Base_T::LeafClusterCount;
        //        using Base_T::RequireClusterMoments;
        
    protected:
            
        static constexpr Int null = static_cast<Int>(0);
        static constexpr Int one  = static_cast<Int>(1);
        
        using Base_T::P_ordering;
        using Base_T::P_inverse_ordering;
        using Base_T::P_score_buffer;
        using Base_T::P_perm_buffer;
        using Base_T::P_serialized;
        using Base_T::P_updated_serialized;
        using Base_T::P_velocities_serialized;
        using Base_T::P_near;
        using Base_T::P_D_near;
        using Base_T::P_far;
        using Base_T::P_D_far;
        using Base_T::thread_P_D_near;
        using Base_T::P_moments;
        using Base_T::P_in;
        using Base_T::P_out;
        
        using Base_T::C_begin;
        using Base_T::C_end;
        using Base_T::C_depth;
        using Base_T::C_right;
        using Base_T::C_left;
        using Base_T::C_desc_count;
        using Base_T::C_next;
        using Base_T::C_serialized;
        using Base_T::C_updated_serialized;
        using Base_T::C_thread_serialized;
        using Base_T::C_far;
        using Base_T::thread_C_D_far;
        using Base_T::C_moments;
        using Base_T::C_in;
        using Base_T::C_out;
        
        using Base_T::leaf_clusters;
        using Base_T::leaf_cluster_ptr;
        using Base_T::leaf_cluster_lookup;
        using Base_T::thread_cluster_counter;
        
        using Base_T::depth;
        using Base_T::settings;
        
        using Base_T::buffer_dim;
        
        //        using Base_T::stack_array;
        //        using Base_T::queue_array;
        
        using Base_T::C_to_P;
        using Base_T::P_to_C;
        using Base_T::hi_pre;
        using Base_T::lo_pre;
        using Base_T::mi_pre;
        using Base_T::hi_post;
        using Base_T::lo_post;
        using Base_T::mi_post;
        
    public:
        
        BoundingVolumeHierarchy()
        :   P_proto (1)
        ,   C_proto (1)
        ,   P_moving(1)
        {
            ptic(className()+" default constructor");
            ptoc(className()+" default constructor");
        }
        
        // To allow polymorphism, we require the user to create instances of the desired types for the primitives and the bounding volumes, so that we can Clone() them.
        BoundingVolumeHierarchy(
            cref<Primitive_T>         P_proto_,
            cref<Tensor2<SReal,Int>>  P_serialized_,
            cref<BoundingVolume_T>    C_proto_,
            cref<Tensor1<Int ,Int>>   P_ordering_,
            cref<SparseBinaryMatrix_T>      Adj,
            cref<BoundingVolumeHierarchySettings> settings_ = BoundingVolumeHierarchySettings()
        )
        :   Base_T( settings_ )
        ,   P_proto      ( ThreadCount() )
        ,   C_proto      ( ThreadCount() )
        ,   P_moving     ( ThreadCount() )
        {
            ptic(className()+"()");
            
            P_serialized = P_serialized_;   // It's a unneccessary copy, but not too expensive.
            
            if( P_ordering_.Dimension(0) == P_serialized.Dimension(0) )
            {
                P_ordering = P_ordering_;
            }
            else
            {
                P_ordering = iota<Int,Int>( P_serialized.Dimension(0) );
            }
            
            ParallelDo(
                [&,this]( const Int thread )
                {
                    P_proto[thread] = P_proto_.Clone();
                    
                    C_proto[thread] = C_proto_.Clone();
                },
                ThreadCount()
            );
            
            this->ComputeClusters();
            
            ptoc(className()+"()");
        }
        
        virtual ~BoundingVolumeHierarchy() override = default;
        
        BoundingVolumeHierarchy(const BoundingVolumeHierarchy & rhs) = default;
        
        BoundingVolumeHierarchy( BoundingVolumeHierarchy && rhs) = default;
        
        const BoundingVolumeHierarchy & operator=( const BoundingVolumeHierarchy & rhs)
        {
            return BoundingVolumeHierarchy ( rhs );
        };
        
        const BoundingVolumeHierarchy & operator=( BoundingVolumeHierarchy && rhs)
        {
            return BoundingVolumeHierarchy ( std::move(rhs) );
        };
        
    protected:
        
        
        // Each thread gets its own primitive prototype to avoid sharing conflicts.
        std::vector<std::shared_ptr<Primitive_T>> P_proto;
        // False sharing prevented by alignment of PrimitiveBase.
        
        // Each thread gets its own bounding volume prototype to avoid sharing conflicts.
        std::vector<std::shared_ptr<BoundingVolume_T>> C_proto;
//        std::vector<std::shared_ptr<BoundingVolume_T>> C_proto_updated;
        // False sharing prevented by alignment of PrimitiveBase.
        
//        mutable std::vector<std::shared_ptr<MovingPrimitive_T>> P_moving;
        
        std::vector<std::vector<Cluster_T *>> tree_rows_ptr;
        
    public:
        
#include "BoundingVolumeHierarchy/Split_Thread.hpp"
#include "BoundingVolumeHierarchy/Serialize_Thread.hpp"
        
    private:
        
        void ComputeClusters()
        {
            ptic(className()+"::ComputeClusters");
            
            // Request some temporary memory for threads.
            
            //        P_ordering         = iota<Int>   ( PrimitiveCount() );
            P_inverse_ordering = Tensor1<Int,Int>( PrimitiveCount() );
            
            // Padding every row to prevent false sharing.
            thread_cluster_counter = Tensor2<Int,Int>( ThreadCount(), CacheLineWidth, 0 );
            
            C_thread_serialized = Tensor3<SReal,Int>( ThreadCount(), 2*PrimitiveCount(), C_proto[0]->Size() );
            
            const Int thread = 0;
            
            ++thread_cluster_counter(thread,0);
            
            auto * root = new Cluster_T( thread, 0, 0, PrimitiveCount(), 0 );
            
            
            ptic(className()+"::ComputeClusters: Initial bounding volume of root node");
            
            C_proto[thread]->SetPointer( C_thread_serialized.data(thread), 0 );
            C_proto[thread]->FromPrimitives( *P_proto[thread], P_serialized.data(), 0, PrimitiveCount(), ThreadCount() );
            
            ptoc(className()+"::ComputeClusters: Initial bounding volume of root node");
            
            
            Split( root );
            
            Serialize( root );
            
            delete root;
            
            // Free memory for threads.
            C_thread_serialized = Tensor3<SReal,Int>();
            
            thread_cluster_counter = Tensor2<Int,Int>();
            
            ptoc(className()+"::ComputeClusters");
        }
        
    public:
        
        constexpr Int AmbDim() const override
        {
            return AMB_DIM;
        }
        
        const BoundingVolume_T & ClusterPrototype() const
        {
            return *C_proto[0];
        }
        
        const Primitive_T & PrimitivePrototype() const
        {
            return *P_proto[0];
        }
      
    protected:
        
        void ComputePrimitiveToClusterMatrix()
        {
            ptic(className()+"::ComputePrimitiveToClusterMatrix");
            
            P_to_C = SparseBinaryMatrix_T(
                ClusterCount(), PrimitiveCount(), PrimitiveCount(), ThreadCount() );
            
            C_to_P = SparseBinaryMatrix_T(
                PrimitiveCount(), ClusterCount(), PrimitiveCount(), ThreadCount() );
            
            ClusterToPrimitiveMatrix().Outer()[PrimitiveCount()] = PrimitiveCount();
            
            {
                mptr<Int> inner__ = ClusterToPrimitiveMatrix().Inner().data();
                
                ParallelDo(
                    [=,this]( const Int i )
                    {
                        const Int leaf  = leaf_clusters[i];
                        const Int begin = C_begin[leaf];
                        const Int end   = C_end  [leaf];

                        for( Int k = begin; k < end; ++k )
                        {
                            inner__[k] = leaf;
                        }
                    },
                    LeafClusterCount(), ThreadCount()
                );
            }
            
            {
                mptr<LInt> i = ClusterToPrimitiveMatrix().Outer().data();
                mptr< Int> j = P_to_C.Inner().data();
                
                const Int primitive_count = PrimitiveCount();
                
                for( Int k = 0; k < primitive_count; ++k )
                {
                    i[k] = k;
                    j[k] = k;
                }
            }
            
            {
                Int cluster_count = ClusterCount();
                
                mptr<LInt> outer__ = P_to_C.Outer().data();
                cptr< Int> left__  = C_left.data();
                cptr< Int> begin__ = C_begin.data();
                cptr< Int> end__   = C_end.data();
                
                for ( Int C = 0; C < cluster_count; ++C )
                {
                    outer__[C+1] =
                        outer__[C]
                        +
                        static_cast<LInt>(left__[C] < 0)
                        *
                        static_cast<LInt>(end__[C] - begin__[C]);
                }
            }
            
            ptoc(className()+"::ComputePrimitiveToClusterMatrix");
        }
        
    public:
        
        
        std::string Stats() const override
        {
            std::stringstream s;
            
            s
            << "\n==== "+className()+" Stats ====" << "\n\n"
            << " AmbDim()                    = " <<  this->AmbDim() << "\n"
            << " ThreadCount()               = " <<  this->ThreadCount() << "\n"
            << "\n"
            << " PrimitiveCount()            = " <<  this->PrimitiveCount() << "\n"
            << " ClusterCount()              = " <<  this->ClusterCount() << "\n"
            << " LeafClusterCount()          = " <<  this->LeafClusterCount() << "\n"
            << " PrimitivePrototype()        = " <<  this->PrimitivePrototype().ClassName() <<"\n"
            << " ClusterPrototype()          = " <<  this->ClusterPrototype().ClassName() <<"\n"
            
            << " SplitThreshold()            = " <<  this->SplitThreshold() << "\n"
            << " TreeDepth()                 = " <<  this->TreeDepth() << "\n"

            << "\n==== "+className()+" Stats ====\n" << std::endl;

            
            return s.str();
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        static std::string className()
        {
            return  "BoundingVolumeHierarchy<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }

    }; // class BoundingVolumeHierarchy
    
} // namespace Repulsor
