#pragma once

#include "BoundingVolumeHierarchy/Settings.hpp"
#include "BoundingVolumeHierarchy/BoundingVolumeHierarchyBase.hpp"

namespace Repulsor
{
    
    // A very general BVH.
    // Using PolytopeBases as primitive.
    // - The number of vertices in the polytope is expected to be the same for all primitives,
    //   but need not be known at compile time.
    // Using arbitrary descendants of BoundingVolumeBase as bounding volume.
    
    template<int AMB_DIM_, typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_>
    class BoundingVolumeHierarchy : public BoundingVolumeHierarchyBase<Real_,Int_,LInt_,SReal_,ExtReal_>
    {
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = LInt_;

        static constexpr Int AMB_DIM = AMB_DIM_;
        
        using Base_T               = BoundingVolumeHierarchyBase<Real,Int,LInt,SReal,ExtReal>;
        
        using SparseBinaryMatrix_T = typename Base_T::SparseBinaryMatrix_T;
        
        using Cluster_T         = Cluster<Int>;
        
        using Primitive_T       = PolytopeBase<AMB_DIM,GJK_Real,Int,SReal>;
        using BoundingVolume_T  = BoundingVolumeBase<AMB_DIM,GJK_Real,Int,SReal>;
        
    public:
        
        using Base_T::SplitThreshold;
        using Base_T::ThreadCount;
        using Base_T::PrimitiveCount;
        using Base_T::ClusterCount;
        using Base_T::LeafClusterCount;
        
    private:
        
        template<typename Int>
        struct Cluster // slim POD container to hold only the data relevant for the construction phase in the tree, before it is serialized
        {
            static_assert(SignedIntQ<Int>,"");
            
        public:
            
            Cluster() = default;

            Cluster( Int thread_, Int l_ID_, Int begin_, Int end_, Int depth_ )
            :   thread(thread_)
            ,   l_ID(l_ID_)
            ,   begin(begin_)
            ,   end(end_)
            ,   depth(depth_)
            ,   max_depth(depth_)
            ,   descendant_count(0)
            ,   descendant_leaf_count(0)
            ,   left(nullptr)
            ,   right(nullptr)
            {}
            
            ~Cluster() = default;
            
            const Int thread = -1;               // thread that created this cluster
            const Int l_ID = -1;                 // local ID of cluster within this thread
            const Int begin = 0;                 // position of first primitive in cluster relative to array ordering
            const Int end = 0;                   // position behind last primitive in cluster relative to array ordering
            const Int depth = 0;                 // depth within the tree -- not absolutely necessary but nice to have for plotting images
            
            Int max_depth = 0;              // maximal depth of all descendants of this cluster
            Int descendant_count = 0;       // number of descendents of cluster, _this cluster included_
            Int descendant_leaf_count = 0;  // number of leaf descendents of cluster
            
            Int g_ID = 0;                   // global ID of cluster
            Int leafs_before_count = 0;     // number of leafs before cluster
            
            Cluster<Int> * left  = nullptr; // left child
            Cluster<Int> * right = nullptr; // right child
            
        }; //Cluster
        
    private:
            
        static constexpr Int null = static_cast<Int>(0);
        static constexpr Int one  = static_cast<Int>(1);
        
        using Base_T::P_ordering;
        using Base_T::P_inverse_ordering;
        using Base_T::P_score_buffer;
        using Base_T::P_perm_buffer;
        using Base_T::P_serialized;
        
        using Base_T::C_begin;
        using Base_T::C_end;
        using Base_T::C_depth;
        using Base_T::C_right;
        using Base_T::C_left;
        using Base_T::C_desc_count;
        using Base_T::C_next;
        using Base_T::C_serialized;
        using Base_T::C_thread_serialized;
        
        using Base_T::leaf_clusters;
        using Base_T::leaf_cluster_ptr;
        using Base_T::leaf_cluster_lookup;
        using Base_T::thread_cluster_counter;
        
        using Base_T::depth;
        using Base_T::settings;
        
//        SparseBinaryMatrix_T Adj;
        SparseBinaryMatrix_T C_to_P;
        SparseBinaryMatrix_T P_to_C;
        
    public:
        
        BoundingVolumeHierarchy()
        :   P_proto (1)
        ,   C_proto (1)
        {
            ptic(className()+" default constructor");
            ptoc(className()+" default constructor");
        }
        
        // To allow polymorphism, we require the user to create instances of the desired types for the primitives and the bounding volumes, so that we can Clone() them.
        BoundingVolumeHierarchy(
            cref<Primitive_T>             P_proto_,
            const Tensor2<SReal,Int>   && P_serialized_,
            cref<BoundingVolume_T>        C_proto_,
            const Tensor1<Int,Int>     && P_ordering_,     // User may provide a preordering.
//            const SparseBinaryMatrix_T && Adj_,            // User has to provided adjacency matrix of the primitives.
            cref<BoundingVolumeHierarchySettings> settings_ = BoundingVolumeHierarchySettings()
        )
        :   Base_T       ( settings_                )
        ,   P_proto      ( ThreadCount()            )
        ,   C_proto      ( ThreadCount()            )
//        ,   Adj          ( std::move(Adj_)          )
        {
            ptic(className()+"()");
            
            
            P_serialized = std::move(P_serialized_);
            
            P_ordering   =  (
                ( P_ordering_.Dimension(0) == P_serialized.Dimension(0) )
                ? std::move( P_ordering_ )
                : iota<Int,Int>( P_serialized.Dimension(0) )
            );
            
//            if( P_serialized.Dimension(0) != Adj.RowCount() )
//            {
//                eprint(className()+": Number of primitives does not match size of adjacency matrix.");
//            }
//
//            if( Adj.RowCount() != Adj.RowColCount() )
//            {
//                eprint(className()+": adjacency matrix must be a square matrix.");
//            }
            
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
        // False sharing prevented by alignment of PrimitiveBase.
    
        
        std::vector<std::vector<Cluster_T *>> tree_rows_ptr;
        
#include "BoundingVolumeHierarchy/Split.hpp"
#include "BoundingVolumeHierarchy/Serialize.hpp"
        
    private:
        
        void ComputeClusters()
        {
            ptic(className()+"::ComputeClusters");
            
            // Request some temporary memory for threads.
            
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
