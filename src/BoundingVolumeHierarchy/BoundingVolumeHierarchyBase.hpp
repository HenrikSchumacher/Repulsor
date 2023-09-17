#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_>
    class BoundingVolumeHierarchyBase : public CachedObject
    {
        ASSERT_FLOAT(Real_   );
        ASSERT_INT  (Int_    );
        ASSERT_INT  (LInt_   );
        ASSERT_FLOAT(SReal_  );
        ASSERT_FLOAT(ExtReal_);
 
        
    public:
        
        using Real                  = Real_;
        using Int                   = Int_;
        using SReal                 = SReal_;
        using ExtReal               = ExtReal_;
        using LInt                  = LInt_;
        
        using SparseBinaryMatrix_T  = Sparse::BinaryMatrixCSR<     Int,LInt>;
        
        explicit BoundingVolumeHierarchyBase(
            const BoundingVolumeHierarchySettings & settings_ = BoundingVolumeHierarchySettings()
        )
        :   settings ( settings_ )
        {}
        
        virtual ~BoundingVolumeHierarchyBase() = default;
        
    public:
        
        BoundingVolumeHierarchyBase( const BoundingVolumeHierarchyBase & rhs ) = default;
        
        BoundingVolumeHierarchyBase( BoundingVolumeHierarchyBase && rhs ) = default;
        
    protected:
        
        static constexpr Int null = static_cast<Int>(0);
        
        const BoundingVolumeHierarchySettings   settings;
        
        // Integer data for the combinatorics of the tree.
        
        Tensor1<Int,Int> P_ordering;           // Reordering of primitives; crucial for communication with outside world
        Tensor1<Int,Int> P_inverse_ordering;   // Inverse ordering of the above; crucial for communication with outside world

        Tensor1<Int,Int> C_begin;
        Tensor1<Int,Int> C_end;
        Tensor1<Int,Int> C_depth;
        Tensor1<Int,Int> C_desc_count;
        Tensor1<Int,Int> C_next;
        Tensor1<Int,Int> C_left;  // list of index of left  child; entry is -1 if no child is present
        Tensor1<Int,Int> C_right; // list of index of right child; entry is -1 if no child is present
        
        Tensor1<Int,Int> leaf_clusters;
        Tensor1<Int,Int> leaf_cluster_lookup;
        Tensor1<Int,Int> leaf_cluster_ptr;

        Tensor1<SReal,Int> P_score_buffer;
        Tensor1<Int,Int>   P_perm_buffer;

        Int depth = -1;
    
        // Container for storing the serialized data of the primitives. Only meant to be accessed by primitive prototypes.
        mutable Tensor2<SReal,Int> P_serialized;
        
        // Container for storing the serialized data of the clusters. Only meant to be accessed by bounding volume prototypes.
        mutable Tensor2<SReal,Int> C_serialized;
        
        Tensor3<SReal,Int> C_thread_serialized;          // False sharing is unlikely as each thread's slice should already be quite large...
        
        Tensor2<Int,Int> thread_cluster_counter;         // TODO: Avoid false sharing!
        
    public:

//##################################################################################################
//##        Get functions
//##################################################################################################
        
        virtual constexpr Int AmbDim() const = 0;
        
        Int ThreadCount() const
        {
            return static_cast<Int>(settings.thread_count);
        }
        
        Int SplitThreshold() const
        {
            return static_cast<Int>(settings.split_threshold);
        }

        Int Depth() const
        {
            return depth;
        }

        Int TreeDepth() const
        {
            return depth;
        }
        
        Int PrimitiveCount() const
        {
            return P_serialized.Dimension(0);
        }

        Int PrimitiveSize() const
        {
            return P_serialized.Dimension(1);
        }

        Int ClusterCount() const
        {
            return C_serialized.Dimension(0);
        }

        Int ClusterSize() const
        {
            return C_serialized.Dimension(1);
        }

        Int LeafClusterCount() const
        {
            return leaf_clusters.Dimension(0);
        }
        
//##################################################################################################
//##        Get functions for containers
//##################################################################################################
        
        
      cref<Tensor1<Int,Int>> LeafClusters() const
        {
            return leaf_clusters;
        }
        
        cref<Tensor1<Int,Int>> LeafClusterLookup() const
        {
            return leaf_cluster_lookup;
        }
        
        cref<Tensor1<Int,Int>> LeafClusterPointers() const
        {
            return leaf_cluster_ptr;
        }
        
        cref<Tensor1<Int,Int>> PrimitiveOrdering() const
        {
            return P_ordering;
        }
        
        cref<Tensor1<Int,Int>> PrimitiveInverseOrdering() const
        {
            return P_inverse_ordering;
        }
        
        cref<Tensor1<Int,Int>> ClusterBegin() const
        {
            return C_begin;
        }

        cref<Tensor1<Int,Int>> ClusterEnd() const
        {
            return C_end;
        }
        
        cref<Tensor1<Int,Int>> ClusterLeft() const
        {
            return C_left;
        }
        
        cref<Tensor1<Int,Int>> ClusterRight() const
        {
            return C_right;
        }

        cref<Tensor1<Int,Int>> ClusterNext() const
        {
            return C_next;
        }
        
        cref<Tensor1<Int,Int>> ClusterDepths() const
        {
            return C_depth;
        }
                
        mref<Tensor2<SReal,Int>> ClusterSerialized() const
        {
            return C_serialized;
        }

        mref<Tensor2<SReal,Int>> PrimitiveSerialized() const
        {
            return P_serialized;
        }
  
    protected:
        
        void ComputePrimitiveToClusterMatrix()
        {
            ptic(className()+"::ComputePrimitiveToClusterMatrix");
            
            SparseBinaryMatrix_T P_to_C (
                ClusterCount(), PrimitiveCount(), PrimitiveCount(), ThreadCount() );
            
            SparseBinaryMatrix_T C_to_P (
                PrimitiveCount(), ClusterCount(), PrimitiveCount(), ThreadCount() );
            
            C_to_P.Outer()[PrimitiveCount()] = PrimitiveCount();
            
            {
                mptr<Int> inner__ = C_to_P.Inner().data();
                
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
                mptr<LInt> i = C_to_P.Outer().data();
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
            
            this->SetCache(
                std::string("PrimitiveToCluster"),
                std::any( std::move( P_to_C ) )
            );
            
            this->SetCache(
                std::string("ClusterToPrimitiveMatrix"),
                std::any( std::move( C_to_P ) )
            );
            
            ptoc(className()+"::ComputePrimitiveToCluster");
        }
        
    public:
        
        
    //        cref<SparseBinaryMatrix_T> PrimitiveAdjacencyMatrix() const
    //        {
    //            return Adj;
    //        }
        
        cref<SparseBinaryMatrix_T> ClusterToPrimitiveMatrix() const
        {
            std::string tag ("ClusterToPrimitiveMatrix");
            
            if( !this->InCacheQ(tag))
            {
                ComputePrimitiveToClusterMatrix();
            }
            
            return std::any_cast<cref<SparseBinaryMatrix_T>>( this->GetCache(tag) );
        }
        
        cref<SparseBinaryMatrix_T> PrimitiveToClusterMatrix() const
        {
            std::string tag ("PrimitiveToClusterMatrix");
            
            if( !this->InCacheQ(tag))
            {
                ComputePrimitiveToClusterMatrix();
            }
            
            return std::any_cast<cref<SparseBinaryMatrix_T>>( this->GetCache(tag) );
        }
      
        
//################################################################################################
//##        General reports
//################################################################################################
        
        virtual std::string Stats() const = 0;
        
    private:
        
        std::string className() const
        {
            return std::string("BoundingVolumeHierarchyBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
    public:
        
        virtual std::string ClassName() const
        {
            return className();
        }

    }; // class BoundingVolumeHierarchyBase
    
} // namespace Repulsor
