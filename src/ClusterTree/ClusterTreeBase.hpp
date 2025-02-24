#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_>
    class ClusterTreeBase : public CachedObject
    {
        static_assert(IntQ< Int_>,"");
        static_assert(IntQ<LInt_>,"");
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
    public:
        
        using Real                  = Real_;
        using Int                   = Int_;
        using SReal                 = SReal_;
        using ExtReal               = ExtReal_;
        using LInt                  = LInt_;
        
        using SparseMatrix_T        = Sparse::MatrixCSR      <Real,Int,LInt>;
        using SparseBinaryMatrix_T  = Sparse::BinaryMatrixCSR<     Int,LInt>;
        
        using DataContainer_T       = Tensor2<Real,Int>;
        using BufferContainer_T     = Tensor1<Real,Int>;

        using DerivativeContainer_T = ThreadTensor3<Real,Int>;
        using Accumulator_T         = ThreadTensor3<Real,Int>;
        
        explicit ClusterTreeBase(
            const ClusterTreeSettings & settings_ = ClusterTreeSettings()
        )
        :   settings ( settings_ )
        {}
        
        virtual ~ClusterTreeBase() = default;
        
    public:
        
        ClusterTreeBase(const ClusterTreeBase & rhs) = default;
        
        ClusterTreeBase(ClusterTreeBase && rhs) = default;
        
    protected:
        
        static constexpr Int null = static_cast<Int>(0);
        
        const ClusterTreeSettings   settings;

        mutable DataContainer_T       P_near;
        mutable DataContainer_T       P_D_near;
        
        mutable DataContainer_T       P_far;
        mutable DataContainer_T       P_D_far;
        
        mutable DataContainer_T       P_moments;
        
        mutable BufferContainer_T     P_in;
        mutable BufferContainer_T     P_out;
        mutable DerivativeContainer_T thread_P_D_near;
        
        mutable Int max_buffer_dim = 0;
        mutable Int     buffer_dim = 0;
        
        mutable DataContainer_T       C_far;
        mutable DataContainer_T       C_moments;
        
        mutable BufferContainer_T     C_in;
        mutable BufferContainer_T     C_out;
        mutable DerivativeContainer_T thread_C_D_far;
        
        mutable Accumulator_T         VF_acc;
        mutable Accumulator_T         NF_acc;
        mutable Accumulator_T         FF_acc;
        
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
    
        
        mutable bool pre_post_initialized = false;
        mutable bool mi_pre_post_initialized = false;
        
        mutable SparseMatrix_T  hi_pre;
        mutable SparseMatrix_T  hi_post;

        mutable SparseMatrix_T  lo_pre;
        mutable SparseMatrix_T  lo_post;
        
        mutable SparseMatrix_T  mi_pre;
        mutable SparseMatrix_T  mi_post;

        mutable SparseBinaryMatrix_T P_to_C;
        mutable SparseBinaryMatrix_T C_to_P;
        
        // Container for storing the serialized data of the primitives. Only meant to be accessed by primitive prototypes.
        mutable Tensor2<SReal,Int> P_serialized;
        mutable Tensor2<SReal,Int> P_velocities_serialized;
        mutable Tensor2<SReal,Int> P_updated_serialized;
        
        // Container for storing the serialized data of the clusters. Only meant to be accessed by bounding volume prototypes.
        mutable Tensor2<SReal,Int> C_serialized;
        mutable Tensor2<SReal,Int> C_updated_serialized;
        
        mutable SReal update_time = Scalar::Zero<SReal>;
        
        // Some temproray shared data that is required for the parallel construction and serialization of the tree.
        
        Tensor3<SReal,Int> C_thread_serialized;          // False sharing is unlikely as each thread's slice should already be quite large...
        
        Tensor2<Int,Int> thread_cluster_counter;         // TODO: Avoid false sharing!
        
    public:

//##################################################################################################
//##        Get functions
//##################################################################################################
        
        virtual constexpr Int AmbDim() const = 0;
        
        virtual           Int NearDim() const = 0;
        
        virtual constexpr Int FarDim() const = 0;
        
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
        
        SReal UpdateTime() const
        {
            return update_time;
        }
        
//##################################################################################################
//##        Get functions for containers
//##################################################################################################
        
        
        Int BufferDimension() const
        {
            return buffer_dim;
        }
        
        Int MaxBufferDimension() const
        {
            return max_buffer_dim;
        }


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
        
        mref<Tensor2<SReal,Int>> ClusterUpdatedSerialized() const
        {
            return C_updated_serialized;
        }

        mref<Tensor2<SReal,Int>> PrimitiveSerialized() const
        {
            return P_serialized;
        }
        
        mref<Tensor2<SReal,Int>> PrimitiveUpdatedSerialized() const
        {
            return P_updated_serialized;
        }
        
        mref<Tensor2<SReal,Int>> PrimitiveVelocitiesSerialized() const
        {
            return P_velocities_serialized;
        }
        
        cref<DataContainer_T> PrimitiveNearFieldData() const
        {
            return P_near;
        }
        
        cref<DataContainer_T> PrimitiveFarFieldData() const
        {
            return P_far;
        }
        
        cref<DataContainer_T> ClusterFarFieldData() const
        {
            return C_far;
        }
        
        virtual cref<Tensor2<Real,Int>> ClusterMoments() const = 0;
        
//#################################################################################
//##        Clustering
//#################################################################################
        
        
//#################################################################################
//##        Matrix multiplication pipeline
//#################################################################################
        
        mref<BufferContainer_T> ClusterInputBuffer() const
        {
            return C_in;
        }
        
        mref<BufferContainer_T> ClusterOutputBuffer() const
        {
            return C_out;
        }
        
        mref<BufferContainer_T> PrimitiveInputBuffer() const
        {
            return P_in;
        }
        
        mref<BufferContainer_T> PrimitiveOutputBuffer() const
        {
            return P_out;
        }
        

        cref<SparseMatrix_T> LowOrderPreProcessor() const
        {
            return lo_pre;
        }
        

        cref<SparseMatrix_T> LowOrderPostProcessor() const
        {
            return lo_post;
        }
        
        
        cref<SparseMatrix_T> HighOrderPreProcessor() const
        {
            return hi_pre;
        }
        
        cref<SparseMatrix_T> HighOrderPostProcessor() const
        {
            return hi_post;
        }
        
        cref<SparseMatrix_T> MixedOrderPreProcessor() const
        {
            return mi_pre;
        }
        
        cref<SparseMatrix_T> MixedOrderPostProcessor() const
        {
            return mi_post;
        }

        
        cref<SparseBinaryMatrix_T> ClusterToPrimitiveMatrix() const
        {
            return C_to_P;
        }
        
        cref<SparseBinaryMatrix_T> PrimitiveToClusterMatrix() const
        {
            return P_to_C;
        }
        
        void RequireBuffers( const Int nrhs ) const
        {
            TOOLS_PTIC(ClassName()+"::RequireBuffers");
            // TODO: parallelize allocation
            if( nrhs > max_buffer_dim )
            {
        //        print("Reallocating buffers to max_buffer_dim = " + std::to_string(nrsh) + "." );
                max_buffer_dim = nrhs;
                
                P_in  = Tensor1<Real,Int> ( PrimitiveCount() * max_buffer_dim );
                
                P_out = Tensor1<Real,Int> ( PrimitiveCount() * max_buffer_dim );

                C_in  = Tensor1<Real,Int> ( ClusterCount()   * max_buffer_dim );
                
                C_out = Tensor1<Real,Int> ( ClusterCount()   * max_buffer_dim );
                
            }
            buffer_dim = nrhs;
            TOOLS_PTOC(ClassName()+"::RequireBuffers");
            
        }; // RequireBuffers

        void CleanseBuffers() const
        {
            TOOLS_PTIC(ClassName()+"::CleanseBuffers");
            
            P_in.SetZero();
            P_out.SetZero();

            C_in.SetZero();
            C_out.SetZero();
        
            TOOLS_PTOC(ClassName()+"::CleanseBuffers");
        }; // CleanseBuffers

        
        void PercolateUp() const
        {
//            TOOLS_PTIC(ClassName()+"::PercolateUp");
            
            switch (settings.tree_perc_alg)
            {
                case TreePercolationAlgorithm::Sequential:
                {
                    TOOLS_PTIC(ClassName()+"::PercolateUp_DFS");
                    PercolateUp_DFS( Scalar::Zero<Int> );
                    TOOLS_PTOC(ClassName()+"::PercolateUp_DFS");
                    break;
                }
                case TreePercolationAlgorithm::Recursive:
                {
                    PercolateUp_Recursive( Scalar::Zero<Int> );
                    break;
                }
                case TreePercolationAlgorithm::Parallel:
                {
                    PercolateUp_Parallel();
                    break;
                }
                default:
                {
                    PercolateUp_DFS( Scalar::Zero<Int> );
                    break;
                }
            }
//            TOOLS_PTOC(ClassName()+"::PercolateUp");
        }; // PercolateUp
        
        void PercolateDown() const
        {
//            TOOLS_PTIC(ClassName()+"::PercolateDown");
            switch (settings.tree_perc_alg)
            {
                case TreePercolationAlgorithm::Sequential:
                {
                    TOOLS_PTIC(ClassName()+"::PercolateDown_DFS");
                    PercolateDown_DFS( Scalar::Zero<Int> );
                    TOOLS_PTOC(ClassName()+"::PercolateDown_DFS");
                    break;
                }
                case TreePercolationAlgorithm::Recursive:
                {
                    PercolateDown_Recursive( Scalar::Zero<Int> );
                    break;
                }
                case TreePercolationAlgorithm::Parallel:
                {
                    PercolateDown_Parallel();
                    break;
                }
                default:
                {
                    PercolateDown_DFS( Scalar::Zero<Int> );
                }
            }
//            TOOLS_PTOC(ClassName()+"::PercolateDown");
        }; // PercolateUp
        
        
        virtual void PercolateUp_DFS( const Int C, const Int max_depth = 64 ) const = 0;
        
        virtual void PercolateDown_DFS( const Int C, const Int max_depth = 64 ) const = 0;
        
        virtual void PercolateUp_Recursive( const Int C ) const = 0;
        
        virtual void PercolateDown_Recursive( const Int C ) const = 0;
        
        virtual void PercolateUp_Parallel() const = 0;
        
        virtual void PercolateDown_Parallel() const = 0;

    public:

        // Reads P_in = X[:,0..ldX],
        // where P_in is the input buffer of the tree.
        virtual void Pre(
            const ExtReal * X, const Int ldX,
            const Int nrsh, const OperatorType op_type
        ) const = 0;
        
        // Writes Y[:,0..ldY] = alpha * P_out + beta * Y[:,0..ldY],
        // where P_out is the output buffer of the tree.
        virtual void Post(
            const ExtReal alpha, const ExtReal beta, ExtReal * Y, const Int ldY,
            const OperatorType op_type
        ) const = 0;
        
        virtual void ClustersToPrimitives( bool add_to = false ) const = 0;
        
        virtual void PrimitivesToClusters( bool add_to = false ) const = 0;
        
        
//#############################################################################################
//##        Neighbors
//#############################################################################################
        
        
    public:
        
        cref<SparseBinaryMatrix_T> PrimitiveAdjacencyMatrix() const
        {
            std::string tag ("PrimitiveAdjacencyMatrix");
            if( !InPersistentCacheQ(tag))
            {
                TOOLS_PTIC(ClassName()+"::" + tag );
                // We compute product of S->lo_pre and  S->lo_post because it is basically the primitive-primitive adjacency matrix, we two primitives are supposed to be adjacencent if they share a common degree of freedom.
                
                SparseBinaryMatrix_T P_adjacency_matrix = lo_pre.DotBinary(lo_post);

                TOOLS_PTOC(ClassName()+"::" + tag );
                
                this->SetPersistentCache(tag, std::move(P_adjacency_matrix) );
            }
            
            return this->template GetPersistentCache<SparseBinaryMatrix_T>(tag);
        }
        
        
//##########################################################################################
//##        Derivatives
//##########################################################################################
        
    public:
        
        cref<DataContainer_T> PrimitiveDNearFieldData() const
        {
            return P_D_near;
        }
        
        cref<DataContainer_T> PrimitiveDFarFieldData() const
        {
            return P_D_far;
        }
        
        mref<DerivativeContainer_T> ThreadPrimitiveDNearFieldData() const
        {
            return thread_P_D_near;
        }
        
        mref<DerivativeContainer_T> ThreadClusterDFarFieldData() const
        {
            return thread_C_D_far;
        }
        
        Accumulator_T & VF_Accumulator() const
        {
            return VF_acc;
        }
        
        Accumulator_T & NF_Accumulator() const
        {
            return NF_acc;
        }
        
        Accumulator_T & FF_Accumulator() const
        {
            return FF_acc;
        }
        
        virtual void CollectNearFieldDerivatives( bool addto = false ) const = 0;
        
        virtual void CollectFarFieldDerivatives ( bool addto = false ) const = 0;
        
        void CleanseDerivativeBuffers() const
        {
            TOOLS_PTIC(ClassName()+"::CleanseDerivativeBuffers");
            
            thread_P_D_near.SetZero();
            thread_C_D_far.SetZero();

            TOOLS_PTOC(ClassName()+"::CleanseDerivativeBuffers");
        }; // CleanseDerivativeBuffers
        
        void CollectDerivatives() const
        {
            TOOLS_PTIC(ClassName()+"::CollectDerivatives");
            
            CollectFarFieldDerivatives ( false );
            CollectNearFieldDerivatives( false );
            
            TOOLS_PTOC(ClassName()+"::CollectDerivatives");
            
        } // CollectDerivatives
        
        
        void CollectPrimitiveDensities() const
        {
            TOOLS_PTIC(ClassName()+"::CollectPrimitiveDensities");
            
            // Collect energies from primitives and clusters and store them in P_out.
            
            // Partially using the pipeline for derivatives. Not fully optimized, but also not performance critical.
            
            const Int far_dim         =  FarDim();
            const Int near_dim        = NearDim();
            const Int cluster_count   = ClusterCount();
            const Int primitive_count = PrimitiveCount();
            const Int thread_count    = ThreadCount();
            
            this->RequireBuffers( static_cast<Int>(1) );
            
            this->CleanseBuffers();
            
            // Collect entries of thread_C_D_far into C_out.
            for( Int thread = 0; thread < thread_count; ++thread )
            {   
                cptr<Real> from = thread_C_D_far.data(thread);
                mptr<Real> to   = C_out.data();
                
                ParallelDo(
                    [=]( const Int i )
                    {
                        to[i] += from[far_dim * i];
                    },
                    cluster_count,
                    ThreadCount()
                );
            }
            
            PercolateDown();
            
            this->ClustersToPrimitives( false );
            
            // Now the simplex far field densities are stored in P_out.
            // Add entries of thread_P_D_near into P_out.
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                cptr<Real> from = thread_P_D_near.data(thread);
                mptr<Real> to   = P_out.data();
                
                ParallelDo(
                    [=]( const Int i )
                    {
                        to[i] += from[near_dim * i];
                    },
                    primitive_count,
                    thread_count
                );
            }
                        
            TOOLS_PTOC(ClassName()+"::CollectPrimitiveDensities");
            
        } // CollectPrimitiveDensities
        
        void CollectPrimitiveDensities( mptr<ExtReal> output, const ExtReal weight, bool addTo = false ) const
        {
            TOOLS_PTIC(ClassName()+"::CollectPrimitiveDensities");
            
            //Collect the energies into P_out.
            CollectPrimitiveDensities();
            
            const Int primitive_count = PrimitiveCount();
            
            // Copy the values to output. We also have to reorder.
            ParallelDo(
                [=,this]( const Int i )
                {
                    const Int j = P_inverse_ordering[i];
                    output[i] = static_cast<ExtReal>(P_out[j]);
                },
                primitive_count,
                ThreadCount()
            );
            
            TOOLS_PTOC(ClassName()+"::CollectPrimitiveDensities");
       
        } // CollectPrimitiveDensities
            
        void CollectVertexDensities( mptr<ExtReal> output, const ExtReal weight, bool addTo = false ) const
        {
            TOOLS_PTIC(ClassName()+"::CollectVertexDensities");
            
            // Partially using the pipeline for derivatives. Not optimized efficient, but also not performance critical.
            this->RequireBuffers( static_cast<Int>(1) );
            
            // Compute dual volume vectors.
            // It's actually not the dual volumes because we do not divide by AMB_DIM+1.
            // But this way we also do not have to divide by AMB_DIM+1 when we average from primitives to simplices.
            fill_buffer( &P_in[0], Scalar::One<Real>, PrimitiveCount() );

            Tensor1<ExtReal,Int> dual_volumes ( lo_post.RowCount() );
        
            lo_post.template Dot<1>(
                Scalar::One<Real>,     P_in.data(),
                Scalar::Zero<ExtReal>, dual_volumes.data(),
                Scalar::One<Int>
            );
        
            //Collect the energies into P_out.
            CollectPrimitiveDensities();
            
            // Distribute the energy densities to energies per vertex. (Note that lo_post also multiplies by the primitives' volumes!)
            lo_post.template Dot<1>(
                Scalar::One<Real>,           P_out.data(),
                static_cast<ExtReal>(addTo), output,
                Scalar::One<Int>
            );

            // Finally, we divide by the dual volumes to obtain the vertex densities.
            ParallelDo(
                [=]( const Int i )
                {
                    output[i] /= dual_volumes[i];
                },
                lo_post.RowCount(),
                ThreadCount()
            );
            
            TOOLS_PTOC(ClassName()+"::CollectVertexDensities");
            
        } // CollectVertexDensities

//#################################################################################
//##        Updates
//#################################################################################

        virtual void SemiStaticUpdate(
            cref<Tensor2<Real,Int>> P_near_,
            cref<Tensor2<Real,Int>> P_far_
        ) const = 0;
        
        virtual void SemiStaticUpdate( cptr<Real> P_near_, cptr<Real> P_far_ ) const = 0;
        
//#################################################################################
//##        Distance queries
//#################################################################################
        
        virtual void FindNearestToPoints(
            cptr<ExtReal> query_points,
            Int           query_point_count,
            mptr<Int>     nearest_simplices,
            mptr<ExtReal> nearest_dist,
            mptr<ExtReal> nearest_points
        ) const = 0;
        
//#################################################################################
//##        Moments
//#################################################################################
        
//        virtual void RequireClusterMoments( const Int moment_degree_ ) const = 0;
        
//#################################################################################
//##        General reports
//#################################################################################
        
        virtual std::string Stats() const = 0;
        
    private:
        std::string className() const
        {
            return std::string("ClusterTreeBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
    public:
        
        virtual std::string ClassName() const
        {
            return className();
        }

    }; // class ClusterTreeBase
    
} // namespace Repulsor
