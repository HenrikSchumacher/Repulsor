#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_>
    class ClusterTreeBase : public CachedObject
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

        
//        mutable Tensor1<Int,Int> stack_array;
//        mutable Tensor1<Int,Int> queue_array;
        
//        mutable Int stack_array [128];
//        mutable Int queue_array [64];

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
        // TODO: Make this private somehow?
        mutable Tensor2<SReal,Int> P_serialized;
        mutable Tensor2<SReal,Int> P_velocities_serialized;
        mutable Tensor2<SReal,Int> P_updated_serialized;
        
        // Container for storing the serialized data of the clusters. Only meant to be accessed by bounding volume prototypes.
        // TODO: Make this private somehow?
        mutable Tensor2<SReal,Int> C_serialized;
        mutable Tensor2<SReal,Int> C_updated_serialized;
        
        mutable SReal update_time = static_cast<SReal>(0);
        
        // Some temproray shared data that is required for the parallel construction and serialization of the tree.
        
        Tensor3<SReal,Int> C_thread_serialized;          // False sharing is unlikely as each thread's slice should already be quite large...
        
        Tensor2<Int,Int> thread_cluster_counter;                          // TODO: Avoid false sharing!
        
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


        const Tensor1<Int,Int> & LeafClusters() const
        {
            return leaf_clusters;
        }
        
        const Tensor1<Int,Int> & LeafClusterLookup() const
        {
            return leaf_cluster_lookup;
        }
        
        const Tensor1<Int,Int> & LeafClusterPointers() const
        {
            return leaf_cluster_ptr;
        }
        
        const Tensor1<Int,Int> & PrimitiveOrdering() const
        {
            return P_ordering;
        }
        
        const Tensor1<Int,Int> & PrimitiveInverseOrdering() const
        {
            return P_inverse_ordering;
        }
        
        const Tensor1<Int,Int> & ClusterBegin() const
        {
            return C_begin;
        }

        const Tensor1<Int,Int> & ClusterEnd() const
        {
            return C_end;
        }
        
        const Tensor1<Int,Int> & ClusterLeft() const
        {
            return C_left;
        }
        
        const Tensor1<Int,Int> & ClusterRight() const
        {
            return C_right;
        }

        const Tensor1<Int,Int> & ClusterNext() const
        {
            return C_next;
        }
        
        const Tensor1<Int,Int> & ClusterDepths() const
        {
            return C_depth;
        }
                
        Tensor2<SReal,Int> & ClusterSerialized() const
        {
            return C_serialized;
        }
        
        Tensor2<SReal,Int> & ClusterUpdatedSerialized() const
        {
            return C_updated_serialized;
        }

        Tensor2<SReal,Int> & PrimitiveSerialized() const
        {
            return P_serialized;
        }
        
        Tensor2<SReal,Int> & PrimitiveUpdatedSerialized() const
        {
            return P_updated_serialized;
        }
        
        Tensor2<SReal,Int> & PrimitiveVelocitiesSerialized() const
        {
            return P_velocities_serialized;
        }
        
        const DataContainer_T & PrimitiveNearFieldData() const
        {
            return P_near;
        }
        
        const DataContainer_T & PrimitiveFarFieldData() const
        {
            return P_far;
        }
        
        const DataContainer_T & ClusterFarFieldData() const
        {
            return C_far;
        }
        
        virtual const Tensor2<Real,Int> & ClusterMoments() const = 0;
        
//##################################################################################################
//##        Clustering
//##################################################################################################
        
        
//##################################################################################################
//##        Matrix multiplication pipeline
//##################################################################################################
        
        BufferContainer_T & ClusterInputBuffer() const
        {
            return C_in;
        }
        
        BufferContainer_T & ClusterOutputBuffer() const
        {
            return C_out;
        }
        
        BufferContainer_T & PrimitiveInputBuffer() const
        {
            return P_in;
        }
        
        BufferContainer_T & PrimitiveOutputBuffer() const
        {
            return P_out;
        }
        

        const SparseMatrix_T & LowOrderPreProcessor() const
        {
            return lo_pre;
        }
        

        const SparseMatrix_T & LowOrderPostProcessor() const
        {
            return lo_post;
        }
        
        
        const SparseMatrix_T & HighOrderPreProcessor() const
        {
            return hi_pre;
        }
        
        const SparseMatrix_T & HighOrderPostProcessor() const
        {
            return hi_post;
        }
        
        const SparseMatrix_T & MixedOrderPreProcessor() const
        {
            return mi_pre;
        }
        
        const SparseMatrix_T & MixedOrderPostProcessor() const
        {
            return mi_post;
        }

        
        const SparseBinaryMatrix_T & ClusterToPrimitiveMatrix() const
        {
            return C_to_P;
        }
        
        const SparseBinaryMatrix_T & PrimitiveToClusterMatrix() const
        {
            return P_to_C;
        }
        
        void RequireBuffers( const Int nrhs ) const
        {
            ptic(ClassName()+"::RequireBuffers");
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
            ptoc(ClassName()+"::RequireBuffers");
            
        }; // RequireBuffers

        void CleanseBuffers() const
        {
            ptic(ClassName()+"::CleanseBuffers");
            
            P_in.SetZero();
            P_out.SetZero();

            C_in.SetZero();
            C_out.SetZero();
        
            ptoc(ClassName()+"::CleanseBuffers");
        }; // CleanseBuffers

        
        void PercolateUp() const
        {
//            ptic(ClassName()+"::PercolateUp");
            
            switch (settings.tree_perc_alg)
            {
                case TreePercolationAlgorithm::Sequential:
                {
                    ptic(ClassName()+"::PercolateUp_DFS");
                    PercolateUp_DFS( Scalar::Zero<Int> );
                    ptoc(ClassName()+"::PercolateUp_DFS");
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
//            ptoc(ClassName()+"::PercolateUp");
        }; // PercolateUp
        
        void PercolateDown() const
        {
//            ptic(ClassName()+"::PercolateDown");
            switch (settings.tree_perc_alg)
            {
                case TreePercolationAlgorithm::Sequential:
                {
                    ptic(ClassName()+"::PercolateDown_DFS");
                    PercolateDown_DFS( Scalar::Zero<Int> );
                    ptoc(ClassName()+"::PercolateDown_DFS");
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
//            ptoc(ClassName()+"::PercolateDown");
        }; // PercolateUp
        
        
        virtual void PercolateUp_DFS( const Int C, const Int max_depth = 64 ) const = 0;
        
        virtual void PercolateDown_DFS( const Int C, const Int max_depth = 64 ) const = 0;
        
        virtual void PercolateUp_Recursive( const Int C ) const = 0;
        
        virtual void PercolateDown_Recursive( const Int C ) const = 0;
        
        virtual void PercolateUp_Parallel() const = 0;
        
        virtual void PercolateDown_Parallel() const = 0;

    public:
        
        virtual void Pre( const ExtReal * input, const Int nrsh, const OperatorType op_type ) const = 0;
        
        virtual void Post( ExtReal * output, const ExtReal alpha, const ExtReal beta, const OperatorType op_type ) const = 0;
        
        virtual void ClustersToPrimitives( bool add_to = false ) const = 0;
        
        virtual void PrimitivesToClusters( bool add_to = false ) const = 0;
        
        
//#############################################################################################
//##        Neighbors
//#############################################################################################
        
        
    public:
        
        const SparseBinaryMatrix_T & PrimitiveAdjacencyMatrix() const
        {
            std::string tag ("PrimitiveAdjacencyMatrix");
            if( !InPersistentCacheQ(tag))
            {
                ptic(ClassName()+"::" + tag );
                // We compute product of S->lo_pre and  S->lo_post because it is basically the primitive-primitive adjacency matrix, we two primitives are supposed to be adjacencent if they share a common degree of freedom.
//
                SparseBinaryMatrix_T P_adjacency_matrix = lo_pre.DotBinary(lo_post);

                ptoc(ClassName()+"::" + tag );
                
                this->SetPersistentCache(tag, std::any( std::move(P_adjacency_matrix) ) );
            }
            
            return std::any_cast<const SparseBinaryMatrix_T & >(this->GetPersistentCache(tag));
        }
        
        
//##################################################################################################
//##        Derivatives
//##################################################################################################
        
    public:
        
        const DataContainer_T & PrimitiveDNearFieldData() const
        {
            return P_D_near;
        }
        
        const DataContainer_T & PrimitiveDFarFieldData() const
        {
            return P_D_far;
        }
        
        DerivativeContainer_T & ThreadPrimitiveDNearFieldData() const
        {
            return thread_P_D_near;
        }
        
        DerivativeContainer_T & ThreadClusterDFarFieldData() const
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
            ptic(ClassName()+"::CleanseDerivativeBuffers");
            
            thread_P_D_near.SetZero();
            thread_C_D_far.SetZero();

            ptoc(ClassName()+"::CleanseDerivativeBuffers");
        }; // CleanseDerivativeBuffers
        
        void CollectDerivatives() const
        {
            ptic(ClassName()+"::CollectDerivatives");
            
            CollectFarFieldDerivatives ( false );
            CollectNearFieldDerivatives( false );
            
            ptoc(ClassName()+"::CollectDerivatives");
            
        } // CollectDerivatives
        
        
        void CollectPrimitiveEnergies() const
        {
            ptic(ClassName()+"::CollectPrimitiveEnergies");
            
            // Collect energies from primitives and clusters and store them in P_out.
            
            // Partially using the pipeline for derivatives. Not really efficient, but also not performance critical.
            
            const Int far_dim         =  FarDim();
            const Int near_dim        = NearDim();
            const Int cluster_count   = ClusterCount();
            const Int primitive_count = PrimitiveCount();
            const Int thread_count    = ThreadCount();
            
            this->RequireBuffers( static_cast<Int>(1) );
            
            this->CleanseBuffers();
            
            // Collect first entry of thread_C_D_far into C_out.
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                ptr<Real> from = thread_C_D_far.data(thread);
                mut<Real> to   = C_out.data();
                
                ParallelDo(
                    [=]( const Int i )
                    {
                        to[i] += from[far_dim * i];
                    },
                    cluster_count,
                    ThreadCount()
                );
            }
            
            //Before percolating the energies down we have to transform them to densities.
            ParallelDo(
                [=]( const Int i )
                {
                    C_out[i] /= C_far[i][0];
                },
                cluster_count,
                ThreadCount()
            );
            
            PercolateDown();
            
            this->ClustersToPrimitives( false );
            
            // Now the simplex far field densities are stored in P_out.
            // We have to convert them back to energies before we add the energies from thread_P_D_near.
            
            ParallelDo(
                [=]( const Int i )
                {
                    P_out[i] *= P_near[i][0];
                },
                primitive_count,
                thread_count
            );
            
            // Add first entries of thread_P_D_near into P_out.
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                ptr<Real> from = thread_P_D_near.data(thread);
                mut<Real> to   = P_out.data();
                
                ParallelDo(
                    [=]( const Int i )
                    {
                        to[i] += from[near_dim * i];
                    },
                    primitive_count,
                    thread_count
                );
            }
                        
            ptoc(ClassName()+"::CollectPrimitiveEnergies");
            
        } // CollectPrimitiveEnergies
        
        void CollectPrimitiveEnergies( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const
        {
            ptic(ClassName()+"::CollectPrimitiveEnergies");
            
            //Collect the energies into P_out.
            CollectPrimitiveEnergies();
            
            const Int primitive_count = PrimitiveCount();
            
            // Copy the values to output. We must not forget reorder!
            ParallelDo(
                [=]( const Int i )
                {
                    const Int j = P_inverse_ordering[i];
                    output[i] = static_cast<ExtReal>(P_out[j]);
                },
                primitive_count,
                ThreadCount()
            );
            
            ptoc(ClassName()+"::CollectPrimitiveEnergies");
       
        } // CollectPrimitiveEnergies
            
        void CollectDensity( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const
        {
            ptic(ClassName()+"::CollectDensity");
            
            // Partially using the pipeline for derivatives. Not really efficient, but also not performance critical.
            this->RequireBuffers( static_cast<Int>(1) );
            
            // Compute dual volume vectors.
            fill_buffer( &P_in[0], static_cast<Real>(1), PrimitiveCount() );

            Tensor1<ExtReal,Int> dual_volumes ( lo_post.RowCount() );
        
            lo_post.template Dot<1>(
                static_cast<Real>(1),    P_in.data(),
                static_cast<ExtReal>(0), dual_volumes.data(),
                static_cast<Int>(1)
            );
        
            //Collect the energies into P_out.
            CollectPrimitiveEnergies();
            
            // Divide by primitive volumes to get densities.
            ParallelDo(
                [=]( const Int j )
                {
                    P_out[j] /= P_near[j][0];
                },
                PrimitiveCount(),
                ThreadCount()
            );
            
            // Distribute the energy densities to energies per vertex. (Note that lo_post also multiplies by the primitives' volumes!)
            lo_post.template Dot<1>(
                static_cast<Real>(1),        P_out.data(),
                static_cast<ExtReal>(addTo), output,
                static_cast<Int>(1)
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
            
            ptoc(ClassName()+"::CollectDensity");
            
        } // CollectDensity

//##################################################################################################
//##        Updates
//##################################################################################################

        virtual void SemiStaticUpdate(
            const Tensor2<Real,Int> & P_near_, const Tensor2<Real,Int> & P_far_ ) const = 0;
        
        virtual void SemiStaticUpdate( ptr<Real> P_near_, ptr<Real> P_far_ ) const = 0;
        
//################################################################################################
//##        Moments
//################################################################################################
        
//        virtual void RequireClusterMoments( const Int moment_degree_ ) const = 0;
        
//################################################################################################
//##        General reports
//################################################################################################
        
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
