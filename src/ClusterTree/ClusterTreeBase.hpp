#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class ClusterTreeBase
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
        
        using DataContainer_T   = Tensor2<Real,Int>;
        using BufferContainer_T = Tensor1<Real,Int>;
        
        // In principle, ThreadTensor3<Real,Int> should have better scaling on multiple socket machines, because I tried to encourages that the thread-local arrays are allocated on local RAM. -- On my tiny Quad Core however, it performs a bit _WORSE_ than Tensor3<Real,Int>.
        using DerivativeContainer_T = ThreadTensor3<Real,Int>;
        using Accumulator_T         = ThreadTensor3<Real,Int>;
        
        explicit ClusterTreeBase(
            const ClusterTreeSettings & settings_ = ClusterTreeSettings()
        )
        :   settings ( settings_ )
        {
//            switch( settings.tree_perc_alg )
//            {
//                case TreePercolationAlgorithm::Sequential:
//                {
//                       print(className()+" using sequential percolation algorithm.");
//                    logprint(className()+" using sequential percolation algorithm.");
//                    break;
//                }
//                case TreePercolationAlgorithm::Recursive:
//                {
//                       print(className()+" using recursive percolation algorithm.");
//                    logprint(className()+" using recursive percolation algorithm.");
//                    break;
//                }
//                case TreePercolationAlgorithm::Tasks:
//                {
//                       print(className()+" using task-based percolation algorithm.");
//                    logprint(className()+" using task-based percolation algorithm.");
//                    break;
//                }
//                case TreePercolationAlgorithm::Parallel:
//                {
//                       print(className()+" using parallel percolation algorithm.");
//                    logprint(className()+" using parallel percolation algorithm.");
//                    break;
//                }
//                default:
//                {
//                       print(className()+" using sequential percolation algorithm.");
//                    logprint(className()+" using sequential percolation algorithm.");
//                }
//            }
        }
        
//        ClusterTreeBase( const ClusterTreeBase & other )
//        :   settings ( other.settings)
//        ,   P_near   ( other.P_near)
//        ,   P_D_near   ( other.P_D_near)
//        ,   P_far   ( other.P_far)
//        ,   P_D_far   ( other.P_D_far)
//        ,   P_moments   ( other.P_moments)
//        ,   P_in   ( other.P_in)
//        {}
        
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
        mutable bool mixed_pre_post_initialized = false;
        
        mutable Sparse::MatrixCSR<Real,Int,Int>  hi_pre;
        mutable Sparse::MatrixCSR<Real,Int,Int>  hi_post;

        mutable Sparse::MatrixCSR<Real,Int,Int>  lo_pre;
        mutable Sparse::MatrixCSR<Real,Int,Int>  lo_post;
        
        mutable Sparse::MatrixCSR<Real,Int,Int>  mixed_pre;
        mutable Sparse::MatrixCSR<Real,Int,Int>  mixed_post;

        mutable Sparse::BinaryMatrixCSR<Int,Int> P_to_C;
        mutable Sparse::BinaryMatrixCSR<Int,Int> C_to_P;
        
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
        
        mutable bool P_adjacency_matrix_initialized = false;

        mutable Sparse::BinaryMatrixCSR<Int,Int> P_adjacency_matrix;
        
        // Some temproray shared data that is required for the parallel construction and serialization of the tree.
        
        Tensor3<SReal,Int> C_thread_serialized;          // False sharing is unlikely as each thread's slice should already be quite large...
        
        Tensor2<Int,Int> thread_cluster_counter;                          // TODO: Avoid false sharing!
        
        mutable bool parallel_perc_roots_initialized = false;
        mutable std::vector<Int> parallel_perc_roots;
        
    public:

//##################################################################################################
//##        Get functions
//##################################################################################################
        
        virtual Int AmbDim() const = 0;
        
        virtual Int NearDim() const = 0;
        
        virtual Int FarDim() const = 0;
        
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
        

        const Sparse::MatrixCSR<Real,Int,Int> & LowOrderPreProcessor() const
        {
            return lo_pre;
        }
        

        const Sparse::MatrixCSR<Real,Int,Int> & LowOrderPostProcessor() const
        {
            return lo_post;
        }
        
        
        const Sparse::MatrixCSR<Real,Int,Int> & HighOrderPreProcessor() const
        {
            return hi_pre;
        }
        
        const Sparse::MatrixCSR<Real,Int,Int> & HighOrderPostProcessor() const
        {
            return hi_post;
        }
        
        const Sparse::MatrixCSR<Real,Int,Int> & MixedOrderPreProcessor() const
        {
            return mixed_pre;
        }
        
        const Sparse::MatrixCSR<Real,Int,Int> & MixedOrderPostProcessor() const
        {
            return mixed_post;
        }

        
        const Sparse::BinaryMatrixCSR<Int,Int> & ClusterToPrimitiveMatrix() const
        {
            return C_to_P;
        }
        
        const Sparse::BinaryMatrixCSR<Int,Int> & PrimitiveToClusterMatrix() const
        {
            return P_to_C;
        }
        
        void RequireBuffers( const Int cols ) const
        {
            ptic(ClassName()+"::RequireBuffers");
            // TODO: parallelize allocation
            if( cols > max_buffer_dim )
            {
        //        print("Reallocating buffers to max_buffer_dim = " + std::to_string(cols) + "." );
                max_buffer_dim = cols;
                
                P_in  = Tensor1<Real,Int> ( PrimitiveCount() * max_buffer_dim );
                
                P_out = Tensor1<Real,Int> ( PrimitiveCount() * max_buffer_dim );

                C_in  = Tensor1<Real,Int> ( ClusterCount()   * max_buffer_dim );
                
                C_out = Tensor1<Real,Int> ( ClusterCount()   * max_buffer_dim );
                
            }
            buffer_dim = cols;
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
        
        void PrimitivesToClusters( bool add_to = false ) const
        {
            ptic(ClassName()+"::PrimitivesToClusters");
            
            P_to_C.Dot(
                static_cast<Real>(1),      P_in.data(),
                static_cast<Real>(add_to), C_in.data(),
                buffer_dim
            );
            
            ptoc(ClassName()+"::PrimitivesToClusters");
        }

        void ClustersToPrimitives( bool add_to = false ) const
        {
            ptic(ClassName()+"::ClustersToPrimitives");
            
            C_to_P.Dot(
                static_cast<Real>(1),      C_out.data(),
                static_cast<Real>(add_to), P_out.data(),
                buffer_dim
            );
            
            ptoc(ClassName()+"::ClustersToPrimitives");
        }

        
        void PercolateUp() const
        {
//            ptic(ClassName()+"::PercolateUp");
            RequireParallelPercolationRoots();
            
            switch (settings.tree_perc_alg)
            {
                case TreePercolationAlgorithm::Sequential:
                {
                    ptic(ClassName()+"::PercolateUp_DFS");
                    PercolateUp_DFS( 0 );
                    ptoc(ClassName()+"::PercolateUp_DFS");
                    break;
                }
                case TreePercolationAlgorithm::Recursive:
                {
                    PercolateUp_Recursive( 0 );
                    break;
                }
                case TreePercolationAlgorithm::Tasks:
                {
                    #pragma omp parallel num_threads( ThreadCount() )
                    {
                        #pragma omp single nowait
                        {
                            PercolateUp_Tasks( 0 );
                        }
                    }
                    break;
                }
                case TreePercolationAlgorithm::Parallel:
                {
                    PercolateUp_Parallel();
                    break;
                }
                default:
                {
                    PercolateUp_DFS( 0 );
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
                    PercolateDown_DFS( 0 );
                    ptoc(ClassName()+"::PercolateDown_DFS");
                    break;
                }
                case TreePercolationAlgorithm::Recursive:
                {
                    PercolateDown_Recursive( 0 );
                    break;
                }
                case TreePercolationAlgorithm::Tasks :
                {
                    #pragma omp parallel num_threads( ThreadCount() )
                    {
                        #pragma omp single nowait
                        {
                            PercolateDown_Tasks( 0 );
                        }
                    }
                    break;
                }
                case TreePercolationAlgorithm::Parallel:
                {
                    PercolateDown_Parallel();
                    break;
                }
                default:
                {
                    PercolateDown_DFS( 0 );
                }
            }
//            ptoc(ClassName()+"::PercolateDown");
        }; // PercolateUp
        

#include "Percolate_DFS.hpp"
#include "Percolate_Recursive.hpp"
#include "Percolate_Tasks.hpp"
#include "Percolate_Parallel.hpp"


        
    public:
        
        void Pre( const ExtReal * input, const Int cols, const OperatorType op_type ) const
        {
            ptic(ClassName()+"::Pre");
            if( pre_post_initialized && mixed_pre_post_initialized )
            {
                Sparse::MatrixCSR<Real,Int,Int> * pre;
                
                switch( op_type )
                {
                    case OperatorType::FractionalOnly:
                    {
                        pre  = &lo_pre ;
                        RequireBuffers( cols );
                        break;
                    }
                    case OperatorType::HighOrder:
                    {
                        pre  = &hi_pre ;
                        RequireBuffers( AmbDim() * cols ); // Beware: The derivative operator increases the number of columns!
                        break;
                    }
                    case OperatorType::LowOrder:
                    {
                        pre  = &lo_pre ;
                        RequireBuffers( cols );
                        break;
                    }
                    case OperatorType::MixedOrder:
                    {
                        pre  = &mixed_pre ;
                        RequireBuffers( (AmbDim()+1) *cols ); // Beware: The mixed preprocessor operator increases the number of columns!
                        break;
                    }
                    default:
                    {
                        eprint("Unknown kernel. Doing no.");
                        ptoc(ClassName()+"::Pre");
                        return;
                    }
                }

                // Caution: Some magic is going on here high order term...
                // Apply diff/averaging operate, reorder and multiply by weights.
                ptic(ClassName()+" pre->Dot");
                pre->Dot(
                    static_cast<Real>(1), input,
                    static_cast<Real>(0), P_in.data(),
                    cols
                );
                ptoc(ClassName()+" pre->Dot");
                
                // Accumulate into leaf clusters.
                PrimitivesToClusters(false);
                
                PercolateUp();
            }
            else
            {
                eprint(ClassName()+"::Pre: Preprocessing matrices are not initialized. Doing nothing.");
            }
            ptoc(ClassName()+"::Pre");
        }; // Pre


        void Post( ExtReal * output, const ExtReal alpha, const ExtReal beta, const OperatorType op_type ) const
        {
            ptic(ClassName()+"::Post");
            
            if( pre_post_initialized )
            {
                Sparse::MatrixCSR<Real,Int,Int> * post;
                
                switch( op_type )
                {
                    case OperatorType::FractionalOnly:
                    {
                        post  = &lo_post;
                        break;
                    }
                    case OperatorType::HighOrder:
                    {
                        post  = &hi_post;
                        break;
                    }
                    case OperatorType::LowOrder:
                    {
                        post  = &lo_post;
                        break;
                    }
                    case OperatorType::MixedOrder:
                    {
                        post  = &mixed_post;
                        break;
                    }
                    default:
                    {
                        eprint("Unknown kernel. Doing nothing.");
                        ptoc(ClassName()+"::Post");
                        return;
                    }
                }
              
                PercolateDown();
                
                ClustersToPrimitives( true );
                                
                // Multiply by weights, restore external ordering, and apply transpose of diff/averaging operator.
                
                ptic(ClassName()+" post->Dot");
                post->Dot(
                    alpha, P_out.data(),
                    beta,  output,
                    ( PrimitiveCount() * buffer_dim ) / post->ColCount()
                );
                ptoc(ClassName()+" post->Dot");
            }
            else
            {
                eprint(ClassName()+"::Post: Postprocessing matrices are not initialized. Doing nothing.");
            }
            ptoc(ClassName()+"::Post");
        }; // Post
        
        
//################################################################################################
//##        Neighbors
//################################################################################################
        
        
    public:
        
        const Sparse::BinaryMatrixCSR<Int,Int> & PrimitiveAdjacencyMatrix() const
        {
            if( !P_adjacency_matrix_initialized )
            {
                ptic(ClassName()+"::PrimitiveAdjacencyMatrix");
                // We compute product of S->lo_pre and  S->lo_post because it is basically the primitive-primitive adjacency matrix, we two primitives are supposed to be adjacencent if they share a common degree of freedom.

                P_adjacency_matrix = lo_pre.DotBinary(lo_post);
                
                P_adjacency_matrix_initialized = true;

                ptoc(ClassName()+"::PrimitiveAdjacencyMatrix");

            } // if( !P_adjacency_matrix_initialized )
            
            return P_adjacency_matrix;
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
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < cluster_count; ++i )
                {
                    to[i] += from[far_dim * i];
                }
            }
            
            //Before percolating the energies down we have to transform them to densities.
            {
                ptr<Real> a = C_far.data();
                mut<Real> e = C_out.data();
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < cluster_count; ++i )
                {
                    e[i] /= a[far_dim * i];
                }
            }
            
            PercolateDown();
            
            ClustersToPrimitives( false );
            
            // Now the simplex far field densities are stored in P_out.
            // We have to convert them back to energies before we add the energies from thread_P_D_near.
            {
                ptr<Real> a = P_near.data();
                mut<Real> e = P_out.data();
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < primitive_count; ++i )
                {
                    e[i] *= a[near_dim * i];
                }
                
            }
            
            // Add first entries of thread_P_D_near into P_out.
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                ptr<Real> from = thread_P_D_near.data(thread);
                mut<Real> to   = P_out.data();
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < primitive_count; ++i )
                {
                    to[i] += from[near_dim * i];
                }
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
            {
                ptr<Int>  o    = P_inverse_ordering.data();
                ptr<Real> from = P_out.data();
            
                #pragma omp parallel for num_threads( ThreadCount() ) schedule( static )
                for( Int i = 0; i < primitive_count; ++i )
                {
                    const Int j = o[i];
                    output[i] = static_cast<ExtReal>(from[j]);
                }
            }
            
            ptoc(ClassName()+"::CollectPrimitiveEnergies");
       
        } // CollectPrimitiveEnergies
            
        void CollectDensity( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const
        {
            ptic(ClassName()+"::CollectDensity");
            
            // Partially using the pipeline for derivatives. Not really efficient, but also not performance critical.
            
            const Int near_dim        = NearDim();
            const Int primitive_count = PrimitiveCount();
            
            this->RequireBuffers( static_cast<Int>(1) );
            
            // Compute dual volume vectors.
            fill_buffer( &P_in[0], primitive_count, static_cast<Real>(1) );

            const Int vertex_count = lo_post.RowCount();

            Tensor1<ExtReal,Int> dual_volumes ( vertex_count );
        
            lo_post.Dot(
                static_cast<Real>(1),    P_in.data(),
                static_cast<ExtReal>(0), dual_volumes.data(),
                1
            );
        
            
            //Collect the energies into P_out.
            CollectPrimitiveEnergies();
            
            // Divide by primitive volumes to get densities.
            {
                ptr<Real> a = P_near.data();
                mut<Real> e = P_out.data();
                
                #pragma omp parallel for num_threads(ThreadCount()) schedule( static )
                for( Int j = 0; j < primitive_count; ++j )
                {
                    e[j] /= a[near_dim * j];
                }
            }
            
            // Distribute the energy densities to energies per vertex. (Note that lo_post also multiplies by the primitives' volumes!)
            lo_post.Dot(
                static_cast<Real>(1),        P_out.data(),
                static_cast<ExtReal>(addTo), output,
                static_cast<Int>(1)
            );

            // Finally, we divide by the dual volumes to obtain the vertex densities.
            {
                ptr<ExtReal> a = dual_volumes.data();
                mut<ExtReal> e = output;
                
                #pragma omp parallel for num_threads(ThreadCount()) schedule( static )
                for( Int i = 0; i < vertex_count; ++i )
                {
                    e[i] /= a[i];
                }
            }
            
            ptoc(ClassName()+"::CollectDensity");
            
        } // CollectDensity

//##################################################################################################
//##        Updates
//##################################################################################################

        virtual void SemiStaticUpdate(
            const Tensor2<Real,Int> & P_near_, const Tensor2<Real,Int> & P_far_ ) const = 0;
        
        virtual void SemiStaticUpdate( ptr<Real> P_near_, ptr<Real> P_far_ ) const = 0;
        
//##################################################################################################
//##        Moments
//##################################################################################################
        
//        virtual void RequireClusterMoments( const Int moment_degree_ ) const = 0;
        
//##################################################################################################
//##        General reports
//##################################################################################################
        
        virtual std::string Stats() const = 0;
        
    private:
        std::string className() const
        {
            return std::string("ClusterTreeBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
    public:
        
        virtual std::string ClassName() const
        {
            return className();
        }

    }; // class ClusterTreeBase
    
} // namespace Repulsor
