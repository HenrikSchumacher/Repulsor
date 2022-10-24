//TODO: sequential computeClusterData?

#pragma once

#define CLASS ClusterTreeBase

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);
 
        
    public:
        
        using DataContainer_T   = Tensor2<Real,Int>;
        using BufferContainer_T = Tensor1<Real,Int>;
        
        // In principle, ThreadTensor3<Real,Int> should have better scaling on multiple socket machines, because I tried to encourages that the thread-local arrays are allocated on local RAM. -- On my tiny Quad Core however, it performs a bit _WORSE_ than Tensor3<Real,Int>.
        using DerivativeContainer_T = ThreadTensor3<Real,Int>;
        using Accumulator_T         = ThreadTensor3<Real,Int>;

        explicit CLASS( const ClusterTreeSettings & settings_ = ClusterTreeSettings() )
        :   settings(settings_)
        {}
        
//        CLASS( const CLASS & other )
//        :   settings ( other.settings)
//        ,   P_near   ( other.P_near)
//        ,   P_D_near   ( other.P_D_near)
//        ,   P_far   ( other.P_far)
//        ,   P_D_far   ( other.P_D_far)
//        ,   P_moments   ( other.P_moments)
//        ,   P_in   ( other.P_in)
//        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
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
        
        mutable SparseMatrixCSR<Real,Int>  hi_pre;
        mutable SparseMatrixCSR<Real,Int>  hi_post;

        mutable SparseMatrixCSR<Real,Int>  lo_pre;
        mutable SparseMatrixCSR<Real,Int>  lo_post;
        
        mutable SparseMatrixCSR<Real,Int>  mixed_pre;
        mutable SparseMatrixCSR<Real,Int>  mixed_post;

        mutable SparseBinaryMatrixCSR<Int> P_to_C;
        mutable SparseBinaryMatrixCSR<Int> C_to_P;
        
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
        mutable SparseBinaryMatrixCSR<Int> P_adjacency_matrix;
        
        // Some temproray shared data that is required for the parallel construction and serialization of the tree.
        
        Tensor3<SReal,Int> C_thread_serialized;          // False sharing is unlikely as each thread's slice should already be quite large...
        
        Tensor2<Int,Int> thread_cluster_counter;                          // TODO: Avoid false sharing!
        
    public:

//##################################################################################################
//##        Get functions
//##################################################################################################
        
        virtual Int AmbDim() const = 0;
        
        virtual Int NearDim() const = 0;
        
        virtual Int FarDim() const = 0;
        
        Int ThreadCount() const
        {
            return settings.thread_count;
        }
        
        Int SplitThreshold() const
        {
            return settings.split_threshold;
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
        
        const SparseMatrixCSR<Real,Int> & LowOrderPreprocessor() const
        {
            return lo_pre;
        }
        
        const SparseMatrixCSR<Real,Int> & LowOrderPostprocessor() const
        {
            return lo_post;
        }
        
        
        const SparseMatrixCSR<Real,Int> & HighOrderPreprocessor() const
        {
            return hi_pre;
        }
        
        const SparseMatrixCSR<Real,Int> & HighOrderPostprocessor() const
        {
            return hi_post;
        }
        
        const SparseMatrixCSR<Real,Int> & MixedOrderPreprocessor() const
        {
            return mixed_pre;
        }
        
        const SparseMatrixCSR<Real,Int> & MixedOrderPostprocessor() const
        {
            return mixed_post;
        }

        
        
        const SparseBinaryMatrixCSR<Int> & ClusterToPrimitiveMatrix() const
        {
            return C_to_P;
        }
        
        const SparseBinaryMatrixCSR<Int> & PrimitiveToClusterMatrix() const
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
            ptic(ClassName()+"::PercolateUp");
            switch (settings.tree_perc_alg)
            {
                case TreePercolationAlgorithm::Sequential:
                {
                    PercolateUp_Sequential( 0 );
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
                            PercolateUp_Tasks( 0, ThreadCount() );
                        }
                    }
                    break;
                }
                default:
                {
                    PercolateUp_Sequential( 0 );
                    break;
                }
            }
            ptoc(ClassName()+"::PercolateUp");
        }; // PercolateUp
        
        void PercolateDown() const
        {
            ptic(ClassName()+"::PercolateDown");
            switch (settings.tree_perc_alg) {
                case TreePercolationAlgorithm::Sequential:
                {
                    PercolateDown_Sequential( 0 );
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
                            PercolateDown_Tasks( 0, ThreadCount() );
                        }
                    }
                    break;
                }
                default:
                {
                    PercolateDown_Sequential( 0 );
                }
            }
            ptoc(ClassName()+"::PercolateDown");
        }; // PercolateUp
        
        
        
        
        
        
        void PercolateUp_Recursive( const Int C ) const
        {
            // C = cluster index
            
            const Int L = C_left [C];
            const Int R = C_right[C];
            
            if( (L >= 0) && (R >= 0) )
            {
                // If not a leaf, compute the values of the children first.
                PercolateUp_Recursive(L);
                PercolateUp_Recursive(R);
                
                // Aftwards, compute the sum of the two children.
                
                Real * restrict const c = C_in.data();
                #pragma omp simd aligned( c : ALIGNMENT )
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    // Overwrite, not add-into. Thus cleansing is not required.
                    c[ buffer_dim * C + k ] = c[ buffer_dim * L + k ] + c[ buffer_dim * R + k ];
                }
            }
            
        }; // PercolateUp_Recursive

        void PercolateDown_Recursive(const Int C)  const
        {
            // C = cluster index
            
            const Int L = C_left [C];
            const Int R = C_right[C];
            
            Real * restrict const c = C_out.data();
            
            if( ( L >= 0 ) && ( R >= 0 ) )
            {
                #pragma omp simd aligned( c : ALIGNMENT )
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    const Real buffer = c[ buffer_dim * C + k ];
                    c[ buffer_dim * L + k ] += buffer;
                    c[ buffer_dim * R + k ] += buffer;
                }
                PercolateDown_Recursive(L);
                PercolateDown_Recursive(R);
            }
        }; // PercolateDown_Recursive
        

        
        // Sequential variant that uses a stack instead of recursion.
        void PercolateUp_Sequential( const Int C_root ) const
        {
//            ptic("PercolateUp_Sequential");
            
            Int stack   [128] = {};
            Int visited [128] = {false};
            
            Int stack_ptr          = static_cast<Int>(0);
            stack[stack_ptr] = C_root;
            
            const Int  * restrict const left  = C_left.data();
            const Int  * restrict const right = C_right.data();
                  Real * restrict const c_in  = C_in.data();
            
            while( (stack_ptr >= 0) && ( stack_ptr < 126 ) )
            {
                // We are at cluster C.
                const Int C = stack[stack_ptr];
                const Int L = left [C];
                const Int R = right[C];
                
                if( !visited[stack_ptr] )
                {
                    visited[stack_ptr] = true;
                    
                    if( (L >= 0) && (R >= 0) )
                    {
                        // If not a leaf, compute the values of the children first.
                        stack[++stack_ptr] = R; // push
                        stack[++stack_ptr] = L; // push
                    }
                }
                else
                {
                    #pragma omp simd aligned( c_in : ALIGNMENT )
                    for( Int k = 0; k < buffer_dim; ++k )
                    {
                        // Overwrite, not add-into. Thus cleansing is not required.
                        c_in[ buffer_dim * C + k ] = c_in[ buffer_dim * L + k ] + c_in[ buffer_dim * R + k ];
                    }
                    visited[stack_ptr] = false;
                    stack_ptr--;  // pop
                }
            }
            
            if( stack_ptr >= 128 )
            {
                eprint(ClassName()+"::PercolateUp_Sequential: stack overflow detected.");
            }
//            ptoc("PercolateUp_Sequential");
        }; // PercolateUp_Sequential
        
        // Sequential variant that uses a stack instead of recursion.
        void PercolateDown_Sequential( const Int C_root ) const
        {
//            ptic("PercolateDown_Sequential");
            
            Int stack [128] = {};

            Int stack_ptr = static_cast<Int>(0);
            stack[stack_ptr] = C_root;
            
            const Int  * restrict const left  = C_left.data();
            const Int  * restrict const right = C_right.data();
                  Real * restrict const c_out = C_out.data();
            
            while( (stack_ptr >= 0) && ( stack_ptr < 126 ) )
            {
                // We are at cluster C.
                const Int C = stack[stack_ptr--]; //pop
                const Int L = left [C];
                const Int R = right[C];
                
                if( (L >= 0) && (R >= 0) )
                {
                    // If not a leaf, compute the values of the children first.
                    #pragma omp simd aligned( c_out : ALIGNMENT )
                    for( Int k = 0; k < buffer_dim; ++k )
                    {
                        const Real buffer = c_out[ buffer_dim * C + k ];
                        c_out[ buffer_dim * L + k ] += buffer;
                        c_out[ buffer_dim * R + k ] += buffer;
                    }
                    
                    stack[++stack_ptr] = R; // push
                    stack[++stack_ptr] = L; // push
                }
            }
            
            if( stack_ptr >= 128 )
            {
                eprint(ClassName()+"::PercolateDown_Sequential: stack overflow detected.");
            }
//            ptoc("PercolateDown_Sequential");
        }; // PercolateDown_Sequential
    
        
        void PercolateUp_Tasks( const Int C, const Int free_thread_count ) const
        {
            // C = cluster index
            
            const Int L = C_left [C];
            const Int R = C_right[C];
            
            if( (L >= 0) && (R >= 0) )
            {
                // If not a leaf, compute the values of the children first.
                #pragma omp task final(free_thread_count<1)  shared( L )
                    PercolateUp_Tasks( L, free_thread_count/2 );
                #pragma omp task final(free_thread_count<1)  shared( R )
                    PercolateUp_Tasks( R, free_thread_count-free_thread_count/2 );
                #pragma omp taskwait
                
                // Aftwards, compute the sum of the two children.
                
                Real * restrict const c = C_in.data();
                #pragma omp simd aligned( c : ALIGNMENT )
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    // Overwrite, not add-into. Thus cleansing is not needed.
                    c[ buffer_dim * C + k ] = c[ buffer_dim * L + k ] + c[ buffer_dim * R + k ];
                }
            }
            
        }; // PercolateUp_Tasks
        
        void PercolateDown_Tasks(const Int C, const Int free_thread_count ) const
        {
            // C = cluster index
            
            const Int L = C_left [C];
            const Int R = C_right[C];
            
            Real * restrict const  c = C_out.data();
            
            if( ( L >= 0 ) && ( R >= 0 ) )
            {
                #pragma omp simd aligned( c : ALIGNMENT )
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    const Real buffer = c[ buffer_dim * C + k ];
                    c[ buffer_dim * L + k ] += buffer;
                    c[ buffer_dim * R + k ] += buffer;
                }
                
                #pragma omp task final(free_thread_count<1)  shared( L )
                    PercolateDown_Tasks( L, free_thread_count/2 );
                #pragma omp task final(free_thread_count<1)  shared( R )
                    PercolateDown_Tasks( R, free_thread_count-free_thread_count/2 );
                #pragma omp taskwait
            }
        }; // PercolateDown_Tasks

        
        //        void PercolateDown_Parallel( const Int max_leaves ) const
        //        {
        //            ptic(ClassName()+"::PercolateDown_Parallel");
        //
        //            DUMP(max_leaves);
        //
        //            std::deque<Int> queue;
        //            queue.push_back(static_cast<Int>(0));
        //
        //            const Int  * restrict const left  = ClusterLeft().data();
        //            const Int  * restrict const right = ClusterRight().data();
        //                  Real * restrict const c     = C_out.data();
        //
        //            while( !queue.empty() && ( queue.size() < max_leaves ) )
        //            {
        //                const Int C = queue.front();
        //                              queue.pop_front();
        //
        //                const Int L = left [C];
        //                const Int R = right[C];
        //
        ////                if( ( L >= 0 ) && ( R >= 0 ) )
        //                if( L >= 0 )
        //                {
        //                    #pragma omp simd aligned( c : ALIGNMENT )
        //                    for( Int k = 0; k < buffer_dim; ++k )
        //                    {
        //                        const Real buffer = c[ buffer_dim * C + k ];
        //                        c[ buffer_dim * L + k ] += buffer;
        //                        c[ buffer_dim * R + k ] += buffer;
        //                    }
        //
        //                    queue.push_back(R);
        //                    queue.push_back(L);
        //                }
        //            }
        //
        ////            ptic("parallel");
        //            #pragma omp parallel for num_threads( ThreadCount() ) schedule( static )
        //            for( Int k = 0; k < queue.size(); ++k )
        //            {
        //                PercolateDown_Sequential(queue[k]);
        //            }
        ////            ptoc("parallel");
        //            ptoc(ClassName()+"::PercolateDown_Parallel");
        //        }
        
        
        
        
        
        void Pre( const ExtReal * input, const Int cols, const KernelType type ) const
        {
            ptic(ClassName()+"::Pre");
            if( pre_post_initialized && mixed_pre_post_initialized )
            {
                SparseMatrixCSR<Real,Int> * pre;
                
                switch (type)
                {
                    case KernelType::FractionalOnly:
                    {
                        pre  = &lo_pre ;
                        RequireBuffers( cols );
                        break;
                    }
                    case KernelType::HighOrder:
                    {
                        pre  = &hi_pre ;
                        RequireBuffers( AmbDim() * cols ); // Beware: The derivative operator increases the number of columns!
                        break;
                    }
                    case KernelType::LowOrder:
                    {
                        pre  = &lo_pre ;
                        RequireBuffers( cols );
                        break;
                    }
                    case KernelType::MixedOrder:
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


        void Post( ExtReal * output, const ExtReal alpha, const ExtReal beta, const KernelType type ) const
        {
            ptic(ClassName()+"::Post");
            
            if( pre_post_initialized )
            {
                SparseMatrixCSR<Real,Int> * post;
                
                switch (type)
                {
                    case KernelType::FractionalOnly:
                    {
                        post  = &lo_post;
                        break;
                    }
                    case KernelType::HighOrder:
                    {
                        post  = &hi_post;
                        break;
                    }
                    case KernelType::LowOrder:
                    {
                        post  = &lo_post;
                        break;
                    }
                    case KernelType::MixedOrder:
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
        
        const SparseBinaryMatrixCSR<Int> & PrimitiveAdjacencyMatrix() const
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
                const Real * restrict const from = thread_C_D_far.data(thread);
                      Real * restrict const to   = C_out.data();
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < cluster_count; ++i )
                {
                    to[i] += from[far_dim * i];
                }
            }
            
            //Before percolating the energies down we have to transform them to densities.
            {
                const Real * restrict const a = C_far.data();
                      Real * restrict const e = C_out.data();
                
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
                const Real * restrict const a = P_near.data();
                      Real * restrict const e = P_out.data();
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < primitive_count; ++i )
                {
                    e[i] *= a[near_dim * i];
                }
                
            }
            
            // Add first entries of thread_P_D_near into P_out.
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                const Real * restrict const from = thread_P_D_near.data(thread);
                      Real * restrict const to   = P_out.data();
                
                #pragma omp parallel for num_threads( thread_count ) schedule( static )
                for( Int i = 0; i < primitive_count; ++i )
                {
                    to[i] += from[near_dim * i];
                }
            }
                        
            ptoc(ClassName()+"::CollectPrimitiveEnergies");
            
        } // CollectPrimitiveEnergies
        
        void CollectPrimitiveEnergies( ExtReal * output, const ExtReal weight, bool addTo = false ) const
        {
            ptic(ClassName()+"::CollectPrimitiveEnergies");
            
            //Collect the energies into P_out.
            CollectPrimitiveEnergies();
            
            const Int primitive_count = PrimitiveCount();
            
            // Copy the values to output. We must not forget reorder!
            {
                const Int  * restrict const o    = P_inverse_ordering.data();
                const Real * restrict const from = P_out.data();
                      Real * restrict const to   = output;
            
                #pragma omp parallel for num_threads( ThreadCount() ) schedule( static )
                for( Int i = 0; i < primitive_count; ++i )
                {
                    const Int j = o[i];
                    to[i] = from[j];
                }
            }
            
            ptoc(ClassName()+"::CollectPrimitiveEnergies");
       
        } // CollectPrimitiveEnergies
            
        void CollectDensity( ExtReal * output, const ExtReal weight, bool addTo = false ) const
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
                const Real * restrict const a = P_near.data();
                      Real * restrict const e = P_out.data();
                
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
                const Real * restrict const a = dual_volumes.data();
                      Real * restrict const e = output;
                
                #pragma omp parallel for simd num_threads(ThreadCount()) aligned( a, e : ALIGNMENT ) schedule( static )
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
            const Tensor2<Real,Int> & P_near_, const Tensor2<Real,Int> & P_far_ ) = 0;
        
        virtual void SemiStaticUpdate(
            const Real * restrict const P_near_, const Real * restrict const P_far_ ) = 0;
        
//##################################################################################################
//##        Moments
//##################################################################################################
        
//        virtual void RequireClusterMoments( const Int moment_degree_ ) const = 0;
        
//##################################################################################################
//##        General reports
//##################################################################################################
        
        virtual std::string Stats() const = 0;
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS) + "<"+TypeName<ExtReal>::Get()+">";
        }

    };
    
} // namespace Repulsor

#undef CLASS
