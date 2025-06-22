#pragma once

#include "ClusterTree/Cluster.hpp"
#include "ClusterTree/ClusterTreeBase.hpp"

namespace Repulsor
{
    template<int AMB_DIM_, typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_>
    class ClusterTree : public ClusterTreeBase<Real_,Int_,LInt_,SReal_,ExtReal_>
    {
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = LInt_;
        
        using Base_T = ClusterTreeBase<Real,Int,LInt,SReal,ExtReal>;
        
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
        
        using Base_T::C_to_P;
        using Base_T::ClusterToPrimitiveMatrix;
        using Base_T::P_to_C;
        using Base_T::PrimitiveToClusterMatrix;
        
        using Base_T::hi_pre;
        using Base_T::lo_pre;
        using Base_T::mi_pre;
        using Base_T::hi_post;
        using Base_T::lo_post;
        using Base_T::mi_post;
        
    public:

        
        // To allow polymorphism, we require the user to create instances of the desired types for the primitives and the bounding volumes, so that we can Clone() them.
        ClusterTree(
            cref<Primitive_T>         P_proto_,
            cref<Tensor2<SReal,Int>>  P_serialized_,
            cref<BoundingVolume_T>    C_proto_,
            cref<Tensor1<Int ,Int>>   P_ordering_,
            cref<Tensor2<Real,Int>>   P_near_, // data used actual interaction computation; assumed to be of size PrimitiveCount() x NearDim(). For a triangle mesh in 3D, we want to feed each triangles i), area ii) barycenter and iii) normal as a 1 + 3 + 3 = 7 vector
            cref<Tensor2<Real,Int>>   P_far_, // data used actual interaction computation; assumed to be of size PrimitiveCount() x FarDim(). For a triangle mesh in 3D, we want to feed each triangles i), area ii) barycenter and iii) orthoprojector onto normal space as a 1 + 3 + 6 = 10 vector
            cref<SparseMatrix_T>      DiffOp,
            cref<SparseMatrix_T>      AvOp,
            cref<ClusterTreeSettings> settings_ = ClusterTreeSettings()
        )
        :   Base_T( settings_ )
        ,   P_proto      ( ToSize_T(ThreadCount()) )
        ,   C_proto      ( ToSize_T(ThreadCount()) )
        ,   P_moving     ( ToSize_T(ThreadCount()) )
        {
            TOOLS_PTIC(className()+"()");
            
            P_serialized = P_serialized_;   // It's a unneccessary copy, but not too expensive.
            
            if( P_ordering_.Dim(0) == P_serialized.Dim(0) )
            {
                P_ordering = P_ordering_;
            }
            else
            {
                P_ordering = iota<Int,Int>( P_serialized.Dim(0) );
            }
            
            ParallelDo(
                [&,this]( const Size_T thread )
                {
                    P_proto[thread] = P_proto_.Clone();
                    
                    C_proto[thread] = C_proto_.Clone();
                },
                ToSize_T(ThreadCount())
            );
            
            this->ComputeClusters();
            
            if( P_near_.Dim(0) != PrimitiveCount() )
            {
                eprint(className()+" : P_near_.Dim(0) != PrimitiveCount()");
                return;
            }
            if( P_far_.Dim(0) != PrimitiveCount() )
            {
                eprint(className()+" : P_far_.Dim(0) != PrimitiveCount()");
                return;
            }
            
            this->AllocateNearFarData( P_near_.Dim(1), P_far_.Dim(1) );
            
            this->ComputePrimitiveData( P_near_.data(), P_far_.data() );
            
            this->ComputeClusterData();
            
            this->ComputePrimitiveToClusterMatrix();
            
            this->ComputePrePost( DiffOp, AvOp );
            
            this->RequireMixedPrePost();
            
//            this->ComputeMixedPrePost( DiffOp, AvOp );
            
            TOOLS_PTOC(className()+"()");
        }
    
        // Default constructor
        ClusterTree()
        :   Base_T()
        ,   P_proto (1)
        ,   C_proto (1)
        ,   P_moving(1)
        {
            TOOLS_PTIMER(timer,className()+" default constructor");
        }
        // Destructor
        virtual ~ClusterTree() override = default;
        // Copy constructor
        ClusterTree( const ClusterTree & other ) = default;
        // Copy assignment operator
        ClusterTree & operator=( const ClusterTree & other ) = default;
        // Move constructor
        ClusterTree( ClusterTree && other ) = default;
        // Move assignment operator
        ClusterTree & operator=( ClusterTree && other ) = default;
        
    protected:
        
        
        // Each thread gets its own primitive prototype to avoid sharing conflicts.
        std::vector<std::shared_ptr<Primitive_T>> P_proto;
        // False sharing prevented by alignment of PrimitiveBase.
        
        // Each thread gets its own bounding volume prototype to avoid sharing conflicts.
        std::vector<std::shared_ptr<BoundingVolume_T>> C_proto;
        std::vector<std::shared_ptr<BoundingVolume_T>> C_proto_updated;
        // False sharing prevented by alignment of PrimitiveBase.
        
        mutable std::vector<std::shared_ptr<MovingPrimitive_T>> P_moving;
        
        // TODO: This is very awkward. Is there a better way?
        std::vector<std::vector<Cluster_T *>> tree_rows_ptr;
        
    public:
        
        //        mutable std::vector<std::shared_ptr<MultipoleMomentsBase<Real,Int>>> M_ker;
        //        mutable Int moment_degree = 0;
        
#include "ClusterTree/Split_Thread.hpp"
#include "ClusterTree/Serialize_Thread.hpp"
#include "ClusterTree/Percolate_DFS.hpp"
#include "ClusterTree/Percolate_Parallel.hpp"
#include "ClusterTree/Percolate_Recursive.hpp"
#include "ClusterTree/ComputeClusterData_Thread.hpp"
        
#include "ClusterTree/DistanceQueries.hpp"

        
    private:
        
        void ComputeClusters()
        {
            TOOLS_PTIC(className()+"::ComputeClusters");
            
            // Request some temporary memory for threads.
            
            P_inverse_ordering = Tensor1<Int,Int>( PrimitiveCount() );
            
            // Padding every row to prevent false sharing.
            thread_cluster_counter = Tensor2<Int,Int>( ThreadCount(), CacheLineWidth, 0 );
            
            C_thread_serialized = Tensor3<SReal,Int>( ThreadCount(), 2*PrimitiveCount(), C_proto[0]->Size() );
            
            const Int thread = 0;
            
            ++thread_cluster_counter(thread,0);
            
            auto * root = new Cluster_T( thread, 0, 0, PrimitiveCount(), 0 );
            
            
            TOOLS_PTIC(className()+"::ComputeClusters: Initial bounding volume of root node");
            
            C_proto[thread]->SetPointer( C_thread_serialized.data(thread), 0 );
            C_proto[thread]->FromPrimitives( *P_proto[thread], P_serialized.data(), 0, PrimitiveCount(), ThreadCount() );
            
            TOOLS_PTOC(className()+"::ComputeClusters: Initial bounding volume of root node");
            
            
            Split( root );
            
            Serialize( root );
            
            delete root;
            
            // Free memory for threads.
            C_thread_serialized = Tensor3<SReal,Int>();
            
            thread_cluster_counter = Tensor2<Int,Int>();
            
            TOOLS_PTOC(className()+"::ComputeClusters");
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
        
        const MovingPrimitive_T & MovingPrimitivePrototype() const
        {
            return *P_moving[0];
        }
        
        Int NearDim() const override
        {
            return P_near.Dim(1);
        }
        
        Int FarDim() const override
        {
            return FAR_DIM;
        }
        
    private:
            
        void AllocateNearFarData( const Int near_dim_, const Int far_dim_ ) // reordering and computing bounding boxes
        {
            TOOLS_PTIC(className()+"::AllocateNearFarData");
            
            P_near   = DataContainer_T( PrimitiveCount(), near_dim_ );
            P_D_near = DataContainer_T( PrimitiveCount(), near_dim_ );
            P_far    = DataContainer_T( PrimitiveCount(), far_dim_  );
            P_D_far  = DataContainer_T( PrimitiveCount(), far_dim_  );
            
            thread_P_D_near = DerivativeContainer_T( ThreadCount(), PrimitiveCount(), near_dim_ );
            
            C_far           = DataContainer_T( ClusterCount(), far_dim_ );
            thread_C_D_far  = DerivativeContainer_T( ThreadCount(), ClusterCount(), far_dim_ );
            
            TOOLS_PTOC(className()+"::AllocateNearFarData");
        }
        
        void ComputePrimitiveData( cptr<Real> P_near_, cptr<Real> P_far_ ) const
        {
            TOOLS_PTIC(className()+"::ComputePrimitiveData");
            
            const Int near_dim = NearDim();
            const Int  far_dim =  FarDim();

            ParallelDo(
                [=,this]( const Int i )
                {
                    const Int j = P_ordering[i];

                    copy_buffer<VarSize,Sequential>( &P_near_[near_dim * j], P_near.data(i), near_dim );

                    copy_buffer<VarSize,Sequential>( &P_far_ [ far_dim * j],  P_far.data(i), far_dim  );
                },
                PrimitiveCount(), ThreadCount()
            );
            
            TOOLS_PTOC(className()+"::ComputePrimitiveData");
        } //ComputePrimitiveData
        
    public:
      
    protected:
        
        void ComputePrimitiveToClusterMatrix()
        {
            TOOLS_PTIC(className()+"::ComputePrimitiveToClusterMatrix");
            
            P_to_C = SparseBinaryMatrix_T(
                ClusterCount(), PrimitiveCount(), PrimitiveCount(), ThreadCount() );
            
            C_to_P = SparseBinaryMatrix_T(
                PrimitiveCount(), ClusterCount(), PrimitiveCount(), ThreadCount() );
            
            C_to_P.Outer()[PrimitiveCount()] = PrimitiveCount();
            
            {
                mptr<Int> inner_ = C_to_P.Inner().data();
                
                ParallelDo(
                    [=,this]( const Int i )
                    {
                        const Int leaf  = leaf_clusters[i];
                        const Int begin = C_begin[leaf];
                        const Int end   = C_end  [leaf];

                        for( Int k = begin; k < end; ++k )
                        {
                            inner_[k] = leaf;
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
                
                mptr<LInt> outer_ = P_to_C.Outer().data();
                cptr< Int> left_  = C_left.data();
                cptr< Int> begin_ = C_begin.data();
                cptr< Int> end_   = C_end.data();
                
                for ( Int C = 0; C < cluster_count; ++C )
                {
                    outer_[C+1] =
                        outer_[C]
                        +
                        static_cast<LInt>(left_[C] < 0)
                        *
                        static_cast<LInt>(end_[C] - begin_[C]);
                }
            }
            
            TOOLS_PTOC(className()+"::ComputePrimitiveToClusterMatrix");
        }
        

        
        void ComputePrePost(
            cref<SparseMatrix_T> DiffOp,
            cref<SparseMatrix_T> AvOp
        )
        {
            if( !this->pre_post_initialized )
            {
                TOOLS_PTIC(className()+"::ComputePrePost");

                
                const Int primitive_count = PrimitiveCount();
                
//                TOOLS_PTIC("hi_pre");
                
                hi_pre = SparseMatrix_T(
                    DiffOp.RowCount(),
                    DiffOp.ColCount(),
                    DiffOp.NonzeroCount(),
                    ThreadCount()
                );
                
                //Assuming that all rows have the same number of entries!
                {
                    
                    hi_pre.Outer().Read( DiffOp.Outer().data() );
                
                    cptr<Int>  ord      = P_ordering.data();
                    cptr<Int>  i_input  = DiffOp.Inner().data();
                    mptr<Int>  i_output = hi_pre.Inner().data();
                    cptr<Real> r_input  = DiffOp.Values().data();
                    mptr<Real> r_output = hi_pre.Values().data();
                    cptr<Real> near     = P_near.data();
                
                    const Int near_dim = NearDim();
                    const Int size = static_cast<Int>(DiffOp.NonzeroCount() / primitive_count);
                    
                    ParallelDo(
                        [=]( const Int i )
                        {
                            const Int from = size * ord[i];
                            const Int to   = size * i;
                            const Real a   = near[ near_dim * i ];

                            for( Int k = 0; k < size; ++k )
                            {
                                i_output[to+k] =     i_input[from+k];
                                r_output[to+k] = a * r_input[from+k];
                            }
                        },
                        PrimitiveCount(), ThreadCount()
                    );
                }
//                TOOLS_PTOC("hi_pre");
                
//                TOOLS_PTIC("hi_post");
                TOOLS_PDUMP(hi_pre.Stats());
                hi_post = hi_pre.Transpose();
//                TOOLS_PTOC("hi_post");
                
//                TOOLS_PTIC("lo_pre");
                
                lo_pre = SparseMatrix_T(
                    AvOp.RowCount(),
                    AvOp.ColCount(),
                    AvOp.NonzeroCount(),
                    ThreadCount()
                );
                
                //Assuming that all rows have the same number of entries!
                {
                    lo_pre.Outer().Read( AvOp.Outer().data() );
                                        
                    cptr<Int>  ord      = P_ordering.data();
                    cptr<Int>  i_input  =   AvOp.Inner().data();
                    mptr<Int>  i_output = lo_pre.Inner().data();
                    cptr<Real> r_input  =   AvOp.Values().data();
                    mptr<Real> r_output = lo_pre.Values().data();
                    cptr<Real> near     = P_near.data();
                    
                    const Int near_dim = NearDim();
                    const Int size = static_cast<Int>(AvOp.NonzeroCount() / primitive_count);
                    
                    ParallelDo(
                        [=]( const Int i )
                        {
                            const Int from = size * ord[i];
                            const Int to   = size * i;
                            const Real a   = near[ near_dim * i ];
    //
                            for( Int k = 0; k < size; ++k )
                            {
                                i_output[to+k] =     i_input[from+k];
                                r_output[to+k] = a * r_input[from+k];
                            }
                        },
                        PrimitiveCount(), ThreadCount()
                    );
                }
//                TOOLS_PTOC("lo_pre");
                
//                TOOLS_PTIC("lo_post");
                TOOLS_PDUMP(lo_pre.Stats());
                lo_post = lo_pre.Transpose();
//                TOOLS_PTOC("lo_post");
                
                this->pre_post_initialized = true;
                
                TOOLS_PTOC(className()+"::ComputePrePost");
            }
            
        } // ComputePrePost
        
        void RequireMixedPrePost()
        {
            if( !this->mi_pre_post_initialized )
            {
                std::string tag = className()+"::RequireMixedPrePost";
                
                TOOLS_PTIC(tag);
                
                const Int primitive_count = PrimitiveCount();
                
//                TOOLS_PTIC(tag + ": mi_pre");
                
                mi_pre = SparseMatrix_T(
                    lo_pre.RowCount()     + hi_pre.RowCount(),
                    lo_pre.ColCount(),
                    lo_pre.NonzeroCount() + hi_pre.NonzeroCount(),
                    ThreadCount()
                );
                
                {
                    // We assume that each row has the same number of nonzeros.
                    TOOLS_DEBUG_ASSERT(
                        lo_pre.Outer(1) - lo_pre.Outer(0) == hi_pre.Outer(1) - hi_pre.Outer(0),
                        ""
                    );
                    
                    const Int row_size = static_cast<Int>( lo_pre.Outer(1) - lo_pre.Outer(0) );
                    
    //                cptr<LInt> lo_outer  = lo_pre.Outer().data();
    //                cptr<LInt> hi_outer  = hi_pre.Outer().data();
                    mptr<LInt> mi_outer  = mi_pre.Outer().data();

                    cptr< Int> lo_inner  = lo_pre.Inner().data();
                    cptr< Int> hi_inner  = hi_pre.Inner().data();
                    mptr< Int> mi_inner  = mi_pre.Inner().data();

                    cptr<Real> lo_values = lo_pre.Values().data();
                    cptr<Real> hi_values = hi_pre.Values().data();
                    mptr<Real> mi_values = mi_pre.Values().data();
                    
                    mi_outer[mi_pre.RowCount()] = mi_pre.RowCount() * row_size;

                    ParallelDo(
                        [=]( const Int i )
                        {
                            const Int row_base = (AMB_DIM+1) * i;
                            
                            {
                                const Int row = row_base;
                                
                                const Int rp  = mi_outer[row] = row * row_size;
                                
                                const Int j   = row_size * i;
                                
                                copy_buffer<VarSize,Sequential>( &lo_inner [j], &mi_inner [rp], row_size );
                                copy_buffer<VarSize,Sequential>( &lo_values[j], &mi_values[rp], row_size );
                                
                            }
                            
                            for( Int k = 1; k < AMB_DIM+1; ++k )
                            {
                                const Int row = row_base + k;
                                
                                const Int rp  = mi_outer[row] = row * row_size;
                                
                                const Int j   = row_size * (AMB_DIM * i + k - 1);
                                
                                copy_buffer<VarSize,Sequential>( &hi_inner [j], &mi_inner [rp], row_size );
                                copy_buffer<VarSize,Sequential>( &hi_values[j], &mi_values[rp], row_size );
                            }
                        },
                        primitive_count, ThreadCount()
                    );
                }
//                TOOLS_PTOC(tag + ": mi_pre");
                
                TOOLS_PDUMP(mi_pre.Stats());
                
//                TOOLS_PTIC(tag + ": mi_post");
                
                mi_post = SparseMatrix_T(
                    lo_post.RowCount(),
                    lo_post.ColCount()     + hi_post.ColCount(),
                    lo_post.NonzeroCount() + hi_post.NonzeroCount(),
                    ThreadCount()
                );
                
                {
                    cptr<LInt> lo_outer  = lo_post.Outer().data();
                    cptr<LInt> hi_outer  = hi_post.Outer().data();
                    mptr<LInt> mi_outer  = mi_post.Outer().data();

                    cptr< Int> lo_inner  = lo_post.Inner().data();
//                    cptr< Int> hi_inner  = hi_post.Inner().data();
                    mptr< Int> mi_inner  = mi_post.Inner().data();

                    cptr<Real> lo_values = lo_post.Values().data();
                    cptr<Real> hi_values = hi_post.Values().data();
                    mptr<Real> mi_values = mi_post.Values().data();
                    
                    ParallelDo(
                        [=]( const Int i )
                        {
                            const LInt lo_begin = lo_outer[i  ];
                            const LInt lo_end   = lo_outer[i+1];
//                            const LInt hi_begin = hi_outer[i]  ;
                            const LInt hi_end   = hi_outer[i+1];
//                            const LInt mi_begin = lo_begin + hi_begin;
                            const LInt mi_end   = lo_end   + hi_end;
                            
                            mi_outer[i+1] = mi_end;
                            
                            for( LInt k = lo_begin; k < lo_end; ++k )
                            {
                                
                                const LInt offset = (AMB_DIM + 1) * k;

                                const Int c = (AMB_DIM + 1) * lo_inner [k];
                                
                                for( Int l = 0; l < (AMB_DIM + 1); ++l )
                                {
                                    mi_inner [offset + l] = c + l;
                                }
                                
                                mi_values[offset] = lo_values[k];
                                copy_buffer<AMB_DIM>(
                                    &hi_values[AMB_DIM * k], &mi_values[offset + 1]
                                );
                            }
                        },
                        lo_post.RowCount(), ThreadCount()
                    );
                }
                
//                TOOLS_PTOC(tag + ": mi_post");

                TOOLS_PDUMP(mi_post.Stats());
                
                this->mi_pre_post_initialized = true;
                
                TOOLS_PTOC(tag);
            }
        }
        
    public:
        
        void PrimitivesToClusters( bool add_to = false ) const override
        {
            TOOLS_PTIC(ClassName()+"::PrimitivesToClusters");
            
            switch( buffer_dim )
            {
                case 1:
                {
                    PrimitiveToClusterMatrix().template Dot<1>(
                        Scalar::One<Real>,         P_in.data(),
                        static_cast<Real>(add_to), C_in.data(),
                        1
                    );
                    break;
                }
                case AMB_DIM:
                {
                    PrimitiveToClusterMatrix().template Dot<AMB_DIM>(
                        Scalar::One<Real>,         P_in.data(),
                        static_cast<Real>(add_to), C_in.data(),
                        AMB_DIM
                    );
                    break;
                }
                case AMB_DIM * AMB_DIM:
                {
                    PrimitiveToClusterMatrix().template Dot<AMB_DIM * AMB_DIM>(
                        Scalar::One<Real>,         P_in.data(),
                        static_cast<Real>(add_to), C_in.data(),
                        AMB_DIM * AMB_DIM
                    );
                    break;
                }
                case (AMB_DIM + 1) * AMB_DIM:
                {
                    PrimitiveToClusterMatrix().template Dot<(AMB_DIM + 1) * AMB_DIM>(
                        Scalar::One<Real>,         P_in.data(),
                        static_cast<Real>(add_to), C_in.data(),
                        (AMB_DIM + 1) * AMB_DIM
                    );
                    break;
                }
                case FAR_DIM:
                {
                    ClusterToPrimitiveMatrix().template Dot<FAR_DIM>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        FAR_DIM
                    );
                    break;
                }
                default:
                {
                    PrimitiveToClusterMatrix().template Dot<0>(
                        Scalar::One<Real>,         P_in.data(),
                        static_cast<Real>(add_to), C_in.data(),
                        buffer_dim
                    );
                    break;
                }
            }
            
            TOOLS_PTOC(ClassName()+"::PrimitivesToClusters");
        }

        void ClustersToPrimitives( bool add_to = false ) const override
        {
            TOOLS_PTIC(ClassName()+"::ClustersToPrimitives");
            
            switch ( buffer_dim )
            {
                case 1:
                {
                    ClusterToPrimitiveMatrix().template Dot<1>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        1
                    );
                    break;
                }
                case AMB_DIM:
                {
                    ClusterToPrimitiveMatrix().template Dot<AMB_DIM>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        AMB_DIM
                    );
                    break;
                }
                case AMB_DIM * AMB_DIM:
                {
                    ClusterToPrimitiveMatrix().template Dot<AMB_DIM * AMB_DIM>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        AMB_DIM * AMB_DIM
                    );
                    break;
                }
                case (AMB_DIM + 1) * AMB_DIM:
                {
                    ClusterToPrimitiveMatrix().template Dot<(AMB_DIM + 1) * AMB_DIM>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        (AMB_DIM + 1) * AMB_DIM
                    );
                    break;
                }
                case FAR_DIM:
                {
                    ClusterToPrimitiveMatrix().template Dot<FAR_DIM>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        FAR_DIM
                    );
                    break;
                }
                default:
                {
                    ClusterToPrimitiveMatrix().template Dot<0>(
                        Scalar::One<Real>,         C_out.data(),
                        static_cast<Real>(add_to), P_out.data(),
                        buffer_dim
                    );
                    break;
                }
            }

            TOOLS_PTOC(ClassName()+"::ClustersToPrimitives");
        }
        
        
        void Pre( 
            const ExtReal * X, const Int ldX,
            const Int nrhs, const OperatorType op_type
        ) const override
        {
            std::string tag = ClassName();
            SparseMatrix_T * pre;
            
            switch( op_type )
            {
                case OperatorType::FractionalOnly:
                {
                    tag += "(lo)";
                    pre  = &lo_pre ;
                    this->RequireBuffers( nrhs );
                    break;
                }
                case OperatorType::HighOrder:
                {
                    tag += "(hi)";
                    pre  = &hi_pre ;
                    this->RequireBuffers( AmbDim() * nrhs ); // Beware: The derivative operator increases the number of columns!
                    break;
                }
                case OperatorType::LowOrder:
                {
                    tag += "(lo)";
                    pre  = &lo_pre ;
                    this->RequireBuffers( nrhs );
                    break;
                }
                case OperatorType::MixedOrder:
                {
                    tag += "(mi)";
                    pre  = &mi_pre ;
                    this->RequireBuffers( (AmbDim()+1) * nrhs ); // Beware: The mixed preprocessor operator increases the number of columns!
                    break;
                }
                default:
                {
                    eprint("Unknown kernel. Doing no.");
                    TOOLS_PTOC(ClassName()+"::Pre");
                    return;
                }
            }
            
            TOOLS_PTIC(tag);

            // Caution: Some magic is going on here high order term...
            // Apply diff/averaging operate, reorder and multiply by weights.
            
            switch( nrhs )
            {
                case 1:
                {
                    pre->template Dot<1>(
                        Scalar::One <Real>, X,           ldX,
                        Scalar::Zero<Real>, P_in.data(), 1,
                        1
                    );
                    break;
                }
                case AMB_DIM:
                {
                    pre->template Dot<AMB_DIM>(
                        Scalar::One <Real>, X,           ldX,
                        Scalar::Zero<Real>, P_in.data(), AMB_DIM,
                        AMB_DIM
                    );
                    break;
                }
                default:
                {
                    pre->template Dot<VarSize>(
                        Scalar::One <Real>, X,           ldX,
                        Scalar::Zero<Real>, P_in.data(), nrhs,
                        nrhs
                    );
                    break;
                }
            }

            
            // Accumulate into leaf clusters.
            PrimitivesToClusters(false);
            
            this->PercolateUp();
            
            TOOLS_PTOC(tag);
            
        }; // Pre


        void Post( 
            const ExtReal alpha,
            const ExtReal beta, ExtReal * Y, const Int ldY,  
            const OperatorType op_type
        ) const override
        {
            SparseMatrix_T * post;
            
            this->PercolateDown();
              
            ClustersToPrimitives( true );
            
            std::string tag = ClassName();
            
            switch( op_type )
            {
                case OperatorType::FractionalOnly:
                {
                    tag += "(lo)";
                    post  = &lo_post;
                    break;
                }
                case OperatorType::HighOrder:
                {
                    tag += "(hi)";
                    post  = &hi_post;
                    break;
                }
                case OperatorType::LowOrder:
                {
                    tag += "(lo)";
                    post  = &lo_post;
                    break;
                }
                case OperatorType::MixedOrder:
                {
                    tag += "(mi)";
                    post  = &mi_post;
                    break;
                }
                default:
                {
                    eprint("Unknown kernel. Doing nothing.");
                    return;
                }
            }
                         
            TOOLS_PTIC(tag);
            
            // Multiply by weights, restore external ordering, and apply transpose of diff/averaging operator.
            
            const Int nrhs = ( PrimitiveCount() * buffer_dim ) / post->ColCount();
            
            switch( nrhs )
            {
                case 1:
                {
                    post->template Dot<1>(
                        alpha, P_out.data(), 1,
                        beta,  Y,            ldY,
                        1
                    );
                    break;
                }
                case AMB_DIM:
                {
                    post->template Dot<AMB_DIM>(
                        alpha, P_out.data(), AMB_DIM,
                        beta,  Y,            ldY,
                        AMB_DIM
                    );
                    break;
                }
                default:
                {
                    post->template Dot<VarSize>(
                        alpha, P_out.data(), nrhs,
                        beta,  Y,            ldY,
                        nrhs
                    );
                }
            }
            
            TOOLS_PTOC(tag);
            
        }; // Post
        
    public:
        
        // Updates only the computational data (primitive/cluster areas, centers of mass and normals).
        // All data related to clustering or multipole acceptance criteria remain unchanged, as well
        // as the preprocessor and postprocessor matrices (that are needed for matrix-vector multiplies of the BCT.)
        
        void SemiStaticUpdate(
            cref<Tensor2<Real,Int>> P_near_,
            cref<Tensor2<Real,Int>> P_far_
        ) const override
        {
            if( P_near_.Dim(0) != PrimitiveCount() )
            {
                eprint(className()+"::SemiStaticUpdate: P_near_.Dim(0) = "+ToString(P_near_.Dim(0))+" != "+ToString(PrimitiveCount())+" = PrimitiveCount(). Skipping update.");
                return;
            }
            if( P_near_.Dim(1) != NearDim() )
            {
                eprint(className()+"::SemiStaticUpdate: P_near_.Dim(1) "+ToString(P_near_.Dim(1))+" != "+ToString(NearDim())+" = NearDim(). Skipping update.");
                return;
            }
            
            if( P_far_.Dim(0) != PrimitiveCount() )
            {
                eprint(className()+"::SemiStaticUpdate: P_far_.Dim(0) "+ToString(P_far_.Dim(0))+" != "+ToString(PrimitiveCount())+" = PrimitiveCount(). Skipping update.");
                return;
            }
            
            if( P_far_.Dim(1) != FarDim() )
            {
                eprint(className()+"::SemiStaticUpdate: P_far_.Dim(1) "+ToString(P_far_.Dim(1))+" != "+ToString(FarDim())+" = FarDim(). Skipping update.");
                return;
            }
            
            SemiStaticUpdate( P_near_.data(), P_far_.data() );
            
        } // SemiStaticUpdate
        
        void SemiStaticUpdate( cptr<Real> P_near_, cptr<Real> P_far_ ) const override
        {
            // Updates only the computational data like primitive/cluster areas, centers of mass and normals. All data related to clustering or multipole acceptance criteria remain are unchanged.
            
            TOOLS_PTIC(className()+"::SemiStaticUpdate");
            
            this->ComputePrimitiveData( P_near_, P_far_ );
            
            this->ComputeClusterData();

            TOOLS_PTOC(className()+"::SemiStaticUpdate");
            
        } // SemiStaticUpdate
        
        void TakeUpdateVectors(
            mref<MovingPrimitive_T> P_moving_,
            mref<Tensor2<SReal,Int>> P_velocities_serialized_,
            const SReal max_time
        ) const
        {
            TOOLS_PTIC(className()+"::TakeUpdateVectors");
            
            this->update_time = max_time;
            // ATTENTION: We assume here that P_velocities_serialized are handed over already in correct order and won't be used anymore and just destroyed by the caller. This allows us to just swap the Tensor2s here, saving us a big copy operation.
            swap( P_velocities_serialized, P_velocities_serialized_ );
            
            // TODO: Potentially superfluous container.
            P_updated_serialized = Tensor2<SReal,Int>( PrimitiveCount(), P_proto[0]->Size() );
            C_updated_serialized = Tensor2<SReal,Int>( ClusterCount(),   C_proto[0]->Size() );
            
            
            cptr<SReal> P_p_ser = P_serialized.data();
            cptr<SReal> P_v_ser = P_velocities_serialized.data();
            
            mptr<SReal> P_up_ser = P_updated_serialized.data();
            mptr<SReal> C_up_ser = C_updated_serialized.data();
            
            TOOLS_PTIC(className()+"::TakeUpdateVectors - Compute the primitives updated with max_time.");
            
            // TODO: Potentially wasteful code.
            ParallelDo(
                [=,this,&P_moving_]( const Int thread )
                {
                    const Int i_begin = JobPointer<Int>(PrimitiveCount(), ThreadCount(), thread     );
                    const Int i_end   = JobPointer<Int>(PrimitiveCount(), ThreadCount(), thread + 1 );
                    
                    P_moving[ToSize_T(thread)] = P_moving_.Clone();
                    MovingPrimitive_T & P_mov = *P_moving[ToSize_T(thread)];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        P_mov.ReadCoordinatesSerialized ( P_p_ser           , i );
                        P_mov.ReadVelocitiesSerialized  ( P_v_ser           , i );
                        P_mov.WriteDeformedSerialized   ( P_up_ser, max_time, i );
                    }
                },
                ThreadCount()
            );
            
            TOOLS_PTOC(className()+"::TakeUpdateVectors - Compute the primitives updated with max_time.");
            
            TOOLS_PTIC(className()+"::TakeUpdateVectors - Compute the AABBs of the updated leaf clusters.");
            
            ParallelDo(
                [=,this]( const Int thread )
                {
                    const Int i_begin = JobPointer<Int>(LeafClusterCount(), ThreadCount(), thread     );
                    const Int i_end   = JobPointer<Int>(LeafClusterCount(), ThreadCount(), thread + 1 );
                    
                    mref<Primitive_T>      P_   = *P_proto[ToSize_T(thread)];
                    mref<BoundingVolume_T> C_bv = *C_proto[ToSize_T(thread)];
                                       
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                       const Int C = leaf_clusters[i];
                       
                       C_bv.SetPointer( C_up_ser, C );
                       C_bv.FromPrimitives( P_, P_up_ser, C_begin[C], C_end[C], 1 );
                    }
                },
                ThreadCount()
            );
            
            TOOLS_PTOC(className()+"::TakeUpdateVectors - Compute the AABBs of the updated leaf clusters.");
            
            TOOLS_PTIC(className()+"::TakeUpdateVectors - Upward pass for AABBs.");
            {
                Int stack [128];
                Int stack_ptr = 0;
                stack[stack_ptr] = 0;

                Tensor1<bool,Int> visited( ClusterCount(), false );

                mref<BoundingVolume_T> C_bv = *C_proto[0];

                while( (0 <= stack_ptr) && (stack_ptr < 126) )
                {
                    const Int C = stack[stack_ptr];
                    const Int L = C_left [C];
                    const Int R = C_right[C];

                    // Using that children are either both interiors or both leaves.
                    if( (L >= 0) /*&& (R >= 0)*/ )
                    {
                        if( visited[C] )
                        {
                            C_bv.SetPointer ( C_up_ser, C );
                            C_bv.Read       ( C_up_ser, L );
                            C_bv.Merge      ( C_up_ser, R );
                            --stack_ptr;
                        }
                        else
                        {
                            stack[++stack_ptr] = R;
                            stack[++stack_ptr] = L;
                            visited[C] = true;
                        }
                    }
                    else
                    {
                        --stack_ptr;
                    }
                    
                } // while
                
                if( stack_ptr >= 126 )
                {
                    eprint( ClassName() + "::TakeUpdateVectors: Stack overflow." );
                }
            }
            TOOLS_PTOC(className()+"::TakeUpdateVectors - Upward pass for AABBs.");
            
            TOOLS_PTOC(className()+"::TakeUpdateVectors");
            
        } // TakeUpdateVectors
        
        void CollectNearFieldDerivatives( bool addto = false ) const override
        {
            // Collect derivative data from the primitives, stored in thread_P_D_near, into P_D_near.
            
            TOOLS_PTIC(className()+"::CollectNearFieldDerivatives");
            
            const Int near_dim = NearDim();
            
            mptr<Real> to      = P_D_near.data();
            cptr<Int>  inv_ord = P_inverse_ordering.data();

            // Add the other slices. (Start with the first slice if we were asked to add into.)
            for( Int thread = 0; thread < ThreadCount(); ++thread )
            {
                cptr<Real> from = thread_P_D_near.data(thread);

                if( thread == 0 && !addto )
                {
                    ParallelDo(
                        [=]( const Int i )
                        {
                            copy_buffer<VarSize,Sequential>(
                                &from[near_dim * inv_ord[i]], &to[near_dim * i], near_dim
                            );
                        },
                        PrimitiveCount(), ThreadCount()
                    );
                }
                else
                {
                    ParallelDo(
                        [=]( const Int i )
                        {
                            add_to_buffer( &from[near_dim * inv_ord[i]], &to[near_dim * i], near_dim );
                        },
                        PrimitiveCount(), ThreadCount()
                    );
                }
            }
            TOOLS_PTOC(className()+"::CollectNearFieldDerivatives");
            
        } // CollectNearFieldDerivatives
        
        
        void CollectFarFieldDerivatives( bool addto = false ) const override
        {
            // Collect derivative data from the clusters, stored in thread_C_D_far, into C_out; then use PercolatedDown and ClustersToPrimitives to collect these into P_D_far.
            
            TOOLS_PTIC(className()+"::CollectFarFieldDerivatives");
            
            const Int far_dim = FarDim();
            
            this->RequireBuffers(FarDim());
            
            thread_C_D_far.AddReduce( C_out.data(), false );
            
            this->PercolateDown();
            
            this->ClustersToPrimitives(false);

            cptr<Real> from    = P_out.data();
            mptr<Real> to      = P_D_far.data();
            cptr<Int>  inv_ord = P_inverse_ordering.data();
            
            // Finally, permute data for the outside world.
            if( addto )
            {
                ParallelDo(
                    [=]( const Int i )
                    {
                        add_to_buffer( &from[far_dim * inv_ord[i]], &to[far_dim * i], far_dim );
                    },
                    PrimitiveCount(), ThreadCount()
                );
            }
            else
            {
                ParallelDo(
                    [=]( const Int i )
                    {
                        copy_buffer<VarSize,Sequential>(
                            &from[far_dim * inv_ord[i]], &to[far_dim * i], far_dim
                        );
                    },
                    PrimitiveCount(), ThreadCount()
                );
            }
            
            TOOLS_PTOC(className()+"::CollectFarFieldDerivatives");
            
        } // CollectFarFieldDerivatives
        
        
        cref<Tensor2<Real,Int>> ClusterMoments() const override
        {
            return C_moments;
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
            return  "ClusterTree<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }

    }; // class ClusterTree
    
} // namespace Repulsor
