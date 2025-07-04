#pragma once

#include "SimplicialMesh/SimplicialMeshBase.hpp"

namespace Repulsor
{
    template<typename Real_, typename Int_, typename ExtReal_, typename ExtInt_>
    class SimplicialRemesherBase;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename ExtReal_, typename ExtInt_>
    class SimplicialRemesher;

    template<
        int DOM_DIM, int AMB_DIM,
        typename Real_, typename Int_ = int, typename LInt_ = Size_T,
        typename SReal_ = Real_, typename ExtReal_ = Real_
    >
    class SimplicialMesh : public SimplicialMeshBase<Real_,Int_,LInt_,SReal_,ExtReal_>
    {   
        static_assert(DOM_DIM >= 0       , "Domain dimension must be nonnegative.");
        static_assert(AMB_DIM >= DOM_DIM , "Ambient dimension must at least as high as the domain dimension.");
        
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
        static_assert(IntQ< Int_>,"");
        static_assert(IntQ<LInt_>,"");
        
    public:
        
        using Base_T = SimplicialMeshBase<Real_,Int_,LInt_,SReal_,ExtReal_>;
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = LInt_;
        
        using TangentVector_T         = typename Base_T::TangentVector_T;
        using CotangentVector_T       = typename Base_T::CotangentVector_T;
        
        using SparseMatrix_T          = typename Base_T::SparseMatrix_T;
        using SparseBinaryMatrix_T    = typename Base_T::SparseBinaryMatrix_T;
        
        using SimplexDataKernel_T     = SimplexDataKernel<DOM_DIM,AMB_DIM,Real,Int>;
        using Primitive_T             = Polytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using MovingPrimitive_T       = MovingPolytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using BoundingVolume_T        = AABB<AMB_DIM,GJK_Real,Int,SReal>;
        
        
        using ClusterTree_T           = ClusterTree      <AMB_DIM,Real,Int,LInt,SReal,ExtReal>;
        using BlockClusterTree_T      = BlockClusterTree <AMB_DIM,Real,Int,LInt,SReal,ExtReal,true>;
        using ObstacleBlockClusterTree_T
                                      = BlockClusterTree <AMB_DIM,Real,Int,LInt,SReal,ExtReal,false>;
        using CollisionTree_T         = CollisionTree    <AMB_DIM,Real,Int,LInt,SReal,ExtReal,true>;
        using ObstacleCollisionTree_T = CollisionTree    <AMB_DIM,Real,Int,LInt,SReal,ExtReal,false>;

        using Remesher_T              = SimplicialRemesher<DOM_DIM,AMB_DIM,Real,Int,Real,Int>;
        
        using RemesherBase_T          = typename Base_T::RemesherBase_T;
        
        using Obstacle_T              = Base_T;
        
        static constexpr Int  FAR_DIM = 1 + AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        static constexpr Int NEAR_DIM = 1 + (DOM_DIM+1) * AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        
        static constexpr Int  SIZE      = DOM_DIM + 1;
        static constexpr Int  HULL_SIZE = AMB_DIM * SIZE;
        static constexpr Real nth       = Inv<Real>( SIZE );
        
        SimplicialMesh(
            cref<Tensor2<Real,Int>> V_coords_,
            // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            cref<Tensor2<Int, Int>> simplices_,
            // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T thread_count_ = 1
        )
        :   SimplicialMesh(
                  V_coords_.data(),
                  V_coords_.Dim(0),
                  false,
                  simplices_.data(),
                  simplices_.Dim(0),
                  false,
                  int_cast<Int>(thread_count_)
            )
        {
            TOOLS_PTIC(className()+"()");
            if( V_coords_.Dim(1) != AMB_DIM )
            {
                eprint(className()+" : V_coords.Dim(1) != AMB_DIM");
                TOOLS_PTOC(className()+"()");
                return;
            }
            if( simplices_.Dim(1) != DOM_DIM+1 )
            {
                eprint(className()+" : simplices_.Dim(1) != DOM_DIM+1");
                TOOLS_PTOC(className()+"()");
                return;
            }
            TOOLS_PTOC(className()+"()");
        }
        
        SimplicialMesh(
            cref<Tensor2<Real,Int>> V_coords_,
            // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            cref<Tensor2<Int, Int>> simplices_,
            // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            cref<Tensor1<Real,Int>> V_charges_,
            // vertex charges; assumed to be of size DOM_DIM+1
            const Size_T thread_count_ = 1
        )
        :   SimplicialMesh(
                  V_coords_.data(),
                  V_coords_.Dim(0),
                  false,
                  simplices_.data(),
                  simplices_.Dim(0),
                  false,
                  V_charges_.data(),
                  int_cast<Int>(thread_count_)
            )
        {
            TOOLS_PTIC(className()+"()");
            if( V_coords_.Dim(1) != AMB_DIM )
            {
                eprint(className()+" : V_coords.Dim(1) != AMB_DIM");
                TOOLS_PTOC(className()+"()");
                return;
            }
            if( simplices_.Dim(1) != DOM_DIM+1 )
            {
                eprint(className()+" : simplices_.Dim(1) != DOM_DIM+1");
                TOOLS_PTOC(className()+"()");
                return;
            }
            if( V_charges_.Dim(0) != V_coords_.Dim(0) )
            {
                eprint(className()+" : V_charges_.Dim(0) != V_coords_.Dim(0) ");
                TOOLS_PTOC(className()+"()");
                return;
            }
            TOOLS_PTOC(className()+"()");
        }

        template<typename ExtReal_2, typename ExtInt>
        SimplicialMesh(
            const ExtReal_2 * V_coords_, // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Size_T vertex_count_,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T simplex_count_,
            const Size_T thread_count_ = 1
        )
        :   SimplicialMesh( V_coords_, vertex_count_, false, simplices_, simplex_count_, false, thread_count_)
        {}
        
        template<typename ExtReal_2, typename ExtInt>
        SimplicialMesh(
            const ExtReal_2 * vertex_coords_, // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Size_T vertex_count_,
            const bool vertex_coords_ColMajorQ,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T simplex_count_,
            const bool simplices_ColMajorQ,
            const Size_T thread_count_ = 1
        )
        :   Base_T    ( int_cast<Int>(thread_count_) )
        ,   V_coords  ( ToTensor2<Real,Int>(
                            vertex_coords_,
                            int_cast<Int>(vertex_count_),
                            int_cast<Int>(AMB_DIM),
                            vertex_coords_ColMajorQ
                        )
                      )
        ,   V_coords_frozen ( V_coords )
        ,   simplices ( ToTensor2<Int,Int>(
                            simplices_,
                            int_cast<Int>(simplex_count_),
                            int_cast<Int>(DOM_DIM+1),
                            simplices_ColMajorQ
                        )
                      )
        ,   V_charges ( int_cast<Int>(vertex_count_ ), Scalar::One<Real> )
        ,   P_charges ( int_cast<Int>(simplex_count_), Scalar::One<Real> )
        ,   details   ( int_cast<Int>(thread_count_ )                    )
        {
            TOOLS_PTIC(className()+" (pointer)");
            TOOLS_PTOC(className()+" (pointer)");
        }

        
        template<typename ExtReal_2, typename ExtReal_3, typename ExtInt>
        SimplicialMesh(
            const ExtReal_2 * vertex_coords_, // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Size_T      vertex_count_,
            const bool        vertex_coords_ColMajorQ,
            const ExtInt    * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T      simplex_count_,
            const bool        simplices_ColMajorQ,
            const ExtReal_3 * vertex_charges_, // vertex charges; assumed to be of size vertex_count_
            const Size_T      thread_count_ = 1
        )
        :   Base_T    ( int_cast<Int>(thread_count_) )
        ,   V_coords  ( ToTensor2<Real,Int>(
                            vertex_coords_,
                            int_cast<Int>(vertex_count_),
                            int_cast<Int>(AMB_DIM),
                            vertex_coords_ColMajorQ
                        )
                      )
        ,   V_coords_frozen ( V_coords )
        ,   simplices ( ToTensor2<Int,Int>(
                            simplices_,
                            int_cast<Int>(simplex_count_),
                            int_cast<Int>(DOM_DIM+1),
                            simplices_ColMajorQ
                        )
                      )
        ,   V_charges ( vertex_charges_, int_cast<Int>(vertex_count_ )                    )
        ,   P_charges (                  int_cast<Int>(simplex_count_), Scalar::One<Real> )
        ,   details   ( int_cast<Int>(thread_count_) )
        {
            TOOLS_PTIC(className()+" (pointer,charges)");
            TOOLS_PTOC(className()+" (pointer,charges)");
        }

        // Default constructor
        SimplicialMesh() = default;
        // Destructor
        virtual ~SimplicialMesh() override = default;
        // Copy constructor
        SimplicialMesh( const SimplicialMesh & other ) = default;
        // Copy assignment operator
        SimplicialMesh & operator=( const SimplicialMesh & other ) = default;
        // Move constructor
        SimplicialMesh( SimplicialMesh && other ) = default;
        // Move assignment operator
        SimplicialMesh & operator=( SimplicialMesh && other ) = default;
        
    public:
        
        using Base_T::block_cluster_tree_settings;
        using Base_T::cluster_tree_settings;
        using Base_T::adaptivity_settings;
        using Base_T::ThreadCount;

    protected:

        
        mutable SReal max_update_step_size = 0;
        
        // Current vertex coordinates; will be updated by SemiStaticUpdate.
        mutable Tensor2<Real,Int> V_coords;
        // Frozen vertex coordinates to be used for the clustering, but not for the energy computation.
        // Won't be updated by SemiStaticUpdate so that energies that depend on the `cluster_tree` will look differentiable.
        Tensor2<Real,Int> V_coords_frozen;
        
        Tensor2<Int, Int> simplices;
        
        Tensor1<Real,Int> V_charges; // vertex charges; currently only used for DOM_DIM = 0;
        
        mutable Tensor1<Real,Int> P_charges; // primitive charges; currently only used for DOM_DIM = 0;
                
        mutable Primitive_T P_proto;
        
        mutable Tensor1<Int,Int> simplex_row_pointers;
        mutable Tensor1<Int,Int> simplex_column_indices;
        
        SimplicialMeshDetails<DOM_DIM,AMB_DIM,Real,Int,LInt> details;
        
        // Two buffers that are supposed to help with solves.
        mutable Tensor2<Real,Int> X_buf;
        mutable Tensor2<Real,Int> Y_buf;

#include "SimplicialMesh/ComputeNearFarDataOps.hpp"
#include "SimplicialMesh/ComputeNearFarData.hpp"
#include "SimplicialMesh/NestedDissectionOrdering.hpp"
#include "SimplicialMesh/ApproximateMinimumDegreeOrdering.hpp"
#include "SimplicialMesh/MetisOrdering.hpp"
#include "SimplicialMesh/H1Metric.hpp"

        
    public:
        
        virtual Int DomDim() const override
        {
            return DOM_DIM;
        }
        
        virtual Int AmbDim() const override
        {
            return AMB_DIM;
        }

        virtual cref<Tensor2<Real,Int>> VertexCoordinates() const override
        {
            return V_coords;
        }
        
        virtual cref<Tensor2<Int,Int>> Simplices() const override
        {
            return simplices;
        }
        
        virtual cref<Tensor1<Real,Int>> VertexCharges() const override
        {
            return V_charges;
        }
        
        virtual cref<Tensor1<Real,Int>> SimplexCharges() const override
        {
            return P_charges;
        }

        virtual constexpr Int FarDim() const override
        {
            return FAR_DIM;
        }
        
        virtual constexpr Int NearDim() const override
        {
            return NEAR_DIM;
        }
        
        virtual Int VertexCount() const override
        {
            return V_coords.Dim(0);
        }
        
        virtual Int SimplexCount() const override
        {
            return simplices.Dim(0);
        }
    
        virtual Int DofCount() const override
        {
            return VertexCount() * AMB_DIM;
        }
        
        virtual cptr<Real> Dofs() const override
        {
            return V_coords.data();
        }
        
        virtual mref<Tensor2<Real,Int>> XBuffer( const Int nrhs ) const override
        {
            X_buf.template RequireSize<false>( VertexCount(), nrhs );
            
            return X_buf;
        }
        
        virtual mref<Tensor2<Real,Int>> YBuffer( const Int nrhs ) const override
        {
            Y_buf.template RequireSize<false>( VertexCount(), nrhs );
            
            return Y_buf;
        }
        
        virtual void SemiStaticUpdate( cptr<ExtReal> V_coords_, const bool transp_ = false ) const override
        {
            TOOLS_PTIC(className()+"::SemiStaticUpdate");
            
            // We read the new coordinates onto V_coords, but not into C_coords_frozen!
            if( transp_ )
            {
                V_coords.ReadTransposed(V_coords_);
            }
            else
            {
                V_coords.Read(V_coords_);
            }
            
            this->ClearCache();
            
            Tensor2<Real,Int> P_near( SimplexCount(), NEAR_DIM );
            Tensor2<Real,Int> P_far ( SimplexCount(), FAR_DIM  );

            ComputeNearFarData( P_near, P_far );
            
            GetClusterTree().SemiStaticUpdate( P_near, P_far );
            
            TOOLS_PTOC(className()+"::SemiStaticUpdate");
        }
       
        
        virtual void LoadUpdateVectors(
            cptr<ExtReal> vecs,
            const ExtReal max_time,
            const bool transp_
        ) const override
        {
            TOOLS_PTIC(className()+"::LoadUpdateVectors");
            
            max_update_step_size = static_cast<SReal>(max_time);
            
            Tensor2<Real,Int> V_updates ( VertexCount(), AMB_DIM );

            if( vecs == nullptr )
            {
                V_updates.Fill(Scalar::Zero<Real>);
            }
            else
            {
                if( transp_ )
                {
                    V_updates.ReadTransposed(vecs);
                }
                else
                {
                    V_updates.Read(vecs);
                }
            }
            // ATTENTION: We reorder already here outside of the cluster tree to save a copy operation of a big Tensor2!
            
            MovingPrimitive_T P_moving;
            
            Tensor2<SReal,Int> P_velocities_serialized ( SimplexCount(), P_moving.VelocitySize(), 0 );
            mptr<SReal> P_v_ser = P_velocities_serialized.data();
            
            cref<Tensor1<Int,Int>> P_ordering = GetClusterTree().PrimitiveOrdering();

            JobPointers<Int> job_ptr ( SimplexCount(), GetClusterTree().ThreadCount() );
            
            ParallelDo(
                [&]( const Int thread )
                {
                    MovingPrimitive_T P_mov;
                    
                    const Int i_begin = job_ptr[thread  ];
                    const Int i_end   = job_ptr[thread+1];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        const Int j = P_ordering[i];
                        P_mov.FromVelocitiesIndexList( V_updates.data(), simplices.data(), j );
                        P_mov.WriteVelocitiesSerialized( P_v_ser, i );
                    }
                },
                job_ptr.ThreadCount()
            );
            
            // P_moving and P_velocities_serialized will be swapped against potentiall empty containers.
            GetClusterTree().TakeUpdateVectors(
                P_moving, P_velocities_serialized, static_cast<SReal>(max_time)
            );
            
            TOOLS_PTOC(className()+"::LoadUpdateVectors");
        }
        
        virtual ExtReal MaximumSafeStepSize(
            cptr<ExtReal> vecs,
            const ExtReal max_time,
            const ExtReal TOL = scalar_cast<ExtReal>(0.0625),
            const bool transp_ = false
        ) override
        {
            TOOLS_PTIC(className()+"::MaximumSafeStepSize");
            
            LoadUpdateVectors( vecs, max_time, transp_ );
            
            ExtReal t = max_time;

            if( this->InPersistentCacheQ("Obstacle") )
            {
                GetObstacle().LoadUpdateVectors( static_cast<ExtReal *>(nullptr), max_time );
                
                t = GetObstacleCollisionTree().MaximumSafeStepSize(t,TOL);
                
                logprint("GetObstacleCollisionTree().MaximumSafeStepSize(t) = "+ToString(t));
            }
                
            t = GetCollisionTree().MaximumSafeStepSize(t,TOL);
            
            logprint("GetCollisionTree().MaximumSafeStepSize(t)         = "+ToString(t));
            
            TOOLS_PTOC(className()+"::MaximumSafeStepSize");
            
            return t;
        }
        
        virtual cref<ClusterTree_T> GetClusterTree() const override
        {
            static std::string tag ( "ClusterTree" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                TOOLS_PTIC(className()+"::GetClusterTree");
                if( (V_coords.Dim(0) > 0) && (simplices.Dim(0) > 0) )
                {
                    auto P_coords      = Tensor2<Real,Int> ( SimplexCount(), AMB_DIM, Scalar::Zero<Real> );
                    auto P_hull_coords = Tensor3<Real,Int> ( SimplexCount(), SIZE, AMB_DIM );
                    auto P_near        = Tensor2<Real,Int> ( SimplexCount(), NEAR_DIM );
                    auto P_far         = Tensor2<Real,Int> ( SimplexCount(), FAR_DIM );

                    Tensor2<SReal,Int> P_serialized ( SimplexCount(), P_proto.Size() );

                    auto DiffOp = SparseMatrix_T(
                        SimplexCount() * AMB_DIM,
                        VertexCount(),
                        SimplexCount() * HULL_SIZE,
                        ThreadCount()
                    );

                    DiffOp.Outer()[SimplexCount() * AMB_DIM] = SimplexCount() * HULL_SIZE;

                    auto AvOp = SparseMatrix_T(
                        SimplexCount(),
                        VertexCount(),
                        SimplexCount() * SIZE,
                        ThreadCount()
                    );

                    AvOp.Outer()[SimplexCount()] = SimplexCount() * SIZE;

                    // What remains is to compute P_coords, P_hull_coords, P_near and P_far and the nonzero values of DiffOp.
                    
                    ComputeNearFarDataOps( P_coords, P_hull_coords, P_near, P_far, DiffOp, AvOp );


                    const JobPointers<Int> job_ptr ( SimplexCount(), ThreadCount() );

//                    TOOLS_PTIC("Creating primitives");
                    ParallelDo(
                        [&]( const Int thread )
                        {
                            Primitive_T P;

                            const Int i_begin = job_ptr[thread];
                            const Int i_end   = job_ptr[thread+1];

                            for( Int i = i_begin; i < i_end; ++i )
                            {
                                P.SetPointer( P_serialized.data(), i );
                                // Beware, we use the frozen coordinates for clustering!
                                P.FromIndexList( V_coords_frozen.data(), simplices.data(), i );
                            }
                        },
                        job_ptr.ThreadCount()
                    );

//                    TOOLS_PTOC("Creating primitives");

//                    TOOLS_PTIC("Initializing cluster prototypes");

                    std::shared_ptr<BoundingVolume_T> C_proto;

                    switch( cluster_tree_settings.bounding_volume_type )
                    {
                        case BoundingVolumeType::AABB_MedianSplit:
                        {
                            logprint("Using AABB_MedianSplit as bounding volume.");
                            C_proto = std::shared_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_MedianSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                            break;
                        }
                        case BoundingVolumeType::AABB_PreorderedSplit:
                        {
                            logprint("Using AABB_PreorderedSplit as bounding volume.");
                            C_proto = std::shared_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_PreorderedSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                            break;
                        }
                        default:
                        {
                            logprint("Using AABB_LongestAxisSplit as bounding volume.");
                            C_proto = std::shared_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_LongestAxisSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                        }
                    }

//                    TOOLS_PTOC("Initializing cluster prototypes");

                    if( cluster_tree_settings.thread_count <= 0 )
                    {
                        cluster_tree_settings.thread_count = ThreadCount();
                    }

                    this->SetPersistentCache( tag,
                        std::make_any<ClusterTree_T>(
                            P_proto, P_serialized, *C_proto, Tensor1<Int,Int>(0),
                            P_near, P_far,
                            DiffOp, AvOp,
                            cluster_tree_settings
                        )
                    );
                }
                else
                {
                    eprint(ClassName()+"(V_coords.Dim(0) <= 0) || (simplices.Dim(0) <= 0)");

                    this->SetPersistentCache( tag, std::make_any<ClusterTree_T>() );
                }

                TOOLS_PTOC(className()+"::GetClusterTree");
            }

            return this->template GetPersistentCache<ClusterTree_T>(tag);
        }

        virtual cref<BlockClusterTree_T> GetBlockClusterTree() const override
        {
            static std::string tag ( "BlockClusterTree" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                TOOLS_PTIC(className()+"::GetBlockClusterTree");
                
                block_cluster_tree_settings.near_field_separation_parameter = adaptivity_settings.theta;
                block_cluster_tree_settings.near_field_intersection_parameter  = adaptivity_settings.intersection_theta;
                
                block_cluster_tree_settings.max_refinement = adaptivity_settings.max_refinement;
                
                this->SetPersistentCache( tag,
                    std::make_any<BlockClusterTree_T>(
                        GetClusterTree(),
                        GetClusterTree(),
                        block_cluster_tree_settings
                    )
                );
                
                TOOLS_PTOC(className()+"::GetBlockClusterTree");
            }
            
            return this->template GetPersistentCache<BlockClusterTree_T>(tag);
        }
        
        virtual cref<CollisionTree_T> GetCollisionTree() const override
        {
            static std::string tag ( "CollisionTree" );
            
            if( !this->InCacheQ( tag ) )
            {
                TOOLS_PTIC(className()+"::GetCollisionTree");

                this->SetCache( tag,
                    std::make_any<CollisionTree_T>( GetClusterTree(), GetClusterTree() )
                );

                TOOLS_PTOC(className()+"::GetCollisionTree");
            }
            
            return this->template GetCache<CollisionTree_T>(tag);
        }
        
        cref<SparseBinaryMatrix_T> DerivativeAssembler() const override
        {
            
            static std::string tag ( "DerivativeAssembler" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                TOOLS_PTIC(className()+"::DerivativeAssembler");
            
                auto A = SparseBinaryMatrix_T(
                    SimplexCount() * SIZE,
                    VertexCount(),
                    SimplexCount() * SIZE,
                    ThreadCount()
                );
                
                A.Outer().iota();
                A.Inner().Read(simplices.data());
                
                this->SetPersistentCache( tag, std::move(A.Transpose()) );
                
                TOOLS_PTOC(className()+"::DerivativeAssembler");
            }
            
            return this->template GetPersistentCache<SparseBinaryMatrix_T>(tag);
            
        } // DerivativeAssembler
        
        void Assemble_ClusterTree_Derivatives(
            const ExtReal alpha, const ExtReal beta, mptr<ExtReal> Y, const Int ldY
        ) const override
        {
            TOOLS_PTIC(className()+"::Assemble_ClusterTree_Derivatives");
            
            Tensor3<Real,Int> buffer ( SimplexCount(), SIZE, AMB_DIM, Scalar::Zero<Real> );
            
            GetClusterTree().CollectDerivatives();

            details.DNearFarToHulls(
                V_coords,
                simplices,
                GetClusterTree().PrimitiveDNearFieldData(),
                GetClusterTree().PrimitiveDFarFieldData(),
                V_charges,
                buffer
            );

            DerivativeAssembler().template Dot<AMB_DIM>(
                alpha, buffer.data(), AMB_DIM,
                beta , Y            , ldY,
                AMB_DIM
            );

            TOOLS_PTOC(className()+"::Assemble_ClusterTree_Derivatives");
        }
        
        void Assemble_ClusterTree_VertexDensities( 
            mptr<ExtReal> output,
            const ExtReal weight,
            bool addTo = false
        ) const override
        {
            TOOLS_PTIC(className()+"::Assemble_ClusterTree_VertexDensities");
            
            GetClusterTree().CollectVertexDensities( output, weight, addTo );
            
            TOOLS_PTOC(className()+"::Assemble_ClusterTree_VertexDensities");
        }
        
        void Assemble_ClusterTree_SimplexDensities( mptr<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            TOOLS_PTIC(className()+"::Assemble_ClusterTree_SimplexDensities");
            
            GetClusterTree().CollectPrimitiveDensities( output, weight, addTo );
            
            TOOLS_PTOC(className()+"::Assemble_ClusterTree_SimplexDensities");
        }

//#####################################################################################
//      Obstacle
//#####################################################################################

        virtual void LoadObstacle( std::unique_ptr<Obstacle_T> obstacle_ ) override
        {
            static std::string tag ( "Obstacle" );
            
            std::shared_ptr<Obstacle_T> obstacle;
            
            // Input obstacle is moved.
            if( obstacle_->AmbDim() != AmbDim() )
            {
                eprint(className()+"::LoadObstacle: Attempted to load obstacle of ambient dimension "+ToString(obstacle_->AmbDim())+" into mesh of ambient dimension "+ToString(AmbDim())+". Setting obstacle to nullptr."
                );
                
                // We have to construct an empy SimplicialMesh because Obstacle_T is abstract.
                obstacle = std::make_shared<SimplicialMesh>();
            
            }
            else
            {
                obstacle = std::move(obstacle_);
            }
            
            this->SetPersistentCache(tag, obstacle );
        }
        
        cref<Obstacle_T> GetObstacle() const override
        {
            static std::string tag ("Obstacle");

            if( !this->InPersistentCacheQ(tag) )
            {
                wprint( ClassName()+"::GetObstacle: Obstacle not initialized.");
                
                // We have to construct an empty SimplicialMesh because Obstacle_T is abstract.
                this->SetPersistentCache(tag, std::make_shared<SimplicialMesh>() );
            }
            
            return *(this->template GetPersistentCache<std::shared_ptr<Obstacle_T>>(tag) );
        }
        
        virtual cref<ClusterTree_T> GetObstacleClusterTree() const override
        {
            return *dynamic_cast<const ClusterTree_T *>( & (GetObstacle().GetClusterTree()) );
        }
        
        virtual cref<ObstacleBlockClusterTree_T> GetObstacleBlockClusterTree() const override
        {
            static std::string tag ( "ObstacleBlockClusterTree" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                TOOLS_PTIC(className()+"::GetObstacleBlockClusterTree");
                
                this->SetPersistentCache( tag,
                    std::make_any<ObstacleBlockClusterTree_T>(
                        GetClusterTree(), GetObstacleClusterTree(), block_cluster_tree_settings
                    )
                );
                
                TOOLS_PTOC(className()+"::GetObstacleBlockClusterTree");

            }
            
            return this->template GetPersistentCache<ObstacleBlockClusterTree_T>(tag);
        }
        
        virtual cref<ObstacleCollisionTree_T> GetObstacleCollisionTree() const override
        {
            static std::string tag ( "ObstacleCollisionTree" );
            
            if( !this->InPersistentCacheQ(tag) )
            {
                TOOLS_PTIC(className()+"::GetObstacleCollisionTree");
                
                this->SetPersistentCache( tag,
                    std::make_any<ObstacleCollisionTree_T>( GetClusterTree(), GetObstacleClusterTree() )
                );

                TOOLS_PTOC(className()+"::GetObstacleCollisionTree");
            }
            
            return this->template GetPersistentCache<ObstacleCollisionTree_T>(tag);
        }

        
//#####################################################################################
//      IO
//#####################################################################################
     
    public:
        
        virtual void WriteToFile( const std::filesystem::path & file ) const override
        {
            TOOLS_PTIC(ClassName()+"::WriteToFile");
            
            print("Writing mesh to file "+ file.string() +".");
            
            std::ofstream s (file);
            
            const Int vertex_count = VertexCount();
            const Int simplex_count = SimplexCount();
            
            cptr<Real> V = V_coords.data();
            cptr<Int>  S = simplices.data();
            
            s << "domain_dimension" << "\t" << DOM_DIM << "\n";
            s << "ambient_dimension" << "\t" << AMB_DIM << "\n";
            s << "vertex_count" << "\t" << vertex_count << "\n";
            s << "simplex_count" << "\t" << simplex_count << "\n";
            
            s << std::scientific << std::uppercase << std::setprecision( std::numeric_limits<Real>::digits10 + 1 );
            
            for( Int i = 0; i < vertex_count; ++i )
            {
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    s << V[AMB_DIM * i + k] << "\t";
                }
                s <<"\n";
            }
            
            for( Int i = 0; i < simplex_count; ++i )
            {
                for( Int k = 0; k < SIZE; ++k )
                {
                    s << S[SIZE * i + k] << "\t";
                }
                s <<"\n";
            }
            
            TOOLS_PTOC(ClassName()+"::WriteToFile");
        }
        
//#####################################################################################
//      Remesher
//#####################################################################################
        
    public:
        
        [[nodiscard]] virtual std::unique_ptr<RemesherBase_T> CreateRemesher() override
        {
            if constexpr( DOM_DIM > 0 )
            {
            TOOLS_PTIC(ClassName()+"::CreateRemesher");
            
                Remesher_T * R = new Remesher_T(
                    VertexCoordinates().data(),
                    VertexCoordinates().Dim(0), false,
                    Simplices().data(),
                    Simplices().Dim(0),         false,
                    VertexCharges().data(),     1,
                    ThreadCount()
                );
                
                TOOLS_PTOC(ClassName()+"::CreateRemesher");
                
                return std::unique_ptr<RemesherBase_T>(R);
                
            }
            else
            {
                wprint(ClassName()+"::CreateRemesher: Remesher does not make sense for meshes with domain dimension 0. Returning nullptr.");
                return std::unique_ptr<RemesherBase_T>(nullptr);
            }
        }
        
//#####################################################################################
//      Standard interface
//#####################################################################################
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    
    private:
  
        static std::string className()
        {
            return "SimplicialMesh<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+TypeName<LInt>+">";
        }
    };
} // namespace Repulsor


#include "SimplicialMesh/SimplicialMesh_Factory.hpp"
