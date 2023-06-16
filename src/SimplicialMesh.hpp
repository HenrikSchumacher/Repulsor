#pragma once

#include "SimplicialMesh/SimplicialMeshBase.hpp"

namespace Repulsor
{
//    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
//    class TangentPoint;
    
    template<typename Real_, typename Int_>
    class SimplicialRemesherBase;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_>
    class SimplicialRemesher;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class SimplicialMesh : public SimplicialMeshBase<Real_,Int_,SReal_,ExtReal_>
    {
    private:
        
        using Base_T = SimplicialMeshBase<Real_,Int_,SReal_,ExtReal_>;
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = typename Base_T::LInt;
        
        using TangentVector_T      = typename Base_T::TangentVector_T;
        using CotangentVector_T    = typename Base_T::CotangentVector_T;
        
        using SparseMatrix_T       = typename Base_T::SparseMatrix_T;
        using SparseBinaryMatrix_T = typename Base_T::SparseBinaryMatrix_T;
        
        using Primitive_T = Polytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using MovingPrimitive_T = MovingPolytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using BoundingVolume_T = AABB<AMB_DIM,GJK_Real,Int,SReal>;
        
        
        using ClusterTree_T           = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        using BlockClusterTree_T      = BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal,true>;
        using ObstacleBlockClusterTree_T
                                      = BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal,false>;
        using CollisionTree_T         = CollisionTree<AMB_DIM,Real,Int,SReal,ExtReal,true>;
        using ObstacleCollisionTree_T = CollisionTree<AMB_DIM,Real,Int,SReal,ExtReal,false>;

        using Remesher_T         = SimplicialRemesher<DOM_DIM,AMB_DIM,Real,Int>;
        
        using RemesherBase_T     = typename Base_T::Remesher_T;
        
        using Obstacle_T         = Base_T;
        
        static constexpr Int  FAR_DIM = 1 + AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        static constexpr Int NEAR_DIM = 1 + (DOM_DIM+1) * AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        
        SimplicialMesh() = default;

        SimplicialMesh(
            const Tensor2<Real,Int> & V_coords_,
            // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Tensor2<Int,Int> & simplices_,
            // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T thread_count_ = 1
        )
        :   SimplicialMesh(
                  V_coords_.data(),
                  V_coords_.Dimension(0),
                  false,
                  simplices_.data(),
                  simplices_.Dimension(0),
                  false,
                  static_cast<Int>(thread_count_)
            )
        {
            ptic(className()+"()");
            if( V_coords_.Dimension(1) != AMB_DIM )
            {
                eprint(className()+" : V_coords.Dimension(1) != AMB_DIM");
                ptoc(className()+"()");
                return;
            }
            if( simplices_.Dimension(1) != DOM_DIM+1 )
            {
                eprint(className()+" : simplices_.Dimension(1) != DOM_DIM+1");
                ptoc(className()+"()");
                return;
            }
            ptoc(className()+"()");
        }
        
#ifdef LTEMPLATE_H
        
        SimplicialMesh(
            const mma::TensorRef<ExtReal> & V_coords_,   // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const mma::TensorRef<long long> & simplices_,   // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const long long thread_count_ = 1
        )
        :   SimplicialMesh(
                V_coords_.data(),
                static_cast<Int>(V_coords_.dimensions()[0]),
                simplices_.data(),
                static_cast<Int>(simplices_.dimensions()[0]),
                static_cast<Int>(thread_count_)
            )
        {
            ptic(className()+"() (from MTensor)");
            if( V_coords_.dimensions()[1] != AMB_DIM )
            {
                eprint(className()+" : V_coords_.dimensions()[1] != AMB_DIM");
                ptoc(className()+"() (from MTensor)");
                return;
            }
            if( simplices_.dimensions()[1] != DOM_DIM+1 )
            {
                eprint(className()+" : simplices_.dimensions()[1] != DOM_DIM+1");
                ptoc(className()+"() (from MTensor)");
                return;
            }
            ptoc(className()+"() (from MTensor)");
        }
        
#endif

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
            const bool vertex_coords_transpose,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T simplex_count_,
            const bool simplices_transpose,
            const Size_T thread_count_ = 1
        )
        :   Base_T    ( static_cast<Int>(thread_count_) )
        ,   V_coords  ( ToTensor2<Real,Int>(
                            vertex_coords_,
                            static_cast<Int>(vertex_count_),
                            static_cast<Int>(AMB_DIM),
                            vertex_coords_transpose
                        )
                      )
        ,   simplices ( ToTensor2<Int,Int>(
                            simplices_,
                            static_cast<Int>(simplex_count_),
                            static_cast<Int>(DOM_DIM+1),
                            simplices_transpose
                        )
                      )
        ,   details   ( static_cast<Int>(thread_count_) )
        {
            ptic(className()+" (pointer)");
            ptoc(className()+" (pointer)");
        }

        
        virtual ~SimplicialMesh() override = default;
        
    public:
        
        using Base_T::block_cluster_tree_settings;
        using Base_T::cluster_tree_settings;
        using Base_T::adaptivity_settings;
        using Base_T::ThreadCount;
        
    protected:

        
        mutable SReal max_update_step_size = 0;
        
        Tensor2<Real,Int> V_coords;
        Tensor2<Int,Int>  simplices;
                
        mutable Primitive_T P_proto;
        
        mutable Tensor1<Int,Int> simplex_row_pointers;
        mutable Tensor1<Int,Int> simplex_column_indices;
        
        SimplicialMeshDetails<DOM_DIM,AMB_DIM,Real_,Int_> details;

        
    public:
        
        virtual Int DomDim() const override
        {
            return DOM_DIM;
        }
        
        virtual Int AmbDim() const override
        {
            return AMB_DIM;
        }

        virtual const Tensor2<Real,Int> & VertexCoordinates() const override
        {
            return V_coords;
        }
        
        virtual const Tensor2<Int,Int> & Simplices() const override
        {
            return simplices;
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
            return V_coords.Dimension(0);
        }
        
        virtual Int SimplexCount() const override
        {
            return simplices.Dimension(0);
        }
    
        virtual Int DofCount() const override
        {
            return VertexCount() * AMB_DIM;
        }
        
        virtual const Real * Dofs() const override
        {
            return V_coords.data();
        }
        
        virtual void SemiStaticUpdate( ptr<ExtReal> new_V_coords_, const bool transp_ = false ) const override
        {
            ptic(className()+"::SemiStaticUpdate");
            
            Tensor2<Real,Int> new_V_coords (VertexCount(),AmbDim());
            
            if( transp_ )
            {
                new_V_coords.ReadTransposed(new_V_coords_);
            }
            else
            {
                new_V_coords.Read(new_V_coords_);
            }
            
            this->ClearCache();
            
            Tensor2<Real,Int> P_near( SimplexCount(), NEAR_DIM );
            Tensor2<Real,Int> P_far ( SimplexCount(), FAR_DIM  );

            details.ComputeNearFarData( new_V_coords, simplices, P_near, P_far );
            
            GetClusterTree().SemiStaticUpdate( P_near, P_far );
            
            ptoc(className()+"::SemiStaticUpdate");
        }
       
        
        virtual void LoadUpdateVectors(
            ptr<ExtReal> vecs,
            const ExtReal max_time,
            const bool transp_
        ) const override
        {
            ptic(className()+"::LoadUpdateVectors");
            
            max_update_step_size = static_cast<SReal>(max_time);
            
            Tensor2<Real,Int> V_updates ( VertexCount(), AMB_DIM );

            if( vecs == nullptr )
            {
                V_updates.Fill(static_cast<Real>(0));
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
            
            const Int thread_count = GetClusterTree().ThreadCount();
            
            MovingPrimitive_T P_moving;
            
            Tensor2<SReal,Int> P_velocities_serialized ( SimplexCount(), P_moving.VelocitySize(), 0 );
            mut<SReal> P_v_ser = P_velocities_serialized.data();
            
            const Tensor1<Int,Int> & P_ordering = GetClusterTree().PrimitiveOrdering();

            JobPointers<Int> job_ptr ( SimplexCount(), thread_count );
            
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
                thread_count
            );
            
            // P_moving and P_velocities_serialized will be swapped against potentiall empty containers.
            GetClusterTree().TakeUpdateVectors(
                P_moving, P_velocities_serialized, static_cast<SReal>(max_time)
            );
            
            ptoc(className()+"::LoadUpdateVectors");
        }
        
        virtual ExtReal MaximumSafeStepSize(
            ptr<ExtReal> vecs,
            const ExtReal max_time,
            const ExtReal TOL = scalar_cast<ExtReal>(0.0625),
            const bool transp_ = false
        ) override
        {
            ptic(className()+"::MaximumSafeStepSize");
            
            LoadUpdateVectors( vecs, max_time, transp_ );
            
            ExtReal t = max_time;

            GetObstacle().LoadUpdateVectors( static_cast<ExtReal *>(nullptr), max_time );

            t = GetObstacleCollisionTree().MaximumSafeStepSize(t,TOL);

            logprint("GetObstacleCollisionTree().MaximumSafeStepSize(t) = "+ToString(t));
            
            t = GetCollisionTree().MaximumSafeStepSize(t,TOL);
            
            logprint("GetCollisionTree().MaximumSafeStepSize(t)         = "+ToString(t));
            
            ptoc(className()+"::MaximumSafeStepSize");
            
            return t;
        }
        
        virtual const ClusterTree_T & GetClusterTree() const override
        {
            static std::string tag ( "ClusterTree" );
            
            if( !this->InCacheQ( tag ) )
            {
                ptic(className()+"::GetClusterTree");
                if( (V_coords.Dimension(0) > 0) && (simplices.Dimension(0) > 0) )
                {
                    ptic("Allocations");
                    auto P_coords      = Tensor2<Real,Int> ( SimplexCount(), AMB_DIM, static_cast<Real>(0) );
                    auto P_hull_coords = Tensor3<Real,Int> ( SimplexCount(), DOM_DIM+1, AMB_DIM );
                    auto P_near        = Tensor2<Real,Int> ( SimplexCount(), NEAR_DIM );
                    auto P_far         = Tensor2<Real,Int> ( SimplexCount(), FAR_DIM );

                    Tensor2<SReal,Int> P_serialized ( SimplexCount(), P_proto.Size() );

                    auto DiffOp = SparseMatrix_T(
                        SimplexCount() * AMB_DIM,
                        VertexCount(),
                        SimplexCount() * AMB_DIM * (DOM_DIM+1),
                        ThreadCount()
                    );

                    DiffOp.Outer()[SimplexCount() * AMB_DIM] = SimplexCount() * AMB_DIM * (DOM_DIM+1);

                    auto AvOp = SparseMatrix_T(
                        SimplexCount(),
                        VertexCount(),
                        SimplexCount() * (DOM_DIM+1),
                        ThreadCount()
                    );

                    AvOp.Outer()[SimplexCount()] = SimplexCount() * (DOM_DIM+1);

                    ptoc("Allocations");

                    // What remains is to compute P_coords, P_hull_coords, P_near and P_far and the nonzero values of DiffOp.
                    details.ComputeNearFarDataOps( V_coords, simplices, P_coords, P_hull_coords, P_near, P_far, DiffOp, AvOp );


                    const JobPointers<Int> job_ptr ( SimplexCount(), ThreadCount() );

                    ptic("Creating primitives");
                    const Int thread_count = job_ptr.ThreadCount();

                    ParallelDo(
                        [&]( const Int thread )
                        {
                            Primitive_T P;

                            const Int i_begin = job_ptr[thread];
                            const Int i_end   = job_ptr[thread+1];

                            for( Int i = i_begin; i < i_end; ++i )
                            {
                                P.SetPointer( P_serialized.data(), i );
                                P.FromIndexList( V_coords.data(), simplices.data(), i );
                            }
                        },
                        thread_count
                    );

                    ptoc("Creating primitives");

                    ptic("Initializing cluster prototypes");

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

                    ptoc("Initializing cluster prototypes");

                    if( cluster_tree_settings.thread_count <= 0 )
                    {
                        cluster_tree_settings.thread_count = ThreadCount();
                    }

                    this->SetCache( tag,
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
                    eprint(ClassName()+"(V_coords.Dimension(0) <= 0) || (simplices.Dimension(0) <= 0)");

                    this->SetCache( tag, std::make_any<ClusterTree_T>() );
                }

                ptoc(className()+"::GetClusterTree");
            }

            return std::any_cast<ClusterTree_T &>( this->GetCache(tag) );
        }

        virtual const BlockClusterTree_T & GetBlockClusterTree() const override
        {
            static std::string tag ( "BlockClusterTree" );
            
            if( !this->InCacheQ( tag ) )
            {
                ptic(className()+"::GetBlockClusterTree");
                
                block_cluster_tree_settings.near_field_separation_parameter = adaptivity_settings.theta;
                block_cluster_tree_settings.near_field_intersection_parameter  = adaptivity_settings.intersection_theta;
                
                block_cluster_tree_settings.max_refinement = adaptivity_settings.max_refinement;
                
                this->SetCache( tag,
                    std::make_any<BlockClusterTree_T>(
                        GetClusterTree(),
                        GetClusterTree(),
                        block_cluster_tree_settings
                    )
                );
                
                ptoc(className()+"::GetBlockClusterTree");
            }
            
            return std::any_cast<BlockClusterTree_T &>( this->GetCache(tag) );
        }
        
        virtual const CollisionTree_T & GetCollisionTree() const override
        {
            static std::string tag ( "CollisionTree" );
            
            if( !this->InCacheQ( tag ) )
            {
                ptic(className()+"::GetCollisionTree");

                this->SetCache( tag,
                    std::make_any<CollisionTree_T>( GetClusterTree(), GetClusterTree() )
                );

                ptoc(className()+"::GetCollisionTree");
            }
            
            return std::any_cast<CollisionTree_T &>( this->GetCache(tag) );
        }
        
        const SparseBinaryMatrix_T & DerivativeAssembler() const override
        {
            
            static std::string tag ( "DerivativeAssembler" );
            
            if( !this->InCacheQ( tag ) )
            {
                ptic(className()+"::DerivativeAssembler");
            
                auto A = SparseBinaryMatrix_T(
                    SimplexCount() * (DOM_DIM+1),
                    VertexCount(),
                    SimplexCount() * (DOM_DIM+1),
                    ThreadCount()
                );
                
                A.Outer().iota();
                A.Inner().Read(simplices.data());
                
                this->SetCache( tag, std::any( std::move(A.Transpose()) ) );
                
                ptoc(className()+"::DerivativeAssembler");
            }
            
            return std::any_cast<SparseBinaryMatrix_T &>( this->GetCache(tag) );
            
        } // DerivativeAssembler
        
        void Assemble_ClusterTree_Derivatives( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_Derivatives");
            
            Tensor3<Real,Int> buffer ( SimplexCount(), DOM_DIM+1, AMB_DIM, static_cast<Real>(0) );
            
            GetClusterTree().CollectDerivatives();
            
            
            details.DNearToHulls( V_coords, simplices, GetClusterTree().PrimitiveDNearFieldData(), buffer, false );
            
            details.DFarToHulls ( V_coords, simplices, GetClusterTree().PrimitiveDFarFieldData(), buffer, true );

            DerivativeAssembler().template Dot<AMB_DIM>(
                static_cast<Real>(weight),   buffer.data(),
                static_cast<ExtReal>(addTo), output,
                AMB_DIM
            );

            ptoc(className()+"::Assemble_ClusterTree_Derivatives");
        }
        
        void Assemble_ClusterTree_Density( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_Density");
            
            GetClusterTree().CollectDensity( output, weight, addTo );
            
            ptoc(className()+"::Assemble_ClusterTree_Density");
        }
        
        void Assemble_ClusterTree_SimplexEnergies( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_SimplexEnergies");
            
            GetClusterTree().CollectPrimitiveEnergies( output, weight, addTo );
            
            ptoc(className()+"::Assemble_ClusterTree_SimplexEnergies");
        }

//##############################################################################################
//      Obstacle
//##############################################################################################

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

            this->SetCache(tag, std::any( obstacle ) );

        }
        
        const Obstacle_T & GetObstacle() const override
        {
            static std::string tag ("Obstacle");

            if( !this->InCacheQ(tag) )
            {
                wprint( ClassName()+"::GetObstacle: Obstacle not initialized.");
                
                // We have to construct an empy SimplicialMesh because Obstacle_T is abstract.
                this->SetCache(tag, std::any(std::make_shared<SimplicialMesh>()) );
            }
            
            return *std::any_cast<std::shared_ptr<Obstacle_T>>( this->GetCache(tag) );
        }
        
        virtual const ClusterTree_T & GetObstacleClusterTree() const override
        {
            return *dynamic_cast<const ClusterTree_T *>( & (GetObstacle().GetClusterTree()) );
        }
        
        virtual const ObstacleBlockClusterTree_T & GetObstacleBlockClusterTree() const override
        {
            static std::string tag ( "ObstacleBlockClusterTree" );
            
            if( !this->InCacheQ( tag ) )
            {
                ptic(className()+"::GetObstacleBlockClusterTree");
                
                this->SetCache( tag,
                    std::make_any<ObstacleBlockClusterTree_T>(
                        GetClusterTree(), GetObstacleClusterTree(), block_cluster_tree_settings
                    )
                );
                
                ptoc(className()+"::GetObstacleBlockClusterTree");

            }
            
            return std::any_cast<ObstacleBlockClusterTree_T &>( this->GetCache(tag) );
        }
        
        virtual const ObstacleCollisionTree_T & GetObstacleCollisionTree() const override
        {
            static std::string tag ( "ObstacleCollisionTree" );
            
            if( !this->InCacheQ(tag) )
            {
                ptic(className()+"::GetObstacleCollisionTree");
                
                this->SetCache( tag,
                    std::make_any<ObstacleCollisionTree_T>( GetClusterTree(), GetObstacleClusterTree() )
                );

                ptoc(className()+"::GetObstacleCollisionTree");
            }
            
            return std::any_cast<ObstacleCollisionTree_T &>( this->GetCache(tag) );
        }

        
//##############################################################################################
//      IO
//##############################################################################################
     
    public:
        
        virtual void WriteToFile( const std::string & file_name ) const override
        {
            ptic(ClassName()+"::WriteToFile");
            
            print("Writing mesh to file "+file_name+".");
            
            std::ofstream s (file_name);
            
            valprint("std::numeric_limits<Real>::digits",std::numeric_limits<Real>::digits10);
            
            s << std::setprecision( std::numeric_limits<Real>::digits );
            
            const Int vertex_count = VertexCount();
            const Int simplex_count = SimplexCount();
            
            constexpr Int simplex_size = DOM_DIM+1;
            
            ptr<Real> V = V_coords.data();
            ptr<Int>  S = simplices.data();
            
            s << "domain_dimension" << "\t" << DOM_DIM << "\n";
            s << "ambient_dimension" << "\t" << AMB_DIM << "\n";
            s << "vertex_count" << "\t" << vertex_count << "\n";
            s << "simplex_count" << "\t" << simplex_count << "\n";
            
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
                for( Int k = 0; k < simplex_size; ++k )
                {
                    s << S[simplex_size * i + k] << "\t";
                }
                s <<"\n";
            }
            
            ptoc(ClassName()+"::WriteToFile");
        }
        
//##############################################################################################
//      Remesher
//##############################################################################################
        
    public:
        
        virtual std::unique_ptr<RemesherBase_T> CreateRemesher() override
        {
            ptic(ClassName()+"::CreateRemesher");
            
            Real * null = nullptr;
            
            Remesher_T * R = new Remesher_T(
                VertexCoordinates().data(), VertexCoordinates().Dimension(0),
                Simplices().data(),         Simplices().Dimension(0),
                null,                       0,
                ThreadCount()
            );
            
            ptoc(ClassName()+"::CreateRemesher");
            
            return std::unique_ptr<RemesherBase_T>(R);
        }
        
//##############################################################################################
//      Standard interface
//##############################################################################################
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    
    private:
  
        static std::string className()
        {
            return "SimplicialMesh<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
    };
} // namespace Repulsor


#include "SimplicialMesh/SimplicialMesh_Factory.hpp"
