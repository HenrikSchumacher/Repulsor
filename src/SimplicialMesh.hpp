#pragma once

#include "SimplicialMesh/SimplicialMeshBase.hpp"

namespace Repulsor
{
//    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
//    class TangentPoint;
    
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class SimplicialRemesherBase;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class SimplicialRemesher;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class SimplicialMesh : public SimplicialMeshBase<Real_,Int_,SReal_,ExtReal_>
    {
    public :
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        
        using Base_T  = SimplicialMeshBase<Real,Int,SReal,ExtReal>;
        
        using TangentVector_T    = typename Base_T::TangentVector_T;
        using CotangentVector_T  = typename Base_T::CotangentVector_T;
        
        using       Primitive_T  =        Polytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using MovingPrimitive_T  =  MovingPolytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using BoundingVolume_T   =           AABB<           AMB_DIM,GJK_Real,Int,SReal>;
        
        
        using ClusterTree_T              =            ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        using BlockClusterTree_T         =       BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal,true>;
        using ObstacleBlockClusterTree_T =       BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal,false>;
        using CollisionTree_T            =          CollisionTree<AMB_DIM,Real,Int,SReal,ExtReal,true>;
        using ObstacleCollisionTree_T    =          CollisionTree<AMB_DIM,Real,Int,SReal,ExtReal,false>;

        using Remesher_T         =     SimplicialRemesher<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>;
        using RemesherBase_T     = SimplicialRemesherBase<                Real,Int,SReal,ExtReal>;
        using Obstacle_T         = Base_T;
        
        SimplicialMesh() = default;

        SimplicialMesh(
            const Tensor2<Real,Int> & V_coords_,
            // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Tensor2<Int,Int> & simplices_,
            // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const long long thread_count_ = 1
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
            const long long vertex_count_,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const long long simplex_count_,
            const long long thread_count_ = 1
        )
        :   SimplicialMesh( V_coords_, vertex_count_, false, simplices_, simplex_count_, false, thread_count_)
        {}
        
        template<typename ExtReal_2, typename ExtInt>
        SimplicialMesh(
            const ExtReal_2 * vertex_coords_, // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const long long vertex_count_,
            const bool vertex_coords_transpose,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const long long simplex_count_,
            const bool simplices_transpose,
            const long long thread_count_ = 1
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
        
        mutable bool cluster_tree_initialized = false;
        mutable std::unique_ptr<ClusterTree_T> cluster_tree
              = std::unique_ptr<ClusterTree_T>( new ClusterTree_T() );
        
        mutable bool block_cluster_tree_initialized = false;
        mutable std::unique_ptr<BlockClusterTree_T> block_cluster_tree
              = std::unique_ptr<BlockClusterTree_T>( new BlockClusterTree_T(*cluster_tree, *cluster_tree) );
        
        mutable bool collision_tree_initialized = false;
        mutable std::unique_ptr<CollisionTree_T> collision_tree
              = std::unique_ptr<CollisionTree_T>( new CollisionTree_T(*cluster_tree, *cluster_tree ) );
        
        mutable SReal max_update_step_size = 0;
        
        Tensor2<Real,Int> V_coords;
        Tensor2<Int,Int>  simplices;
        
        mutable bool derivative_assembler_initialized = false;
        mutable SparseBinaryMatrixCSR<Int,Int> derivative_assembler;
        
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

        virtual Int FarDim() const override
        {
            return 1 + AMB_DIM + (AMB_DIM * (AMB_DIM + 1))/2;
        }
        
        virtual Int NearDim() const override
        {
            return 1 + (DOM_DIM+1)*AMB_DIM + (AMB_DIM * (AMB_DIM + 1))/2;
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
        
        virtual void SemiStaticUpdate( const ExtReal * restrict const new_V_coords_ ) override
        {
            ptic(className()+"::SemiStaticUpdate");
            
            V_coords.Read(new_V_coords_);
            
            this->ClearCache();
            
            Tensor2<Real,Int> P_near( SimplexCount(), NearDim() );
            Tensor2<Real,Int> P_far ( SimplexCount(), FarDim()  );

            details.ComputeNearFarData( V_coords, simplices, P_near, P_far );
            
            cluster_tree->SemiStaticUpdate( P_near, P_far );
            
            ptoc(className()+"::SemiStaticUpdate");
        }
       
        
        virtual void LoadUpdateVectors( const ExtReal * const vecs, const ExtReal max_time ) override
        {
            ptic(className()+"::LoadUpdateVectors");
            
            max_update_step_size = static_cast<SReal>(max_time);
            collision_tree_initialized = false;
            
            auto V_updates = Tensor2<Real,Int>( VertexCount(), AMB_DIM );
            if( vecs == nullptr )
            {
                V_updates.Fill(static_cast<Real>(0));
            }
            else
            {
                V_updates.Read(vecs);
            }
            // ATTENTION: We reorder already here outside of the cluster tree to save a copy operation of a big Tensor2!
            
            const Int thread_count = GetClusterTree().ThreadCount();
            
            MovingPrimitive_T P_moving;
            
            Tensor2<SReal,Int> P_velocities_serialized ( SimplexCount(), P_moving.VelocitySize(), 0 );
            SReal * restrict const P_v_ser = P_velocities_serialized.data();
            
            const Tensor1<Int,Int> P_ordering = GetClusterTree().PrimitiveOrdering();

            JobPointers<Int> job_ptr ( SimplexCount(), thread_count );
            
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
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
            }
            
            // P_moving and P_velocities_serialized will be swapped against potentiall empty containers.
            GetClusterTree().TakeUpdateVectors(
                P_moving, P_velocities_serialized, static_cast<SReal>(max_time)
            );
            
            ptoc(className()+"::LoadUpdateVectors");
        }
        
        virtual ExtReal MaximumSafeStepSize( const ExtReal * restrict const vecs, const ExtReal max_time ) override
        {
            ptic(className()+"::MaximumSafeStepSize");
            
            LoadUpdateVectors( vecs, max_time );
            
            ExtReal t = max_time;

            if( obstacle_initialized )
            {
                obstacle->LoadUpdateVectors( static_cast<ExtReal *>(nullptr), max_time );

                t = GetObstacleCollisionTree().MaximumSafeStepSize(t);

                logprint("GetObstacleCollisionTree().MaximumSafeStepSize(t) = "+ToString(t));
            }
            
            t = GetCollisionTree().MaximumSafeStepSize(t);
            
            logprint("GetCollisionTree().MaximumSafeStepSize(t)         = "+ToString(t));
            
            ptoc(className()+"::MaximumSafeStepSize");
            
            return t;
        }
        
        virtual const ClusterTree_T & GetClusterTree() const override
        {
            if( !cluster_tree_initialized )
            {
                ptic(className()+"::GetClusterTree");
                if( (V_coords.Dimension(0) > 0) && (simplices.Dimension(0) > 0) )
                {
                    ptic("Allocations");
                    auto P_coords      = Tensor2<Real,Int> ( SimplexCount(), AMB_DIM, static_cast<Real>(0) );
                    auto P_hull_coords = Tensor3<Real,Int> ( SimplexCount(), DOM_DIM+1, AMB_DIM );
                    auto P_near        = Tensor2<Real,Int> ( SimplexCount(), NearDim() );
                    auto P_far         = Tensor2<Real,Int> ( SimplexCount(), FarDim() );
                    
                    Tensor2<SReal,Int> P_serialized ( SimplexCount(), P_proto.Size() );
                    
                    auto DiffOp = SparseMatrixCSR<Real,Int,Int>(
                        SimplexCount() * AMB_DIM,
                        VertexCount(),
                        SimplexCount() * AMB_DIM * (DOM_DIM+1),
                        ThreadCount()
                    );
        
                    DiffOp.Outer()[SimplexCount() * AMB_DIM] = SimplexCount() * AMB_DIM * (DOM_DIM+1);
                    
                    auto AvOp = SparseMatrixCSR<Real,Int,Int>(
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
                    const Int thread_count = job_ptr.Size()-1;
                    
                    #pragma omp parallel for num_threads(thread_count)
                    for( Int thread = 0; thread < thread_count; ++thread )
                    {
                        Primitive_T P;
                        
                        Int i_begin = job_ptr[thread];
                        Int i_end   = job_ptr[thread+1];

                        for( Int i = i_begin; i < i_end; ++i )
                        {
                            P.SetPointer( P_serialized.data(), i );
                            P.FromIndexList( V_coords.data(), simplices.data(), i );
                        }
                    }
                    ptoc("Creating primitives");
                    
                    ptic("Initializing cluster prototypes");
                    
                    std::unique_ptr<BoundingVolume_T> C_proto;

                    switch( cluster_tree_settings.bounding_volume_type )
                    {
                        case BoundingVolumeType::AABB_MedianSplit:
                        {
                            logprint("Using AABB_MedianSplit as bounding volume.");
                            C_proto = std::unique_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_MedianSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                            break;
                        }
                        case BoundingVolumeType::AABB_PreorderedSplit:
                        {
                            logprint("Using AABB_PreorderedSplit as bounding volume.");
                            C_proto = std::unique_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_PreorderedSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                            break;
                        }
                        default:
                        {
                            logprint("Using AABB_LongestAxisSplit as bounding volume.");
                            C_proto = std::unique_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_LongestAxisSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                        }
                    }

                    ptoc("Initializing cluster prototypes");
                    
                    if( cluster_tree_settings.thread_count <= 0 )
                    {
                        cluster_tree_settings.thread_count = ThreadCount();
                    }
                    
                    cluster_tree = std::make_unique<ClusterTree_T>(
                        P_proto, P_serialized, *C_proto, Tensor1<Int,Int>(0),
                        P_near, P_far,
                        DiffOp, AvOp,
                        cluster_tree_settings
                    );
                    
                    cluster_tree_initialized = true;
                }
                else
                {
                    cluster_tree = std::make_unique<ClusterTree_T>();
                }
                
                ptoc(className()+"::GetClusterTree");
            }
            return *cluster_tree;
        }

        virtual const BlockClusterTree_T & GetBlockClusterTree() const override
        {
            if( !block_cluster_tree_initialized )
            {
                ptic(className()+"::GetBlockClusterTree");
                
                block_cluster_tree_settings.near_field_separation_parameter = adaptivity_settings.theta;
                block_cluster_tree_settings.near_field_intersection_parameter  = adaptivity_settings.intersection_theta;
                
                block_cluster_tree = std::make_unique<BlockClusterTree_T>(
                    GetClusterTree(),
                    GetClusterTree(),
                    block_cluster_tree_settings
                );
                
                block_cluster_tree_initialized = true;
                
                ptoc(className()+"::GetBlockClusterTree");
            }
            return *block_cluster_tree;
        }
        
        virtual const CollisionTree_T & GetCollisionTree() const override
        {
            if( !collision_tree_initialized )
            {
                ptic(className()+"::GetCollisionTree");
                
                collision_tree = std::make_unique<CollisionTree_T>( GetClusterTree(), GetClusterTree() );

                collision_tree_initialized = true;

                ptoc(className()+"::GetCollisionTree");
            }
            return *collision_tree;
        }
        
        
        const SparseBinaryMatrixCSR<Int,Int> & DerivativeAssembler() const override
        {
            if( !derivative_assembler_initialized )
            {
                ptic(className()+"::DerivativeAssembler");
                
                auto A = SparseBinaryMatrixCSR<Int,Int>(
                    SimplexCount() * (DOM_DIM+1),
                    VertexCount(),
                    SimplexCount() * (DOM_DIM+1),
                    ThreadCount()
                );
                
                A.Outer().iota();
                A.Inner().Read(simplices.data());
                
                derivative_assembler = A.Transpose();
                
                derivative_assembler_initialized = true;
                
                ptoc(className()+"::DerivativeAssembler");
            }
            return derivative_assembler;
            
        } // DerivativeAssembler
        
        void Assemble_ClusterTree_Derivatives( ExtReal * output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_Derivatives");
            
            Tensor3<Real,Int> buffer ( SimplexCount(), DOM_DIM+1, AMB_DIM, static_cast<Real>(0) );
            
            GetClusterTree().CollectDerivatives();
            
            
            details.DNearToHulls( V_coords, simplices, GetClusterTree().PrimitiveDNearFieldData(), buffer, false );
            
            details.DFarToHulls ( V_coords, simplices, GetClusterTree().PrimitiveDFarFieldData(), buffer, true );

            DerivativeAssembler().Dot(
                static_cast<Real>(weight),   buffer.data(),
                static_cast<ExtReal>(addTo), output,
                AMB_DIM
            );
             
            ptoc(className()+"::Assemble_ClusterTree_Derivatives");
        }
        
        void Assemble_ClusterTree_Density( ExtReal * output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_Density");
            
            GetClusterTree().CollectDensity( output, weight, addTo );
            
            ptoc(className()+"::Assemble_ClusterTree_Density");
        }
        
        void Assemble_ClusterTree_SimplexEnergies( ExtReal * output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_SimplexEnergies");
            
            GetClusterTree().CollectPrimitiveEnergies( output, weight, addTo );
            
            ptoc(className()+"::Assemble_ClusterTree_SimplexEnergies");
        }

//##############################################################################################
//      Obstacle
//##############################################################################################
    
        
        mutable bool obstacle_initialized = false;
        mutable std::unique_ptr<Obstacle_T> obstacle = nullptr;
        
        mutable bool obstacle_block_cluster_tree_initialized = false;
        mutable std::unique_ptr<ObstacleBlockClusterTree_T> obstacle_block_cluster_tree
              = std::unique_ptr<ObstacleBlockClusterTree_T>(
                    new ObstacleBlockClusterTree_T(*cluster_tree, *cluster_tree) );
        
        
        mutable bool obstacle_collision_tree_initialized = false;
        mutable std::unique_ptr<ObstacleCollisionTree_T> obstacle_collision_tree
              = std::unique_ptr<ObstacleCollisionTree_T>(
                    new ObstacleCollisionTree_T(*cluster_tree, *cluster_tree) );
        
        bool ObstacleInitialized() const override
        {
            return obstacle != nullptr;
        }
        
        virtual void LoadObstacle( std::unique_ptr<Obstacle_T> obstacle_ ) override
        {
            // Input obstacle is moved.
            if( obstacle_->AmbDim() != AmbDim() )
            {
                eprint(className()+"::LoadObstacle: Attempted to load obstacle of ambient dimension "+ToString(obstacle_->AmbDim())+" into mesh of ambient dimension "+ToString(AmbDim())+". Setting obstacle to nullptr."
                );
                obstacle = nullptr;
                return;
            }
            
            obstacle = std::move(obstacle_);
            
            obstacle_block_cluster_tree = nullptr;
            obstacle_block_cluster_tree_initialized = false;
            
            obstacle->cluster_tree_settings = cluster_tree_settings;
            
            obstacle_initialized = true;
        }
        
        const Obstacle_T & GetObstacle() const override
        {
            if( !obstacle_initialized )
            {
                wprint(className()+"::GetObstacle: Obstacle not initialized.");
                obstacle = std::unique_ptr<Base_T>(new SimplicialMesh());
            }
            return *obstacle;
        }
        
        virtual const ClusterTree_T & GetObstacleClusterTree() const override
        {
            return *dynamic_cast<const ClusterTree_T *>( & (GetObstacle().GetClusterTree()) );
        }
        
        virtual const ObstacleBlockClusterTree_T & GetObstacleBlockClusterTree() const override
        {
            if( !obstacle_block_cluster_tree_initialized )
            {
                ptic(className()+"::GetObstacleBlockClusterTree");
                if( obstacle_initialized )
                {
                    obstacle_block_cluster_tree = std::make_unique<ObstacleBlockClusterTree_T>(
                        GetClusterTree(), GetObstacleClusterTree(), block_cluster_tree_settings
                    );
                    obstacle_block_cluster_tree_initialized = true;
                }
                else
                {
                    wprint(className()+"::GetObstacleBlockClusterTree: Obstacle not initialized.");
                }
                ptoc(className()+"::GetObstacleBlockClusterTree");
                
            }
            return *obstacle_block_cluster_tree;
        }
        
        virtual const ObstacleCollisionTree_T & GetObstacleCollisionTree() const override
        {
            if( !obstacle_collision_tree_initialized )
            {
                ptic(className()+"::GetObstacleCollisionTree");
                
                if( obstacle_initialized )
                {
                    obstacle_collision_tree = std::make_unique<ObstacleCollisionTree_T>(
                            GetClusterTree(), GetObstacleClusterTree()
                    );

                    obstacle_collision_tree_initialized = true;
                }
                else
                {
                    print(ClassName()+"GetObstacleCollisionTree: Obstacle not initialized.");
                }

                ptoc(className()+"::GetObstacleCollisionTree");
            }
            return *obstacle_collision_tree;
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
            
            const Real * restrict const V = V_coords.data();
            const Int  * restrict const S = simplices.data();
            
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
            
            Remesher_T * R = new Remesher_T(*this);
            
            ptoc(ClassName()+"::CreateRemesher");
            
            return std::unique_ptr<RemesherBase_T>(R);
        }
        
//##############################################################################################
//      Standard interface
//##############################################################################################
 
    public:
        
//        virtual Real TangentPointEnergy_New() const override
//        {
//            ptic(ClassName()+"::TangentPointEnergy_New");
//            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false,false> T (
//                GetBlockClusterTree(),
//                GetTangentPointExponents().first,
//                GetTangentPointExponents().second
//            );
//            
//            const Real result = T.Compute();
//            
//            ptoc(ClassName()+"::TangentPointEnergy_New");
//            
//            return result;
//        }
//        
//
//        virtual Real TangentPointEnergy_New_Differential( ExtReal * output, bool addTo = false ) const  override
//        {
//            ptic(ClassName()+"::TangentPointEnergy_New_Differential");
//            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,true,false,false> T (
//                GetBlockClusterTree(),
//                GetTangentPointExponents ().first,
//                GetTangentPointExponents ().second
//            );
//            
//            wprint(ClassName()+"::TangentPointEnergy_New_Differential is not fully implemented.");
//            
//            const Real result = T.Compute();
//            
//            ptoc(ClassName()+"::TangentPointEnergy_New_Differential");
//            
//            return result;
//        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    
    private:
  
        static std::string className()
        {
            return "SimplicialMesh<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
    };
} // namespace Repulsor


#include "SimplicialMesh/SimplicialMesh_Factory.hpp"
