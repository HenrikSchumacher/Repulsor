#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_,typename ExtReal_,typename ExtInt_>
    class SimplicialRemesherBase;
    
    template<
        typename Real_, typename Int_, typename LInt_, 
        typename SReal_ = Real_, typename ExtReal_ = Real_
    >
    class SimplicialMeshBase : public CachedObject
    {
        static_assert(IntQ< Int_>,"");
        static_assert(IntQ<LInt_>,"");
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = LInt_;
        using SparseMatrix_T                 = Sparse::MatrixCSR      <Real,Int,LInt>;
        using SparseBinaryMatrix_T           = Sparse::BinaryMatrixCSR<     Int,LInt>;
        
        using ClusterTreeBase_T              =        ClusterTreeBase<Real,Int,LInt,SReal,ExtReal>;
        using BlockClusterTreeBase_T         =   BlockClusterTreeBase<Real,Int,LInt,SReal,ExtReal,true>;
        using CollisionTreeBase_T            =      CollisionTreeBase<Real,Int,LInt,SReal,ExtReal,true>;
        using ObstacleBlockClusterTreeBase_T =   BlockClusterTreeBase<Real,Int,LInt,SReal,ExtReal,false>;
        using ObstacleCollisionTreeBase_T    =      CollisionTreeBase<Real,Int,LInt,SReal,ExtReal,false>;
        using RemesherBase_T                 = SimplicialRemesherBase<Real,Int,Real,Int>;
        
        using TangentVector_T                = Tensor2<ExtReal,Int>;
        using CotangentVector_T              = Tensor2<ExtReal,Int>;
        
        // Default constructor
        SimplicialMeshBase() = default;
        // Destructor
        virtual ~SimplicialMeshBase() = default;
        // Copy constructor
        SimplicialMeshBase( const SimplicialMeshBase & other ) = default;
        // Copy assignment operator
        SimplicialMeshBase & operator=( const SimplicialMeshBase & other ) = default;
        // Move constructor
        SimplicialMeshBase( SimplicialMeshBase && other ) = default;
        // Move assignment operator
        SimplicialMeshBase & operator=( SimplicialMeshBase && other ) = default;
        
        
        explicit SimplicialMeshBase( const Int thread_count_ )
        :   thread_count( Max(Int(1),thread_count_) )
        {};
        
        mutable      ClusterTreeSettings       cluster_tree_settings;
        mutable BlockClusterTreeSettings block_cluster_tree_settings;
        mutable       AdaptivitySettings         adaptivity_settings;
        
    protected:
        
        const Int thread_count = 1;
            
    public:
        
        virtual Int DomDim() const = 0;
        
        virtual Int AmbDim() const = 0;

        virtual cref<Tensor2<Real,Int>> VertexCoordinates() const = 0;
        
        virtual cref<Tensor2<Int,Int>> Simplices() const = 0;

        virtual cref<Tensor1<Real,Int>> VertexCharges() const = 0;
        
        virtual cref<Tensor1<Real,Int>> SimplexCharges() const = 0;
        
        virtual Int FarDim() const = 0;
        
        virtual Int NearDim() const = 0;
        
        virtual Int VertexCount() const = 0;
        
        virtual Int SimplexCount() const = 0;
        
        virtual Int DofCount() const = 0;
        
        virtual cptr<Real> Dofs() const = 0;
        
        virtual Int ThreadCount() const
        {
            return thread_count;
        }
        
        virtual mref<Tensor2<Real,Int>> XBuffer( const Int nrhs ) const = 0;
        
        virtual mref<Tensor2<Real,Int>> YBuffer( const Int nrhs ) const = 0;
        
        virtual void SemiStaticUpdate( cptr<ExtReal> V_coords_, const bool transp_ = false ) const = 0;
        
        virtual cref<ClusterTreeBase_T> GetClusterTree() const = 0;
        
        virtual cref<BlockClusterTreeBase_T> GetBlockClusterTree() const = 0;
        
        virtual cref<CollisionTreeBase_T> GetCollisionTree() const = 0;
                
        virtual cref<SparseBinaryMatrix_T> DerivativeAssembler() const = 0;
        
        virtual void Assemble_ClusterTree_Derivatives(
            const ExtReal alpha, const ExtReal beta, mptr<ExtReal> Y, const Int ldY
        ) const = 0;

        virtual void Assemble_ClusterTree_SimplexDensities(
            mptr<ExtReal> output,
            const ExtReal weight,
            bool addTo = false
        ) const = 0;
        
        virtual void Assemble_ClusterTree_VertexDensities(
            mptr<ExtReal> output,
            const ExtReal weight,
            bool addTo = false
        ) const = 0;
        
        virtual void LoadUpdateVectors( 
            cptr<ExtReal> vecs, const ExtReal max_time, const bool transp_ = false
        ) const = 0;

        virtual ExtReal MaximumSafeStepSize(
            cptr<ExtReal> vecs,
            const ExtReal max_time,
            const ExtReal TOL = scalar_cast<ExtReal>(0.0625),
            const bool transp_ = false
        ) = 0;
        
        
        virtual cref<SparseMatrix_T> H1Metric() const = 0;
        
        virtual Real H1StiffnessWeight() const = 0;

        virtual void SetH1StiffnessWeight( const Real weight ) = 0;

        virtual Real H1MassWeight() const  = 0;

        virtual void SetH1MassWeight( const Real weight )  = 0;
        
        virtual cref<SparseMatrix_T> StiffnessMatrix() const = 0;
        
        virtual cref<SparseMatrix_T> MassMatrix() const = 0;
        
        virtual cref<Permutation<Int>> NestedDissectionOrdering( 
            const Int local_thread_count = 1
        ) const = 0;
        
        virtual cref<Permutation<Int>> ApproximateMinimumDegreeOrdering() const = 0;
        
        virtual cref<Permutation<Int>> MetisOrdering() const = 0;
        
        
        virtual mref<Sparse::CholeskyDecomposition<Real,Int,LInt>> H1Solver() const = 0;
        
        virtual Int H1SolverID() const = 0;

        virtual void SetH1SolverID( const Int ID ) const = 0;
        
        
        virtual mref<Sparse::CholeskyDecomposition<Real,Int,LInt>> MassSolver() const = 0;
        
        virtual Int MassSolverID() const = 0;

        virtual void SetMassSolverID( const Int ID ) const = 0;
        
//#######################################################################################
//      Obstacle
//#######################################################################################
        
    // TODO: Make list of obstacles.
        
    public:
        
        virtual void  LoadObstacle( std::unique_ptr<SimplicialMeshBase> obstacle ) = 0;

        virtual cref<SimplicialMeshBase> GetObstacle() const = 0;
        
//        virtual bool  ObstacleInitialized() const = 0;
        
        virtual cref<ClusterTreeBase_T> GetObstacleClusterTree() const = 0;
        
        virtual cref<ObstacleBlockClusterTreeBase_T> GetObstacleBlockClusterTree() const = 0;
        
        virtual cref<ObstacleCollisionTreeBase_T> GetObstacleCollisionTree() const = 0;
            
//#################################################################################
//      IO
//#################################################################################
     
    public:
        
        virtual void WriteToFile( const std::filesystem::path & file ) const = 0;
        
//#################################################################################
//      Remesher
//#################################################################################
        
    public:
        
        virtual std::unique_ptr<RemesherBase_T> CreateRemesher() = 0;
        
//#################################################################################
//      Standard interface
//#################################################################################
        
    public:

        static std::string className()
        {
            return std::string("SimplicialMeshBase")+"<"+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    }; // class SimplicialMeshBase
    
} // namespace Repulsor
