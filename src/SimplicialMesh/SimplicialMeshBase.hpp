#pragma once

#define CLASS SimplicialMeshBase


//#define ENABLE_ENERGIES
//#define ENABLE_METRICS

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class EnergyBase;
    
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class SimplicialRemesherBase;
    
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        ASSERT_FLOAT(SReal);
        ASSERT_FLOAT(ExtReal);
        
        
    public:
        
        using ClusterTree_T              =        ClusterTreeBase<Real,Int,SReal,ExtReal>;
        using BlockClusterTree_T         =   BlockClusterTreeBase<Real,Int,SReal,ExtReal,true>;
        using CollisionTree_T            =      CollisionTreeBase<Real,Int,SReal,ExtReal,true>;
        using ObstacleBlockClusterTree_T =   BlockClusterTreeBase<Real,Int,SReal,ExtReal,false>;
        using ObstacleCollisionTree_T    =      CollisionTreeBase<Real,Int,SReal,ExtReal,false>;
        using Remesher_T                 = SimplicialRemesherBase<Real,Int,SReal,ExtReal>;
        
        CLASS () {}
        
        explicit CLASS( const Int thread_count_ )
        :   thread_count( std::max( static_cast<Int>(1), thread_count_) )
        {};
        
        virtual ~CLASS() = default;

        mutable      ClusterTreeSettings       cluster_tree_settings;
        mutable BlockClusterTreeSettings block_cluster_tree_settings;
        mutable       AdaptivitySettings         adaptivity_settings;
        
    protected:
        
        const Int thread_count = 1;
        
        mutable std::unordered_map<std::string,std::any> cache;
        mutable std::unordered_map<std::string,std::any> persistent_cache;
            
    public:
        
        virtual Int DomDim() const = 0;
        
        virtual Int AmbDim() const = 0;

        virtual const Tensor2<Real,Int> & VertexCoordinates() const = 0;
        
        virtual const Tensor2<Int,Int> & Simplices() const = 0;

        virtual Int FarDim() const = 0;
        
        virtual Int NearDim() const = 0;
        
        virtual Int VertexCount() const = 0;
        
        virtual Int SimplexCount() const = 0;
        
        virtual Int DofCount() const = 0;
        
        virtual const Real * Dofs() const = 0;
        
        virtual Int ThreadCount() const
        {
            return thread_count;
        }
        
        virtual void SemiStaticUpdate( const ExtReal * restrict const V_coords_ ) = 0;
        
        virtual const ClusterTree_T & GetClusterTree() const = 0;
        
        virtual const BlockClusterTree_T & GetBlockClusterTree() const = 0;
        
        virtual const CollisionTree_T & GetCollisionTree() const = 0;
                
        virtual const SparseBinaryMatrixCSR<Int> & DerivativeAssembler() const = 0;
        
        virtual void Assemble_ClusterTree_Derivatives(
            ExtReal * output,
            const ExtReal weight,
            bool addTo = false
        ) const = 0;

        virtual void Assemble_ClusterTree_SimplexEnergies(
            ExtReal * output,
            const ExtReal weight,
            bool addTo = false
        ) const = 0;
        
        virtual void Assemble_ClusterTree_Density(
            ExtReal * output,
            const ExtReal weight,
            bool addTo = false
        ) const = 0;
        
        virtual void LoadUpdateVectors( const ExtReal * restrict const vecs, const ExtReal max_time ) = 0;

        virtual ExtReal MaximumSafeStepSize( const ExtReal * restrict const vecs, const ExtReal max_time ) = 0;
        
//#######################################################################################
//      Obstacle
//#######################################################################################
        
    public:
        
        virtual void LoadObstacle( std::unique_ptr<CLASS> obstacle_ ) = 0;

        virtual const CLASS & GetObstacle() const = 0;
        
        virtual bool ObstacleInitialized() const = 0;
        
        virtual const ClusterTree_T & GetObstacleClusterTree() const = 0;
        
        virtual const ObstacleBlockClusterTree_T & GetObstacleBlockClusterTree() const = 0;
        
        virtual const ObstacleCollisionTree_T & GetObstacleCollisionTree() const = 0;
        
//##############################################################################################
//      Tangent-point
//##############################################################################################

    public:
        
        virtual std::pair<Real,Real> GetTangentPointExponents() const = 0;
        
        virtual void SetTangentPointExponents( const Real alpha, const Real beta ) const = 0;
        
        virtual ExtReal GetTangentPointWeight() const = 0;
        
        virtual void SetTangentPointWeight( const ExtReal weight ) const = 0;
        
//##############################################################################################
//      TangentPointEnergy
//##############################################################################################

#ifdef ENABLE_ENERGIES
        
        virtual ExtReal TangentPointEnergy() const = 0;

        virtual ExtReal TangentPointEnergy_Differential( ExtReal * output, bool addTo = false ) const  = 0;
        
        virtual ExtReal TangentPointEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const  = 0;
        
        virtual ExtReal TangentPointEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const  = 0;
        
        virtual void TangentPointEnergy_Density( ExtReal * output, bool addTo = false ) const  = 0;
        
        virtual void TangentPointEnergy_Density( Tensor1<ExtReal,Int> & output, bool addTo = false ) const  = 0;
        
        virtual void TangentPointEnergy_SimplexEnergies( ExtReal * output, bool addTo = false ) const  = 0;
        
        virtual void TangentPointEnergy_SimplexEnergies( Tensor1<ExtReal,Int> & output, bool addTo = false ) const  = 0;
        
#endif
        
//##############################################################################################
//      TangentPointObstacleEnergy
//##############################################################################################

#ifdef ENABLE_ENERGIES
        
        virtual ExtReal TangentPointObstacleEnergy() const = 0;

        virtual ExtReal TangentPointObstacleEnergy_Differential( ExtReal * output, bool addTo = false ) const  = 0;
        
        virtual ExtReal TangentPointObstacleEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const = 0;
        
        virtual ExtReal TangentPointObstacleEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const = 0;
        
        virtual void TangentPointObstacleEnergy_Density( ExtReal * output, bool addTo = false ) const  = 0;
        
        virtual void TangentPointObstacleEnergy_Density( Tensor1<ExtReal,Int> & output, bool addTo = false ) const = 0;
   
        virtual void TangentPointObstacleEnergy_SimplexEnergies( ExtReal * output, bool addTo = false ) const  = 0;
        
        virtual void TangentPointObstacleEnergy_SimplexEnergies( Tensor1<ExtReal,Int> & output, bool addTo = false ) const  = 0;
        
#endif
        
        
//##############################################################################################
//      Trivial energy (for debugging purposes)
//##############################################################################################
        
#ifdef ENABLE_ENERGIES
        
        virtual ExtReal GetTrivialEnergyWeight() const = 0;

        virtual void SetTrivialEnergyWeight( const ExtReal weight ) const = 0;

        virtual ExtReal TrivialEnergy() const = 0;

        virtual ExtReal TrivialEnergy_Differential( ExtReal * output, bool addTo = false ) const = 0;

        virtual ExtReal TrivialEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const = 0;
        
        virtual ExtReal TrivialEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const = 0;
        
#endif
        
//##############################################################################################
//      TrivialObstacleEnergy (for debugging purposes)
//##############################################################################################
         
#ifdef ENABLE_ENERGIES
        
        virtual ExtReal TrivialObstacleEnergy() const = 0;

        virtual ExtReal TrivialObstacleEnergy_Differential( ExtReal * output, bool addTo = false ) const = 0;
        
        virtual ExtReal TrivialObstacleEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const = 0;
        
        virtual ExtReal TrivialObstacleEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const = 0;

#endif
        
//##############################################################################################
//      TangentPointMetric
//##############################################################################################

#ifdef ENABLE_METRICS
        
        virtual void TangentPointMetric_Multiply(
            const ExtReal alpha, const ExtReal * U,
            const ExtReal  beta,       ExtReal * V,
            Int cols
        ) const  = 0;
        
        virtual void TangentPointMetric_Multiply(
            const ExtReal alpha, const Tensor1<ExtReal,Int> & U,
            const ExtReal  beta,       Tensor1<ExtReal,Int> & V
        ) const = 0;
        
        virtual void TangentPointMetric_Multiply(
            const ExtReal alpha, const Tensor2<ExtReal,Int> & U,
            const ExtReal  beta,       Tensor2<ExtReal,Int> & V
        ) const = 0;
        
        virtual void TangentPointMetric_Multiply(
            const ExtReal alpha, const ExtReal * U,
            const ExtReal  beta,       ExtReal * V,
            Int cols,
            KernelType kernel
        ) const = 0;
        
        virtual void TangentPointMetric_Multiply(
            const ExtReal alpha, const Tensor1<ExtReal,Int> & U,
            const ExtReal  beta,       Tensor1<ExtReal,Int> & V,
            KernelType kernel
        ) const = 0;
        
        virtual void TangentPointMetric_Multiply(
            const ExtReal alpha, const Tensor2<ExtReal,Int> & U,
            const ExtReal  beta,       Tensor2<ExtReal,Int> & V,
            KernelType kernel
        ) const = 0;
        
        
        virtual const Tensor1<Real,Int> & TangentPointMetric_Values(
            const bool farQ,
            const KernelType kernel
        ) const = 0;
        
        virtual void TangentPointMetric_ApplyKernel(
            const bool farQ,
            const KernelType kernel
        ) const = 0;
        
#endif
        
////##############################################################################################
////      TangentPointSingularMetric
////##############################################################################################
//
//        
//        virtual void TangentPointSingularMetric_Multiply(
//            const ExtReal alpha, const ExtReal * U,
//            const ExtReal  beta,       ExtReal * V,
//            Int cols
//        ) const  = 0;
//        
//        virtual void TangentPointSingularMetric_Multiply(
//            const ExtReal alpha, const Tensor1<ExtReal,Int> & U,
//            const ExtReal  beta,       Tensor1<ExtReal,Int> & V
//        ) const = 0;
//        
//        virtual void TangentPointSingularMetric_Multiply(
//            const ExtReal alpha, const Tensor2<ExtReal,Int> & U,
//            const ExtReal  beta,       Tensor2<ExtReal,Int> & V
//        ) const = 0;
//        
//        virtual const Tensor1<Real,Int> & TangentPointSingularMetric_Values(
//            const bool farQ
//        ) const = 0;
//        
//        virtual void TangentPointSingularMetric_ApplyKernel(
//            const bool farQ
//        ) const = 0;
        
//##############################################################################################
//      IO
//##############################################################################################
     
    public:
        
        virtual void WriteToFile( const std::string & file_name ) const = 0;
        
//##############################################################################################
//      Remesher
//##############################################################################################
        
    public:
        
        virtual std::unique_ptr<Remesher_T> CreateRemesher() = 0;
//        
//        virtual Real TangentPointEnergy_New() const = 0;
//        
//        virtual Real TangentPointEnergy_New_Differential( ExtReal * output, bool addTo = false ) const = 0;
        
//##############################################################################################
//      Standard interface
//##############################################################################################
        
    public:
        
        bool IsCached( const std::string & s ) const
        {
            #pragma omp critical (cache)
            {
                return static_cast<bool>( cache.count(s) );
            }
        }
        
        std::any & GetCache( const std::string & s ) const
        {
            #pragma omp critical (cache)
            {
                try{
                  return cache.at(s);
                }
                catch( const std::out_of_range & e )
                {
                    eprint(ClassName()+"GetCache: Key \""+s+"\" not found!.");
                    throw; //an internal catch block forwards the exception to its external level
                }
            }
        }
        
        void SetCache( const std::string & s, std::any & thing ) const
        {
            #pragma omp critical (cache)
            {
                cache[s] = thing;
            }
        }
        
        void ClearCache() const
        {
            cache = std::unordered_map<std::string,std::any>();
        }
        
                                                          
        virtual CLASS & DownCast() = 0;

        virtual const CLASS & DownCast() const = 0;
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
    };
    
} // namespace Repulsor

#undef CLASS
