#include <iostream>

//#define NDEBUG

/// enable profiler
#define TOOLS_ENABLE_PROFILER

#ifdef __APPLE__
/// Use these while on a mac. Don't forget to issue the compiler flag `-framework Accelerate`.
///
    #include "../submodules/Tensors/Accelerate.hpp"
#else
/// This should work for OpenBLAS.
    #include "../submodules/Tensors/OpenBLAS.hpp"
#endif

#include "../Repulsor.hpp"

using Int     = int;
using LInt    = std::size_t;
using Real    = double;

using namespace Repulsor;
using namespace Tensors;
using namespace Tools;


int main( void )
{
    print("");
    print("###############################################################");
    print("###      Example program for Repulsor::SimplicialMesh      ####");
    print("###############################################################");
    print("");
    
    /// Set up profiler to write to `~/Tools_Profile.tsv` and `~/Tools_Log.txt`
    Profiler::Clear();

    int thread_count = 8;
    
    
    /// Setting up the scalar types to be used by the meshes:
    using Mesh_T = SimplicialMeshBase<Real,Int,LInt>;
    
    /// Meaning of template parameters:
    ///
    /// `SimplicialMeshBase<Real, Int, SReal, ExtReal>`;
    ///
    /// `Real`      Floating point type used for computations regarding energy and metric.
    ///             Typically, you want to use `Real = double`.
    ///
    /// `Int`       Signed integer type used internally, in particular as indices
    ///             for sparse arrays.
    ///             `Int = int`       corresponds to MKL's  LP64 mode
    ///             `Int = long long` corresponds to MKL's ILP64 mode
    ///             I recommend using `Int = int`.
    ///
    /// `LInt`      Integer type used for the row pointers of a sparse array.
    ///             I recommend using `LInt = std::size_t`.

    /// `SimplicialMesh_Factory<Mesh_T,n_min,n_max,m_min,m_max>`
    /// is a factory that can create meshes with base type `Mesh_T` and domain dimension from `n_min` to `n_max` and with ambient dimension from `m_min` to `m_max`.
    /// The main purpose of this factory is to tell the compiler which template instantiations of `SimplicialMesh` it shall generate. So this factory is basically about saving compile time.
    /// The narrower these intervals are, the less code has to compiled.
    /// The wieder these intervals are, the more flexible the runtime polumorphism will be.
    /// Note that with obstacles you might want to use meshes of various dimensions,
    /// so you have to use an according factory.

    /// Creating a factory type whose instances can create instances of `SimplicialMesh` with domain dimension in the range `2,...,2` and ambient dimension in the range `3,...,3`.
    
    using Mesh_Factory_T = SimplicialMesh_Factory<Mesh_T,2,2,3,3>;
    
    /// We need also factories for the energy and the metric.
    using Energy_T       = EnergyBase<Mesh_T>;
    using Metric_T       = MetricBase<Mesh_T>;

    
    
    /// Instantiate factory.
    Mesh_Factory_T mesh_factory;

#include "../Meshes/TwoSpheres.hpp"
    
    tic("Initializing mesh");
    
    /// Use`mesh_factory` to initialize a mesh from arrays that are already in memory.
    std::unique_ptr<Mesh_T> M_ptr = mesh_factory.Make(
        &vertex_coordinates[0][0],  vertex_count, amb_dim,   false,
        &simplices[0][0],          simplex_count, dom_dim+1, false,
        thread_count
    );

    auto & M = *M_ptr;  // I don't like pointers. Give me a reference.

    toc("Initializing mesh");

    /// Setup for bounding volume hierarchy and block cluster tree.
    /// These are the default values, so you will typically not have to bother about these.
    
    /// How many simplices to allow at most in the leave nodes of the bounding volume hierarchy.
    M.cluster_tree_settings.split_threshold                        =  2;
    
    /// Set the separation parameter for the multipole acceptance criterion.
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    
    /// Set the separation parameter for the multipole acceptance criterion in adaptive refinement.
    M.adaptivity_settings.theta                                    = 10.0;

    tic("Creating ClusterTree");
    /// It is not necessary to call this separately.
    /// It will automatically be called by all routines that require it.
    M.GetClusterTree();
    toc("Creating ClusterTree");

    tic("Creating BlockClusterTree");
    /// It is not necessary to call this separately.
    /// It will automatically be called by all routines that require it.
    M.GetBlockClusterTree();
    toc("Creating BlockClusterTree");

    dump(M.GetClusterTree().ThreadCount());

    dump(M.GetClusterTree().NearDim());

    dump(M.GetBlockClusterTree().PrimitiveIntersectionCount());

    print("");

    /// Instantiate factories for multipole-accelerated energy and metric, resp.
    TangentPointEnergy_Factory<Mesh_T,2,2,3,3> TPE_factory;
    TangentPointMetric_Factory<Mesh_T,2,2,3,3> TPM_factory;
    
    TangentPointEnergy0_Factory<Mesh_T,2,2,3,3> TPE0_factory;
    TangentPointMetric0_Factory<Mesh_T,2,2,3,3> TPM0_factory;
    
    const Real q = 6;
    const Real p = 12;

    std::unique_ptr<Energy_T> tpe_ptr = TPE_factory.Make( dom_dim, amb_dim, q, p );
    std::unique_ptr<Metric_T> tpm_ptr = TPM_factory.Make( dom_dim, amb_dim, q, p );
    
    std::unique_ptr<Energy_T> tpe0_ptr = TPE0_factory.Make( dom_dim, amb_dim, q, p );
    std::unique_ptr<Metric_T> tpm0_ptr = TPM0_factory.Make( dom_dim, amb_dim, q, p );

    /// Get references.
    const auto & tpe  = *tpe_ptr;
    const auto & tpm  = *tpm_ptr;
    const auto & tpe0 = *tpe0_ptr;
    const auto & tpm0 = *tpm0_ptr;

    double en;
    double en0;

    tic("tpe.Energy(M)");
    en = tpe.Value(M);
    toc("tpe.Energy(M)");
    
    tic("tpe0.Energy(M)");
    en0 = tpe0.Value(M);
    toc("tpe0.Energy(M)");

    dump(en);
    dump(en0);

    
    /// Require differential of tangent-point energy and load it into the variable `diff`.
    /// `Mesh_T::CotangentVector_T` is `Tensor2<Real,Int>` in this case. ]
    /// It is a simple container class for heap-allocated matrices.
    Mesh_T::CotangentVector_T diff;
    tic("tpe.Differential(M)");
    diff = tpe.Differential(M);
    toc("tpe.Differential(M)");
    
    /// Alternatively, you can write into any pointer-array like `diff.data()`.
    tic("tpe.Differential(M)");
    tpe.Differential(M, diff.data() );
    toc("tpe.Differential(M)");

    print("");
    print("");

    /// Let's test the metrics.
    /// `Mesh_T::TangentVector_T` and `Mesh_T::CotangentVector_T` are both equal
    /// to `Tensor2<Real,Int>` in this example.
    Mesh_T::TangentVector_T   X ( M.VertexCount(), M.AmbDim() );
    Mesh_T::CotangentVector_T Y ( M.VertexCount(), M.AmbDim() );

    /// Load some random data into X.
    X.Random();

    const Real alpha = 1.0;
    const Real beta  = 0.0;
    
    ///Performs generalized matrix-matrix product Y = alpha * A * X + beta * Y, where A is the tangent-point metric.

    tic("tpm.MetricValues(M)");
    /// It is not necessary to call this separately.
    /// It will automatically be called by all routines that require it.
    tpm.MetricValues(M);
    toc("tpm.MetricValues(M)");

    tic("tpm.MultiplyMetric");
    /// Compute `Y = alpha * A * x + beta * Y`, where `A` is the Gram matrix if the metric.
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm.MultiplyMetric");

    tic("tpm.MultiplyMetric");
    /// Compute `Y = alpha * A * x + beta * Y`, where `A` is the Gram matrix if the metric.
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm.MultiplyMetric");
    
    print("");

    /// Now the same for the 0-variantes.
    tic("tpm0.MetricValues(M)");
    /// It is not necessary to call this separately.
    /// It will automatically be called by all routines that require it.
    tpm0.MetricValues(M);
    toc("tpm0.MetricValues(M)");

    tic("tpm0.MultiplyMetric");
    /// Compute `Y = alpha * A * x + beta * Y`, where `A` is the Gram matrix if the metric.
    tpm0.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm0.MultiplyMetric");

    tic("tpm0.MultiplyMetric");
    /// Compute `Y = alpha * A * x + beta * Y`, where `A` is the Gram matrix if the metric.
    tpm0.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm0.MultiplyMetric");
    
    
    print("");
    print("");

    /// Initialize obstacle mesh and load it into M.
    tic("Load obstacle");
    M.LoadObstacle(
       mesh_factory.Make(
           &obstacle_vertex_coordinates[0][0], obstacle_vertex_count,  amb_dim,   false,
           &obstacle_simplices[0][0],          obstacle_simplex_count, dom_dim+1, false,
           thread_count
       )
    );
    toc("Load obstacle");
    print("");

    tic("Create obstacle trees");
    /// It is not necessary to call this separately.
    /// It will automatically be called by all routines that require it.
    M.GetObstacleBlockClusterTree();
    toc("Create obstacle trees");
    print("");

    /// Compute tangent-point energy between mesh and obstacle.
    TangentPointObstacleEnergy_Factory<Mesh_T,2,2,2,2,3,3> TPOE_factory;

    std::unique_ptr<Energy_T> tpo_ptr = TPOE_factory.Make( dom_dim, dom_dim, amb_dim, q, p );
    const auto & tpo = *tpo_ptr;

    tic("Compute tangent-point energy between mesh and obstacle");
    en = tpo.Value(M);
    toc("Compute tangent-point energy between mesh and obstacle");
    dump(en);
    print("");

    ///Compute derivative also of this energy.
    tic("Compute derivative of obstacle energy");
    diff = tpo.Differential(M);
    toc("Compute derivative of obstacle energy");

    print("");
    print("");


    dump(tpe0.Value(M));
    
    return 0;
}
