#include <iostream>

/// This is a simple program meant to make you familiar with how to
/// - generate a simplicial mesh from data _in memory_;
/// - compute the tangent-point energy (and derivative) of such a mesh;
/// - perform matrix-vector operations with the tangent-point metric;
/// - compute the gradient of the tangent-point energy w.r.t. that metric; and
/// - find a collision-free step size for a displacement field along the manifold.

/// Moreover, it also shows you how to
/// - load a an obstacle mesh into a simplicial mesh and
/// - compute the tangent-point obstacle energy (and its derivative).



/// Enable the shipped profiler.
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

using namespace Tools;

using Int     = int;
using LInt    = std::size_t;
using Real    = double;


int main(void)
{
    print("");
    print("###############################################################");
    print("###      Example program for Repulsor::SimplicialMesh      ####");
    print("###############################################################");
    print("");
    
    /// Set up profiler to write to `~/Tools_Profile.tsv` and `~/Tools_Log.txt`
    Profiler::Clear();

    int thread_count = 1;
    
    
    /// Setting up the scalar types to be used by the meshes:
    using Mesh_T = Repulsor::SimplicialMeshBase<Real,Int,LInt>;
    
    /// Meaning of template parameters of `SimplicialMeshBase<Real,Int,LInt>`:
    ///
    /// `Real`      Floating point type used for computations regarding energy and
    ///             metric. Typically, you want to use `Real = double`.
    ///
    /// `Int`       Signed integer type used internally, in particular as indices
    ///             for sparse arrays.
    ///             `Int = int`       corresponds to MKL's  LP64 mode
    ///             `Int = long long` corresponds to MKL's ILP64 mode
    ///             I recommend using `Int = int` when you can. This requires a bit less memory for the sparse matrices and this speeds up matrix-vector multiplication.
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

    
    
    /// Instantiate factory that can create instances of `SimplicialMesh<2,3,Real,Int,LInt>`.
    
    Repulsor::SimplicialMesh_Factory<Mesh_T,2,2,3,3> mesh_factory;
    
    /// We need also base types for the energies and the metrics.
    using Energy_T = Repulsor::EnergyBase<Mesh_T>;
    using Metric_T = Repulsor::MetricBase<Mesh_T>;
    
    /// Instantiate factories for multipole-accelerated energy and metric, resp.
    Repulsor::TangentPointEnergy0_Factory<Mesh_T,2,2,3,3> TPE_factory;
    Repulsor::TangentPointMetric0_Factory<Mesh_T,2,2,3,3> TPM_factory;


#include "../Meshes/TwoSpheres.hpp"
    
    /// `vertex_coordinates`, `vertex_count`, `simplices`, `simplex_count`, `amb_dim`, and `dom_dim` have been defined in `TwoSpheres.hpp`.
    
    tic("Initializing mesh from memory");
    
    /// Use`mesh_factory` to initialize a mesh from arrays that are already in memory.
    std::unique_ptr<Mesh_T> M_ptr = mesh_factory.Make(
        &vertex_coordinates[0][0],  vertex_count, amb_dim,   false,
        &simplices[0][0],          simplex_count, dom_dim+1, false,
        thread_count
    );
    /// I don't like pointers. Give me a reference.
    auto & M = *M_ptr;

    toc("Initializing mesh from memory");
    
    
    TOOLS_DUMP(M.VertexCount());
    TOOLS_DUMP(M.SimplexCount());
    TOOLS_DUMP(M.ThreadCount());

    /// Setup for bounding volume hierarchy and block cluster tree.
    /// These are the default values, so you will typically not have to bother about these.
    
    /// How many simplices to allow at most in the leave nodes of the bounding volume hierarchy.
    M.cluster_tree_settings.split_threshold                        =  2;
    
    /// Set the separation parameter for the multipole acceptance criterion.
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    
    /// Set the separation parameter for the multipole acceptance criterion in adaptive refinement.
    M.adaptivity_settings.theta                                    = 10.0;

    print("");

    
    /// These turned out to be good choices for the parameters of the tangent-point energy for 2-dimensional surfaces.
    const Real q = 6;
    const Real p = 12;
    
    std::unique_ptr<Energy_T> tpe_ptr = TPE_factory.Make( dom_dim, amb_dim, q, p );
    std::unique_ptr<Metric_T> tpm_ptr = TPM_factory.Make( dom_dim, amb_dim, q, p );

    /// Get references.
    const auto & tpe  = *tpe_ptr;
    const auto & tpm  = *tpm_ptr;

    Real en;

    tic("Compute tangent-point energy");
    en = tpe.Value(M);
    toc("Compute tangent-point energy");

    TOOLS_DUMP(en);

    
    /// Require differential of tangent-point energy and load it into the variable `diff`.
    /// Here `Mesh_T::CotangentVector_T` is a simple container class for heap-allocated matrices. (It is a  `Tensors::Tensor2<Real,Int>`.)
    Mesh_T::CotangentVector_T diff;
    tic("Compute differential of energy");
    diff = tpe.Differential(M);
    toc("Compute differential of energy");
    
    /// Alternatively, you can write into any pointer-array like `diff.data()`.
    
    tic("Compute differential of energy (pointer version)");
    tpe.Differential(M,diff.data());
    toc("Compute differential of energy (pointer version)");
    
    tic("Do both at once");
    en = tpe.Differential(M,diff.data());
    toc("Do both at once");
    
    print("");

    /// Let's do something with the metric.

    /// First create a tangent and a cotangent vector.
    /// `Mesh_T::TangentVector_T` and `Mesh_T::CotangentVector_T` are both equal
    /// to `Tensor2<Real,Int>` in this example.
    Mesh_T::TangentVector_T   X ( M.VertexCount(), M.AmbDim() );
    Mesh_T::CotangentVector_T Y ( M.VertexCount(), M.AmbDim() );

    /// Load some random data into U.
    X.Random();

    const Real alpha = 1.0;
    const Real beta  = 0.0;

    tic("Generalized multiplication with metric");
    /// Compute `Y = alpha * A * x + beta * Y`, where `A` is the Gram matrix if the metric.
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("Generalized multiplication with metric");

    /// This will get substantially faster for the second time because `A` is precomputed now.
    tic("Generalized multiplication with metric");
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("Generalized multiplication with metric");
    
    tic("Generalized multiplication with metric");
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("Generalized multiplication with metric");

    print("");

    /// This is how you can get the gradient of the tangent-point energy:
    const Int  max_iter  = 100;
    const Real relative_tolerance = 0.00001;
    
    Mesh_T::TangentVector_T   gradient ( M.VertexCount(), M.AmbDim() );
    
    tic("Solving for gradient");
    /// We have to solve M.AmbDim() equations at once.
    tpm.Solve(M, 
        Real(1), diff.data(),     M.AmbDim(),
        Real(0), gradient.data(), M.AmbDim(),
        M.AmbDim(), max_iter, relative_tolerance
    );
    toc("Solving for gradient");
    
    TOOLS_DUMP(tpm.Solver_IterationCount());
    TOOLS_DUMP(tpm.Solver_RelativeResiduals());
    
    
    Mesh_T::TangentVector_T downward_gradient = gradient;
    
    downward_gradient*= -1.;
    
    print("");
    
    /// Given a displacement field like `downward_gradient` along the manifold `M`, you can use `MaximumSafeStepSize` to compute an (almost) maximal collision-free step size in the interval [0,1].
    tic("MaximumSafeStepSize");
    const Real t = M.MaximumSafeStepSize(downward_gradient.data(), 1.);
    
    toc("MaximumSafeStepSize");
    
    valprint("Collision-free step size",t);
    
    
    print("");
    
    tic("SemiStaticUpdate");
    /// X = M.VertexCoordinates() + t * downward_gradient.
    
    X = downward_gradient;
    X *= t;
    X += M.VertexCoordinates();

    
    /// Update mesh without rebuilding the bounding volume hierarchy.
    /// Use this for line search!
    M.SemiStaticUpdate( X.data() );
    
    toc("SemiStaticUpdate");
    
    print("");
    
    /// Let's see whether the energy has decreased...
    
    Real en_new;

    tic("Compute tangent-point energy of updated mesh");
    en_new = tpe.Value(M);
    toc("Compute tangent-point energy of updated mesh");

    TOOLS_DUMP(en_new);
    
    TOOLS_DUMP(en_new < en);
    
    print("");
    print("");

    /// Initialize obstacle mesh and load it into M.
    ///
    /// `obstacle_vertex_coordinates`, `obstacle_simplices`, `obstacle_vertex_count`, `obstacle_simplex_count` have been defined in the file `TwoSpheres.hpp`.
    ///
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

    TOOLS_DUMP(M.GetObstacle().VertexCount());
    TOOLS_DUMP(M.GetObstacle().SimplexCount());
    
    print("");
    
    /// Compute tangent-point energy between mesh and obstacle.
    
    /// Instantiate a factory that for creating instacle energies between
    /// two instances of `SimplicialMesh<2,3,Real,Int,LInt>`.
    Repulsor::TangentPointObstacleEnergy_Factory<Mesh_T,2,2,2,2,3,3> TPOE_factory;

    std::unique_ptr<Energy_T> tpo_ptr = TPOE_factory.Make( dom_dim, dom_dim, amb_dim, q, p );
    const auto & tpo = *tpo_ptr;

    tic("Compute tangent-point energy between mesh and obstacle");
    en = tpo.Value(M);
    toc("Compute tangent-point energy between mesh and obstacle");
    
    TOOLS_DUMP(en);
    
    print("");

    ///Compute derivative also of this energy.
    tic("Compute derivative of obstacle energy");
    diff = tpo.Differential(M);
    toc("Compute derivative of obstacle energy");

    print("");
    print("");
    
    return 0;
}
