#include <iostream>

#include <sys/types.h>
#include <pwd.h>

// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
#define INT     int32_t
#define LINT    std::size_t
//#define INT     int64_t
#define EXTINT  int64_t
#define REAL    double
#define SREAL   float
#define EXTREAL double

//#define REMESHER_VERBATIM

#define TOOLS_ENABLE_PROFILER // enable profiler

//#include <Accelerate/Accelerate.h>
#include <cblas.h>
#include <lapacke.h>

#include "../Repulsor.hpp"
#include "../Tensors/MyBLAS.hpp"
#include "../Tensors/Sparse.hpp"



//#include "../Tensors/ConjugateGradient.hpp"
//#include "../Tensors/GMRES.hpp"
//#include "../Tensors/Sparse.hpp"

using namespace Repulsor;
using namespace Tensors;
using namespace Tools;


int main(int argc, const char * argv[])
{
    const char * homedir = getenv("HOME");

    if( homedir == nullptr)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string path ( homedir );
    
    Profiler::Clear( path );

    
    int thread_count = 8;
//    int thread_count = 1;    
    
#include "Meshes.hpp"
    
    // Meaning of template parameters:
    
    // SimplicialMeshBase<Real, Int, SReal, ExtReal>;
    //
    // Real    = floating point type used for computations regarding energy and metric
    // Int     = signed integer type used internally, in particular as indices for sparse arrays.
    //           Int = int       corresponds to MKL's  LP64 mode
    //           Int = long long corresponds to MKL's ILP64 mode
    // SReal   = ("short real") floating point type for storage of clustering information
    // ExtReal = ("external real") floating point type that is used by the outer world, e.g. for submitting the mesh's vertex coordinates and vectors to multiply against the metric.
    

    using Mesh_T   = SimplicialMeshBase<REAL,INT,SREAL,EXTREAL>;
    using Energy_T = EnergyBase<REAL,INT,SREAL,EXTREAL>;
    using Metric_T = MetricBase<REAL,INT,SREAL,EXTREAL>;

    // Create a factory that can make instances of SimplicialMesh with domain dimension in the range 2,...,2 and ambient dimension in the range 3,...,3.
    // (The main purpose of this factory is to tell the compiler which template instantiations of SimplicialMesh it shall generate. Note that with obstacles one might want to use meshes of various dimensions. So this factory basically about saving compile time.)
    SimplicialMesh_Factory<Mesh_T,2,2,3,3> mesh_factory;


    // This factory allows run-time polymorphism. But know that all specializations of SimplicialMesh in theses ranges have to be compiled! That may lead to horrific compile times, thus the ability to restrict the ranges in the factory.

    // Initialize mesh by the mesh_factory to allow runtime polymorphism.
    tic("Initializing mesh");
    std::unique_ptr<Mesh_T> M_ptr = mesh_factory.Make(
        &vertex_coordinates[0][0],  vertex_count, amb_dim,   false,
        &simplices[0][0],          simplex_count, dom_dim+1, false,
        thread_count
    );

    auto & M = *M_ptr;  // I don't like pointers. Give me a reference.

    dump(M.ThreadCount());

    toc("Initializing mesh");

    // Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;

    tic("Creating ClusterTree");
    M.GetClusterTree();           // Not necessary. Will automatically called by al routines that require it.
    toc("Creating ClusterTree");

    tic("Creating BlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically called by al routines that require it.
    toc("Creating BlockClusterTree");

    valprint("M.GetClusterTree().ThreadCount()",M.GetClusterTree().ThreadCount());

    valprint("M.GetClusterTree().NearDim()",M.GetClusterTree().NearDim());

    valprint("number of detected intersections",M.GetBlockClusterTree().PrimitiveIntersectionCount());

    print("");

    const REAL q = 6;
    const REAL p = 12;

    // Again some factories.
    TangentPointEnergy_Factory<Mesh_T,2,2,3,3> TPE_factory;
    TangentPointMetric_Factory<Mesh_T,2,2,3,3> TPM_factory;

    std::unique_ptr<Energy_T> tpe_ptr = TPE_factory.Make( dom_dim, amb_dim, q, p );

    const auto & tpe = *tpe_ptr;

    std::unique_ptr<Metric_T> tpm_ptr = TPM_factory.Make( dom_dim, amb_dim, q, p );

    const auto & tpm = *tpm_ptr;

    double en;
    // Mesh_T::CotangentVector_T is Tensor2<REAL,INT> in this case. It is a simple container class for heap-allocated matrices.
    Mesh_T::CotangentVector_T diff;

    tic("tpe.Energy(M)");
    en = tpe.Value(M);
    toc("tpe.Energy(M)");

    dump(en);

    tic("tpe.Differential(M)");
    diff = tpe.Differential(M);
    toc("tpe.Differential(M)");

    print("");

    TangentPointEnergy_AllPairs_Factory<Mesh_T,2,2,3,3> TPE_AllPairs_factory;

    std::unique_ptr<Energy_T> tpe_slow_ptr = TPE_AllPairs_factory.Make( dom_dim, amb_dim, q, p );

    const auto & tpe_slow = *tpe_slow_ptr;

    tic("tpe_slow.Energy(M)");
    en = tpe_slow.Value(M);
    toc("tpe_slow.Energy(M)");

    dump(en);
    print("");

    // Mesh_T::TangentVector_T and Mesh_T::CotangentVector_T are both Tensor2<REAL,INT> in this example.
    Mesh_T::TangentVector_T   X ( M.VertexCount(), M.AmbDim() );
    Mesh_T::CotangentVector_T Y ( M.VertexCount(), M.AmbDim() );

    // Load some random data into U.
    X.Random();

    const REAL alpha = 1.0;
    const REAL beta  = 0.0;
    //Performs generalized matrix-matrix product Y = alpha * A * X + beta * Y, where A is the tangent-point metric.

    tic("tpm.MetricValues(M)");
    tpm.MetricValues(M);         // Not necessary. Will automatically called by al routines that require it.
    toc("tpm.MetricValues(M)");

    tic("tpm.MultiplyMetric");
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm.MultiplyMetric");

    tic("tpm.MultiplyMetric");
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm.MultiplyMetric");

    dump(X.MaxNorm());
    dump(Y.MaxNorm());

    print("");


    // Initialize mesh by the factory to allow runtime polymorphism.
    tic("Initialize obstacle mesh");
    std::unique_ptr<Mesh_T> Q_ptr = mesh_factory.Make(
        &obstacle_vertex_coordinates[0][0], obstacle_vertex_count,  amb_dim,   false,
        &obstacle_simplices[0][0],          obstacle_simplex_count, dom_dim+1, false,
        thread_count
    );
    toc("Initialize obstacle mesh");
    print("");

    // Load obstacle into mesh.
    tic("Load obstacle");
    M.LoadObstacle( std::move(Q_ptr) );
    toc("Load obstacle");
    print("");

    tic("Create obstacle trees");
    M.GetObstacleBlockClusterTree();    // Not necessary. Will automatically called by al routines that require it.
    toc("Create obstacle trees");
    print("");

    // Compute tangent-point energy between mesh and obstacle.

    TangentPointObstacleEnergy_Factory<Mesh_T,2,2,2,2,3,3> TPOE_factory;

    std::unique_ptr<Energy_T> tpo_ptr = TPOE_factory.Make( dom_dim, dom_dim, amb_dim, q, p );

    const auto & tpo = *tpo_ptr;

    tic("Compute tangent-point energy between mesh and obstacle");
    en = tpo.Value(M);
    toc("Compute tangent-point energy between mesh and obstacle");
    dump(en);
    print("");
//
    //Compute derivative also of this energy.

    tic("Compute derivative of obstacle energy");
    diff = tpo.Differential(M);
    toc("Compute derivative of obstacle energy");

    print("");
    print("");


    TangentPointEnergy0_Factory<Mesh_T,2,2,3,3> TPE0_factory;
    TangentPointMetric0_Factory<Mesh_T,2,2,3,3> TPM0_factory;
    PseudoLaplacian_Factory    <Mesh_T,2,2,3,3> Prec_factory;

    const REAL s = (p - 2)/q;

    dump(s);
    dump(2-s);

    std::unique_ptr<Energy_T> tpe0_ptr = TPE0_factory.Make( dom_dim, amb_dim, q, p );
    std::unique_ptr<Metric_T> tpm0_ptr = TPM0_factory.Make( dom_dim, amb_dim, q, p );
    std::unique_ptr<Metric_T> prec_ptr = Prec_factory.Make( dom_dim, amb_dim, 2-s  );

    const auto & tpe0 = *tpe0_ptr;
    const auto & tpm0 = *tpm0_ptr;
    const auto & prec = *prec_ptr;

    dump(tpe0.Value(M));

    tic("tpm0.MetricValues(M)");
    tpm0.MetricValues(M);
    toc("tpm0.MetricValues(M)");

    tic("tpm0.MultiplyMetric");
    tpm0.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm0.MultiplyMetric");

    tic("tpm0.MultiplyMetric");
    tpm0.MultiplyMetric(M, alpha, X, beta, Y);
    toc("tpm0.MultiplyMetric");

    dump(X.MaxNorm());
    dump(Y.MaxNorm());

    tic("prec.MetricValues(M)");
    prec.MetricValues(M);
    toc("prec.MetricValues(M)");

    tic("prec.MultiplyMetric");
    prec.MultiplyMetric(M, alpha, X, beta, Y);
    toc("prec.MultiplyMetric");

    tic("prec.MultiplyMetric");
    prec.MultiplyMetric(M, alpha, X, beta, Y);
    toc("prec.MultiplyMetric");

    dump(X.MaxNorm());
    dump(Y.MaxNorm());
    
    
    // A = M.StiffnessMatrix() + M.MassMatrix()
    auto A = M.H1Metric(1,1);
    
    auto perm = M.NestedDissectionOrdering();
    
//    print( perm.ToString() );
    
    Sparse::CholeskyDecomposition<REAL, INT, LINT> S (
        A.Outer().data(), A.Inner().data(), perm.data(), A.RowCount(),
        A.ThreadCount(), static_cast<INT>(4)
    );
    
    S.SymbolicFactorization();
    
    S.NumericFactorization( A.Values().data(), Scalar::Zero<REAL> );
    
//    print( S.GetPermutation().GetPermutation().ToString() );

    Tensor1<REAL,INT> b ( M.VertexCount(), 1 );
    Tensor1<REAL,INT> x ( M.VertexCount(), 0 );
    
    S.Solve( b.data(), x.data() );
    
//    print( x.ToString() );
    
    print("");
    print("Testing remesher.");

    std::unique_ptr<Mesh_T::Remesher_T> R = M_ptr->CreateRemesher();

    std::vector<INT> edges {1,2,3};

    R->CollapseEdges( edges );

    R->SplitEdges( edges );

//    R->UnifyEdgeLengths(<#const Real collapse_threshold#>, <#const Real split_threshold#>);
    
    
    return 0;
}
