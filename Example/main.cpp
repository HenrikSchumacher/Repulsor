#include <iostream>

#include <sys/types.h>
#include <pwd.h>

// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
#define INT long long
#define REAL double

#define TOOLS_ENABLE_PROFILER // enable profiler

#define ENABLE_POINTCLOUDS  0
#define ENABLE_CURVES       0
#define ENABLE_SURFACES     1

#define ENABLE_1D           0
#define ENABLE_2D           0
#define ENABLE_3D           1
#define ENABLE_4D           0

#include "../Repulsor.hpp"

int main(int argc, const char * argv[])
{
    
    using namespace Repulsor;
    using namespace Tensors;
    using namespace Tools;
    
    const char * homedir = getenv("HOME");

    if( homedir == NULL)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string path ( homedir );
    
    Profiler::Clear( path );

    
    int thread_count = 8;
    omp_set_num_threads(thread_count);
    
    
#include "Meshes.hpp"
    
    // Meaning of template parameters:
    
    // Make_SimplicialMesh<Real, Int, SReal, ExtReal, ExtInt>;
    //
    // Real    = floating point type used for computations regarding energy and metric
    // Int     = signed integer type used internally, in particular as indices for sparse arrays.
    //           Int = int       corresponds to MKL's  LP64 mode
    //           Int = long long corresponds to MKL's ILP64 mode
    // SReal   = ("short real") floating point type for storage of clustering information
    // ExtReal = ("external real") floating point type that is used by the outer world, e.g. for submitting the mesh's vertex coordinates and vectors to multiply against the metric.
    // ExtInt = ("external real") integer type that is used by the outer world, e.g. for submitting the mesh's vertex coordinates and vectors to multiply against the metric.
    
    
    using Mesh_T = SimplicialMeshBase<REAL,INT,REAL,REAL>;
    
    // Initialize mesh by the factory Make_SimplicialMesh to allow runtime polymorphism.
    tic("Initializing mesh");
    std::unique_ptr<Mesh_T> M_ptr = Make_SimplicialMesh<REAL,INT,REAL,REAL>(
        &vertex_coordinates[0][0],  vertex_count, amb_dim,
        &simplices[0][0],          simplex_count, dom_dim+1,
        thread_count
    );
    
    auto & M = *M_ptr;
    
    dump(M.ThreadCount());
    
//    tic("Initializing mesh");
//    SimplicialMesh<dom_dim,amb_dim,REAL,INT,REAL,REAL> M (
//        &vertex_coordinates[0][0],  vertex_count,
//        &simplices[0][0],          simplex_count,
//        thread_count
//    );

    toc("Initializing mesh");
    
//    // Alternatively, you can do this, independent on which dimensions are enabled for Make_SimplicialMesh. However, dom_dim, and amb_dim have to be known at compile time:
//
//    auto M = std::make_unique<SimplicialMesh<dom_dim,amb_dim,REAL,INT,REAL,REAL>>(
//        &vertex_coordinates[0][0], vertex_count,
//        &simplices[0][0],          simplex_count,
//        thread_pool.ThreadCount()
//    );
    
    // Some quite decent settings for 2-dimensional surfaces.
    
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.5;
    M.adaptivity_settings.theta                                    = 10.0;

    tic("Creating ClusterTree");
    M.GetClusterTree();
    toc("Creating ClusterTree");
    
    tic("Creating BlockClusterTree");
    M.GetBlockClusterTree();
    toc("Creating BlockClusterTree");

//    print(M.GetBlockClusterTree().Stats());
    
    valprint("M.GetClusterTree().ThreadCount()",M.GetClusterTree().ThreadCount());

    valprint("M.GetClusterTree().NearDim()",M.GetClusterTree().NearDim());
//    print(M.GetBlockClusterTree().Stats());

    valprint("number of detected intersections",M.GetBlockClusterTree().PrimitiveIntersectionCount());

    print("");

    const double q = 6;
    const double p = 12;
//    const double weight = 1;
    
    std::unique_ptr<EnergyBase<REAL,INT,REAL,REAL>>
        tpe_ptr = Make_TangentPointEnergy<REAL,INT,REAL,REAL> (dom_dim,amb_dim,q,p);
    
    auto & tpe = *tpe_ptr;
    
    std::unique_ptr<MetricBase<REAL,INT,REAL,REAL>>
        tpm_ptr = Make_TangentPointMetric<REAL,INT,REAL,REAL> (dom_dim,amb_dim,q,p);

    auto & tpm = *tpm_ptr;
    
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

   
    // Mesh_T::TangentVector_T and Mesh_T::CotangentVector_T are both Tensor2<REAL,INT> in this example.
    Mesh_T::TangentVector_T   X ( M.VertexCount(), M.AmbDim() );
    Mesh_T::CotangentVector_T Y ( M.VertexCount(), M.AmbDim() );

    // Load some random data into U.
    X.Random();

    const REAL alpha = 1.0;
    const REAL beta  = 0.0;
    //Performs generalized matrix-matrix product Y = alpha * A * X + beta * Y, where A is the tangent-point metric.
    
    tic("tpm.MetricValues(M)");
    tpm.MetricValues(M);
    toc("tpm.MetricValues(M)");
    
    tic("Matrix multiplication");
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("Matrix multiplication");

    tic("Matrix multiplication");
    tpm.MultiplyMetric(M, alpha, X, beta, Y);
    toc("Matrix multiplication");
   


    // Initialize mesh by the factory Make_SimplicialMesh to allow runtime polymorphism.
    tic("Initialize obstacle mesh");
    std::unique_ptr<Mesh_T> Q = Make_SimplicialMesh<REAL,INT,REAL,REAL>(
        &obstacle_vertex_coordinates[0][0],  obstacle_vertex_count, amb_dim,
        &obstacle_simplices[0][0],          obstacle_simplex_count, dom_dim+1,
        thread_count
      );
    toc("Initialize obstacle mesh");
    print("");

//    // Alternatively, you can do this, independent on which dimensions are enabled for Make_SimplicialMesh. However, dom_dim, and amb_dim have to be known at compile time:
////    Q = std::unique_ptr<SimplicialMeshBase<double,int,double,double>>(
////        new SimplicialMesh<dom_dim,amb_dim,double,int,double,double>(
////            &obstacle_vertex_coordinates[0][0],  obstacle_vertex_count,
////            &obstacle_simplices[0][0],          obstacle_simplex_count,
////            thread_pool.ThreadCount()
////        )
////    );

    // Load obstacle into mesh.
    tic("Load obstacle");
    M.LoadObstacle( std::move(Q) );
    toc("Load obstacle");
    print("");

    tic("Create obstacle trees");
    M.GetObstacleBlockClusterTree();
    toc("Create obstacle trees");
    print("");

//    print(M.GetObstacleBlockClusterTree().Stats() );

//    //Compute tangent-point energy between mesh and obstacle.
//
//    tic("Compute tangent-point energy between mesh and obstacle");
//    double tpo = M.TangentPointObstacleEnergy();
//    toc("Compute tangent-point energy between mesh and obstacle");
//    print("");
//
//    valprint("tangent-point obstacle energy", tpo);
//    print("");
//
//    //Compute derivative also of this energy.
//
//    add_to = true;
//
//    tic("Compute derivative of obstacle energy");
//    M.TangentPointObstacleEnergy_Differential(diff,add_to);        // Add into diff.
//    toc("Compute derivative of obstacle energy");
//    print("");
//    // Some vector-valued functions on vertices, represented by matrices.
//
//    const INT cols = M.AmbDim();
//
//    // Tensor2 is a simple container class for matrices.
//    Tensor2<REAL,INT> U ( M.VertexCount(), cols );
//    Tensor2<REAL,INT> V ( M.VertexCount(), cols );
//
//    // Load some random data into U.
//    U.Random();
//
//    const REAL a = 1.0;
//    const REAL b = 0.0;
//
//    print("");
//    print("");
//
//    valprint("M.ThreadCount()",M.ThreadCount());
//
//    //Performs generalized matrix-matrix product V = a * A * U + b * Y, where A is the tangent-point metric. Passing matrices by pointers (U.data() and V.data() are of type double*).
//    tic("Matrix multiplication");
//    M.TangentPointMetric_Multiply( a, U.data(), b, V.data(), cols, KernelType::LowOrder       );
//    M.TangentPointMetric_Multiply( a, U.data(), b, V.data(), cols, KernelType::HighOrder      );
//    M.TangentPointMetric_Multiply( a, U.data(), b, V.data(), cols, KernelType::FractionalOnly );
//        // This performs multiplication with high order _and low order term and returns the sum.
//    M.TangentPointMetric_Multiply( a, U.data(), b, V.data(), cols );
//    toc("Matrix multiplication");
//
//
//    print("");
//    print("Btw.: M.TangentPointMetric_Multiply has some overhead when you call it the first time. So next time it will be faster.");
//
//    print("");
//    tic("Matrix multiplication");
//    // You can also pass Tensor1 and Tensor2 objects directly.
//    M.TangentPointMetric_Multiply( a, U, b, V, KernelType::LowOrder       );
//    M.TangentPointMetric_Multiply( a, U, b, V, KernelType::HighOrder      );
//    M.TangentPointMetric_Multiply( a, U, b, V, KernelType::FractionalOnly );
//
//    // This performs multiplication with high order _and low order term and returns the sum.
//    M.TangentPointMetric_Multiply( a, U, b, V );
//    toc("Matrix multiplication");

    print("");
    print("");

    return 0;
}
