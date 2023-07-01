#include <iostream>

#include <sys/types.h>
#include <pwd.h>

// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.

//#define NDEBUG
#define TOOLS_DEBUG

#define TOOLS_ENABLE_PROFILER // enable profiler

#define LAPACK_DISABLE_NAN_CHECK
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
//#include <cblas.h>
//#include <lapacke.h>

#include "../Repulsor.hpp"
#include "../Tensors/MyBLAS.hpp"
#include "../Tensors/Sparse.hpp"


using namespace Repulsor;
using namespace Tensors;
using namespace Tools;

using Int     = Int64;
using LInt    = Int64;
using ExtInt  = Int64;

using Real    = Real64;
using SReal   = Real64;
using ExtReal = Real64;


int main(int argc, const char * argv[])
{
    const char * homedir = getenv("HOME");

    if( homedir == nullptr)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string home ( homedir );
    
    Profiler::Clear( home );

    std::string path = home + "/github/BAEMM/Meshes/";
    std::string name = "TorusMesh_00038400T.txt";
//    std::string name = "Spot_00005856T.txt";
    
    
    using Mesh_T   = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
    using Energy_T = EnergyBase<Mesh_T>;
    using Metric_T = MetricBase<Mesh_T>;

    const Int thread_count = 1;
    const Real q = 6;
    const Real p = 12;
    
    SimplicialMesh_Factory<Mesh_T,2,2,3,3> mesh_factory;
    TangentPointEnergy_Factory<Mesh_T,2,2,3,3> TPE_factory;
//    TangentPointMetric_Factory<Mesh_T,2,2,3,3> TPM_factory;
    
    std::unique_ptr<Energy_T> tpe_ptr = TPE_factory.Make(2,3,q,p);
    const auto & tpe = *tpe_ptr;

//    std::unique_ptr<Metric_T> tpm_ptr = TPM_factory.Make(2,3,q,p);
//    const auto & tpm = *tpm_ptr;
    
    tic("Initializing mesh");
    std::unique_ptr<Mesh_T> M_ptr = mesh_factory.Make_FromFile( path + name, thread_count );
    
    auto & M = *M_ptr;
    
    dump(M.ThreadCount());

    toc("Initializing mesh");
    
    dump( min_buffer( M.Simplices().data(), (M.DomDim()+1) * M.SimplexCount() ) );
    dump( max_buffer( M.Simplices().data(), (M.DomDim()+1) * M.SimplexCount() ) );
    
    print("");
    
    // Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;

    tic("GetClusterTree");
    M.GetClusterTree();           // Not necessary. Will automatically called by all routines that require it.
    toc("GetClusterTree");

    tic("GetBlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically called by all routines that require it.
    toc("GetBlockClusterTree");
//
    
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
    
    tic("MaximumSafeStepSize");
    Tensor2<Real,Int> V ( M.VertexCount(), M.AmbDim() );
    V.Random();
//    V.Fill(1.);
    const Real t = M.MaximumSafeStepSize(V.data(), 1.);
    toc("MaximumSafeStepSize");
    dump(t);
    
    print("");
    print("###############################################");
    print("");
    tic("SemiStaticUpdate");
    V *= t;
    V += M.VertexCoordinates();
    
    M.SemiStaticUpdate(V.data());
    toc("SemiStaticUpdate");
    print("");
    
    
    
    tic("GetClusterTree");
    M.GetClusterTree();           // Not necessary. Will automatically called by all routines that require it.
    toc("GetClusterTree");

    tic("GetBlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically called by all routines that require it.
    toc("GetBlockClusterTree");
    
    tic("tpe.Energy(M)");
    en = tpe.Value(M);
    toc("tpe.Energy(M)");

    dump(en);

    tic("tpe.Differential(M)");
    diff = tpe.Differential(M);
    toc("tpe.Differential(M)");
    

    return 0;
}
