#include <iostream>

#include <sys/types.h>
#include <pwd.h>

//#define NDEBUG
#define TOOLS_DEBUG

#define TOOLS_ENABLE_PROFILER // enable profiler

#define LAPACK_DISABLE_NAN_CHECK

/// Use these while on a mac. Don't forget to issue the
/// compiler flag `-framework Accelerate`.
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>


/// Use these instead of Accelerate under Windows or Linux, e.g. together with OpenBLAS, Intel oneMKL, AMD AOCL-BLAS. Of course, your path variables or compiler flags should hint the compiler to these files. And you also have to link the according libraries.
//#include <cblas.h>
//#include <lapack.h>

#include "../Repulsor.hpp"

using namespace Repulsor;
using namespace Tensors;
using namespace Tools;


// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
using Int     = Int64;
using LInt    = Int64;
using ExtInt  = Int64;

using Real    = Real64;
using SReal   = Real64;
using ExtReal = Real64;


int main(int argc, const char * argv[])
{
    print("");
    print("###############################################################");
    print("###              Test program preconditioner               ####");
    print("###############################################################");
    print("");
    
    
    const char * homedir = getenv("HOME");

    if( homedir == nullptr)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string home ( homedir );
    
    Profiler::Clear( home );

    
    int thread_count = 8;

    std::string path = home + "/github/BAEMM/Meshes/";
//    std::string name = "TorusMesh_00038400T.txt";
//    std::string name = "Spot_00005856T.txt";
    std::string name = "Spot_00093696T.txt";
    
    
    using MeshBase_T     = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
    using Mesh_T         = SimplicialMesh<2,3,Real,Int,LInt,SReal,ExtReal>;
    using RemesherBase_T = SimplicialRemesherBase<Real,Int,ExtReal,ExtInt>;
    
    SimplicialMesh_Factory<MeshBase_T,2,2,3,3> mesh_factory;
    
    SimplicialRemesher_Factory<RemesherBase_T,2,2,3,3> remesher_factory;
    
    
    tic("Initializing mesh");
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile( path + name, thread_count );
    
    auto & M = *M_ptr;

    dump(M.ThreadCount());

    toc("Initializing mesh");

    print("");

    
    print("");
    print("Testing remesher.");
    print("");
    
    
    std::unique_ptr<RemesherBase_T> R = remesher_factory.Make(
        M.VertexCoordinates().data(), M.VertexCount(),  M.AmbDim(), false,
        M.Simplices().data(),         M.SimplexCount(), M.AmbDim(), false,
        nullptr,                                        0,          false,
        thread_count
    );
    
    Tensor1<Real,Int> squared_edge_lengths = R->SquaredEdgeLengths();
    
    Sort( squared_edge_lengths.begin(), squared_edge_lengths.end() );
    
    const Real lower_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * Real(0.1))] );
    const Real upper_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * Real(0.9))]  );
    
    valprint("Minimum edge length            ", std::sqrt(squared_edge_lengths.First()) );
    valprint("Lower bound used for collapsing", lower_bound );
    valprint("Upper bound used for splitting ", upper_bound );
    valprint("Maximum edge length            ", std::sqrt(squared_edge_lengths.Last()) );
    
    const Int  collapse_split_iter = 4;

    print("");
    
    
    //Apply collapse_split_iter rounds of edge collapse and edge split.
    tic("UnifyEdgeLengths");
    R->UnifyEdgeLengths( lower_bound, upper_bound, collapse_split_iter );
    toc("UnifyEdgeLengths");
    
    print("");
    
    const Int flip_iter = 10;
    
    print("");
    
    //Apply flip_max_iter rounds of edge flips to improve triangle shapes.
    tic("DelaunayFlip");
    R->DelaunayFlip( flip_iter );
    toc("DelaunayFlip");
    
    const Int smooth_iter = 10;
    
    print("");
    
    // Apply smooth_iter iterations of tangential smoothening to further improve mesh regularity.
    tic("TangentialSmoothing");
    R->TangentialSmoothing( smooth_iter );
    toc("TangentialSmoothing");
    
    
    print("");
    
    R->SelfCheck();
    
    print("");
    
    tic("Compress");
    R->Compress();
    toc("Compress");
    
    print("");
    
    R->SelfCheck();
    
    print("");
    
    squared_edge_lengths = R->SquaredEdgeLengths();
    
    valprint("Minimum edge length            ", std::sqrt(squared_edge_lengths.Min()) );
    valprint("Maximum edge length            ", std::sqrt(squared_edge_lengths.Max()) );

    print("");
    
    // Creating a new mesh object.
    
    std::unique_ptr<MeshBase_T> N_ptr = mesh_factory.Make(
        R->VertexCoordinates().data(), R->VertexCount(),  R->AmbDim(),   false,
        R->Simplices().data(),         R->SimplexCount(), R->DomDim()+1, false,
        thread_count
    );
    
    dump(N_ptr->VertexCount());
    dump(N_ptr->SimplexCount());
    
    return 0;
}
