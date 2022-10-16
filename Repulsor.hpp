#ifndef REPULSOR_HPP

    #define REPULSOR_HPP

    #include <tuple>
    #include <unordered_map>

    //#include "MyMath.hpp"

    #include "Tensors/Tensors.hpp"
    #include "GJK/GJK.hpp"
    //#include "Collision.hpp"

    namespace Repulsor {
        
        using namespace Tools;
        using namespace Tensors;
        using namespace GJK;
    }


    #include "src/Enums.hpp"
    #include "src/Settings.hpp"
    //#include "src/MultipoleMoments.hpp"
    #include "src/SimplexHierarchy.hpp"
    #include "src/ClusterTree.hpp"
    #include "src/Traversor.hpp"

    #include "src/BlockClusterTree.hpp"
    #include "src/CollisionTree.hpp"

    #include "src/SimplicialMesh/SimplicialMeshDetails.hpp"
    #include "src/SimplicialMesh.hpp"

    #include "src/FMM.hpp"

    // TODO: Finalize this!
    #include "src/SimplicialRemesher/SimplicialRemesherBase.hpp"
    #include "src/SimplicialRemesher/SimplicialRemesher.hpp"


//    // toggle whether primitive data should be copied by kernels.
//    #define NearField_S_Copy
//    #define NearField_T_Copy
//    #include "src/Kernels/NearFieldKernelBase.hpp"
//
//    // toggle whether cluster data should be copied by kernels.
//    #define FarField_S_Copy
//    #define FarField_T_Copy
//    #include "src/Kernels/FarFieldKernelBase_FMM.hpp"
//
//    #include "src/Energies.hpp"
//    #include "src/Metrics.hpp"

#endif
