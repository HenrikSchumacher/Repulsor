#pragma once

#ifndef REPULSOR_HPP

    #define REPULSOR_HPP

    #include <tuple>
    #include <unordered_map>

    #include "submodules/Tensors/Tensors.hpp"
    #include "submodules/Tensors/Sparse.hpp"
    #include "submodules/Tensors/ConjugateGradient.hpp"

    namespace Repulsor {
        
        using namespace Tools;
        using namespace Tensors;
        
        using GJK_Real = double;
    }

//    #include "src/GJK_Old.hpp"
    #include "src/GJK.hpp"
    #include "src/CollisionFinder.hpp"

    #include "src/Enums.hpp"
    #include "src/Settings.hpp"

//    #include "src/SpaceFillingCurve.hpp"

    //#include "src/MultipoleMoments.hpp"
    #include "src/SimplexHierarchy.hpp"
    #include "src/ClusterTree.hpp"
    #include "src/ClusterTreePairTraversor.hpp"

    #include "src/BlockClusterTree.hpp"
    #include "src/CollisionTree.hpp"


    #include "src/SimplexDataKernel.hpp"
    #include "src/SimplicialMesh/SimplicialMeshDetails.hpp"
    #include "src/SimplicialMesh.hpp"


    // TODO: Finalize this!
    #include "src/SimplicialRemesher/SimplicialRemesherBase.hpp"
    #include "src/SimplicialRemesher.hpp"
    #include "src/SimplicialRemesher/SimplicialRemesher_Factory.hpp"

//    // toggle whether primitive data should be copied by kernels.
//#define NearField_S_Copy
//#define NearField_T_Copy

//    // toggle whether cluster data should be copied by kernels.
//#define FarField_S_Copy
//#define FarField_T_Copy

    #include "src/FMM.hpp"
    #include "src/Energies.hpp"
    #include "src/Metrics.hpp"

#endif
