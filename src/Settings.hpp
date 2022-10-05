#pragma once

namespace Repulsor
{
    struct ClusterTreeSettings
    {
        // An auxilliary data structure in order to feed and hold options for ClusterTree.
        long long split_threshold   = 1;
        long long threads_available = 1;
        long long thread_count      = 1;
        
//        TreePercolationAlgorithm tree_perc_alg = TreePercolationAlgorithm::Tasks;
//        TreePercolationAlgorithm tree_perc_alg = TreePercolationAlgorithm::Sequential;
        TreePercolationAlgorithm tree_perc_alg = TreePercolationAlgorithm::Recursive;

        BoundingVolumeType bounding_volume_type = BoundingVolumeType::AABB;
        
        ClusterTreeSettings() = default;
        
        ~ClusterTreeSettings() = default;
        
    }; // ClusterTreeSettings
    
    struct BlockClusterTreeSettings
    {
        // handling symmetry properties of the BCT
        bool exploit_symmetry = true;

        bool near_upper_triangular = true;
        bool  far_upper_triangular = true;
        
//        bool near_upper_triangular = false;
//        bool  far_upper_triangular = false;
        

        // If exploit_symmetry == false, S == T is assumed and only roughly half the block clusters are generated during the split pass performed by RequireBlockClusters.
        // If upper_triangular == true and if exploit_symmetry == true, only the upper triangle of the interaction matrices will be generated. --> RequireBlockClusters will be faster.
        // If exploit_symmetry == true and upper_triangular == false then the block cluster twins are generated _at the end_ of the splitting pass by RequireBlockClusters.
        // CAUTION: Currently, we have no faster matrix-vector multiplication for upper triangular matrices, so upper_triangular == false is the default.
        
        long long threads_available = 1;
        
        double  far_field_separation_parameter = 0.5;
        double near_field_separation_parameter = 10.;
        double near_field_intersection_parameter = 10000000000.;
        
        BlockSplitMethod block_split_method = BlockSplitMethod::Parallel;
        
        BlockClusterTreeSettings() = default;
        ~BlockClusterTreeSettings() = default;
    };
    
    struct AdaptivitySettings
    {
        long long  max_level = 30;
        long long  min_level = 0;
        double theta = 10.;
        double intersection_theta = 10000000000.;
        
        AdaptivitySettings() = default;
        
        ~AdaptivitySettings() = default;
    };

} // namespace Repulsor
