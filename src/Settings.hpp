#pragma once

namespace Repulsor
{
    struct ClusterTreeSettings
    {
        // An auxilliary data structure in order to feed and hold options for ClusterTree.
        int split_threshold     = 1;
        int thread_count        = 1;
        int parallel_perc_depth = 5;
        
//        TreePercolationAlgorithm tree_perc_alg = TreePercolationAlgorithm::Tasks;
//        TreePercolationAlgorithm tree_perc_alg = TreePercolationAlgorithm::Sequential;
        TreePercolationAlgorithm tree_perc_alg = TreePercolationAlgorithm::Parallel;

        BoundingVolumeType bounding_volume_type = BoundingVolumeType::AABB_LongestAxisSplit;
        
    }; // ClusterTreeSettings
    
    struct BlockClusterTreeSettings
    {
        using Int = int;
        
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
        
        double  far_field_separation_parameter = 0.5;
        double near_field_separation_parameter = 10.;
        double near_field_intersection_parameter = 10000000000.;
        
        // Maximal number of subdivisions to apply for simplices in the very-near field./
        int max_refinement = 30;
        
        BlockSplitMethod block_split_method = BlockSplitMethod::Parallel;
    };
    
    struct AdaptivitySettings
    {
        int    max_refinement = 30;
        double theta = 10.;
        double intersection_theta = 10000000000.;
    };

} // namespace Repulsor
