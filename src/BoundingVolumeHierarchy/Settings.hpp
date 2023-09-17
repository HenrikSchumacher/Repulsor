#pragma once

namespace Repulsor
{
    struct BoundingVolumeHierarchySettings
    {
        // An auxilliary data structure in order to feed and hold options for ClusterTree.
        long long split_threshold     = 1;
        long long thread_count        = 1;
        long long parallel_perc_depth = 5;
        
        BoundingVolumeType bounding_volume_type = BoundingVolumeType::AABB_LongestAxisSplit;
        
    }; // BoundaryVolumeHierarchySettings
    
} // namespace Repulsor

