#pragma once

namespace Repulsor
{
    enum class BoundingVolumeType
    {
        AABB,
        AABB_LongestAxisSplit,
        AABB_MedianSplit,
        AABB_PreorderedSplit
    };
    
    enum class SplitMethod
    {
        Parallel,
        Sequential,
        Recursive
    };
    
    enum class BlockSplitMethod
    {
        Parallel,
        Sequential,
        Recursive
    };
    
    enum class TreePercolationAlgorithm
    {
        Tasks,
        Sequential,
        Recursive,
        BruteForce
    };

    enum class InteractionKernels
    {
        Default,
        MinDist,
        Degenerate
    };
    
    enum class KernelType
    {
        FractionalOnly,
        HighOrder,
        MixedOrder,
        LowOrder,
        SquaredDistance
    };
    
    static std::map<KernelType, std::string> KernelTypeName
    {
        {KernelType::FractionalOnly, "FractionalOnly"},
        {KernelType::HighOrder, "HighOrder"},
        {KernelType::MixedOrder, "MixedOrder"},
        {KernelType::LowOrder, "LowOrder"},
        {KernelType::SquaredDistance, "SquaredDistance"}
    };
    
} // namespace Repulsor
