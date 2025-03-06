#pragma once

namespace Repulsor
{
    enum class BoundingVolumeType
    {
//        AABB,
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
        Sequential,
        Recursive,
        Parallel
    };

    enum class InteractionKernels
    {
        Default,
        MinDist,
        Degenerate
    };
    
    enum class OperatorType
    {
        FractionalOnly,
        HighOrder,
        MixedOrder,
        LowOrder,
        SquaredDistance
    };
  
    enum class InOut : int
    {
        In,
        Out
    };
    
    namespace FMM
    {
        enum class Type : int
        {
            IN = -1,
            VF =  0,
            NF =  1,
            FF =  2
        };
    }

    static std::map<OperatorType, std::string> OperatorTypeName
    {
        {OperatorType::FractionalOnly, "FractionalOnly"},
        {OperatorType::HighOrder, "HighOrder"},
        {OperatorType::MixedOrder, "MixedOrder"},
        {OperatorType::LowOrder, "LowOrder"},
        {OperatorType::SquaredDistance, "SquaredDistance"}
    };
    
} // namespace Repulsor
