#pragma once

namespace Repulsor
{
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
