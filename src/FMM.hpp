#pragma once

namespace Repulsor
{
    enum FMM_Type : int
    {
        VF = 0,
        NF = 1,
        FF = 2
    };

} // namespace Repulsor


#include "FMM/FMM_Configurator.hpp"
#include "FMM/FMM_Kernel.hpp"
#include "FMM/FMM_Kernel_FF.hpp"
#include "FMM/FMM_Kernel_NF.hpp"
#include "FMM/FMM_Kernel_VF.hpp"
#include "FMM/FMM_Traversor.hpp"
