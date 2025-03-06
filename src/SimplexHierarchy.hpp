#pragma once


namespace Repulsor
{

    template<int DOM_DIM_, int AMB_DIM_, typename Real, typename Int, typename SReal>
    class alignas(ObjectAlignment) SimplexHierarchy
    {
        
        static constexpr Int DOM_DIM = DOM_DIM_;
        static constexpr Int AMB_DIM = AMB_DIM_;
        
#include "SimplexHierarchy/SimplexHierarchy_Details.hpp"
        
        constexpr Int ChildCount() const
        {
            return 0;
        }
        
        void ToChild( const Int k )
        {
            (void)k;
            print(ClassName()+"::ToChild not implemented.");
        }
        
        void ToParent()
        {
            print(ClassName()+"::ToParent not implemented.");
        }
        
    }; // SimplexHierarchy
    
} // Repulsor

#include "SimplexHierarchy/SimplexHierarchy_0D.hpp"
#include "SimplexHierarchy/SimplexHierarchy_1D.hpp"
#include "SimplexHierarchy/SimplexHierarchy_2D.hpp"
