#pragma once

namespace Repulsor
{
    
#define DOM_DIM 0
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class alignas( OBJECT_ALIGNMENT ) SimplexHierarchy<DOM_DIM,AMB_DIM,Real,Int,SReal>
    {
        
#include "SimplexHierarchy_Details.hpp"
        
        constexpr Int ChildCount() const
        {
            return 1;
        }
        
        void ToChild( const Int k )
        {
            current_simplex_computed = false;
            
            former_child_id = child_id;
            child_id = k;
            column = ChildCount() * column + k;
            ++level;
            
            // Just don't subdivide.
        }
        
        void ToParent()
        {
            if( level > 0)
            {
                current_simplex_computed = false;

//                Int k = child_id;
                
                column = (column-child_id) / ChildCount();
                former_child_id = child_id;
                child_id = column % ChildCount();
                --level;
                
            }
            else
            {
                wprint("Level is equal to 0 already.");
            }
        }
        
    }; // SimplexHierarchy
    
#undef DOM_DIM
    
} // namespace Repulsor
