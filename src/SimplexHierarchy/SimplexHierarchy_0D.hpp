#pragma once

namespace Repulsor
{
    
    template<int AMB_DIM_, typename Real, typename Int, typename SReal>
    class alignas(ObjectAlignment) SimplexHierarchy<0,AMB_DIM_,Real,Int,SReal>
    {
    public:
        
        static constexpr Int DOM_DIM = 0;
        static constexpr Int AMB_DIM = AMB_DIM_;
        
#include "SimplexHierarchy_Details.hpp"
        
        static constexpr Int ChildCount()
        {
            return 1;
        }
        
        static constexpr Level_T MaxLevel()
        {
            return 0;
        }
        
        void ToChild( const Child_T k )
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
    
} // namespace Repulsor
