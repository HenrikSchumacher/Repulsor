#pragma once

namespace Repulsor
{
    
#define DOM_DIM 1
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class alignas(ObjectAlignment) SimplexHierarchy<DOM_DIM,AMB_DIM,Real,Int,SReal>
    {
        
#include "SimplexHierarchy_Details.hpp"
        
        static constexpr Int ChildCount()
        {
            return 2;
        }
        
        static constexpr Level_T MaxLevel()
        {
            return static_cast<Level_T>( std::floor(
                std::log( static_cast<double>(std::numeric_limits<Column_T>::max())))
                /
                std::log(static_cast<double>(ChildCount())
            ) );
        }
        
        void ToChild( const Child_T k )
        {
            current_simplex_computed = false;
            
            former_child_id = child_id;
            child_id = k;
            column = ChildCount() * column + k;
            ++level;
            
            scale  *= half;
            weight *= half;
            
            if(k == 0)
            {
                corners[1][0] = corners[0][0] + half * (corners[1][0] - corners[0][0]);
                corners[1][1] = corners[0][1] + half * (corners[1][1] - corners[0][1]);
                
                center[0]= corners[0][0] + half * (center[0] - corners[0][0]);
                center[1]= corners[0][1] + half * (center[1] - corners[0][1]);
                
//                corners[1][0] = center[0];
//                corners[1][1] = center[1];
//
//                center[0] = half * ( corners[0][0] + center[0] );
//                center[1] = half * ( corners[0][1] + center[1] );
                
            }
            else // if( k == 1 )
            {
                corners[0][0] = corners[1][0] + half * (corners[0][0] - corners[1][0]);
                corners[0][1] = corners[1][1] + half * (corners[0][1] - corners[1][1]);
                
                center[0] = corners[1][0] + half * (center[0] - corners[1][0]);
                center[1] = corners[1][1] + half * (center[1] - corners[1][1]);
                
//                corners[0][0] = center[0];
//                corners[0][1] = center[1];
//
//                center[0] = half * ( corners[1][0] + center[0] );
//                center[1] = half * ( corners[1][1] + center[1] );
                

                
            }
        }
        
        void ToParent()
        {
            if( level > 0)
            {
                current_simplex_computed = false;
                
                const Child_T k = child_id;
                
                column = (column-child_id) / ChildCount();
                former_child_id = child_id;
                child_id = column % ChildCount();
                --level;
                
                scale  *= two;
                weight *= two;
                
                if( k == 0 )
                {
                    corners[1][0] = corners[0][0] + two * (corners[1][0] - corners[0][0]);
                    corners[1][1] = corners[0][1] + two * (corners[1][1] - corners[0][1]);
                    
                    center[0] = corners[0][0] + two * (center[0] - corners[0][0]);
                    center[1] = corners[0][1] + two * (center[1] - corners[0][1]);
                    
//                    center[0] = corners[1][0];
//                    center[1] = corners[1][1];
//
//                    corners[1][0] = corners[0][0] + two * (corners[1][0] - corners[0][0]);
//                    corners[1][1] = corners[0][1] + two * (corners[1][1] - corners[0][1]);
                }
                else // if( k == 1 )
                {
                    
                    corners[0][0] = corners[1][0] + two * (corners[0][0] - corners[1][0]);
                    corners[0][1] = corners[1][1] + two * (corners[0][1] - corners[1][1]);
                
                    center[0] = corners[1][0] + two * (center[0] - corners[1][0]);
                    center[1] = corners[1][1] + two * (center[1] - corners[1][1]);
                    
//                    center[0] = corners[0][0];
//                    center[1] = corners[0][1];
//                    
//                    corners[0][0] = corners[1][0] + two * (corners[0][0] - corners[1][0]);
//                    corners[0][1] = corners[1][1] + two * (corners[0][1] - corners[1][1]);
                }
            }
            else
            {
                wprint("Level is equal to 0 already.");
            }
        }
        
    }; // SimplexHierarchy
    
#undef DOM_DIM
    
} // namespace Repulsor
