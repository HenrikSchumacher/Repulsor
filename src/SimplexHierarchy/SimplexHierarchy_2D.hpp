#pragma once

namespace Repulsor
{
    template<int AMB_DIM_, typename Real, typename Int, typename SReal>
    class alignas(ObjectAlignment) SimplexHierarchy<2,AMB_DIM_,Real,Int,SReal>
    {
    public:
        
        static constexpr Int DOM_DIM = 2;
        static constexpr Int AMB_DIM = AMB_DIM_;
        
#include "SimplexHierarchy_Details.hpp"
        
    public:
        
        static constexpr Child_T ChildCount()
        {
            return 4;
        }
        
        static constexpr Level_T MaxLevel()
        {
            return static_cast<Level_T>( std::floor(
                std::log( static_cast<double>(std::numeric_limits<Column_T>::max())))
                /
                std::log(static_cast<double>(ChildCount())
            ) );
        }
        
    public:
        
        void ToChild( const Child_T k )
        {
            current_simplex_computed = false;
            
            former_child_id = child_id;
            child_id = k;
            column = ChildCount() * column + k;
            ++level;
            
            scale  *= Scalar::Half<SReal>;
            weight *= Scalar::Quarter<SReal>;
            
            if( (0<=k) && (k<=2) )
            {
                for( Int i = 0; i < VertexCount(); ++i )
                {
                    if( i != k )
                    {
                        for( Int j = 0; j < VertexCount(); ++j )
                        {
                            corners[i][j] = corners[k][j] + half * (corners[i][j] - corners[k][j]);
                        }
                    }
                    else
                    {
                        for( Int j = 0; j < VertexCount(); ++j )
                        {
                            center[j]     = corners[k][j] + half * (center[j]     - corners[k][j]);
                        }
                    }
                }
            }
            else if( child_id == 3 )
            {
                for( Int i = 0; i < VertexCount(); ++ i )
                {
                    for( Int j = 0; j < VertexCount(); ++j )
                    {
                        corners[i][j] = center[j] - half * (corners[i][j] - center[j]);
                    }
                }
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
                
                scale  *= Scalar::Two<SReal>;
                weight *= Scalar::Four<SReal>;
                
                if( (0<=k) && (k<=2) )
                {
                    for( Int i = 0; i < VertexCount(); ++i )
                    {
                        if( i != k )
                        {
                            for( Int j = 0; j < VertexCount(); ++j )
                            {
                                corners[i][j] = corners[k][j] + two * (corners[i][j] - corners[k][j]);
                            }
                        }
                        else
                        {
                            for( Int j = 0; j < VertexCount(); ++j )
                            {
                                center[j]     = corners[k][j] + two * (center[j]     - corners[k][j]);
                            }
                        }
                    }
                }
                else if( k == 3 )
                {
                    for( Int i = 0; i < VertexCount(); ++ i )
                    {
                        for( Int j = 0; j < VertexCount(); ++j )
                        {
                            corners[i][j] = center[j] - two * (corners[i][j] - center[j]);
                        }
                    }
                }
            }
            else
            {
                wprint("Level is equal to 0 already.");
            }
        }
        
    }; // SimplexHierarchy
    
} // namespace Repulsor
