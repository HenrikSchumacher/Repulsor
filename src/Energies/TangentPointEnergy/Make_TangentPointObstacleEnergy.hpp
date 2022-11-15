
//##############################################################################################
//      Factory-like code
//##############################################################################################

namespace Repulsor
{
    
    template<int DOM_DIM_S, int DOM_DIM_T, typename Real, typename Int, typename SReal, typename ExtReal>
    std::unique_ptr<ROOT> CONCAT(Make_,CLASS)(
        const int amb_dim,
        const Real q,
        const Real p
    )
    {
        switch(amb_dim)
        {
#if ENABLE_1D == 1
            case 1:
            {
                return std::unique_ptr<ROOT>(
                    new CLASS<DOM_DIM_S,DOM_DIM_T,1,Real,Int,SReal,ExtReal>(q,p)
                );
            }
#endif
#if ENABLE_2D == 1
            case 2:
            {
                return std::unique_ptr<ROOT>(
                    new CLASS<DOM_DIM_S,DOM_DIM_T,2,Real,Int,SReal,ExtReal>(q,p)
                );
            }
#endif
#if ENABLE_3D == 1
            case 3:
            {
                return std::unique_ptr<ROOT>(
                    new CLASS<DOM_DIM_S,DOM_DIM_T,3,Real,Int,SReal,ExtReal>(q,p)
                );
            }
#endif
#if ENABLE_4D == 1
            case 4:
            {
                return std::unique_ptr<ROOT>(
                    new CLASS<DOM_DIM_S,DOM_DIM_T,4,Real,Int,SReal,ExtReal>(q,p)
                );
            }
#endif
            default:
            {
                eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported. Returning nullptr.");
                return std::unique_ptr<ROOT>(nullptr);
            }
                
        }
    }
    
    
    template<int DOM_DIM_S, typename Real, typename Int, typename SReal, typename ExtReal>
    std::unique_ptr<ROOT> CONCAT(Make_,CLASS)(
        const int dom_dim_T,
        const int amb_dim,
        const Real q,
        const Real p
    )
    {
        switch(dom_dim_T)
        {
#if ENABLE_CURVES == 1
            case 1:
            {
                return CONCAT(Make_,CLASS)<DOM_DIM_S,1,Real,Int,SReal,ExtReal>(amb_dim,q,p);
            }
#endif
#if ENABLE_SURFACES == 1
            case 2:
            {
                return CONCAT(Make_,CLASS)<DOM_DIM_S,2,Real,Int,SReal,ExtReal>(amb_dim,q,p);
            }
#endif
#if ENABLE_POINTCLOUDS == 1
            case 0:
            {
                return CONCAT(Make_,CLASS)<DOM_DIM_S,0,Real,Int,SReal,ExtReal>(amb_dim,q,p);
            }
#endif
            default:
            {
                eprint("Make_"+TO_STD_STRING(CLASS)+": domain dimension "+ToString(dom_dim_T)+" not supported. Returning nullptr.");
                return std::unique_ptr<ROOT>(nullptr);
            }
                
        }
    }
    
    
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    std::unique_ptr<ROOT> CONCAT(Make_,CLASS)(
        const int dom_dim_S,
        const int dom_dim_T,
        const int amb_dim,
        const Real q,
        const Real p
    )
    {
        switch(dom_dim_S)
        {
#if ENABLE_CURVES == 1
            case 1:
            {
                return  CONCAT(Make_,CLASS)<1,Real,Int,SReal,ExtReal>(dom_dim_T,amb_dim,q,p);
            }
#endif
#if ENABLE_SURFACES == 1
            case 2:
            {
                return  CONCAT(Make_,CLASS)<2,Real,Int,SReal,ExtReal>(dom_dim_T,amb_dim,q,p);
            }
#endif
#if ENABLE_POINTCLOUDS == 1
            case 0:
            {
                return CONCAT(Make_,CLASS)<0,Real,Int,SReal,ExtReal>(dom_dim_T,amb_dim,q,p);
            }
#endif
            default:
            {
                eprint("Make_"+TO_STD_STRING(CLASS)+": domain dimension "+ToString(dom_dim_S)+" not supported. Returning nullptr.");
                return std::unique_ptr<ROOT>( nullptr );
            }
        }
    }
    
    
    
} // namespace Repulsor


#undef ROOT

