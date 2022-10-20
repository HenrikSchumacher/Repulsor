
//##############################################################################################
//      Factory-like code
//##############################################################################################

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    std::unique_ptr<ROOT> CONCAT(Make_,Area)(
        const int dom_dim,
        const int amb_dim,
        const ExtReal weight
    )
    {
        ptic("Make_"+TO_STD_STRING(CLASS)+" ("+ToString(dom_dim)+","+ToString(amb_dim)+")");
        
        ROOT * r = nullptr;
        
        switch(dom_dim)
        {
#if ENABLE_CURVES == 1
            case 1:
            {
                switch(amb_dim)
                {
#if ENABLE_2D == 1
                    case 2:
                    {
                        r = new CLASS<1,2,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
#if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<1,3,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
#if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<1,4,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
                    default:
                    {
                        eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported for domain dimension "+ToString(dom_dim)+". Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
#endif
#if ENABLE_SURFACES == 1
            case 2:
            {
                switch(amb_dim)
                {
#if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<2,3,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
#if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<2,4,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
                    default:
                    {
                        eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported for domain dimension "+ToString(dom_dim)+". Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
#endif
#if ENABLE_POINTCLOUDS == 1
            case 0:
            {
                switch(amb_dim)
                {
#if ENABLE_1D == 1
                    case 1:
                    {
                        r = new CLASS<0,1,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
#if ENABLE_2D == 1
                    case 2:
                    {
                        r = new CLASS<0,2,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
#if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<0,3,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
#if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<0,4,Real,Int,SReal,ExtReal>(weight);
                        break;
                    }
#endif
                    default:
                    {
                        eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported for domain dimension "+ToString(dom_dim)+". Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
#endif
            default:
            {
                eprint("Make_"+TO_STD_STRING(CLASS)+": domain dimension "+ToString(dom_dim)+" not supported. Returning nullptr.");
                r = nullptr;
                break;
            }
                
        }
        
        ptoc("Make_"+TO_STD_STRING(CLASS)+" ("+ToString(dom_dim)+","+ToString(amb_dim)+")");
        
        return std::unique_ptr<ROOT>(r);
    }
    
} // namespace Repulsor


#undef ROOT
