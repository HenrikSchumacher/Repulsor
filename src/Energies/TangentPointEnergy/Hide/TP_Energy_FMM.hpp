#pragma once

#define CLASS TP_Energy_FMM
#define BASE  Energy_FMM<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        CLASS(
            const Real alpha_,
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1)
        )
        :   BASE(
                 *Make_TP_Energy__NFK<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>( alpha_, beta_ ),
                 *Make_TP_Energy__FFK_FMM<AMB_DIM,DEGREE,Real,Int>( alpha_, beta_ ),
            weight
            )
        ,   alpha(alpha_)
        ,   beta(beta_)
        {}

        virtual ~CLASS() override = default;
                
    protected:
        
        const Real alpha;
        const Real beta;
        
    public:
        
        Real GetAlpha() const
        {
            return alpha;
        }
        
        Real GetBeta() const
        {
            return beta;
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
    
    
    
    template<int Degree, typename Real, typename Int, typename SReal, typename ExtReal>
    std::unique_ptr<EnergyBase<Real,Int,SReal,ExtReal>> CONCAT(Make_,CLASS) (
        const Int dom_dim, const Int amb_dim,
        const Real alpha, const Real beta, const Real weight
    )
    {
        EnergyBase<Real,Int,SReal,ExtReal> * r = nullptr;

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
                        r = new CLASS<1,2,Degree,Real,Int,SReal,ExtReal>(alpha, beta, weight);
                        break;
                    }
    #endif
    #if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<1,3,Degree,Real,Int,SReal,ExtReal>(alpha, beta, weight);
                        break;
                    }
    #endif
    #if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<1,4,Degree,Real,Int,SReal,ExtReal>(alpha, beta, weight);
                        break;
                    }
    #endif
                    default:
                    {
                        eprint("Make"+TO_STD_STRING(CLASS)+": "+TO_STD_STRING(CLASS) +"<"+ToString(dom_dim)+","+ToString(amb_dim)+","+ToString(Degree)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"> not implemented. Returning nullptr.");
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
                switch( amb_dim)
                {
    #if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<2,3,Degree,Real,Int,SReal,ExtReal>(
                            alpha, beta, weight
                        );
                        break;
                    }
    #endif
    #if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<2,4,Degree,Real,Int,SReal,ExtReal>(
                            alpha, beta, weight
                        );
                        break;
                    }
    #endif
                    default:
                    {
                        eprint("Make"+TO_STD_STRING(CLASS)+": "+TO_STD_STRING(CLASS) +"<"+ToString(dom_dim)+","+ToString(amb_dim)+","+ToString(Degree)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"> not implemented. Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
#endif
#if ENABLE_POINTCLOUNDS == 1
            case 0:
            {
                switch( amb_dim)
                {
    #if ENABLE_1D == 1
                    case 1:
                    {
                        r = new CLASS<0,1,Degree,Real,Int,SReal,ExtReal>(alpha, beta, weight);
                        break;
                    }
    #endif
    #if ENABLE_2D == 1
                    case 2:
                    {
                        r = new CLASS<0,2,Degree,Real,Int,SReal,ExtReal>(alpha, beta, weight);
                        break;
                    }
    #endif
    #if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<0,3,Degree,Real,Int,SReal,ExtReal>(alpha, beta, weight);
                        break;
                    }
    #endif
    #if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<0,4,Degree,Real,Int,SReal,ExtReal>( alpha, beta, weight);
                        break;
                    }
    #endif
                    default:
                    {
                        eprint("Make"+TO_STD_STRING(CLASS)+": "+TO_STD_STRING(CLASS) +"<"+ToString(dom_dim)+","+ToString(amb_dim)+","+ToString(Degree)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"> not implemented. Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
#endif
            default:
            {
                eprint("Make"+TO_STD_STRING(CLASS)+": "+TO_STD_STRING(CLASS) +"<"+ToString(dom_dim)+","+ToString(amb_dim)+","+ToString(Degree)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"> not implemented. Returning nullptr.");
                r = nullptr;
                break;
            }

        }
        
        return std::unique_ptr<EnergyBase<Real,Int,SReal,ExtReal>>(r);
    }

}// namespace Repulsor

#undef BASE
#undef CLASS
