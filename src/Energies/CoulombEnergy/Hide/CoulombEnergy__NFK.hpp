#pragma once

#define CLASS CoulombEnergy__NFK
#define BASE  Energy__NFK<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename T2, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_near;
        using BASE::S_D_near;
        
        using BASE::T_near;
        using BASE::T_D_near;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::p;
        
        using BASE::b;
        using BASE::y;
        using BASE::q;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::S_scale;
        using BASE::T_scale;
        
    public:
        
        explicit CLASS( const T2 betahalf_ )
        :
            BASE(),
            beta(static_cast<Real>(2)*betahalf_),
            betahalf(betahalf_),
            minus_betahalf(-betahalf),
            minus_betahalf_minus_1(-betahalf-1)
        {}

        // Copy constructor
        CLASS( const CLASS & other )
        :
            BASE( other ),
            beta(other.beta),
            betahalf(other.betahalf),
            minus_betahalf(other.minus_betahalf),
            minus_betahalf_minus_1(other.minus_betahalf_minus_1)
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        const Real beta;
        const T2 betahalf;
        const T2 minus_betahalf;
        const T2 minus_betahalf_minus_1;
                
    public:
        
        virtual Real Energy() override
        {
            Real v    [AMB_DIM] = {};
            
            Real r2 = static_cast<Real>(0);
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            const Real rMinusBetaMinus2 = MyMath::pow<Real,T2>( r2, minus_betahalf_minus_1 );

            const Real E = rMinusBetaMinus2 * r2;
            
            return a * E * b;
        }
        
        virtual Real DEnergy() override
        {
            Real v    [AMB_DIM] = {};
            Real dEdv [AMB_DIM] = {};
            
            Real r2 = static_cast<Real>(0);
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            const Real rMinusBetaMinus2 = MyMath::pow<Real,T2>( r2, minus_betahalf_minus_1 );

            const Real E = rMinusBetaMinus2 * r2;

            const Real H = - beta * rMinusBetaMinus2;
            
            Real dEdvx = static_cast<Real>(0);
            Real dEdvy = static_cast<Real>(0);
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                dEdv[i] = H * v[i];
                dEdvx += dEdv[i] * x[i];
                dEdvy += dEdv[i] * y[i];
                
                for( Int ii = 0; ii < DOM_DIM1+1; ++ii )
                {
                    DX[1+AMB_DIM*ii+i] -= dEdv[i] * S_scale * b;
                }
                
                for( Int ii = 0; ii < DOM_DIM2+1; ++ii )
                {
                    DY[1+AMB_DIM*ii+i] += dEdv[i] * T_scale * a;
                }
            }
        
            DX[0] += b * ( E + dEdvx );
            DY[0] += a * ( E - dEdvy );
            
            return a * E * b;
        }
        
    public:
        
        virtual std::string Stats() const override
        {
            return ClassName()+": beta = "+ToString(beta);
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<T2>::Get()+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+">";
        }
  
    };
    
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    std::unique_ptr<BASE> CONCAT(Make_,CLASS)( const Real beta )
    {
        double betahalf_intpart;
        bool betahalf_int = (std::modf( beta/static_cast<double>(2), &betahalf_intpart) == 0.);
        
        BASE * r;
        
        if( betahalf_int )
        {
            r = new CLASS<DOM_DIM1,DOM_DIM2,AMB_DIM,Int,Real,Int,SReal>(static_cast<Real>(betahalf_intpart));
        }
        else
        {
            r = new CLASS<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Real,Int,SReal>(beta/2);
        }
        
        return std::unique_ptr<BASE>(r);
    }
    
    
} // namespace Repulsor

#undef CLASS
#undef BASE
