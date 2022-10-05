#pragma once

#define CLASS CoulombEnergy__FFK_FMM
#define BASE  Energy__FFK_FMM<AMB_DIM,DEGREE,Real,Int>

namespace Repulsion
{
    template<int AMB_DIM, int DEGREE, typename T2, typename Real, typename Int>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_far;
        using BASE::S_D_far;
        
        using BASE::T_far;
        using BASE::T_D_far;
        
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
        
        REPULSION__ADD_CLONE_CODE(CLASS)
        
    protected:
        
        const Real beta;
        const T2 betahalf;
        const T2 minus_betahalf;
        const T2 minus_betahalf_minus_1;
        
    public:
        
        virtual Real Energy() override
        {
            Real v [AMB_DIM] = {};
            
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
            Real v [AMB_DIM] = {};
            
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
                const Real dEdv_i = H * v[i];
                
                dEdvx += dEdv_i * x[i];
                dEdvy += dEdv_i * y[i];
                
                DX[1+i] -= b * dEdv_i;
                DY[1+i] += a * dEdv_i;
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
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<T2>::Get()+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
    template<int AMB_DIM, int DEGREE, typename Real, typename Int>
    std::unique_ptr<BASE> CONCAT(Make_,CLASS)( const Real beta )
    {
        double betahalf_intpart;
        bool betahalf_int = (std::modf( beta/static_cast<double>(2), &betahalf_intpart) == 0.);
        
        BASE * r;
        
        if( betahalf_int )
        {
            r = new CLASS<AMB_DIM,DEGREE,Int,Real,Int>(static_cast<Real>(betahalf_intpart));
        }
        else
        {
            r = new CLASS<AMB_DIM,DEGREE,Real,Real,Int>(beta/2);
        }
        
        return std::unique_ptr<BASE>(r);
    }
    
} // namespace Repulsion

#undef CLASS
#undef BASE
