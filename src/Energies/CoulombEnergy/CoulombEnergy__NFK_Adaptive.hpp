#pragma once

#define CLASS CoulombEnergy__NFK_Adaptive
#define BASE  Energy__NFK_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>


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
        using BASE::x_buffer;
        using BASE::p;
        
        using BASE::b;
        using BASE::y;
        using BASE::y_buffer;
        using BASE::q;
        
        using BASE::settings;
        using BASE::S;
        using BASE::T;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
    public:
        
        explicit CLASS(
            const T2 betahalf_,
            const AdaptivitySettings & settings_ = AdaptivitySettings()
        )
        :
            BASE(settings_),
            beta(static_cast<Real>(2)*betahalf_),
            betahalf(betahalf_),
            minus_betahalf(-betahalf),
            minus_betahalf_minus_1(-betahalf-1)
        {}

        // Copy constructor
        CLASS( const CLASS & other )
        :
            BASE(other),
            beta(other.beta),
            betahalf(other.betahalf),
            minus_betahalf(-other.betahalf),
            minus_betahalf_minus_1(-other.betahalf-1)
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        using BASE::lambda;
        using BASE::mu;
        
        const Real beta;
        const T2 betahalf;
        const T2 minus_betahalf;
        const T2 minus_betahalf_minus_1;
    
    public:
        
        virtual Real energy() const override
        {
            Real v    [AMB_DIM] = {};
            Real r2 = static_cast<Real>(0);
            
            const Real w = S.Weight() * T.Weight();
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                x[k] = lambda[0] * x_buffer[AMB_DIM * 0 + k];
                y[k] =     mu[0] * y_buffer[AMB_DIM * 0 + k];
                
                for( Int kk = 1; kk < DOM_DIM1+1; ++kk )
                {
                    x[k] += lambda[kk] * x_buffer[AMB_DIM * kk + k];
                }
                
                for( Int kk = 1; kk < DOM_DIM2+1; ++kk )
                {
                    y[k] +=     mu[kk] * y_buffer[AMB_DIM * kk + k];
                }
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            const Real rMinusBetaMinus2 = MyMath::pow<Real,T2>( r2, minus_betahalf_minus_1 );

            const Real E = rMinusBetaMinus2 * r2;
            
            return w * a * E * b;
        }
    
        
        virtual Real denergy() const override
        {
            Real v    [AMB_DIM] = {};
            Real dEdv [AMB_DIM] = {};
            
            Real r2 = static_cast<Real>(0);
            
            const Real w = S.Weight() * T.Weight();
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                x[k] = lambda[0] * x_buffer[AMB_DIM * 0 + k];
                y[k] =     mu[0] * y_buffer[AMB_DIM * 0 + k];
                
                for( Int kk = 1; kk < DOM_DIM1+1; ++kk )
                {
                    x[k] += lambda[kk] * x_buffer[AMB_DIM * kk + k];
                }
                
                for( Int kk = 1; kk < DOM_DIM2+1; ++kk )
                {
                    y[k] +=     mu[kk] * y_buffer[AMB_DIM * kk + k];
                }
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            
//            this->center_stream << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
//            this->center_stream << y[0] << "\t" << y[1] << "\t" << y[2] << "\n";
//
//            const Real * restrict const X = this->S.SimplexSerialized()+1+AMB_DIM;
//
//            this->corner_stream
//            << X[0] << "\t" << X[1] << "\t" << X[2] << "\t"
//            << X[3] << "\t" << X[4] << "\t" << X[5] << "\t"
//            << X[6] << "\t" << X[7] << "\t" << X[8] << "\n";
//
//            const Real * restrict const Y = this->T.SimplexSerialized()+1+AMB_DIM;
//
//            this->corner_stream
//            << Y[0] << "\t" << Y[1] << "\t" << Y[2] << "\t"
//            << Y[3] << "\t" << Y[4] << "\t" << Y[5] << "\t"
//            << Y[6] << "\t" << Y[7] << "\t" << Y[8] << "\n";
            
            const Real rMinusBetaMinus2 = MyMath::pow<Real,T2>( r2, minus_betahalf_minus_1 );

            const Real E = rMinusBetaMinus2 * r2;

            const Real H = - beta * rMinusBetaMinus2;
            
            const Real wa = w * a;
            const Real wb = w * b;
            
            Real dEdvx = static_cast<Real>(0);
            Real dEdvy = static_cast<Real>(0);

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                dEdv[i] = H * v[i];
                dEdvx += dEdv[i] * x[i];
                dEdvy += dEdv[i] * y[i];
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int ii = 0; ii < DOM_DIM1+1; ++ii )
                {
                    DX[1+AMB_DIM*ii+i] -= wb * dEdv[i] * lambda[ii];
                }
                
                for( Int ii = 0; ii < DOM_DIM2+1; ++ii )
                {
                    DY[1+AMB_DIM*ii+i] += wa * dEdv[i] *     mu[ii];
                }
            }
            
            DX[0] += wb * ( E + dEdvx );
            DY[0] += wa * ( E - dEdvy );
            
            return wa * E * b;
        }
        
    public:
        
        virtual bool IsRepulsive() const override
        {
            return ( - beta + DOM_DIM1 + DOM_DIM2 < static_cast<Real>(0) );
        }
        
        virtual std::string Stats() const override
        {
            return ClassName()+": beta = "+ToString(beta)+", theta = "+ToString(settings.theta);
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<T2>::Get()+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+">";
        }
  
    };
    
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    std::unique_ptr<BASE> CONCAT(Make_,CLASS)(
        const Real beta,
        const AdaptivitySettings & settings = AdaptivitySettings()
    )
    {
        double betahalf_intpart;
        bool betahalf_int = (std::modf( beta/static_cast<double>(2), &betahalf_intpart) == 0.);
        
        BASE * r;
        
        if( betahalf_int )
        {
            r = new CLASS<DOM_DIM1,DOM_DIM2,AMB_DIM,Int,Real,Int,SReal>(static_cast<Real>(betahalf_intpart),settings);
        }
        else
        {
            r = new CLASS<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Real,Int,SReal>(beta/2,settings);
        }
        
        return std::unique_ptr<BASE>(r);
    }

} // namespace Repulsr

#undef CLASS
#undef BASE
