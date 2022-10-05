#pragma once

#define CLASS TP_SingularMetric__NFK_Adaptive
#define BASE  Metric__NFK_Adaptive<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        static constexpr Int nnz_per_block = AMB_DIM+2;
        
        using BASE::S_near;
        using BASE::S_D_near;
        using BASE::S_serialized;
        
        using BASE::T_near;
        using BASE::T_D_near;
        using BASE::T_serialized;
        
        using BASE::a;
        using BASE::x_buffer;
        using BASE::p;
        
        using BASE::b;
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
        
        using BASE::lambda;
        using BASE::mu;

    public:
        
        using BASE::NearDimS;
        using BASE::NearDimT;
        using BASE::CoordDimS;
        using BASE::CoordDimT;
        using BASE::ProjectorDim;
        
    public:
        
        CLASS( const Real alpha_, const Real beta_,
               const AdaptivitySettings & settings_ = AdaptivitySettings()
        )
        :   BASE( settings_ )
        ,   alpha(alpha_)
        ,   beta(beta_)
        ,   alpha_half_minus_one(static_cast<Real>(0.5)*alpha_-static_cast<Real>(1.0))
        ,   minus_beta_half(static_cast<Real>(-0.5)*beta_)
        ,   exp_s( (beta_ - DOM_DIM) / alpha_ )
        ,   hi_exponent( static_cast<Real>(-0.5) * (static_cast<Real>(2) * (exp_s - static_cast<Real>(1)) + DOM_DIM) ) // multiplying by 0.5 because we use r * r instead of r for saving a sqrt
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   BASE        ( other )
        ,   alpha       ( other.alpha )
        ,   beta        ( other.beta )
        ,   alpha_half_minus_one ( other.alpha )
        ,   minus_beta_half      ( other.beta )
        ,   exp_s       ( other.exp_s )
        ,   hi_exponent ( other.hi_exponent )
        ,   hi_values   ( other.hi_values )
        {}
        
        // Move constructor
        CLASS( CLASS && other ) noexcept
        :   BASE        ( other )
        ,   alpha       ( other.alpha )
        ,   beta        ( other.beta )
        ,   alpha_half_minus_one ( other.alpha )
        ,   minus_beta_half      ( other.beta )
        ,   exp_s       ( other.exp_s )
        ,   hi_exponent ( other.hi_exponent )
        ,   hi_values   ( other.hi_values )
        {
            other.hi_values = nullptr;
        }
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
            
        const Real alpha;
        const Real beta;
        const Real alpha_half_minus_one;
        const Real minus_beta_half;
        
        const Real exp_s; // differentiability of the energy space
        const Real hi_exponent; // The only exponent we have to use for pow to compute matrix entries. All other exponents have been optimized away.
        
        Real * restrict hi_values = nullptr;
        
        Real hi_values_loc [nnz_per_block];
        
    public:
        
        virtual void AllocateValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & near_values, const Int nnz ) override
        {
            near_values[KernelType::MixedOrder]  = Tensor1<Real,Int> (nnz);
            
            hi_values = near_values[KernelType::MixedOrder].data();
        }
        
        virtual void LoadValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & near_values ) override
        {
            hi_values = near_values[KernelType::LowOrder].data();
        }
        
        virtual void StartEntry() override
        {
            fill_buffer(&hi_values_loc[0], nnz_per_block, static_cast<Real>(0) );
        }
        
        virtual void FinishEntry( const Int pos ) override
        {
            copy_buffer(&hi_values_loc[0], &hi_values[nnz_per_block*pos], nnz_per_block);
        }

        virtual void metric() override
        {
            Real x [AMB_DIM] = {};
            Real y [AMB_DIM] = {};
            Real v [AMB_DIM] = {};
            
            const Real delta = static_cast<Real>(S_ID == T_ID);
            
            Real r2 = delta;
            Real rCosPhi2 = static_cast<Real>(0);
            Real rCosPsi2 = static_cast<Real>(0);
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                x[i] = lambda[0] * x_buffer[AMB_DIM*0 +i];
                y[i] = mu    [0] * y_buffer[AMB_DIM*0 +i];

                for( Int ii = 1; ii < DOM_DIM+1; ++ii )
                {
                    x[i] += lambda[ii] * x_buffer[AMB_DIM*ii +i];
                }
                
                for( Int ii = 1; ii < DOM_DIM+1; ++ii )
                {
                    y[i] += mu    [ii] * y_buffer[AMB_DIM*ii +i];
                }
            }

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                Real Pv_i = static_cast<Real>(0);
                Real Qv_i = static_cast<Real>(0);
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    const Int k = lin_k[i][j];
                    Pv_i += p[k] * v[j];
                    Qv_i += q[k] * v[j];
                }
                rCosPhi2 += v[i] * Pv_i;
                rCosPsi2 += v[i] * Qv_i;
            }
            
            // Ordering of hi_values: { K_ij, K_ji, v[0],...,v[AMB_DIM-1] }
            const Real w = S.Weight() * T.Weight();
            
            const Real factor = w * (static_cast<Real>(1) -delta) * MyMath::pow(r2,minus_beta_half);
            
            hi_values_loc[0] += MyMath::pow(rCosPhi2, alpha_half_minus_one) * factor;
            hi_values_loc[1] += MyMath::pow(rCosPsi2, alpha_half_minus_one) * factor;
            
            copy_buffer( &v[0], &hi_values_loc[2], AMB_DIM );
        }
        
        virtual void TransposeBlock( const Int from, const Int to ) override
        {
            const Int from_ = nnz_per_block * from;
            const Int to_   = nnz_per_block * to;
            
            hi_values[to_  ] = hi_values[from_+1];
            hi_values[to_+1] = hi_values[from_  ];
            
            for( Int k = 2; k < AMB_DIM+2; ++k )
            {
                hi_values[to_+k] = - hi_values[from_+k];
            }
        }
        
    public:
        
        virtual std::string Stats() const override
        {
            return ClassName()+": alpha = "+ToString(alpha)+", beta = "+ToString(beta);
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef BASE
#undef CLASS

