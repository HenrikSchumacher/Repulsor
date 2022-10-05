#pragma once

#define CLASS TP_SingularMetric__FFK_FMM
#define BASE  Metric__FFK_FMM<AMB_DIM,DEGREE,Real,Int>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int>
    class CLASS : public BASE
    {
    protected:
        
        static constexpr Int nnz_per_block = AMB_DIM + 2;
        
        using BASE::S_far;
        using BASE::T_far;
        
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

        using BASE::LoadS;
        using BASE::LoadT;
        
    public:
        
        CLASS( const Real alpha_, const Real beta_ )
        :   BASE()
        ,   alpha(alpha_)
        ,   beta(beta_)
        ,   alpha_half_minus_one(static_cast<Real>(0.5)*alpha_-static_cast<Real>(1.0))
        ,   minus_beta_half(static_cast<Real>(-0.5)*beta_)
        ,   exp_s( (beta_ - DOM_DIM) / alpha_ )
        ,   hi_exponent( static_cast<Real>(-0.5) * (static_cast<Real>(2.0) * (exp_s - static_cast<Real>(1)) + DOM_DIM) ) // multiplying by 0.5 because we use r * r instead of r for saving a sqrt
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
        ,   hi_values ( other.hi_values )
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
        
        const Real exp_s;
        const Real hi_exponent;
        
        Real * restrict hi_values = nullptr;
        
    public:

        virtual void AllocateValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & far_values, const Int nnz ) override
        {
            far_values[KernelType::MixedOrder]  = Tensor1<Real,Int> (nnz * nnz_per_block);
            
            hi_values = far_values[KernelType::MixedOrder].data();
        }
        
        virtual void LoadValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & far_values ) override
        {
            hi_values = far_values[KernelType::MixedOrder].data();
        }
        
        
        virtual void Metric( const Int pos ) override
        {
            Real v [AMB_DIM] = {};
            
            Real r2       = static_cast<Real>(0);
            Real rCosPhi2 = static_cast<Real>(0);
            Real rCosPsi2 = static_cast<Real>(0);
            
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
            
            const Real factor = MyMath::pow( r2, minus_beta_half );
            
            const Int loc = nnz_per_block * pos ;
            // Ordering of hi_values: { K_ij, K_ji, v[0],...,v[AMB_DIM-1] }
            
            hi_values[loc + 0] = MyMath::pow(rCosPhi2, alpha_half_minus_one) * factor;
            hi_values[loc + 1] = MyMath::pow(rCosPsi2, alpha_half_minus_one) * factor;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                hi_values[loc + 2 + k] = v[k];
            }
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
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef CLASS
#undef BASE
