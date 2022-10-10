#pragma once

#define CLASS TP_Metric__NFK
#define BASE  Metric__NFK<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>


namespace Repulsor
{
    
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_near;
        using BASE::T_near;
        
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

        using BASE::NearDimS;
        using BASE::NearDimT;
        
        using BASE::LoadS;
        using BASE::LoadT;
        
    public:
        
        static constexpr Int BLOCK_SIZE = AMB_DIM * (DOM_DIM+1);
        static constexpr Int VALUE_SIZE = 2;
        
        CLASS( const Real alpha_, const Real beta_ )
        :
            BASE(),
            alpha(alpha_),
            beta(beta_),
            exp_s( (beta_ - DOM_DIM) / alpha_ ),
            hi_exponent( static_cast<Real>(-0.5) * (static_cast<Real>(2) * (exp_s - static_cast<Real>(1)) + DOM_DIM) ) // multiplying by 0.5 because we use r * r instead of r for saving a sqrt
        {}

        // Copy constructor
        CLASS( const CLASS & other )
        :
            BASE        ( other ),
            alpha       ( other.alpha ),
            beta        ( other.beta ),
            exp_s       ( other.exp_s ),
            hi_exponent ( other.hi_exponent )
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        const Real alpha;
        const Real beta;
        
        const Real exp_s; // differentiability of the energy space
        const Real hi_exponent; // The only exponent we have to use for pow to compute matrix entries. All other exponents have been optimized away.
        
        Real * restrict fr_values = nullptr;
        Real * restrict hi_values = nullptr;
        Real * restrict lo_values = nullptr;
        
        Real * restrict X = nullptr;
        Real * restrict Y = nullptr;
        
        Real z [BLOCK_SIZE] = {};
                
    public:

        virtual void AllocateValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & near_values, const Int nnz ) override
        {
            near_values[KernelType::FractionalOnly] = Tensor1<Real,Int> (nnz);
            near_values[KernelType::HighOrder]      = Tensor1<Real,Int> (nnz);
            near_values[KernelType::LowOrder]       = Tensor1<Real,Int> (nnz);
            
            fr_values = near_values[KernelType::FractionalOnly].data();
            hi_values = near_values[KernelType::HighOrder].data();
            lo_values = near_values[KernelType::LowOrder].data();
        }
        
        virtual void LoadValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & near_values ) override
        {
            fr_values = near_values[KernelType::FractionalOnly].data();
            hi_values = near_values[KernelType::HighOrder].data();
            lo_values = near_values[KernelType::LowOrder].data();
        }
        
        virtual void ComputeBlock( const Int pos ) override
        {
            Real v [AMB_DIM] = {};
            
            const Real delta = static_cast<Real>(S_ID == T_ID);
            
            Real r2 = delta;
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

            const Real r4 = r2 * r2;
            const Real mul = r4 * MyMath::pow<Real,Int>( r2, static_cast<Int>(DOM_DIM-1) );

            // The following line makes up approx 2/3 of this function's runtime! This is why we avoid pow as much as possible and replace it with MyMath::pow.;
            // I got it down to this single call to pow. We might want to generate a lookup table for it...;
            // The factor of (-2.) is here, because we assemble the _metric_, not the kernel.;
            const Real power = static_cast<Real>(-2) * MyMath::pow(r2, hi_exponent);

            const Real factor = (static_cast<Real>(1) - delta);
            
            const Real hi = factor* power;
            hi_values[pos] = hi;

            const Real fr = static_cast<Real>(4) * factor / (power * mul);
            fr_values[pos] = fr;

            const Real lo = static_cast<Real>(0.5) * factor * (rCosPhi2 + rCosPsi2) / r4 * power;
            lo_values[pos] = lo;
        }
        
        virtual void TransposeBlock( const Int from, const Int to ) override
        {
            hi_values[to] = hi_values[from];
            fr_values[to] = fr_values[from];
            lo_values[to] = lo_values[from];
        }
        
        virtual void ApplyBlock( const Real alpha, const Int pos, const Int j ) override
        {
            
        }
                
        virtual void ClearVector() override
        {
            
        }
        
        virtual void WriteVector( const Int i ) = 0;
        
        virtual void LoadOutputBuffer( const Real * restrict const Y ) = 0;
        
        virtual void LoadInputBuffer ( const Real * restrict const X ) = 0;
        
    public:
        
        virtual std::string Stats() const override
        {
            return ClassName()+": alpha = "+ToString(alpha)+", beta = "+ToString(beta);
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef CLASS
#undef BASE
