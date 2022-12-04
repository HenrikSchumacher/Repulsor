#pragma once

#define BASE FMM_Kernel_NF<                              \
    S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,is_symmetric_,  \
    false,false,true                                     \
>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool high_order_
    >
    class PseudoLaplacian_Kernel_NF : public BASE
    {
    private:
        
        using Base_T = BASE;
        
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        using Values_T           = typename Base_T::Values_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;

        
        using Base_T::AMB_DIM;
        using Base_T::PROJ_DIM;
        using Base_T::S_DOM_DIM;
        using Base_T::T_DOM_DIM;
        using Base_T::S_COORD_DIM;
        using Base_T::T_COORD_DIM;
        using Base_T::T_DATA_DIM;
        using Base_T::S_DATA_DIM;
        using Base_T::energy_flag;
        using Base_T::diff_flag;
        using Base_T::metric_flag;
        
        static constexpr Real high_order = high_order_;
        static constexpr Int  ROWS       = high_order ? AMB_DIM : 1;
        static constexpr Int  COLS       = high_order ? AMB_DIM : 1;
        static constexpr Int  BLOCK_NNZ  = 1;
        static constexpr Int  DIAG_NNZ   = 1;
        
        using Base_T::S;
        using Base_T::T;
//        using Base_T::S_Tree;
//        using Base_T::T_Tree;
        
        using Base_T::zero;
        using Base_T::one;
        using Base_T::two;
        using Base_T::is_symmetric;
        
        
#include "../../FMM/FMM_Kernel_Common.hpp"
        
    public:
        
        PseudoLaplacian_Kernel_NF() = delete;
        
        PseudoLaplacian_Kernel_NF( Configurator_T & conf, const Real s_ )
        :   Base_T ( conf                            )
        ,   s      ( s_                              )
        ,   s_exp  ( -S_DOM_DIM/two - (s-high_order) )
        {}
        
        PseudoLaplacian_Kernel_NF( PseudoLaplacian_Kernel_NF & other )
        :   Base_T (other       )
        ,   s      (other.s     )
        ,   s_exp  (other.s_exp )
        {}
        
        ~PseudoLaplacian_Kernel_NF() = default;
        
    protected:
        
        using Base_T::metric_data;

        using Base_T::S_diag;
        using Base_T::T_diag;
        
        using Base_T::DX;
        using Base_T::DY;
        
        using Base_T::a;
        using Base_T::x;
        using Base_T::P;
        
        using Base_T::b;
        using Base_T::y;
        using Base_T::Q;
        
        using Base_T::lin_k;
        
        using Base_T::S_scale;
        using Base_T::T_scale;
        
        const Real s;
        const Real s_exp;
        
        Real ii_block [BLOCK_NNZ]  = {};
        Real jj_block [BLOCK_NNZ]  = {};
        
    public:
        
        
        force_inline Real Compute( const LInt k_global )
        {
            Real v [AMB_DIM ] = {};
            
            Real r2 = zero;
            
            for( Int l = 0; l < AMB_DIM; ++l )
            {
                v[l] = y[l] - x[l];
                r2  += v[l] * v[l];
            }

            const Real val = - MyMath::pow(r2, s_exp);
            
            const Real b_over_a = b/a;
            const Real a_over_b = a/b;
            
            metric_data[BLOCK_NNZ * k_global] = val;
            ii_block[0] -=   b_over_a * val;
            jj_block[0]  = - a_over_b * val;
            
            return 0;
        }
        
        
        force_inline void LoadS( const Int i_global )
        {
            this->loadS( i_global );
        }
        
        force_inline void WriteS( const Int i_global )
        {
            this->writeS( i_global );
            
            S_diag[DIAG_NNZ * i_global] += ii_block[0];
        }
        
        force_inline void LoadT( const Int j_global )
        {
            this->loadT( j_global );
        }

        force_inline void WriteT( const Int j_global )
        {
            T_diag[DIAG_NNZ * j_global] += jj_block[0];
        }
        
        std::string className() const
        {
            return "PseudoLaplacian_Kernel_NF<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
            + this->GetS().ClassName() + ","
            + ToString(is_symmetric) + ","
            + ToString(high_order)
            + ">";
        }
    };

} // namespace Repulsor

#undef BASE

