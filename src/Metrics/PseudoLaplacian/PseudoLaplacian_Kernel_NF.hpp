#pragma once

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool symmetricQ_,
        bool high_order_
    >
    class PseudoLaplacian_Kernel_NF : public FMM_Kernel_NF<
        S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,symmetricQ_,
        false,false,true,false
    >
    {
    private:
        
        using Base_T = FMM_Kernel_NF<
            S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,symmetricQ_,
            false,false,true,false
        >;
        
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        using Values_T           = typename ValueContainer_T::Values_T;

        
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
        using Base_T::symmetricQ;
        
        
#include "../../FMM/FMM_Kernel_Common.hpp"
#include "../../FMM/FMM_Kernel_NF_Common.hpp"
        
    public:
        
        PseudoLaplacian_Kernel_NF( mref<Configurator_T> conf, const Int thread_, const Real s_ )
        :   Base_T ( conf, thread_                   )
        ,   s      ( s_                              )
        ,   s_exp  ( -S_DOM_DIM/two - (s-high_order) )
        {}
        
        PseudoLaplacian_Kernel_NF( mref<PseudoLaplacian_Kernel_NF> other, const Int thread_ )
        :   Base_T (other, thread_  )
        ,   s      (other.s         )
        ,   s_exp  (other.s_exp     )
        {}
        
        // Default constructor
        PseudoLaplacian_Kernel_NF() = delete;
        // Destructor
        virtual ~PseudoLaplacian_Kernel_NF() override = default;
        // Copy constructor
        PseudoLaplacian_Kernel_NF( const PseudoLaplacian_Kernel_NF & other ) = default;
        // Copy assignment operator
        PseudoLaplacian_Kernel_NF & operator=( const PseudoLaplacian_Kernel_NF & other ) = default;
        // Move constructor
        PseudoLaplacian_Kernel_NF( PseudoLaplacian_Kernel_NF && other ) = default;
        // Move assignment operator
        PseudoLaplacian_Kernel_NF & operator=( PseudoLaplacian_Kernel_NF && other ) = default;
        
    protected:
        
        using Base_T::metric_data;

        using Base_T::S_diag;
        using Base_T::T_diag;
        
        using Base_T::a;
        using Base_T::x;
        using Base_T::P;
        
        using Base_T::b;
        using Base_T::y;
        using Base_T::Q;
        
        using Base_T::S_scale;
        using Base_T::T_scale;
        
        const Real s;
        const Real s_exp;
        
        Real ii_block;
        Real jj_block;
        
    public:
        
        
        TOOLS_FORCE_INLINE Real Compute( const LInt k_global )
        {
            Real v [AMB_DIM ] = {};
            
            Real r2 = zero;
            
            for( Int l = 0; l < AMB_DIM; ++l )
            {
                v[l] = y[l] - x[l];
                r2  += v[l] * v[l];
            }

            const Real val = - Power(r2, s_exp);
            
            const Real b_over_a = b/a;
            const Real a_over_b = a/b;
            
            metric_data[BLOCK_NNZ * k_global] = val;
            ii_block -=   b_over_a * val;
            jj_block  = - a_over_b * val;
            
            return Scalar::Zero<Real>;
        }
        
        
        TOOLS_FORCE_INLINE void LoadS( const Int i_global )
        {
            this->loadS( i_global );
            
            ii_block = Scalar::Zero<Real>;
        }
        
        TOOLS_FORCE_INLINE void WriteS( const Int i_global )
        {
            this->writeS( i_global );
            
            S_diag[DIAG_NNZ * i_global] += ii_block;
        }
        
        TOOLS_FORCE_INLINE void LoadT( const Int j_global )
        {
            this->loadT( j_global );
            
            jj_block = Scalar::Zero<Real>;
        }

        TOOLS_FORCE_INLINE void WriteT( const Int j_global )
        {
            T_diag[DIAG_NNZ * j_global] += jj_block;
        }
        
        std::string className() const
        {
            return "PseudoLaplacian_Kernel_NF<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
//            + this->GetS().ClassName() + ","
            + "...,"
            + ToString(symmetricQ) + ","
            + ToString(high_order)
            + ">";
        }
    };

} // namespace Repulsor
