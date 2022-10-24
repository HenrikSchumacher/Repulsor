#pragma once

#define CLASS TP_Kernel_VF
#define BASE  FMM_Kernel_VF<S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_, typename T1, typename T2,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        using typename BASE::Configurator_T;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::S_DOM_DIM;
        using BASE::T_DOM_DIM;
        using BASE::S_COORD_DIM;
        using BASE::T_COORD_DIM;
        using BASE::T_DATA_DIM;
        using BASE::S_DATA_DIM;
        using BASE::energy_flag;
        using BASE::diff_flag;
        using BASE::metric_flag;
        
        static constexpr Int ROWS      = 1 + AMB_DIM;
        static constexpr Int COLS      = 1 + AMB_DIM;
        static constexpr Int BLOCK_NNZ = 1 + 2 * AMB_DIM;
        static constexpr Int DIAG_NNZ  = ROWS * COLS;
        
        using BASE::S;
        using BASE::T;
        
        using BASE::zero;
        using BASE::one;
        using BASE::two;
        
        static constexpr bool is_symmetric = is_symmetric_;
        
    public:
        
        CLASS() = delete;
        
        CLASS(
            Configurator_T & conf,
            const Real theta_, const Int max_level_, const T1 q_half_, const T2 p_half_
        )
        :   BASE                 (conf, theta_, max_level_ )
        ,   q                    (two*q_half_)
        ,   q_half               (q_half_    )
        ,   q_half_minus_1       (q_half-1   )
        ,   p                    (two*p_half_)
        ,   p_half               (p_half_    )
        ,   minus_p_half         (-p_half    )
        ,   minus_p_half_minus_1 (-p_half-1  )
        {}
        
        CLASS( const CLASS & other )
        :   BASE                 (other                     )
        ,   q                    (other.q                   )
        ,   q_half               (other.q_half              )
        ,   q_half_minus_1       (other.q_half_minus_1      )
        ,   p                    (other.p                   )
        ,   p_half               (other.p_half              )
        ,   minus_p_half         (other.minus_p_half        )
        ,   minus_p_half_minus_1 (other.minus_p_half_minus_1)
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        using BASE::metric_data;
        
        using BASE::S_data;
        using BASE::S_D_data;
        
        using BASE::T_data;
        using BASE::T_D_data;
        
        using BASE::S_diag;
        using BASE::T_diag;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::x_buffer;
        using BASE::P;
        
        using BASE::b;
        using BASE::y;
        using BASE::y_buffer;
        using BASE::Q;
        
        using BASE::i_global;
        using BASE::j_global;
        using BASE::k_global;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::S_Tree;
        using BASE::T_Tree;
        
        using BASE::lambda;
        using BASE::mu;
        
        const Real q;
        const T1   q_half;
        const T1   q_half_minus_1;
        
        const Real p;
        const T2   p_half;
        const T2   minus_p_half;
        const T2   minus_p_half_minus_1;
        
        
        Real ij_block [BLOCK_NNZ];
        Real ii_block [ROWS][COLS];
        Real jj_block [ROWS][COLS];
        
    protected:
        
        virtual Real compute() override
        {
            Real x    [AMB_DIM]  = {};
            Real y    [AMB_DIM]  = {};
            Real v    [AMB_DIM]  = {};
            Real Pv   [AMB_DIM]  = {};
            Real Qv   [AMB_DIM]  = {};
            Real dEdv [AMB_DIM]  = {};
            Real V    [PROJ_DIM] = {};
            
            Real r2        = zero;
            Real rCosPhi_2 = zero;
            Real rCosPsi_2 = zero;
            
            const Real w = S_Tree.Weight() * T_Tree.Weight();
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                x[i] = lambda[0] * x_buffer[AMB_DIM*0 +i];
                y[i] = mu    [0] * y_buffer[AMB_DIM*0 +i];

                for( Int ii = 1; ii < S_DOM_DIM+1; ++ii )
                {
                    x[i] += lambda[ii] * x_buffer[AMB_DIM*ii +i];
                }
                
                for( Int ii = 1; ii < T_DOM_DIM+1; ++ii )
                {
                    y[i] += mu    [ii] * y_buffer[AMB_DIM*ii +i];
                }
                
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    const Int k = lin_k[i][j];
                    Pv[i] += P[k] * v[j];
                    Qv[i] += Q[k] * v[j];
                    if( j >= i )
                    {
                        V[k] = ( static_cast<Real>(1) + static_cast<Real>(i!=j) )*v[i]*v[j];
                    }
                }
                rCosPhi_2 += v[i] * Pv[i];
                rCosPsi_2 += v[i] * Qv[i];
            }
            
            // TODO: It probably suffices to multiply w only once against r_minus_p_minus_2, because the latter is multiplied to almost everything.
            
            // |P*(y-x)|^{q-2}
            const Real rCosPhi_q_minus_2 = MyMath::pow<Real,T1>( fabs(rCosPhi_2), q_half_minus_1);
            // |Q*(y-x)|^{q-2}
            const Real rCosPsi_q_minus_2 = MyMath::pow<Real,T1>( fabs(rCosPsi_2), q_half_minus_1);
            // r^{2-p}
            const Real r_minus_p_minus_2 = MyMath::pow<Real,T2>( r2, minus_p_half_minus_1 );
            // |y-x|^-p
            const Real r_minus_p = r_minus_p_minus_2 * r2;
            // |P*(y-x)|^q
            const Real rCosPhi_q = rCosPhi_q_minus_2 * rCosPhi_2;
            // |Q*(y-x)|^q
            const Real rCosPsi_q = rCosPsi_q_minus_2 * rCosPsi_2;
            
            const Real Num = ( rCosPhi_q + rCosPsi_q );

            const Real E = Num * r_minus_p;

            if constexpr ( diff_flag || metric_flag )
            {
                // Needed for both the differential and the metric.
                
                const Real factor = q_half * r_minus_p;         // = q / |y-x|^p
                
                const Real K_xy = factor * rCosPhi_q_minus_2;   // = |P*(y-x)|^(q-2) / |y-x|^p
                const Real K_yx = factor * rCosPsi_q_minus_2;   // = |Q*(y-x)|^(q-2) / |y-x|^p

                if constexpr ( diff_flag )
                {
                    const Real H = - p * r_minus_p_minus_2 * Num;
                    
                    Real dEdvx = zero;
                    Real dEdvy = zero;
                    
                    const Real wa = w * a;
                    const Real wb = w * b;
                    
                    for( Int i = 0; i < AMB_DIM; ++i )
                    {
                        dEdv[i] = two * ( K_xy * Pv[i] + K_yx * Qv[i] ) + H * v[i];
                        dEdvx  += dEdv[i] * x[i];
                        dEdvy  += dEdv[i] * y[i];
                        
                        for( Int ii = 0; ii < S_DOM_DIM+1; ++ii )
                        {
                            DX[1+AMB_DIM*ii+i] -= lambda[ii] * wb * dEdv[i];
                        }
                        
                        for( Int ii = 0; ii < T_DOM_DIM+1; ++ii )
                        {
                            DY[1+AMB_DIM*ii+i] +=     mu[ii] * wa * dEdv[i];
                        }
                    }
                    
                    DX[0] += wb * ( E - factor * rCosPhi_q + dEdvx );
                    DY[0] += wa * ( E - factor * rCosPsi_q - dEdvy );
                    
                    const Real w_b_K_xy = wb * K_xy;
                    const Real w_a_K_yx = wa * K_yx;
                    
                    for( Int k = 0; k < PROJ_DIM; ++k )
                    {
                        DX[1+S_COORD_DIM+k] += w_b_K_xy * V[k];
                        DY[1+T_COORD_DIM+k] += w_a_K_yx * V[k];
                    }
                }
                
                if constexpr ( metric_flag )
                {
                    const Real w_K_xy = w * K_xy;
                    const Real w_K_yx = w * K_yx;
// ij_block
//    /                                                                 \
//    |   - K_xy - K_yx    K_yx * v[0]    K_yx * v[1]    K_yx * v[2]    |
//    |                                                                 |
//    |   - K_xy * v[0]         0              0              0         |
//    |                                                                 |
//    |   - K_xy * v[1]         0              0              0         |
//    |                                                                 |
//    |   - K_xy * v[2]         0              0              0         |
//    \                                                                 /
//
//    This are 1 + 2 * AMB_DIM nonzero values.
//    It is tempting to compress this into 2 + AMB_DIM values.
//    But it did not appear to be faster

                    
// ii block
//    /                                                                                 \
//    |     K_xy + K_yx      K_xy * v[0]         K_xy * v[1]         K_xy * v[2]        |
//    |                                                                                 |
//    |     K_xy * v[0]   K_xy * v[0] * v[0]  K_xy * v[0] * v[1]  K_xy * v[0] * v[2]    |
//    |                                                                                 |
//    |     K_xy * v[1]   K_xy * v[1] * v[0]  K_xy * v[1] * v[1]  K_xy * v[1] * v[2]    |
//    |                                                                                 |
//    |     K_xy * v[2]   K_xy * v[2] * v[0]  K_xy * v[2] * v[1]  K_xy * v[2] * v[2]    |
//    \                                                                                 /
                    
// jj block
//
//    /                                                                                 \
//    |     K_xy + K_yx     - K_yx * v[0]       - K_yx * v[1]       - K_yx * v[2]       |
//    |                                                                                 |
//    |   - K_yx * v[0]   K_yx * v[0] * v[0]  K_yx * v[0] * v[1]  K_yx * v[0] * v[2]    |
//    |                                                                                 |
//    |   - K_yx * v[1]   K_yx * v[1] * v[0]  K_yx * v[1] * v[1]  K_yx * v[1] * v[2]    |
//    |                                                                                 |
//    |   - K_yx * v[2]   K_yx * v[2] * v[0]  K_yx * v[2] * v[1]  K_yx * v[2] * v[2]    |
//    \                                                                                 /
        
                    
                    const Real w_K_sym = w_K_xy + w_K_yx;
                    
                    ij_block[0]    -= w_K_sym;
                    ii_block[0][0] += w_K_sym;
                    jj_block[0][0] += w_K_sym;
                    
                    for( Int l = 0; l < AMB_DIM; ++l )
                    {
                        const Real w_K_xy_v = w_K_xy * v[l];
                        const Real w_K_yx_v = w_K_yx * v[l];
                        
                        ij_block[1+l]         += w_K_yx_v;
                        ij_block[1+AMB_DIM+l] -= w_K_xy_v;
                        
                        // store only upper triangle
                        ii_block[0][1+l] +=   w_K_xy_v;
                        jj_block[0][1+l] -=   w_K_yx_v;
                    }
                    
                    // store only upper triangle
                    for( Int i = 1; i < ROWS; ++i )
                    {
                        for( Int j = i; j < COLS; ++j )
                        {
                            const Real vv = v[i-1] * v[j-1];
                            
                            ii_block[i][j] += K_xy * vv;
                            jj_block[i][j] += K_yx * vv;
                        }
                    }
                }
                
            }
            
            if constexpr ( energy_flag )
            {
                return  a * E * b;
            }
            else
            {
                return zero;
            }
        }
        
    protected:
    
        
        virtual void loadS() override
        {
            if constexpr ( metric_flag )
            {
                zerofy_buffer( &ii_block[0][0], DIAG_NNZ );
            }
        }
        
        virtual void writeS() override
        {
            if constexpr ( metric_flag )
            {
                copy_buffer( &ii_block[0][0], &S_diag[DIAG_NNZ * i_global], DIAG_NNZ );
            }
        }
        
        virtual void loadT() override
        {
            if constexpr ( metric_flag )
            {
                zerofy_buffer( &ij_block[0], BLOCK_NNZ );
                
                zerofy_buffer( &jj_block[0][0], DIAG_NNZ );
            }
        }
        
        virtual void writeT() override
        {
            if constexpr ( metric_flag )
            {
                copy_buffer( &ij_block[0], &metric_data[BLOCK_NNZ * k_global], BLOCK_NNZ );
                
                add_to_buffer( &jj_block[0][0], &T_diag[DIAG_NNZ * j_global], DIAG_NNZ );
            }
        }
        
    public:
        
        virtual Int NonzeroCount() const override
        {
            return BLOCK_NNZ;
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
            + S.ClassName() + ","
            + TypeName<T1>::Get() + ","
            + TypeName<T2>::Get() + ","
            + ToString(is_symmetric) + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) +
            ">";
        }
    };

} // namespace Repulsor

#undef BASE
#undef CLASS


