#pragma once

#define CLASS Area_Kernel_VF
#define BASE  FMM_Kernel_VF<S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,hess_flag_,metric_flag_>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool hess_flag_, bool metric_flag_
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
        using BASE::hess_flag;
        using BASE::metric_flag;

        static constexpr Int METRIC_NNZ = 1 + 2 * AMB_DIM;
        static constexpr Int PREC_NNZ   = 1;
        
        using BASE::S;
        using BASE::T;
//        using BASE::S_Tree;
//        using BASE::T_Tree;
        
        using BASE::zero;
        using BASE::one;
        using BASE::two;
        
        static constexpr bool is_symmetric = is_symmetric_;
        
    public:
        
        CLASS() = delete;
        
        CLASS(
            Configurator_T & conf,
            const Real theta_, const Int max_level_
        )
        :   BASE                 (conf, theta_, max_level_ )
        {}
        
        CLASS( const CLASS & other )
        :   BASE                 (other                     )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        using BASE::metric_values;
        using BASE::prec_values;
        
        using BASE::S_data;
        using BASE::S_D_data;
        
        using BASE::T_data;
        using BASE::T_D_data;
        
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
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::S_Tree;
        using BASE::T_Tree;
        using BASE::S_scale;
        using BASE::T_scale;
        
        using BASE::lambda;
        using BASE::mu;
        
    public:
        
        virtual force_inline Real compute( const Int block_ID ) override
        {
            const Real w = S.Weight() * T.Weight();

            if constexpr ( diff_flag || metric_flag )
            {
                if constexpr ( diff_flag )
                {
                    DX[0] += static_cast<Real>(2) * b * w;
                    DY[0] += a * static_cast<Real>(2) * w;
                }
                
                if constexpr ( metric_flag )
                {
                    Real * restrict const m_vals = &metric_values[ METRIC_NNZ * block_ID ];
                }
            }
            
            if constexpr ( energy_flag )
            {
                return a * static_cast<Real>(2) * b;
            }
            else
            {
                return zero;
            }
        }
        
    public:
        
        virtual Int MetricNonzeroCount() const override
        {
            return METRIC_NNZ;
        }

        virtual Int PreconditionerNonzeroCount() const override
        {
            return PREC_NNZ;
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
            + ToString(is_symmetric) + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(hess_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }
    };

} // namespace Repulsor

#undef BASE
#undef CLASS


