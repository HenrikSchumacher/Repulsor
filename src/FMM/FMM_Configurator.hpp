#pragma once

namespace Repulsor
{
    template<typename Real_, typename LInt_>
    struct MetricValueContainer
    {
        using Real = Real_;
        using LInt = LInt_;
        
        using Values_T = typename Tensors::Tensor2<Real,LInt>;
        
        Values_T FF;
        Values_T NF;
        Values_T VF;
        
        Values_T FF_diag;
        Values_T NF_diag;
        Values_T VF_diag;
    };
    
    template<typename ClusterTree_T_>
    class alignas( ObjectAlignment ) FMM_Configurator
    {
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename ClusterTree_T::Real;
        using SReal              = typename ClusterTree_T::SReal;
        using ExtReal            = typename ClusterTree_T::ExtReal;
        using Int                = typename ClusterTree_T::Int;
        using LInt               = typename ClusterTree_T::LInt;
        
        
        using ValueContainer_T   = MetricValueContainer<Real,LInt>;
        
        
    public:

        FMM_Configurator() = delete;
        
        FMM_Configurator(
            cref<ClusterTree_T> S_,
            cref<ClusterTree_T> T_,
            mref<ValueContainer_T> metric_values_
        )
        :   S             ( S_             )
        ,   T             ( T_             )
        ,   metric_values ( metric_values_ )
        {}
        
        // Copy constructor
        FMM_Configurator( cref<FMM_Configurator> other )
        :   S             ( other.S             )
        ,   T             ( other.T             )
        ,   metric_values ( other.metric_values )
        {}
        
        ~FMM_Configurator() = default;

    protected:

        cref<ClusterTree_T> S;
        cref<ClusterTree_T> T;
        
        mref<ValueContainer_T> metric_values;

    public:
        
        cref<ClusterTree_T> GetS() const
        {
            return S;
        }
        
        cref<ClusterTree_T> GetT() const
        {
            return T;
        }
        
        mref<ValueContainer_T> MetricValues()
        {
            return metric_values;
        }
        
    public:
        
        std::string ClassName() const
        {
            return "FMM_Configurator<"+S.ClassName()+">";
        }

    };

} // namespace Repulsor

