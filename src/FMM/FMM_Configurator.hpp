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

        FMM_Configurator(
            cref<ClusterTree_T> S_,
            cref<ClusterTree_T> T_,
            mref<ValueContainer_T> metric_values_
        )
        :   S             ( S_             )
        ,   T             ( T_             )
        ,   metric_values ( metric_values_ )
        {}
        
        // Default constructor
        FMM_Configurator() = delete;
        // Destructor
        virtual ~FMM_Configurator() = default;
        // Copy constructor
        FMM_Configurator( const FMM_Configurator & other ) = default;
        // Copy assignment operator
        FMM_Configurator & operator=( const FMM_Configurator & other ) = default;
        // Move constructor
        FMM_Configurator( FMM_Configurator && other ) = default;
        // Move assignment operator
        FMM_Configurator & operator=( FMM_Configurator && other ) = default;

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

