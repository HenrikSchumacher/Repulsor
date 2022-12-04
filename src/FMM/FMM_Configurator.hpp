#pragma once

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class alignas( OBJECT_ALIGNMENT ) FMM_Configurator
    {
    public:
        
        using ClusterTree_T       = ClusterTree_T_;
        
        using Real                = typename ClusterTree_T::Real;
        using SReal               = typename ClusterTree_T::SReal;
        using ExtReal             = typename ClusterTree_T::ExtReal;
        using Int                 = typename ClusterTree_T::Int;
        using LInt                = size_t;
        
        using Values_T           = Tensor2<Real,LInt>;
        using ValueContainer_T   = std::unordered_map<std::string,Values_T>;

        
    public:

        FMM_Configurator() = delete;
        
        FMM_Configurator( const ClusterTree_T & S_, const ClusterTree_T & T_, ValueContainer_T & metric_values_ )
        :   S             ( S_             )
        ,   T             ( T_             )
        ,   metric_values ( metric_values_ )
        {}
        
        // Copy constructor
        FMM_Configurator( const FMM_Configurator & other )
        :   S             ( other.S             )
        ,   T             ( other.T             )
        ,   metric_values ( other.metric_values )
        {}
        
        ~FMM_Configurator() = default;

    protected:

        const ClusterTree_T & S;
        const ClusterTree_T & T;
        
        ValueContainer_T & metric_values;

    public:
        
        const ClusterTree_T & GetS() const
        {
            return S;
        }
        
        const ClusterTree_T & GetT() const
        {
            return T;
        }
        
        ValueContainer_T & MetricValues()
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

