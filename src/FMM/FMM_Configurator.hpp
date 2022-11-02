#pragma once

#define CLASS FMM_Configurator

namespace Repulsor
{
    template< typename BlockClusterTree_T_ >
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
    public:
        
        using BlockClusterTree_T  = BlockClusterTree_T_;
        
        using ClusterTree_T       = typename BlockClusterTree_T::ClusterTree_T;
        using Values_T            = typename BlockClusterTree_T::Values_T;
        using ValueContainer_T    = typename BlockClusterTree_T::ValueContainer_T;
        
        using Real                = typename BlockClusterTree_T::Real;
        using SReal               = typename BlockClusterTree_T::SReal;
        using ExtReal             = typename BlockClusterTree_T::ExtReal;
        using Int                 = typename BlockClusterTree_T::Int;
        using LInt                = typename BlockClusterTree_T::LInt;
        
        
    public:

        CLASS() = delete;
        
        CLASS(
            const BlockClusterTree_T & bct_,
            Values_T & metric_values_
        )
        :   bct           ( bct_                )
        ,   metric_values ( metric_values_      )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   bct           ( other.bct           )
        ,   metric_values ( other.metric_values )
        {}
        
        ~CLASS() = default;

    protected:

        const BlockClusterTree_T & bct;
        
        Values_T & metric_values;

    public:

        const BlockClusterTree_T & GetBlockClusterTree() const
        {
            return bct;
        }
        
        const ClusterTree_T & GetS() const
        {
            return bct.GetS();
        }
        
        const ClusterTree_T & GetT() const
        {
            return bct.GetT();
        }
        
        Values_T & MetricValues() const
        {
            return metric_values;
        }
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+bct.ClassName()+">";
        }

    };

} // namespace Repulsor

#undef CLASS

