#pragma once

#define CLASS FMM_Configurator

namespace Repulsor
{
    template<
        typename ClusterTree_T_
    >
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
    public:
        
        using ClusterTree_T  = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
        using ValueContainer_T = Tensor2<Real,Int>;
        
    protected:
        
        
    public:

        CLASS() = delete;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_,
            ValueContainer_T & metric_values_,
            ValueContainer_T & prec_values_
        )
        :   S             ( S_                  )
        ,   T             ( T_                  )
        ,   metric_values ( metric_values_      )
        ,   prec_values   ( prec_values_        )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   S             ( other.S             )
        ,   T             ( other.T             )
        ,   metric_values ( other.metric_values )
        ,   prec_values   ( other.prec_values   )
        {}
        
        virtual ~CLASS() = default;

    protected:
        
        const ClusterTree_T & S;
        const ClusterTree_T & T;
        
        ValueContainer_T & metric_values;
        ValueContainer_T & prec_values;
        
    public:
        
        const ClusterTree_T & GetS() const
        {
            return S;
        }
        
        const ClusterTree_T & GetT() const
        {
            return T;
        }
        
        ValueContainer_T & MetricValues() const
        {
            return metric_values;
        }
        
        ValueContainer_T & PreconditionerValues() const
        {
            return prec_values;
        }
        
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+S.ClassName()+">";
        }

    };

} // namespace Repulsor

#undef CLASS

