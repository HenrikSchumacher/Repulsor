#pragma once

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool symmetricQ_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class alignas( OBJECT_ALIGNMENT ) FMM_Kernel
    {
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename ClusterTree_T::Real;
        using SReal              = typename ClusterTree_T::SReal;
        using ExtReal            = typename ClusterTree_T::ExtReal;
        using Int                = typename ClusterTree_T::Int;
        using LInt               = typename ClusterTree_T::LInt;

        using Configurator_T     = FMM_Configurator<ClusterTree_T>;
        using ValueContainer_T   = typename Configurator_T::ValueContainer_T;
        using Values_T           = typename ValueContainer_T::Values_T;
        
        static constexpr bool symmetricQ = symmetricQ_;
        static constexpr bool energy_flag  = energy_flag_;
        static constexpr bool diff_flag    = diff_flag_;
        static constexpr bool metric_flag  = metric_flag_;
        
        static constexpr Int AMB_DIM   = ClusterTree_T::AMB_DIM;
        static constexpr Int PROJ_DIM  = (AMB_DIM*(AMB_DIM+1))/2;
        
    protected:
        
        static constexpr Real zero = static_cast<Real>(0);
        static constexpr Real one  = static_cast<Real>(1);
        static constexpr Real two  = static_cast<Real>(2);
        
        static constexpr Real symmetry_factor = one / (one + !static_cast<Real>(symmetricQ) );
                                            
        const Int thread = 0;
        
    public:

        FMM_Kernel() = delete;
        
        explicit FMM_Kernel( Configurator_T & conf, const Int thread_ )
        :   thread        ( thread_              )
        ,   S             ( conf.GetS()          )
        ,   T             ( conf.GetT()          )
        ,   metric_values ( conf.MetricValues()  )   // In configure mode, kernels needs
        {
            debug_print(std::string( "Initializing " + this-> ClassName() + " from Configurator_T on thread " + ToString(thread)) );
        }
        
        FMM_Kernel( const FMM_Kernel & other, const Int thread_ )
        :   thread        ( thread_              )
        ,   S             ( other.S              )
        ,   T             ( other.T              )
        ,   metric_values ( other.metric_values  )
        {
            debug_print(std::string( "Initializing " + this->ClassName() + " from "+ClassName()+" on thread " + ToString(thread)) );
        }
        
        
        Int Thread() const
        {
            return thread;
        }
        
        std::string ThreadString() const
        {
            return "(" + ToString(thread) +")";
        }
        
        ~FMM_Kernel() = default;

    public:
        
        const ClusterTree_T & GetS() const
        {
            return S;
        }
        
        const ClusterTree_T & GetT() const
        {
            return T;
        }
        
        Int ThreadCount() const
        {
            return std::min( S.ThreadCount(), T.ThreadCount() );
        }
        
    protected:
        
        const ClusterTree_T & S;
        const ClusterTree_T & T;
        
        ValueContainer_T & metric_values;
        
    public:
        
        std::string ClassName() const
        {
//            return "FMM_Kernel<"+S.ClassName()+">";
            return "FMM_Kernel<...>";
        }

    };

} // namespace Repulsor
