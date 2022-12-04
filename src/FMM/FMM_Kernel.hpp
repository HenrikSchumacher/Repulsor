#pragma once

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool is_symmetric_,
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
        using LInt               = size_t;

        using Configurator_T     = FMM_Configurator<ClusterTree_T>;
        using Values_T           = typename Configurator_T::Values_T;
        using ValueContainer_T   = typename Configurator_T::ValueContainer_T;
        
        static constexpr bool is_symmetric = is_symmetric_;
        static constexpr bool energy_flag  = energy_flag_;
        static constexpr bool diff_flag    = diff_flag_;
        static constexpr bool metric_flag  = metric_flag_;
        
        static constexpr Int AMB_DIM   = ClusterTree_T::AMB_DIM;
        static constexpr Int PROJ_DIM  = (AMB_DIM*(AMB_DIM+1))/2;
        
    protected:
        
        static constexpr Real zero = static_cast<Real>(0);
        static constexpr Real one  = static_cast<Real>(1);
        static constexpr Real two  = static_cast<Real>(2);
        
        static constexpr Real symmetry_factor = one / (one + !static_cast<Real>(is_symmetric) );
                                                        
        Int tri_i [PROJ_DIM] = {};
        Int tri_j [PROJ_DIM] = {};
        Int lin_k [AMB_DIM][AMB_DIM] = {};
        
        void Init()
        {
            Int k = 0;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                lin_k[i][i] = k;
                tri_i[k]    = i;
                tri_j[k]    = i;
                ++k;
                for( Int j = i+1; j < AMB_DIM; ++j )
                {
                    tri_i[k]    = i;
                    tri_j[k]    = j;
                    lin_k[i][j] = lin_k[j][i] = k;
                    ++k;
                }
            }
        }
        
        void CheckInit()
        {
            print(ClassName()+"::CheckInit");
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    print("{ "+ToString(i)+"," +ToString(j)+" } -> " + ToString(lin_k[i][j]) );
                }
            }

            for( Int k = 0; k < PROJ_DIM; ++k )
            {
                print( ToString(k) +" -> { "+ToString(tri_i[k])+"," +ToString(tri_j[k])+" }");
            }
        }
        
        
    public:

        FMM_Kernel() = delete;
        
        explicit FMM_Kernel( Configurator_T & conf )
        :   S             ( conf.GetS() )
        ,   T             ( conf.GetT() )
        ,   metric_values ( conf.MetricValues()        )   // In configure mode, kernels needs
        {
            Init();
        }
        
        FMM_Kernel( const FMM_Kernel & other )
        :   S             ( other.S              )
        ,   T             ( other.T              )
        ,   metric_values ( other.metric_values  )
        {
            Init();
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
            return "FMM_Kernel<"+S.ClassName()+">";
        }

    };

} // namespace Repulsor
