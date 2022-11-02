#pragma once

#define CLASS FMM_Kernel

namespace Repulsor
{
    template<
        typename BlockClusterTree_T_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
    public:
        
        using BlockClusterTree_T = BlockClusterTree_T_;
        
        using ClusterTree_T      = typename BlockClusterTree_T::ClusterTree_T;
        using Values_T           = typename BlockClusterTree_T::Values_T;
        using ValueContainer_T   = typename BlockClusterTree_T::ValueContainer_T;
        
        using Real               = typename BlockClusterTree_T::Real;
        using SReal              = typename BlockClusterTree_T::SReal;
        using ExtReal            = typename BlockClusterTree_T::ExtReal;
        using Int                = typename BlockClusterTree_T::Int;
        using LInt               = typename BlockClusterTree_T::LInt;
        
        using Configurator_T     = FMM_Configurator<BlockClusterTree_T>;
        
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
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

        CLASS() = delete;
        
        CLASS( Configurator_T & conf )
        :   bct           ( conf.GetBlockClusterTree() )
        ,   metric_values ( conf.MetricValues()        )   // In configure mode, kernels needs
        {
            Init();
        }
        
        CLASS( const CLASS & other )
        :   bct           ( other.bct            )
        ,   metric_values ( other.metric_values  )
        // In compute mode the pointers are needed!
        ,   metric_data   ( metric_values.data() )
        {
            Init();
        }
        
        ~CLASS() = default;

    public:
        const ClusterTree_T & GetS() const
        {
            return bct.GetS();
        }
        
        const ClusterTree_T & GetT() const
        {
            return bct.GetT();
        }
        
    protected:
        
        const BlockClusterTree_T & bct;
        
        Values_T & metric_values;
        
        Real * restrict const metric_data = nullptr;
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+bct.ClassName()+">";
        }

    };

} // namespace Repulsor

#undef CLASS
