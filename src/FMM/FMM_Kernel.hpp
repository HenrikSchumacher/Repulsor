#pragma once

#define CLASS FMM_Kernel

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag, bool diff_flag, bool hess_flag, bool metric_flag
    >
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
        
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
        static constexpr bool is_symmetric = is_symmetric_;
        
        static constexpr Int AMB_DIM   = ClusterTree_T::AMB_DIM;
        static constexpr Int PROJ_DIM  = (AMB_DIM*(AMB_DIM+1))/2;
        
    protected:
        
        static constexpr Real zero = static_cast<Real>(0);
        static constexpr Real one  = static_cast<Real>(1);
        static constexpr Real two  = static_cast<Real>(2);
        
        static constexpr Real symmetry_factor = one / (one + static_cast<Real>(is_symmetric) );
                                                        
        Int tri_i [PROJ_DIM] = {};
        Int tri_j [PROJ_DIM] = {};
        Int lin_k [AMB_DIM][AMB_DIM] = {};
        
        void Init()
        {
            Int k = 0;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                lin_k[i][i] = k;
                tri_i[k] = i;
                tri_j[k] = i;
                ++k;
                for( Int j = i+1; j < AMB_DIM; ++j )
                {
                    tri_i[k] = i;
                    tri_j[k] = j;
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
//
//        CLASS()
//        {
//            Init();
//        };
        
        CLASS( const ClusterTree_T & S_, const ClusterTree_T & T_ )
        :   S( S_ )
        ,   T( T_ )
        {
            Init();
        }
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   S ( other.S )
        ,   T ( other.T )
        {
            Init();
        }
        
        virtual ~CLASS() = default;

    public:
        
        virtual void LoadS( const Int i ) = 0;
        
        virtual void LoadT( const Int j ) = 0;

//        void PrefetchS( const Int i ) const
//        {}
        
        virtual void PrefetchT( const Int j ) const = 0;
        
        virtual Real compute() = 0;
                                                       
        virtual Real Compute() = 0;
        
        virtual void WriteS() = 0;
        
        virtual void WriteT() = 0;
        
    protected:
        
        const ClusterTree_T & S;
        const ClusterTree_T & T;
        
        Int S_ID = -1;
        Int T_ID = -1;
        
    public:
        
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
