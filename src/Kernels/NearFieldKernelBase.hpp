#pragma once

#define CLASS NearFieldKernelBase

// This sets up the basic I/O routines for all near field kernels, both for energies and metrics.

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        
    public:
        
        static constexpr Int PROJ_DIM   = (AMB_DIM*(AMB_DIM+1))/2;
        static constexpr Int COORD_DIMS = (DOM_DIM1+1)*AMB_DIM;
        static constexpr Int COORD_DIMT = (DOM_DIM2+1)*AMB_DIM;
        static constexpr Int NEAR_DIMS  = 1 + COORD_DIMS + PROJ_DIM;
        static constexpr Int NEAR_DIMT  = 1 + COORD_DIMT + PROJ_DIM;
        
        using S_TREE_T = SimplexHierarchy<DOM_DIM1,AMB_DIM,GJK_Real,Int,SReal>;
        using T_TREE_T = SimplexHierarchy<DOM_DIM2,AMB_DIM,GJK_Real,Int,SReal>;
        
        using GJK_T    = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
        CLASS()
        {
            Init();
        }
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   S_near        ( other.S_near )
        ,   S_D_near      ( other.S_D_near )
        ,   S_serialized  ( other.S_serialized  )
        ,   T_near        ( other.T_near )
        ,   T_D_near      ( other.T_D_near )
        ,   T_serialized  ( other.T_serialized  )
        {
            copy_buffer( &other.tri_i[0],    &tri_i[0],    PROJ_DIM        );
            copy_buffer( &other.tri_j[0],    &tri_j[0],    PROJ_DIM        );
            copy_buffer( &other.lin_k[0][0], &lin_k[0][0], AMB_DIM*AMB_DIM );
        }
        
        virtual ~CLASS() = default;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(CLASS)
        
    public:
        
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
        

    protected:
        
        const  Real * restrict S_near       = nullptr;
               Real * restrict S_D_near     = nullptr;
        const SReal * restrict S_serialized = nullptr;

        
        const  Real * restrict T_near       = nullptr;
               Real * restrict T_D_near     = nullptr;
        const SReal * restrict T_serialized = nullptr;
        
        mutable Real a = static_cast<Real>(0);
        mutable Real b = static_cast<Real>(0);
        
        
#ifdef NearField_S_Copy
        mutable Real x [AMB_DIM] = {};
        mutable Real p [PROJ_DIM] = {};
        mutable Real x_buffer [NEAR_DIMS] = {};
#else
        mutable Real x [AMB_DIM] = {};
        mutable Real const * restrict p = nullptr;
        mutable Real const * restrict x_buffer = nullptr;
#endif
        
#ifdef NearField_T_Copy
        mutable Real y [AMB_DIM] = {};
        mutable Real q [PROJ_DIM] = {};
        mutable Real y_buffer [NEAR_DIMT] = {};
#else
        mutable Real y [AMB_DIM] = {};
        mutable Real const * restrict q = nullptr;
        mutable Real const * restrict y_buffer = nullptr;
#endif
        
        mutable Real DX [NEAR_DIMS] = {};
        mutable Real DY [NEAR_DIMT] = {};
        
        Int tri_i [PROJ_DIM] = {};
        Int tri_j [PROJ_DIM] = {};
        Int lin_k [AMB_DIM][AMB_DIM] = {};
        
        mutable Int S_ID = -1;
        mutable Int T_ID = -1;
        
        static const constexpr Real S_scale = static_cast<Real>(1)/static_cast<Real>(DOM_DIM1+1);
        static const constexpr Real T_scale = static_cast<Real>(1)/static_cast<Real>(DOM_DIM2+1);

        mutable S_TREE_T S;
        mutable T_TREE_T T;
        
        const SReal * restrict const lambda = S.Center();
        const SReal * restrict const mu     = T.Center();
        
        mutable GJK_T gjk;
        
    public:
        
        bool LoadNearField(
            const Tensor2<Real,Int> & S_near_,
            const Tensor2<Real,Int> & T_near_
        )
        {
            bool check = true;
            if( (S_near_.Dimension(1) != NEAR_DIMS) || (T_near_.Dimension(1) != NEAR_DIMT) )
            {
                eprint(ClassName()+"::LoadNearField : (S_near_.Dimension(1) != NearDimS()) || (T_near_.Dimension(1) != NearDimT()).");
                check = false;
            }
            else
            {
                S_near = S_near_.data();
                T_near = T_near_.data();
            }
            return check;
        }

        
        bool LoadDNearField(
            Tensor2<Real,Int> & S_D_near_,
            Tensor2<Real,Int> & T_D_near_
        )
        {
            bool check = true;
            if( (S_D_near_.Dimension(1) != NEAR_DIMS) || (T_D_near_.Dimension(1) != NearDimT()) )
            {
                eprint(ClassName()+"::LoadDNearField : (S_D_near_.Dimension(1) != NEAR_DIMS) || (T_D_near_.Dimension(1) != NearDimT()).");
                check = false;
            }
            else
            {
                S_D_near = S_D_near_.data();
                T_D_near = T_D_near_.data();
            }
            return check;
        }
        
        bool LoadPrimitiveSerializedData(
            const Tensor2<SReal,Int> & S_serialized_,
            const Tensor2<SReal,Int> & T_serialized_
        )
        {
            bool check = true;
            
            if( (S_serialized_.Dimension(1) != S.SimplexPrototype().Size()) || (T_serialized_.Dimension(1) != T.SimplexPrototype().Size()) )
            {
                check = false;
                eprint(ClassName()+"::LoadPrimitiveSerializedData.");
            }
            else
            {
                S_serialized = S_serialized_.data();
                T_serialized = T_serialized_.data();
            }
            
            return check;
        }
        
        static constexpr Int AmbDim()
        {
            return AMB_DIM;
        }
        
        static constexpr Int DomDimS()
        {
            return DOM_DIM1;
        }
        
        static constexpr Int DomDimT()
        {
            return DOM_DIM2;
        }
        
        static constexpr Int CoordDimS()
        {
            return COORD_DIMS;
        }
        
        static constexpr Int CoordDimT()
        {
            return COORD_DIMT;
        }
        
        static constexpr Int ProjectorDim()
        {
            return PROJ_DIM;
        }
        
        static constexpr Int NearDimS()
        {
            return NEAR_DIMS;
        }
        
        static constexpr Int NearDimT()
        {
            return NEAR_DIMT;
        }
        
        virtual void LoadS( const Int i )
        {
            S_ID = i;
            const Real * const restrict X = &S_near[NEAR_DIMS * S_ID];
            
            a = X[0];
            
#ifdef NearField_S_Copy
            copy_buffer( &X[1+COORD_DIMS], &p[0], PROJ_DIM );
#else
            p = &X[1+COORD_DIMS];
#endif
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                x[k] = S_scale * X[1 + AMB_DIM * 0 + k];
                
                for( Int kk = 1; kk < DOM_DIM1+1; ++kk )
                {
                    x[k] += S_scale * X[1 + AMB_DIM * kk + k];
                }
            }
        }
        
        virtual void LoadT( const Int j )
        {
            T_ID = j;
            const Real * const restrict Y = &T_near[NEAR_DIMT * T_ID];
            
            b = Y[0];
            
#ifdef NearField_T_Copy
            copy_buffer( &Y[1+COORD_DIMT], &q[0], PROJ_DIM );
#else
            q = &Y[1+COORD_DIMT];
#endif
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                y[k] = T_scale * Y[1 + AMB_DIM * 0 + k];

                for( Int kk = 1; kk < DOM_DIM2+1; ++kk )
                {
                    y[k] += T_scale * Y[1 + AMB_DIM * kk + k];
                }
            }
        }
        
        virtual void CleanseDBufferS()
        {
            zerofy_buffer( &DX[0], NEAR_DIMS );
        }
        
        virtual void WriteDBufferS()
        {
            Real * restrict const to = &S_D_near[NEAR_DIMS * S_ID];

            for( Int k = 0; k < NEAR_DIMS; ++k )
            {
                to[k] += DX[k];
            }
        }
        
        virtual void CleanseDBufferT()
        {
            zerofy_buffer( &DY[0], NEAR_DIMT );
        }
        
        virtual void WriteDBufferT()
        {
            Real * restrict const to = &T_D_near[NEAR_DIMT * T_ID];

            for( Int k = 0; k < NEAR_DIMT; ++k )
            {
                to[k] += DY[k];
            }
        }
        
        
    public:
        
        virtual std::string Stats() const
        {
            return this->ClassName();
        }
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef CLASS


