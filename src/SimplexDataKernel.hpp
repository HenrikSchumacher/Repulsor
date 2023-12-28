#pragma once

namespace Repulsor
{
    
    template<int DOM_DIM_, int AMB_DIM_, typename Real_, typename Int_>
    class SimplexDataKernel
    {
    public:
        
        using Real = Real_;
        using Int  = Int_;
        
        static constexpr Int DOM_DIM    = DOM_DIM_;
        static constexpr Int AMB_DIM    = AMB_DIM_;
        static constexpr Int SIZE       = DOM_DIM + 1;
        static constexpr Int PROJ_DIM   = (AMB_DIM*(AMB_DIM+1))/2;
        static constexpr Int  HULL_SIZE = AMB_DIM * SIZE;
        static constexpr Real nth       = Inv<Real>( SIZE );
        
    private:
        
        
        cref<Tensor2<Real,Int>> V_coords;
        cref<Tensor2<Int ,Int>> simplices;
        cref<Tensor1<Real,Int>> V_charges;
        
        Tiny::Vector<SIZE,Int,Int>   simplex;
        Tiny::Vector<SIZE,Int,Int> s_simplex;
        
        SortNet<SIZE,Int> sort;
        
        Tiny::Vector<AMB_DIM,Real,Int> center;
        
        Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
        Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
        Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
        Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
        
        Tiny::SelfAdjointMatrix<DOM_DIM,Real,Int> g_inv;
        
        Real a;
        Real charge;
        
    public:
        
//        SimplexDataKernel( const Mesh_T & M_ )
//        :   M ( M_ )
//        {
//            if constexpr ( DOM_DIM == 0 )
//            {
//                P.SetZero();
//                for( Int k = 0; k < AMB_DIM; ++k )
//                {
//                    P[k][k] = Scalar::One<Real>;
//                }
//            }
//        }
        
        SimplexDataKernel(
            cref<Tensor2<Real,Int>> V_coords_,
            cref<Tensor2<Int ,Int>> simplices_,
            cref<Tensor1<Real,Int>> V_charges_
        )
        :   V_coords ( V_coords_  )
        ,   simplices( simplices_ )
        ,   V_charges( V_charges_ )
        {
            if constexpr ( DOM_DIM == 0 )
            {
                P.SetZero();
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    P[k][k] = Scalar::One<Real>;
                }
            }
        }
        
        ~SimplexDataKernel() = default;
        
        void ReadSimplex( const Int i )
        {
            simplex.Read( simplices.data(i) );
            s_simplex.Read( simplex.data()      );
            
            // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
            sort( &s_simplex[0] );
            
            for( Int l = 0; l < SIZE; ++l )
            {
                copy_buffer<AMB_DIM>( V_coords.data(simplex[l]), hull[l] );
            }
            
            copy_buffer<AMB_DIM>(hull[0],center.data());
            
            charge = V_charges[simplex[0]];
            for( Int l = 1; l < SIZE; ++l )
            {
                add_to_buffer<AMB_DIM>( hull[l], center.data() );
                
                charge += V_charges[simplex[l]];
            }
            center *= nth;
            charge *= nth;
            
            a = charge * StandardSimplexVolume<Real>(DOM_DIM);
            
            if constexpr ( DOM_DIM > 0 )
            {
                for( Int l = 0; l < DOM_DIM; ++l )
                {
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        df[k][l] =
                        V_coords[s_simplex[l+1]][k]
                        -
                        V_coords[s_simplex[0  ]][k];
                    }
                }
                
                df.Transpose( dfdagger );
                
                // Compute inverse metric.
                
                // g = df^T * df.
                // At the moment dfdagger is just the transpose of df.
                Tiny::gemm<
                    Op::Id,Op::Id,DOM_DIM,DOM_DIM,AMB_DIM,
                    Scalar::Flag::Plus,Scalar::Flag::Zero
                >(
                    Scalar::One<Real>,  dfdagger.data(), AMB_DIM,
                                        df.data(),       DOM_DIM,
                    Scalar::Zero<Real>, g_inv.data(),    DOM_DIM
                );
                
                // Factorize g in place.
                g_inv.Cholesky();
                
                for( Int l = 0; l < DOM_DIM; ++l )
                {
                    a *= g_inv[l][l];
                }
                
                
                // Compute pseudo-inverse of df
                
                //  dfdagger = g^{-1} * df^T
                g_inv.CholeskySolve( dfdagger );
                
                
                // Compute normal projektor.
                
                // P = id - df * dfdagger
                Tiny::gemm<
                    Op::Id,Op::Id,AMB_DIM,AMB_DIM,DOM_DIM,
                    Scalar::Flag::Minus,Scalar::Flag::Zero
                >(
                    -Scalar::One<Real>,  df.data(),       DOM_DIM,
                                         dfdagger.data(), AMB_DIM,
                    Scalar::Zero<Real>,  P.data(),        AMB_DIM
                );
                
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    P[k][k] += Scalar::One<Real>;
                }
            }
        }
        
        void WriteNear( mptr<Real> near ) const
        {
            near[0] = a;
            
            hull.Write( &near[1] );
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < PROJ_DIM; ++k )
            {
                near[1 + AMB_DIM * SIZE + k] = P[ tri_i<AMB_DIM>(k) ][ tri_j<AMB_DIM>(k) ];
            }
        }
        
        void WriteFar( mptr<Real> far ) const
        {
            far[0] = a;
            
            center.Write( &far[1] );
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < PROJ_DIM; ++k )
            {
                far[1 + AMB_DIM + k] = P[ tri_i<AMB_DIM>(k) ][ tri_j<AMB_DIM>(k) ];
            }
        }
        
        void WriteSortedSimplex( mptr<Int> simplex_out ) const
        {
            s_simplex.Write(simplex_out);
        }
        
        
        
        void WriteDiffOp( mptr<Real> diff )
        {
            // Compute derivative operator (AMB_DIM x SIZE matrix).
            
            if constexpr ( DOM_DIM == 0 )
            {
                zerofy_buffer<DOM_DIM>(diff);
            }
            else
            {
                
                Tiny::Matrix<AMB_DIM,SIZE,Real,Int> D_f;
                
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    D_f[k][0] = 0;
                    
                    for( Int l = 0; l < DOM_DIM; ++l )
                    {
                        D_f[k][0  ] -= dfdagger[l][k];
                        D_f[k][l+1]  = dfdagger[l][k];
                    }
                }
                
                D_f.Write(diff);
            }
        }
        
        void WriteErrorQuadric( mptr<Real> Q_out ) const
        {
            Tiny::Matrix<AMB_DIM+1,AMB_DIM+1,Real,Int> Q;
            
            Tiny::Vector<AMB_DIM,Real,Int> b;
            
            Dot<Overwrite>( P, center, b );
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                Q[AMB_DIM][i] = -b[i];
                Q[i][AMB_DIM] = -b[i];
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    Q[i][j] = P[i][j];
                }
            }
            
            Q[AMB_DIM][AMB_DIM] = Dot(b,b);
            
                
            Q.Write( Q_out );
        }
        
        
        const Real Volume() const
        {
            return a;
        }
        
        const Tiny::Matrix<SIZE,AMB_DIM,Real,Int> & Hull() const
        {
            return hull;
        }
        
        const Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> & NormalProjector() const
        {
            return P;
        }
        
        const Tiny::Vector<AMB_DIM,Real,Int> & Center() const
        {
            return center;
        }
        
        const Tiny::Vector<AMB_DIM,Int,Int> & Simplex() const
        {
            return simplex;
        }
        
        
    }; // class SimplexDataKernel
    
}
