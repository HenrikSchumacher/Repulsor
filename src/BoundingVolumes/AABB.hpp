#pragma once

#define CLASS AABB
#define BASE  BoundingVolumeBase<AMB_DIM,Real,Int,SReal>

namespace Repulsor
{
    // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type SReal by calling the member SetPointer.
    
    // DATA LAYOUT
    // serialized_data[0] = squared radius
    // serialized_data[1],...,serialized_data[AMB_DIM] = center
    // serialized_data[AMB_DIM + 1],...,serialized_data[AMB_DIM + AMB_DIM] = half the edge lengths.
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        Real id_matrix[AMB_DIM][AMB_DIM];
    
        const int bits = static_cast<int>(64-1) / static_cast<int>(AMB_DIM);
        const SReal power  = static_cast<SReal>(0.9999) * std::pow( Scalar::Two<SReal>, bits );
        const SReal offset = Scalar::Half<SReal> * power;
        
    protected:
        
        void Initialize()
        {
            for( Int k1 = 0; k1 < AMB_DIM; ++k1 )
            {
                for( Int k2 = 0; k2 < AMB_DIM; ++k2 )
                {
                    id_matrix[k1][k2] = static_cast<SReal>(k1==k2);
                }
            }
        }
    
    public:
        
        CLASS() : BASE()
        {
            Initialize();
        }

        // Copy constructor
        CLASS( const CLASS & other ) : BASE( other )
        {
            Initialize();
        }
        
        // Move constructor
        CLASS( CLASS && other ) noexcept : BASE( other )
        {
            Initialize();
        }
        
        virtual ~CLASS() override = default;
        
        
        static constexpr Int SIZE = 1 + AMB_DIM + AMB_DIM;
        
        virtual constexpr Int Size() const override
        {
            return SIZE;
        }
        
//    protected:
//
//        mutable SReal self_buffer [1 + AMB_DIM + AMB_DIM];
 
    public:
        
#include "../Primitives/Primitive_Common.hpp"
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(CLASS)
        
        
    public:
        
        // array p is suppose to represent a matrix of size N x AMB_DIM
        void FromPointCloud( cptr<SReal> coords_in, const Int N ) const override
        {
//            tic(ClassName()+"::FromPointCloud");
            mref<SReal> r2 = serialized_data[0];
            
            // Abusing serialized_data temporily as working space.
            mptr<SReal> box_min = serialized_data + 1;
            mptr<SReal> box_max = serialized_data + 1 + AMB_DIM;
            
            cptr<SReal> coords = coords_in;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                box_min[k] = std::numeric_limits<SReal>::max();
                box_max[k] = std::numeric_limits<SReal>::lowest();
            }
            
            // No parallelization here because AABB objects are supposed to be created per thread anyways.
            // A desperate attempt to ask for simdization. The data layout of coords is wrong, so AVX will be useful only for AMB_DIM >=4. Nonetheless SSE instructs could be used for AMB_DIM = 2 and AMB_DIM = 3...
            
            
            for( Int i = 0; i < N; ++i )
            {
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    const SReal x = static_cast<SReal>(coords[ AMB_DIM * i + k ]);
                    box_min[k] = Min( box_min[k], x );
                    box_max[k] = Max( box_max[k], x );
                }
            }
            
            r2 = Scalar::Zero<SReal>;

            for( Int k = 0; k < AMB_DIM; ++k )
            {
                const SReal diff = Scalar::Half<SReal> * (box_max[k] - box_min[k]);
                r2 += diff * diff;
                
                // adding half the edge length to obtain the k-th coordinate of the center
                box_min[k] += diff;
                // storing half the edge length in the designated storage.
                box_max[k]  = diff;
            }
//            toc(ClassName()+"::FromPointCloud");
        }
        

        // array p is supposed to represent a matrix of size N x AMB_DIM
        void FromPrimitives(
            mref<PolytopeBase<AMB_DIM,Real,Int,SReal>> P,  // primitive prototype
            mptr<SReal> P_serialized,                      // serialized data of primitives
            const Int begin,                               // which _P_rimitives are in question
            const Int end,                                 // which _P_rimitives are in question
            Int thread_count = 1                           // how many threads to utilize
        ) const
        {
            using Vector_T = std::array<SReal,AMB_DIM>;
            
            Vector_T lower;
            Vector_T upper;
         
            lower.fill( Scalar::Max<SReal> );
            upper.fill( Scalar::Min<SReal> );
            
            if( thread_count <= 1 )
            {
                for( Int i = begin; i < end; ++i )
                {
                    P.SetPointer( P_serialized, i );
                    P.BoxMinMax( &lower[0], &upper[0] );
                }
            }
            else
            {
                std::vector<Vector_T> thread_lower ( thread_count );
                std::vector<Vector_T> thread_upper ( thread_count );
                
                ParallelDo(
                    [&]( const Int thread )
                    {
                        std::shared_ptr<PolytopeBase<AMB_DIM,Real,Int,SReal>> Q = P.Clone();
                        
                        Vector_T L;
                        Vector_T U;
                        
                        L.fill( Scalar::Max<SReal> );
                        U.fill( Scalar::Min<SReal> );
                        
                        const Int i_begin = begin + JobPointer( end - begin, thread_count, thread    );
                        const Int i_end   = begin + JobPointer( end - begin, thread_count, thread +1 );
                        
                        for( Int i = i_begin; i < i_end; ++i )
                        {
                            Q->SetPointer( P_serialized, i );
                            Q->BoxMinMax( &L[0], &U[0] );
                        }
                        
                        thread_lower[thread] = L;
                        thread_upper[thread] = U;
                        
                    },
                    thread_count
                );
                
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        lower[k] = Min( lower[k], thread_lower[thread][k] );
                        upper[k] = Max( upper[k], thread_upper[thread][k] );
                    }
                }
            }
            
            SReal r2 = Scalar::Zero<SReal>;

            for( Int k = 0; k < AMB_DIM; ++k )
            {
                const SReal diff = Scalar::Half<SReal> * (upper[k] - lower[k]);
                r2 += diff * diff;

                // adding half the edge length to obtain the k-th coordinate of the center
                lower[k] += diff;
                // storing half the edge length in the designated storage.
                upper[k]  = diff;
            }
            
            serialized_data[0] = r2;
            copy_buffer<AMB_DIM>( &lower[0], &serialized_data[1          ] );
            copy_buffer<AMB_DIM>( &upper[0], &serialized_data[1 + AMB_DIM] );
        }
        
        // array p is supposed to represent a matrix of size N x AMB_DIM
        virtual void FromPrimitives(
            mref<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> P,    // primitive prototype
            mptr<SReal> P_serialized,                               // serialized data of primitives
            const Int begin,                                        // which _P_rimitives are in question
            const Int end,                                          // which _P_rimitives are in question
            Int thread_count = 1                                    // how many threads to utilize
        ) const override
        {
            using Vector_T = std::array<SReal,AMB_DIM>;
            
            Vector_T lower;
            Vector_T upper;

            lower.fill( Scalar::Max<SReal> );
            upper.fill( Scalar::Min<SReal> );
            
            if( thread_count <= 1 )
            {
                for( Int i = begin; i < end; ++i )
                {
                    P.SetPointer( P_serialized, i );

                    for( Int j = 0; j < AMB_DIM; ++j )
                    {
                        Real min_val;
                        Real max_val;

                        P.MinMaxSupportValue( &id_matrix[j][0], min_val, max_val );

                        lower[j] = Min( lower[j], static_cast<SReal>(min_val) );
                        upper[j] = Max( upper[j], static_cast<SReal>(max_val) );
                    }
                }
            }
            else
            {
                std::vector<Vector_T> thread_lower ( thread_count );
                std::vector<Vector_T> thread_upper ( thread_count );

                ParallelDo(
                    [&]( const Int thread )
                    {
                        const Int i_begin = begin + JobPointer( end - begin, thread_count, thread    );
                        const Int i_end   = begin + JobPointer( end - begin, thread_count, thread +1 );

//                        print( ToString(thread) + " -> { " + ToString(i_begin) + ", " + ToString(i_end) + " }" );

                        std::shared_ptr<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> Q = P.Clone();

                        Vector_T L;
                        Vector_T U;
                        
                        L.fill( Scalar::Max<SReal> );
                        U.fill( Scalar::Min<SReal> );

                        for( Int i = i_begin; i < i_end; ++i )
                        {
                            Q->SetPointer( P_serialized, i );

                            for( Int j = 0; j < AMB_DIM; ++j )
                            {
                                Real min_val;
                                Real max_val;

                                Q->MinMaxSupportValue( &id_matrix[j][0], min_val, max_val );

                                L[j] = Min( L[j], static_cast<SReal>(min_val) );
                                U[j] = Max( U[j], static_cast<SReal>(max_val) );
                            }
                        }

                        thread_lower[thread] = L;
                        thread_upper[thread] = U;
                    },
                    thread_count
                );

                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        lower[k] = Min( lower[k], thread_lower[thread][k] );
                        upper[k] = Max( upper[k], thread_upper[thread][k] );
                    }
                }
            }

            SReal r2 = Scalar::Zero<SReal>;

            for( Int k = 0; k < AMB_DIM; ++k )
            {
                const SReal diff = Scalar::Half<SReal> * (upper[k] - lower[k]);
                r2 += diff * diff;
                
                // adding half the edge length to obtain the k-th coordinate of the center
                lower[k] += diff;
                // storing half the edge length in the designated storage.
                upper[k]  = diff;
            }
            
            serialized_data[0] = r2;
            copy_buffer<AMB_DIM>( &lower[0], &serialized_data[1          ] );
            copy_buffer<AMB_DIM>( &upper[0], &serialized_data[1 + AMB_DIM] );
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> v, mptr<Real> s ) const override
        {
            cptr<SReal> x = serialized_data + 1;
            cptr<SReal> L = serialized_data + 1 + AMB_DIM;
            
            Real R2 = Scalar::Zero<Real>;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                const Real x_k = static_cast<Real>(x[k]);
                const Real L_k = static_cast<Real>(L[k]);
                
                const Real R1  = v[k] * L_k;
                
                R2  += v[k] * x_k + Abs(R1);
                s[k] = x_k + Sign(R1) * L_k;
            }

            return R2;
        }


        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> v, mptr<Real> s ) const override
        {
            cptr<SReal> x = serialized_data + 1;
            cptr<SReal> L = serialized_data + 1 + AMB_DIM;
            
            Real R2 = Scalar::Zero<Real>;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                const Real x_k = static_cast<Real>(x[k]);
                const Real L_k = static_cast<Real>(L[k]);
                
                const Real R1  = v[k] * L_k;
                
                R2  += v[k] * x_k - Abs(R1);
                s[k] = x_k - Sign(R1) * L_k;
            }

            return R2;
        }
                
        
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            // Could be implemented more efficiently, but since this routine is unlikely to be used...
            min_val = MinSupportVector( dir, &this->Real_buffer[0] );
            max_val = MaxSupportVector( dir, &this->Real_buffer[AMB_DIM] );
        }
        
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
        
        inline friend Real AABB_SquaredDistance( cref<CLASS> P, cref<CLASS> Q )
        {
            cptr<SReal> P_x = P.serialized_data+1;              // center of box P
            cptr<SReal> P_L = P.serialized_data+1+AMB_DIM;      // edge half-lengths of box P
            
            cptr<SReal> Q_x = Q.serialized_data+1;              // center of box Q
            cptr<SReal> Q_L = Q.serialized_data+1+AMB_DIM;      // edge half-lengths of box Q
            
            Real d2 = Scalar::Zero<Real>;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                Real x = static_cast<Real>(
                    Max(
                        Scalar::Zero<SReal>,
                        Max( P_x[k]-P_L[k], Q_x[k]-Q_L[k] )
                        -
                        Min( P_x[k]+P_L[k], Q_x[k]+Q_L[k] )
                    )
                );
                d2 += x * x;
            }
            return d2;
        }
        
        void Merge( mptr<SReal> C_Serialized, const Int i = 0 ) const
        {
            mptr<SReal> p = C_Serialized + SIZE * i;
            
            if( serialized_data != p )
            {
                
                mptr<SReal> x1 = serialized_data + 1;
                mptr<SReal> L1 = serialized_data + 1 + AMB_DIM;
                
                mptr<SReal> x2 = p + 1;
                mptr<SReal> L2 = p + 1 + AMB_DIM;
                
                SReal r2 = Scalar::Zero<SReal>;
                
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    const SReal box_min = Min( x1[k] - L1[k], x2[k] - L2[k] );
                    const SReal box_max = Max( x1[k] + L1[k], x2[k] + L2[k] );
                    
                    x1[k] = Scalar::Half<SReal> * ( box_max + box_min );
                    L1[k] = Scalar::Half<SReal> * ( box_max - box_min );
                    
                    r2 += L1[k] * L1[k];
                }
                serialized_data[0] = r2;
            }
        }
        
    }; // CLASS

} // namespace Repulsor

#undef CLASS
#undef BASE
    


