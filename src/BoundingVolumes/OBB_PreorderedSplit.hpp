#pragma once

#define BASE  OBB<AMB_DIM,Real,Int,SReal>
#define CLASS OBB_PreorderedSplit

namespace Repulsor
{
    
    // The OBB is the image of Cuboid[ {-L[0],...,-L[AMB_DIM-1]}, {L[0],...,L[AMB_DIM-1]} ] under the ORTHOGONAL mapping x \mapsto rotation * x + center.
    
    
    // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
    
    // DATA LAYOUT
    // serialized_data[0] = squared radius
    // serialized_data[1],...,serialized_data[AMB_DIM] = center
    // serialized_data[1+AMB_DIM],...,serialized_data[AMB_DIM+AMB_DIM] = L (vector of half the edge lengths)
    // serialized_data[1+AMB_DIM+AMB_DIM],...,serialized_data[AMB_DIM + AMB_DIM + AMB_DIM x AMB_DIM] = rotation^T. BEWARE THE TRANSPOSITION!!!!!!!!!!!!!
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    public:
        
        CLASS() : BASE() {}

        // Copy constructor
        CLASS( const CLASS & other ) : BASE( other ) {}
        
        // Move constructor
        CLASS( CLASS && other ) noexcept : BASE( other ) {}
        
        __ADD_CLONE_CODE__(CLASS)
        
        virtual ~CLASS() override = default;
        
        static constexpr Int SIZE = 1 + AMB_DIM + AMB_DIM + AMB_DIM * AMB_DIM;
        
        virtual constexpr Int Size() const override
        {
            return SIZE;
        }
        
//protected:
//        
//        using BASE::self_buffer;
 
public:

#include "../Primitives/Primitive_Common.hpp"
        
        
        // Array coords_in is supposes to represent a matrix of size n x AMB_DIM
        void FromPointCloud( const Real * const coords_in, const Int n ) const override
        {
            if( n == 0 )
            {
                eprint(ClassName()+"FromPointCloud : 0 members in point could.");
                return;
            }
            
            // Zero bounding volume's data.
            zerofy_buffer<SIZE>(serialized_data);
            
            // Abusing serialized_data temporily as working space.
            mref<SReal> r2 = serialized_data[0];
            mptr<SReal> average    = serialized_data + 1;
            mptr<SReal> covariance = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            // Compute average of the points
            for( Int i = 0; i < n; ++i )
            {
                cptr<Real> p = coords_in + AMB_DIM * i;

                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    average[k] += p[k];
                }
            }

            const Real n_inv = Inv<Real>(n);
            
            scale_buffer<AMB_DIM>(n_inv,average);
            
            // Compute covariance matrix.
            for( Int i = 0; i < n; ++i )
            {
                cptr<Real> p = coords_in + AMB_DIM * i;

                for( Int k1 = 0; k1 < AMB_DIM; ++k1 )
                {
                    const Real delta = (p[k1] - average[k1]);

                    for( Int k2 = k1; k2 < AMB_DIM; ++k2 )
                    {
                        covariance[AMB_DIM * k1 + k2] += delta * (p[k2] - average[k2]);
                    }
                }
            }

            for( Int k1 = 0; k1 < AMB_DIM; ++k1 )
            {
                for( Int k2 = k1; k2 < AMB_DIM; ++k2 )
                {
                    covariance[AMB_DIM * k1 + k2] *= n_inv;
                }
            }
            
            // Abusing serialized_data temporily as working space.
            mptr<SReal> box_min = serialized_data + 1;
            mptr<SReal> box_max = serialized_data + 1 + AMB_DIM;
            
            (void)SymmetricEigenSolve( covariance, box_min );
            
            // Now "covariance" stores the eigenbasis.
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                box_min[k] = std::numeric_limits<SReal>::max();
                box_max[k] = std::numeric_limits<SReal>::lowest();
            }

            // TODO: Check whether reversing this loop would be helpful.
            // TODO: Can this loop be parallelized?
            for( Int i = 0; i < n; ++i )
            {
                cptr<Real> p = coords_in + AMB_DIM * i;

                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    cptr<Real> vec = covariance + AMB_DIM * j;

                    Real x = p[0] * vec[0];

                    for( Int k = 1; k < AMB_DIM; ++k )
                    {
                        x += p[k] * vec[k];
                    }

                    box_min[j] = Min( box_min[j], x );
                    box_max[j] = Max( box_max[j], x );
                }
            }

            for( Int k = 0; k < AMB_DIM; ++k )
            {
                Real diff = static_cast<Real>( Scalar::Half<SReal> * (box_max[k] - box_min[k]) );
                r2 += diff * diff;

                // adding half the edge length to obtain the k-th coordinate of the center
                this->buffer[k] = box_min[k] + diff;
        
                // storing half the edge length in the designated storage.
                box_max[k]  = diff;
            }
            
            mptr<Real> center    = serialized_data + 1;
            mptr<Real> rotationT = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            for( Int j = 0; j < AMB_DIM; ++j )
            {
                center[j] = rotationT[j] * this->buffer[0];

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    center[j] += this->buffer[k] * rotationT[AMB_DIM * k + j];
                }
            }
        } // FromPointCloud

        
        // array p is supposed to represent a matrix of size N x AMB_DIM
        virtual void FromPrimitives(
            PrimitiveSerialized<AMB_DIM,Real,Int,SReal> & P,      // primitive prototype
            mptr<SReal> P_serialized,                  // serialized data of primitives
            const Int begin,                           // which _P_rimitives are in question
            const Int end,                             // which _P_rimitives are in question
            Int thread_count = 1                       // how many threads to utilize
        ) const override
        {
//            ptic(ClassName()+"::FromPrimitives (PrimitiveSerialized)");
            if( begin >= end )
            {
                eprint(ClassName()+"::FromPrimitives : begin = "+ToString(begin)+" >= "+ToString(end)+" = end");
                return;
            }
            
            // Zero bounding volume's data.
            zerofy_buffer<SIZE>(serialized_data);
            
            // Abusing serialized_data temporily as working space.
            mref<SReal> r2 = serialized_data[0];
            mptr<SReal> average    = serialized_data + 1;
            mptr<SReal> covariance = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            const Int P_Size = P.Size();
            
            // Compute average of the InterPoints of all primitives.
            for( Int i = begin; i < end; ++i )
            {
                cptr<SReal> p = P_serialized + 1 + P_Size * i;

                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    average[k] += p[k];
                }
            }

            const SReal n_inv = Inv<Real>(end-begin);

            scale_buffer<AMB_DIM>(n_inv,average);
            
            // Compute covariance matrix.
            for( Int i = begin; i < end; ++i )
            {
                cptr<SReal> p = P_serialized + 1 + P_Size * i;

                for( Int k1 = 0; k1 < AMB_DIM; ++k1 )
                {
                    const SReal delta = (p[k1] - average[k1]);

                    for( Int k2 = k1; k2 < AMB_DIM; ++k2 )
                    {
                        covariance[AMB_DIM * k1 + k2] += delta * (p[k2] - average[k2]);
                    }
                }
            }

            for( Int k1 = 0; k1 < AMB_DIM; ++k1 )
            {
                for( Int k2 = k1; k2 < AMB_DIM; ++k2 )
                {
                    covariance[AMB_DIM * k1 + k2] *= n_inv;
                }
            }
            
            // Abusing serialized_data temporily as working space.
            mptr<SReal> box_min = serialized_data + 1;
            mptr<SReal> box_max = serialized_data + 1 + AMB_DIM;
        
            (void)SymmetricEigenSolve( covariance, box_min );
            
            // Now "covariance" stores the eigenbasis.
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                box_min[k] = std::numeric_limits<SReal>::max();
                box_max[k] = std::numeric_limits<SReal>::lowest();
            }
            
            // TODO: Check whether reversing this loop would be helpful.
            // TODO: Can this loop be parallelized?
            for( Int i = begin; i < end; ++i )
            {
                P.SetPointer( P_serialized, i );
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    SReal min_val;
                    SReal max_val;
                    P.MinMaxSupportValue( covariance + AMB_DIM * j, min_val, max_val );
                    box_min[j] = Min( box_min[j], min_val );
                    box_max[j] = Max( box_max[j], max_val );
                }
            }
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                Real diff = static_cast<Real>( Scalar::Half<SReal> * (box_max[k] - box_min[k]) );
                r2 += diff * diff;

                // adding half the edge length to obtain the k-th coordinate of the center (within the transformed coordinates)
                this->buffer[k] = box_min[k] + diff;
        
                // storing half the edge length in the designated storage.
                box_max[k]  = diff;
            }
            
            // Rotate this->buffer so that it falls onto the true center of the bounding box.
            
            mptr<SReal> center    = serialized_data + 1;
            mptr<SReal> rotationT = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            for( Int j = 0; j < AMB_DIM; ++j )
            {
                center[j] = rotationT[j] * this->buffer[0];

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    center[j] += this->buffer[k] * rotationT[AMB_DIM * k + j];
                }
            }
//            ptoc(ClassName()+"::FromPrimitives (PrimitiveSerialized)");
        }
        
        virtual Int Split(
            PrimitiveSerialized<AMB_DIM,Real,Int,SReal> & P,                          // primitive prototype; to be "mapped" over P_serialized, thus not const.
            SReal * const P_serialized, const Int begin, const Int end,    // which _P_rimitives are in question
            Int   * const P_ordering,                                        // to keep track of the permutation of the primitives
            SReal * const C_serialized, const Int C_ID,                     // where to get   the bounding volume info for _C_urrent bounding volume
            SReal * const L_serialized, const Int L_ID,                     // where to store the bounding volume info for _L_eft  child (if successful!)
            SReal * const R_serialized, const Int R_ID,                     // where to store the bounding volume info for _R_ight child (if successful!)
            SReal *       score,                                             // some scratch buffer for one scalar per primitive
            Int   *       perm,                                              // some scratch buffer for one Int per primitive (for storing local permutation)
            Int   *       inv_perm,                                          // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count = 1                                           // how many threads to utilize
        ) override
        {
//            ptic(ClassName()+"::Split");
            Int split_index = begin + ((end-begin)/2);
            
            // Prevent further computations with fewer than two primitives.
            if( ( begin+1 < split_index ) && ( split_index < end-1 ) )
            {
                // Compute bounding volume of left child.
                SetPointer( L_serialized, L_ID );
                FromPrimitives( P, P_serialized,   begin, split_index,   thread_count );
                // Compute bounding volume of right child.
                SetPointer( R_serialized, R_ID );
                FromPrimitives( P, P_serialized,   split_index, end,   thread_count );
            }
            else
            {
                split_index = -1;
            }
            // ... otherwise we assume that the bounding volume hierarchy / cluster tree won't do the split.
//            ptoc(ClassName()+"::Split");
            
            return split_index;
        
        } // Split
        
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<SReal> supp ) const override
        {
            cptr<SReal> x = serialized_data + 1;
            cptr<SReal> L = serialized_data + 1 + AMB_DIM;
            cptr<SReal> A = serialized_data + 1 + AMB_DIM + AMB_DIM;
            cptr<Real> v = dir;
            mptr<Real> s = supp;

            Real R1;
            Real R2;
            Real R3 = Scalar::Zero<Real>;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                s[i] = static_cast<Real>(x[i]);
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                // Multiply v with i-th row.
                R1 = static_cast<Real>(A[ AMB_DIM * i]) * v[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += v[j] * static_cast<Real>(A[ AMB_DIM * i + j ]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? static_cast<Real>(L[i]) : -static_cast<Real>(L[i]);
                
                R3 += static_cast<Real>(x[i]) * v[i] + R1 * R2;
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    s[j] +=  static_cast<Real>(A[ AMB_DIM * i + j ]) * R2;
                }
            }

            return R3;
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            cptr<SReal> x = serialized_data + 1;
            cptr<SReal> L = serialized_data + 1 + AMB_DIM;
            cptr<SReal> A = serialized_data + 1 + AMB_DIM + AMB_DIM;
            cptr< Real> v = dir;
            mptr< Real> s = supp;

            Real R1;
            Real R2;
            Real R3 = Scalar::Zero<Real>;

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                s[i] = static_cast<Real>(x[i]);
            }

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                // Multiply v with i-th row.
                R1 = static_cast<Real>(A[AMB_DIM *i]) * v[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += v[j] * static_cast<Real>(A[ AMB_DIM * i + j ]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? -static_cast<Real>(L[i]) : static_cast<Real>(L[i]);
                
                R3 += static_cast<Real>(x[i]) * v[i] + R1 * R2;
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    s[j] += static_cast<Real>(A[ AMB_DIM * i + j ]) * R2;
                }
            }

            return R3;
        }
                
        virtual void MinMaxSupportValue( const Real * const dir, Real & min_val, Real & max_val ) const override
        {
            min_val = MinSupportVector( dir, &this->Real_buffer[0] );
            max_val = MaxSupportVector( dir, &this->Real_buffer[AMB_DIM] );
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    }; // CLASS
    
} // namespace Repulsor


#undef CLASS
#undef BASE
