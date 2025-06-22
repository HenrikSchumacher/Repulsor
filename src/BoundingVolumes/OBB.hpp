#pragma once

namespace Repulsor
{
    // The OBB is the image of Cuboid[ {-L[0],...,-L[AMB_DIM-1]}, {L[0],...,L[AMB_DIM-1]} ] under the ORTHOGONAL mapping x \mapsto rotation * x + center.
    
    
    // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
    
    // DATA LAYOUT
    // serialized_data[0] = squared radius
    // serialized_data[1],...,serialized_data[AMB_DIM] = center
    // serialized_data[1+AMB_DIM],...,serialized_data[AMB_DIM+AMB_DIM] = L (vector of half the edge lengths)
    // serialized_data[1+AMB_DIM+AMB_DIM],...,serialized_data[AMB_DIM + AMB_DIM + AMB_DIM x AMB_DIM] = rotation^T. BEWARE THE TRANSPOSITION!!!
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class OBB : public BoundingVolumeBase<AMB_DIM,Real,Int,SReal>
    {
    public:
        
        using Base_T = BoundingVolumeBase<AMB_DIM,Real,Int,SReal>;
        
        // Default constructor
        OBB()
        : Base_T()
        {}
        // Destructor
        virtual ~OBB() override = default;
        // Copy constructor
        OBB( const OBB & other ) = default;
        // Copy assignment operator
        OBB & operator=( const OBB & other ) = default;
        // Move constructor
        OBB( OBB && other ) = default;
        // Move assignment operator
        OBB & operator=( OBB && other ) = default;
        
        
        static constexpr Int SIZE = 1 + AMB_DIM + AMB_DIM + AMB_DIM * AMB_DIM;
        
        virtual constexpr Int Size() const override
        {
            return SIZE;
        }
        
    public:
        
        [[nodiscard]] std::shared_ptr<OBB> Clone () const
        {
            return std::shared_ptr<OBB>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual OBB * CloneImplementation() const override = 0;
        
//protected:
//
//        mutable Real self_buffer [SIZE];
 
public:

#include "../Primitives/Primitive_Common.hpp"
        
        
        // Array coords_in is supposes to represent a matrix of size n x AMB_DIM
        void FromPointCloud( cptr<SReal> coords_in, const Int n ) const override
        {
            if( n == 0 )
            {
                eprint(ClassName()+"::FromPointCloud : 0 members in point could.");
                return;
            }
            
            // Zero bounding volume's data.
            zerofy_buffer<SIZE>(serialized_data);
            
            // Abusing serialized_data temporily as working space.
            mref<SReal> r2         = serialized_data[0];
            mptr<SReal> average    = serialized_data + 1;
            mptr<SReal> covariance = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            // Compute average of the points
            for( Int i = 0; i < n; ++i )
            {
                cptr<SReal> p = coords_in + AMB_DIM * i;

                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    average[k] += p[k];
                }
            }

            const SReal n_inv = Inv<SReal>(n);
            
            scale_buffer<AMB_DIM>(n_inv,average);
            
            // Compute covariance matrix.
            for( Int i = 0; i < n; ++i )
            {
                cptr<SReal> p = coords_in + AMB_DIM * i;

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
            
            int info = SymmetricEigenSolve( covariance, box_min );
            
            if( info )
            {
                print(ClassName()+"::Split:  SymmetricEigenSolve returned info = "+ToString(info)+".");
            }
            
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
                cptr<SReal> p = coords_in + AMB_DIM * i;

                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    cptr<SReal> vec = covariance + AMB_DIM * j;

                    SReal x = p[0] * vec[0];

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
                SReal diff = Scalar::Half<SReal> * (box_max[k] - box_min[k]);
                r2 += diff * diff;

                // adding half the edge length to obtain the k-th coordinate of the center
                this->buffer[k] = box_min[k] + diff;
        
                // storing half the edge length in the designated storage.
                box_max[k]  = diff;
            }
            
            mptr<SReal> center    = serialized_data + 1;
            mptr<SReal> rotationT = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            for( Int j = 0; j < AMB_DIM; ++j )
            {
                center[j] = rotationT[j] * this->SReal_buffer[0];

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    center[j] += this->SReal_buffer[k] * rotationT[AMB_DIM * k + j];
                }
            }
        } // FromPointCloud

        
        // array p is supposed to represent a matrix of size N x AMB_DIM
        virtual void FromPrimitives(
            mref<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> P,      // primitive prototype
            mptr<SReal> P_serialized,                  // serialized data of primitives
            const Int begin,                           // which _P_rimitives are in question
            const Int end,                             // which _P_rimitives are in question
            Int thread_count = 1                       // how many threads to utilize
        ) const override
        {
//            TOOLS_PTIC(ClassName()+"::FromPrimitives (PrimitiveSerialized)");
            if( begin >= end )
            {
                eprint(ClassName()+"::FromPrimitives : begin = "+ToString(begin)+" >= "+ToString(end)+" = end");
                return;
            }
            
            
            
            // Zero bounding volume's data.
            zerofy_buffer<SIZE>(serialized_data);
            
            // Abusing serialized_data temporily as working space.
            mref<SReal> r2         = serialized_data[0];
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

            const SReal n_inv = Inv<SReal>(end-begin);
            
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
                SReal diff = Scalar::Half<SReal> * (box_max[k] - box_min[k]);
                r2 += diff * diff;

                // adding half the edge length to obtain the k-th coordinate of the center (within the transformed coordinates)
                this->SReal_buffer[k] = box_min[k] + diff;
        
                // storing half the edge length in the designated storage.
                box_max[k]  = diff;
            }
            
            // Rotate this->SReal_buffer so that it falls onto the true center of the bounding box.
            
            mptr<SReal> center    = serialized_data + 1;
            mptr<SReal> rotationT = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            for( Int j = 0; j < AMB_DIM; ++j )
            {
                center[j] = rotationT[j] * this->SReal_buffer[0];

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    center[j] += this->SReal_buffer[k] * rotationT[AMB_DIM * k + j];
                }
            }
//            TOOLS_PTOC(ClassName()+"::FromPrimitives (PrimitiveSerialized)");
        }
        
        virtual Int Split(
            mref<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> P,                          // primitive prototype; to be "mapped" over P_serialized, thus not const.
            mptr<SReal> P_serialized, const Int begin, const Int end,    // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                        // to keep track of the permutation of the primitives
            mptr<SReal> C_serialized, const Int C_ID,                     // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_serialized, const Int L_ID,                     // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_serialized, const Int R_ID,                     // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                             // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                              // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                          // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count = 1                                           // how many threads to utilize
        ) override
        {
//            TOOLS_PTIC(ClassName()+"::Split");
//            valprint("begin",begin);
//            valprint("end",end);
            
            const Int P_Size = P.Size();
            
            SetPointer( C_serialized, C_ID );
            
            cptr<SReal> center    = serialized_data + 1;
            cptr<SReal> L         = serialized_data + 1 + AMB_DIM;
            cptr<SReal> rotationT = serialized_data + 1 + AMB_DIM + AMB_DIM;
            
            // Find the longest axis `split_dir` of primitives's bounding box.
            Int split_dir = 0;
            
            SReal L_max = L[0];
            
            for( Int k = 1; k < AMB_DIM; ++k )
            {
                if( L[k] > L_max )
                {
                    L_max = L[k];
                    split_dir = k;
                }
            }
            
            if( L_max <= Scalar::Zero<SReal> )
            {
                eprint(ClassName()+"::Split: longest axis has length <=0.");
                return -1;
            }
            
            // TODO: It can happen then the InteriorPoints of all primitives lie on one side of mid because mid is the center of the bounding box of the primitives --- and not the center of the bounding box of the InteriorPoints!!

            // Swap primitives according to their interior points in direction split_dir;
            Int split_index = begin;
            for( Int i = begin; i < end; ++i )
            {
                cptr<SReal> p = P_serialized + 1 + P_Size * i;

                SReal x = rotationT[ AMB_DIM * split_dir ] * (p[0] - center[0]);

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    x += rotationT[ AMB_DIM * split_dir + k ] * (p[k] - center[k]);
                }

                if( x < Scalar::Zero<SReal> )
                {
                    std::swap( P_ordering[split_index], P_ordering[i] );

                    P.Swap( P_serialized, split_index, P_serialized, i );

                    ++split_index;
                }
            }
            
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
//            TOOLS_PTOC(ClassName()+"::Split");
            
            return split_index;
            
        } // Split
        
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<SReal> dir, mptr<SReal> supp ) const override
        {
            cptr<SReal> x = serialized_data + 1;
            cptr<SReal> L = serialized_data + 1 + AMB_DIM;
            cptr<SReal> A = serialized_data + 1 + AMB_DIM + AMB_DIM;

            Real R1;
            Real R2;
            Real R3 = Scalar::Zero<Real>;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                supp[i] = x[i];
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                // Multiply v with i-th row.
                R1 = static_cast<Real>(A[ AMB_DIM * i]) * dir[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += dir[j] * static_cast<Real>(A[ AMB_DIM * i + j ]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? static_cast<Real>(L[i]) : -static_cast<Real>(L[i]);
                
                R3 += static_cast<Real>(x[i]) * dir[i] + R1 * R2;
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    supp[j] +=  static_cast<Real>(A[ AMB_DIM * i + j ]) * R2;
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

            Real R1;
            Real R2;
            Real R3 = Scalar::Zero<Real>;
            
            copy_buffer<AMB_DIM>( x, supp );
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                // Multiply v with i-th row.
                R1 = static_cast<Real>(A[AMB_DIM *i]) * dir[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += dir[j] * static_cast<Real>(A[ AMB_DIM * i + j ]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? -static_cast<Real>(L[i]) : static_cast<Real>(L[i]);
                
                R3 += static_cast<Real>(x[i]) * dir[i] + R1 * R2;

                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    supp[j] += static_cast<Real>(A[ AMB_DIM * i + j ]) * R2;
                }
            }

            return R3;
        }
                
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            min_val = MinSupportVector( dir, &this->Real_buffer[0] );
            max_val = MaxSupportVector( dir, &this->Real_buffer[AMB_DIM] );
        }
        
        virtual std::string ClassName() const override
        {
            return std::string("OBB")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+">";
        }
    
    protected:
        
        static lapack_int SymmetricEigenSolve( double * const matrix, double * const eigenvalues )
        {
            return LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'L', AMB_DIM, matrix, AMB_DIM, eigenvalues );
        }
        
        static lapack_int SymmetricEigenSolve(  float * const matrix,  float * const eigenvalues )
        {
            return LAPACKE_ssyev( LAPACK_COL_MAJOR, 'V', 'L', AMB_DIM, matrix, AMB_DIM, eigenvalues );
        }
        
    }; // OBB
    
} // namespace Repulsor
