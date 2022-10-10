#pragma once

#define CLASS SparseKernelMatrixCSR

namespace Repulsor
{
    template<typename Kernel_T>
    class CLASS : public BASE
    {
    public:
        
        using T   = typename Kernel_T::T;
        using Int = typename Kernel_T::Int;
        
        using SparsityPattern_T = SparseCSR<Int>;
        
        CLASS() = default;
        
        CLASS(
            const SparsityPattern_T & pattern_
        )
        :   pattern ( pattern_ )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   pattern ( pattern_ )
        {}

        virtual ~CLASS() override = default;
        
    protected:
        
        const SparsityPattern_T & pattern;
        Kernel_T                  kernel;
        const Tensor1<T,Int>    & values;
        
    public:
        
        
        const Tensor1<T,Int> & NonzeroValues() const
        {
            return values;
        }
        
        Int RowCount() const
        {
            return pattern.RowCount() * kernel.RowCount();
        }
        
        Int ColCount() const
        {
            return pattern.ColCount() * kernel.ColCount();
        }
        
        Int NonzeroCount() const
        {
            return pattern.NonzeroCount() * kernel.NonzeroCount();
        }

    public:
        
//##############################################################################################
//      Symmetrization
//##############################################################################################
        
        void FillLowerTriangleFromUpperTriangle() const
        {
            ptic(ClassName()+"::FillLowerTriangleFromUpperTriangle");
            
            if( pattern.WellFormed() && (pattern.RowCount()>= pattern.ColCount()) )
            {
                const Int * restrict const diag   = pattern.Diag().data();
                const Int * restrict const outer  = pattern.Outer().data();
                const Int * restrict const inner  = pattern.Inner().data();
                
                const auto & job_ptr = pattern.LowerTriangularJobPtr();
                
                const Int thread_count = job_ptr.Size()-1;
                
                #pragma omp parallel for num_threads( thread_count )
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    auto Kernel_T ker ( values.data() );
                    
                    const Int i_begin = job_ptr[thread];
                    const Int i_end   = job_ptr[thread+1];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        const Int k_begin = outer[i];
                        const Int k_end   =  diag[i];
                        
                        for( Int k = k_begin; k < k_end; ++k )
                        {
                            const Int j = inner[k];
                            
                            Int L =  diag[j];
                            Int R = outer[j+1]-1;
                            
                            while( L < R )
                            {
                                const Int M   = R - (R-L)/static_cast<Int>(2);
                                const Int col = inner[M];

                                if( col > i )
                                {
                                    R = M-1;
                                }
                                else
                                {
                                    L = M;
                                }
                            }
                            
                            ker.TransposeBlock(L, k);
                            
                        } // for( Int k = k_begin; k < k_end; ++k )

                    } // for( Int i = i_begin; i < i_end; ++i )
                    
                } // #pragma omp parallel
            }
            
            ptoc(ClassName()+"::FillLowerTriangleFromUpperTriangle");
        }
        
//##############################################################################################
//      Matrix multiplication
//##############################################################################################

        
    virtual void Dot(
        const T                      alpha,
        const T_in  * restrict const X,
        const T_out                  beta,
              T_out * restrict const Y
    ) const
    {
        ptoc(ClassName()+"::Dot" );
        
        const Int  * restrict const rp = pattern.Outer().data();
        const Int  * restrict const ci = pattern.Inner().data();
        
        const auto & job_ptr = pattern.JobPtr();
        
        const Int thread_count = job_ptr.Size()-1;
        
        #pragma omp parallel for num_threads( thread_count )
        for( Int thread = 0; thread < thread_count; ++thread )
        {
            // Initialize local kernel and feed it all information that is going to be constant along its life time.
            Kernel_T ker ( values, alpha, X, beta, Y );
            
            const block_rows = ker.RowCount();
            const block_cols = ker.ColCount();
            const nnz        = ker.NonzeroCount();
            
            // Kernel is supposed the following rows of pattern:
            const Int i_begin = job_ptr[thread  ];
            const Int i_end   = job_ptr[thread+1];
            
            for( Int i = i_begin; i < i_end; ++i )
            {
                // These are the corresponding nonzero blocks in i-th row.
                const Int k_begin = rp[i  ];
                const Int k_end   = rp[i+1];
                
                if( k_end > k_begin )
                {
                    // Clear the local vector chunk of the kernel.
                    ker.CleanseVector();
                    
                    // Perform all but the last calculation in row with prefetch.
                    for( int k = k_begin; k < k_end-1; ++k )
                    {
                        const Int j = ci[k];
                        
                        __builtin_prefetch( &X[cols * ci[k+1]] );
                        
                        __builtin_prefetch( &a[ nnz * (k+1) ] );
                        
                        // Let the kernel apply to the k-th block to the j-th chunk of the input.
                        // The result is stored in the kernel's local vector chunk X.
                        ker.ApplyBlock( k, j );
                    }
                    
                    // Perform last calculation in row without prefetch.
                    {
                        const Int k = k_end-1;
                        
                        const Int j = ci[k];
                        
                        // Should we include a guard here against prefecting from forbidden space?
                        __builtin_prefetch( &ci[k_end] );
                        
                        // Let the kernel apply to the k-th block to the j-th chunk of the input X.
                        // The result is stored in the kernel's local vector chunk.
                        ker.ApplyBlock( k, j );
                    }
                    
                    // Incorporate the kernel's local vector chunk into the i-th chunk if the output Y.
                    ker.WriteVector(i);
                }
                else
                {
                    // Just zerogy the i-th chunk if the output Y.
                    zerofy_buffer( &Y[rows * i ], rows);
                }
                
                // Incoporate the local vector chunk into the i-th chunk of the output.
                
            }
        }
        
        ptoc(ClassName()+"::Dot" );
    }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+kernel.ClassName()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
