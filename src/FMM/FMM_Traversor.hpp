#pragma once

// TODO: Handle symmetric/asymmetric cases correctly!

namespace Repulsor
{
    
    // FMM_Traversor is a class to map kernel over all nonzero positions of a sparsity pattern (SparsityPatternCSR<Int>).
    // This can be used, for example, to compute nonlocal energies, their derivatives, and the nonzero values of hierarchical matrices.
    // Mostly intended to be used with the results of BlockClusterTree's GetMatrix<...>() routine the verynear/near/far block matrices.
    
    template<class Pattern_T, class Kernel_T>
    class FMM_Traversor
    {
    public:
        
        using Real     = typename Kernel_T::Real;
        using SReal    = typename Kernel_T::SReal;
        using ExtReal  = typename Kernel_T::ExtReal;
        using Int      = typename Kernel_T::Int;
        using LInt     = typename Kernel_T::LInt ;
                
        using Values_T          = Tensor2<Real,LInt>;
        using ValueContainer_T  = std::unordered_map<std::string,Values_T>;
        
        static constexpr bool symmetricQ = Kernel_T::symmetricQ;
        
    public:
        
        FMM_Traversor( cref<Pattern_T> pattern_, mref<Kernel_T> kernel_ )
        :   pattern ( pattern_ )
        ,   kernel  ( kernel_  )
        {}
        
        // Default constructor
        FMM_Traversor() = delete;
        // Destructor
        virtual ~FMM_Traversor() = default;
        // Copy constructor
        FMM_Traversor( const FMM_Traversor & other ) = default;
        // Copy assignment operator
        FMM_Traversor & operator=( const FMM_Traversor & other ) = default;
        // Move constructor
        FMM_Traversor( FMM_Traversor && other ) = default;
        // Move assignment operator
        FMM_Traversor & operator=( FMM_Traversor && other ) = default;
        
    protected:
    
        cref<Pattern_T> pattern;
        
        Kernel_T kernel;
        
    public:

        TOOLS_FORCE_FLATTENING Real Compute()
        {
            TOOLS_PTIC(ClassName()+"::Compute");

            if( pattern.NonzeroCount() <= 0 )
            {
                TOOLS_PTOC(ClassName()+"::Compute");
                return static_cast<Real> (0);
            }

            cref<JobPointers<Int>> job_ptr = (
                symmetricQ ?
                pattern.UpperTriangularJobPtr() :
                pattern.JobPtr()
            );

            // Make sure that pattern.Diag() is created before the parallel region.
            if constexpr ( symmetricQ )
            {
                (void)pattern.Diag();
            }

            kernel.Allocate( pattern.NonzeroCount() );
            
            // DEBUG THIS!
            Real global_sum = ParallelDoReduce(
                [=,this,&job_ptr]( const Int thread )
                {
                    Kernel_T ker ( kernel, thread );

                    Real local_sum (0);
                    
                    cptr<LInt> diag = (symmetricQ ? pattern.Diag().data() : nullptr);
                    cptr<LInt> rp   = pattern.Outer().data();
                    cptr<Int>  ci   = pattern.Inner().data();
                    
                    // Kernel is supposed the following rows of pattern:
                    const Int i_begin = job_ptr[thread  ];
                    const Int i_end   = job_ptr[thread+1];
                    
                    const Time start_time = Clock::now();
                    
//                    print("Loop on " + ToString(thread) + " begins");
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        // These are the corresponding nonzero blocks in i-th row.
                        const LInt k_begin = ( symmetricQ ? diag[i] : rp[i] );
                        const LInt k_end   = rp[i+1];
                        
                        if( k_end > k_begin )
                        {
                            // Clear the local vector chunk of the kernel.
                            ker.LoadS(i);
                            
                            // Perform all but the last calculation in row with prefetch.
                            for( LInt k = k_begin; k < k_end-1; ++k )
                            {
                                const Int j = ci[k];
                                
                                ker.LoadT(j);
                                
                                ker.Prefetch(ci[k+1]);

                                local_sum += ker.Compute(k);

                                ker.WriteT(j);
                            }
                            
                            // Perform last calculation in row without prefetch.
                            {
                                const LInt k = k_end-1;
                                
                                const Int j = ci[k];
                                
                                ker.LoadT(j);
                                
                                local_sum += ker.Compute(k);
                                
                                ker.WriteT(j);
                            }
                            
                            // Incorporate the kernel's local vector chunk into the i-th chunk if the output Y.
                            
                            ker.WriteS(i);
                        }
                        
                        // Incorporate the local vector chunk into the i-th chunk of the output.
                        
                    }
                    
                    const Time stop_time = Clock::now();
                    
                    ker.PrintReport( Tools::Duration(start_time,stop_time) );
                    
                    return local_sum; 
                },
                AddReducer<Real, Real>(),
                Scalar::Zero<Real>,
                job_ptr.ThreadCount()
            );
            
            TOOLS_PTOC(ClassName()+"::Compute");

            return global_sum;
        }
        
    public:
        
        std::string ClassName() const
        {
            return "FMM_Traversor<"+kernel.ClassName() + ","+TypeName<LInt>+">";
        }
        
    }; // class FMM_Traversor
    
}// namespace Repulsor

