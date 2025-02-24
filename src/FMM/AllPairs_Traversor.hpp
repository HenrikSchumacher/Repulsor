#pragma once

// TODO: Handle symmetric/asymmetric cases correctly!

namespace Repulsor
{
    
    // FMM_Traversor is a class to map kernel over all positions of an m x n matrix.
    // This can be used, for example, to compute nonlocal energies, their derivatives, and the nonzero values of hierarchical matrices.
    // Mostly intended to be used for comparison against FMM_Traversor.
    
    template<class Kernel_T>
    class AllPairs_Traversor
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
        
        AllPairs_Traversor() = default;
        
        template<typename I0, typename I1>
        AllPairs_Traversor( const I0 m_, const I1 n_, mref<Kernel_T> kernel_ )
        :   m       ( m_      )
        ,   n       ( n_      )
        ,   kernel  ( kernel_ )
        {
            static_assert(IntQ<I0>,"");
            static_assert(IntQ<I1>,"");
        }
        
        // Copy constructor
        AllPairs_Traversor( cref<AllPairs_Traversor> other )
        :   m       ( other.m      )
        ,   n       ( other.n      )
        ,   kernel  ( other.kernel )
        {}

        ~AllPairs_Traversor() = default;
        
    protected:
    
        const Int m = 0;
        const Int n = 0;
        
        Kernel_T kernel;
        
    public:

        TOOLS_FORCE_FLATTENING Real Compute()
        {
            TOOLS_PTIC(ClassName()+"::Compute");
            
            const Int thread_count = kernel.ThreadCount();
            
            // TODO: Symmetrization
            JobPointers<Int> job_ptr;
            
            if constexpr ( symmetricQ )
            {
                Tensor1<LInt, Int> costs (m+1);
                costs[0] = 0;
                for( Int i = 0; i < m; ++i )
                {
                    costs[i+1] = costs[i] + m-1-i;
                }
            
                job_ptr = JobPointers<Int>( m, costs.data(), thread_count, false );
            }
            else
            {
                job_ptr = JobPointers<Int>( m, thread_count );
            }
            
            kernel.Allocate( m * n );

            Real global_sum = Scalar::Zero<Real>;

            global_sum = ParallelDoReduce(
                [&,this]( const Int thread ) -> Real
                {
                    // Initialize local kernel and feed it all the information that is going to be constant along its life time.
                    
                    Kernel_T ker ( kernel );

                    Real local_sum = Scalar::Zero<Real>;
                    
                    // Kernel is supposed the following rows of pattern:
                    const Int i_begin = job_ptr[thread  ];
                    const Int i_end   = job_ptr[thread+1];

                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        // These are the corresponding nonzero blocks in i-th row.
                        const Int j_begin = ( symmetricQ ? i+1 : 0 );
                        const Int j_end   = n;
                        
                        // Clear the local vector chunk of the kernel.
                        ker.LoadS(i);
                        
                        const LInt offset = n*i;
                        
                        for( Int j = j_begin; j < j_end; ++j )
                        {
                            ker.LoadT(j);
                            
                            local_sum += ker.Compute(offset+j);

                            ker.WriteT(j);
                        }
                        
                        // Incorporate the kernel's local vector chunk into the i-th chunk if the output Y.
                        
                        ker.WriteS(i);
                        
                        // Incoporate the local vector chunk into the i-th chunk of the output.
                        
                    }
                    
                    return local_sum;
                },
                AddReducer<Real, Real>(),
                Scalar::Zero<Real>,
                thread_count
            );

            TOOLS_PTOC(ClassName()+"::Compute");

            return global_sum;
        }
        
    public:
        
        std::string ClassName() const
        {
            return "AllPairs_Traversor<"+kernel.ClassName()+TypeName<LInt>+">";
        }
        
    };
    
}// namespace Repulsor

