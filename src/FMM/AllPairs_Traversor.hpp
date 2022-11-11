#pragma once

#define CLASS AllPairs_Traversor

// TODO: Handle symmetric/asymmetric cases correctly!

namespace Repulsor
{
    
    // FMM_Traversor is a class to map kernel over all positions of an m x n matrix.
    // This can be used, for example, to compute nonlocal energies, their derivatives, and the nonzero values of hierarchical matrices.
    // Mostly intended to be used for comparison against FMM_Traversor.
    
    template<class Kernel_T>
    class CLASS
    {
    public:
        
        using Real     = typename Kernel_T::Real;
        using SReal    = typename Kernel_T::SReal;
        using ExtReal  = typename Kernel_T::ExtReal;
        using Int      = typename Kernel_T::Int;
        using LInt     = typename Kernel_T::LInt ;
                
        using Values_T          = Tensor2<Real,LInt>;
        using ValueContainer_T  = std::unordered_map<std::string,Values_T>;
        
        static constexpr bool is_symmetric = Kernel_T::is_symmetric;
        
        CLASS() = default;
        
        template<typename I0, typename I1>
        CLASS( const I0 m_, const I1 n_, const Kernel_T & kernel_ )
        :   m       ( m_      )
        ,   n       ( n_      )
        ,   kernel  ( kernel_ )
        {
            ASSERT_INT(I0);
            ASSERT_INT(I1);
        }
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   m       ( other.m      )
        ,   n       ( other.n      )
        ,   kernel  ( other.kernel )
        {}

        ~CLASS() = default;
        
    protected:
    
        const LInt m = 0;
        const LInt n = 0;
        
        Kernel_T kernel;
        
    public:

        __attribute__((flatten)) Real Compute()
        {
            ptic(ClassName()+"::Compute");
            
            const Int thread_count = kernel.ThreadCount();
            
            // TODO: Symmetrization
            JobPointers<LInt> job_ptr;
            
            if constexpr ( is_symmetric )
            {
                Tensor1<LInt, LInt> costs (m+1);
                costs[0] = 0;
                for( LInt i = 0; i < m; ++i )
                {
                    costs[i+1] = costs[i] + m-1-i;
                }
            
                job_ptr = JobPointers<LInt>( m, costs.data(), thread_count, false );
            }
            else
            {
                job_ptr = JobPointers<LInt>( m, thread_count );
            }
            
            print(job_ptr.ToString());
            
            kernel.Allocate( m * n );

            Real global_sum = static_cast<Real>(0);

            if( thread_count > 1 )
            {
                #pragma omp parallel for num_threads( thread_count ) reduction( + : global_sum )
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    // Initialize local kernel and feed it all the information that is going to be constant along its life time.
                    
                    Kernel_T ker ( kernel );

                    
                    Real local_sum = static_cast<Real>(0);
                    
                    // Kernel is supposed the following rows of pattern:
                    const LInt i_begin = job_ptr[thread  ];
                    const LInt i_end   = job_ptr[thread+1];

                    #pragma clang loop unroll_count(4)
                    for( LInt i = i_begin; i < i_end; ++i )
                    {
                        // These are the corresponding nonzero blocks in i-th row.
                        const LInt j_begin = COND( is_symmetric, i+1, 0 );
                        const LInt j_end   = n;
                        
                        // Clear the local vector chunk of the kernel.
                        ker.LoadS(i);
                        
                        const LInt offset = n*i;
                        
                        #pragma clang loop unroll_count(4)
                        for( LInt j = j_begin; j < j_end; ++j )
                        {
                            ker.LoadT(j);
                            
                            local_sum += ker.Compute(offset+j);

                            ker.WriteT(j);
                        }
                        
                        // Incorporate the kernel's local vector chunk into the i-th chunk if the output Y.
                        
                        ker.WriteS(i);
                        
                        // Incoporate the local vector chunk into the i-th chunk of the output.
                        
                    }
                    
                    global_sum += local_sum;
                }
            }
            else
            {
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    // Initialize local kernel and feed it all the information that is going to be constant along its life time.
                    
                    Kernel_T ker ( kernel );

                    Real local_sum = static_cast<Real>(0);

                    #pragma clang loop unroll_count(4)
                    for( LInt i = 0; i < m; ++i )
                    {
                        // These are the corresponding nonzero blocks in i-th row.
                        const LInt j_begin = COND( is_symmetric, i+1, 0 );
                        const LInt j_end   = n;
                        
                        // Clear the local vector chunk of the kernel.
                        ker.LoadS(i);
                        
                        const LInt offset = n*i;
                        
                        // Perform all but the last calculation in row with prefetch.
                        #pragma clang loop unroll_count(4)
                        for( LInt j = j_begin; j < j_end; ++j )
                        {
                            ker.LoadT(j);
                            
                            local_sum += ker.Compute(offset+j);

                            ker.WriteT(j);
                        }
                        
                        // Incorporate the kernel's local vector chunk into the i-th chunk if the output Y.
                        
                        ker.WriteS(i);
                        
                        // Incoporate the local vector chunk into the i-th chunk of the output.
                        
                    }
                    
                    global_sum += local_sum;
                }
            }

            ptoc(ClassName()+"::Compute");

            return global_sum;
        }
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+kernel.ClassName()+TypeName<LInt>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef CLASS

