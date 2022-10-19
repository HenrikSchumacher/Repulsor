#pragma once

#define CLASS FMM_Traversor

// TODO: Handle symmetric/asymmetric cases correctly!

namespace Repulsor
{    
    template<typename Kernel_T>
    class CLASS
    {
    public:
        
        using Real    = typename Kernel_T::Real;
        using Int     = typename Kernel_T::Int;
        using SReal   = typename Kernel_T::SReal;
        using ExtReal = typename Kernel_T::ExtReal;
        
        using ValueContainer_T  = Tensor2<Real,Int>;
        
        using SparsityPattern_T = SparsityPatternCSR<Int>;
        
        static constexpr bool is_symmetric = Kernel_T::is_symmetric;
        
        CLASS() = default;
        
        explicit CLASS(
            const SparsityPattern_T & pattern_,
            const Kernel_T & kernel_
        )
        :   pattern ( pattern_ )
        ,   kernel  ( kernel_  )
        {
        }
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   pattern ( other.pattern )
        ,   kernel  ( other.kernel )
        {}

        ~CLASS() = default;
        
    protected:
        
        const SparsityPattern_T   & pattern;
        
        Kernel_T kernel;
        
    public:
        

//##############################################################################################
//      Matrix multiplication
//##############################################################################################
        
        Real Compute() const
        {
            ptic(ClassName()+"::Compute");

            if( pattern.NonzeroCount() <= 0 )
            {
                ptoc(ClassName()+"::Compute");
                return static_cast<Real> (0);
            }
            
            const auto & job_ptr = COND(
                is_symmetric,
                pattern.UpperTriangularJobPtr(),
                pattern.JobPtr()
            );
            
            // Make sure that pattern.Diag() is created before the parallel region.
            if constexpr ( is_symmetric )
            {
                (void)pattern.Diag();
            }
            
            if constexpr ( Kernel_T::metric_flag )
            {
                kernel.GetS().CleanseDerivativeBuffers();
                
                if( !is_symmetric )
                {
                    kernel.GetT().CleanseDerivativeBuffers();
                }
            }
            
            const Int thread_count = job_ptr.Size()-1;
            
            Real global_sum = static_cast<Real>(0);
            
            #pragma omp parallel for num_threads( thread_count ) reduction( + : global_sum )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                // Initialize local kernel and feed it all the information that is going to be constant along its life time.

                Kernel_T ker ( kernel );
                
                Real local_sum = static_cast<Real>(0);
                
                const Int * restrict const diag = COND(is_symmetric,
                    pattern.Diag().data(),
                    nullptr
                );
                
                const Int * restrict const rp = pattern.Outer().data();
                const Int * restrict const ci = pattern.Inner().data();
                
                // Kernel is supposed the following rows of pattern:
                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    // These are the corresponding nonzero blocks in i-th row.
                    const Int k_begin = COND( is_symmetric, diag[i]+1, rp[i] );
                    const Int k_end   = rp[i+1];
                    
                    if( k_end > k_begin )
                    {
                        // Clear the local vector chunk of the kernel.
                        ker.LoadS(i);
                        
                        // Perform all but the last calculation in row with prefetch.
                        for( Int k = k_begin; k < k_end-1; ++k )
                        {
                            const Int j = ci[k];

                            ker.LoadT(j);
                            
                            ker.PrefetchT(ci[k+1]);

                            local_sum += ker.Compute(k);
                            
                            ker.WriteT();
                        }
                        
                        // Perform last calculation in row without prefetch.
                        {
                            const Int k = k_end-1;
                            
                            const Int j = ci[k];
                            
                            ker.LoadT(j);
                            
                            local_sum += ker.Compute(k);
                            
                            ker.WriteT();
                        }
                        
                        // Incorporate the kernel's local vector chunk into the i-th chunk if the output Y.
                        
                        ker.WriteS();
                    }
                    
                    // Incoporate the local vector chunk into the i-th chunk of the output.
                    
                }
                
                global_sum += local_sum;
            }
        
            ptoc(ClassName()+"::Compute");

            return global_sum;
        }
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+kernel.ClassName()+">";
        }
        
    };
    
}// namespace Repulsor

#undef CLASS

