public:

    void RequireParallelPercolationRoots() const
    {
        if( ! parallel_perc_roots_initialized )
        {
            ptic(ClassName()+"::RequireParallelPercolationRoots");

            const Int max_depth = static_cast<Int>(settings.parallel_perc_depth);
            
            Tensor1<Int,Int> stack  ( 2*max_depth+2 );
            Tensor1<Int,Int> depths ( 2*max_depth+2 );
            
            Int stack_ptr     = null;
            stack [stack_ptr] = null;
            depths[stack_ptr] = null;
            
            while( stack_ptr >= null )
            {
                // We are at cluster C.
                const Int d = depths[stack_ptr];
                const Int C = stack [stack_ptr--]; //pop
                const Int L = C_left[C];
                const Int R = C_right[C];
                
                if( ( d < max_depth ) && (L >= null) && (R >= null) )
                {
                    // push
                    ++stack_ptr;
                    stack [stack_ptr] = R;
                    depths[stack_ptr] = d+1;
                    
                    // push
                    ++stack_ptr;
                    stack [stack_ptr] = L;
                    depths[stack_ptr] = d+1;
                }
                else
                {
                    parallel_perc_roots.push_back(C);
                }
            }
            
            parallel_perc_roots_initialized = true;
            
            ptoc(ClassName()+"::RequireParallelPercolationRoots");
        }
    }

    void PercolateUp_Parallel() const
    {
        ptic(ClassName()+"::PercolateUp_Parallel");
        
        Int n_roots = static_cast<Int>(parallel_perc_roots.size());
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int i = 0; i < n_roots; ++i )
        {
            PercolateUp_DFS(parallel_perc_roots[i]);
        }
        
        PercolateUp_DFS( null, static_cast<Int>(settings.parallel_perc_depth) );

        ptoc(ClassName()+"::PercolateUp_Parallel");
    } // PercolateUp_Parallel

    void PercolateDown_Parallel() const
    {
        ptic(ClassName()+"::PercolateDown_Parallel");
        
        Int n_roots = static_cast<Int>(parallel_perc_roots.size());
        
        PercolateDown_DFS( null, static_cast<Int>(settings.parallel_perc_depth) );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int i = 0; i < n_roots; ++i )
        {
            PercolateDown_DFS(parallel_perc_roots[i]);
        }

        ptoc(ClassName()+"::PercolateDown_Parallel");
        
    } // PercolateDown_Parallel
