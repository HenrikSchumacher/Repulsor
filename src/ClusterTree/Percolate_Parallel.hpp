protected:

    using PercolationRoots_T = std::vector<Int>;

    const PercolationRoots_T & ParallelPercolationRoots() const
    {
        static std::string tag ( "ParallelPercolationRoots" );
        
        if( !this->InPersistentCacheQ( tag ) )
        {
            ptic(ClassName()+"::" + tag );
            
            const Int max_depth = static_cast<Int>(settings.parallel_perc_depth);
            
            PercolationRoots_T parallel_perc_roots;
            
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
            
            this->SetPersistentCache( tag, std::any( std::move(parallel_perc_roots) ) );
            
            ptoc( ClassName()+"::" + tag );
        }
        
        return std::any_cast<const PercolationRoots_T &>( this->GetPersistentCache( tag ) );
    }

public:

    void PercolateUp_Parallel() const override
    {
        ptic(ClassName()+"::PercolateUp_Parallel");
        
        const PercolationRoots_T & parallel_perc_roots = ParallelPercolationRoots();
        
        // TODO: Add thread balancing.
        ParallelDo(
            [this,&parallel_perc_roots]( const Int i )
            {
                PercolateUp_DFS(parallel_perc_roots[i]);
            },
            static_cast<Int>(parallel_perc_roots.size()),
            ThreadCount()
        );
        
        PercolateUp_DFS( null, static_cast<Int>(settings.parallel_perc_depth) );

        ptoc(ClassName()+"::PercolateUp_Parallel");
    } // PercolateUp_Parallel

    void PercolateDown_Parallel() const override
    {
        ptic(ClassName()+"::PercolateDown_Parallel");
        
        const PercolationRoots_T & parallel_perc_roots = ParallelPercolationRoots();
        
        PercolateDown_DFS( null, static_cast<Int>(settings.parallel_perc_depth) );
        
        // TODO: Add thread balancing.
        ParallelDo(
            [this,&parallel_perc_roots]( const Int i )
            {
                PercolateDown_DFS(parallel_perc_roots[i]);
            },
            static_cast<Int>(parallel_perc_roots.size()),
            ThreadCount()
        );

        ptoc(ClassName()+"::PercolateDown_Parallel");
        
    } // PercolateDown_Parallel

