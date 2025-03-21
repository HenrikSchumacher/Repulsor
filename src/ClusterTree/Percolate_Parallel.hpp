public:
    using TopTreeRows_T = std::vector<std::vector<Int>>;

    class Forest_T
    {
        TopTreeRows_T top_tree_rows;
        
        Tensor1<Int,Int> tree_begin;
        Tensor1<Int,Int> tree_end;
        
        JobPointers<Int> job_ptr;
        
    public:
        
        Forest_T() noexcept = default;
        
        Forest_T(
            TopTreeRows_T && top_tree_rows_,
            cref<Tensor1<Int,Int>> C_next,
            const Int thread_count_
        ) noexcept
        :   top_tree_rows ( std::move(top_tree_rows_) )
        ,   tree_begin ( static_cast<Int>(top_tree_rows.back().size()) )
        ,   tree_end   ( static_cast<Int>(top_tree_rows.back().size()) )
        {
            Tensor1<Int,Int> tree_costs ( TreeCount() + 1, static_cast<Int>(0) );
            
            for( Int i = 0; i < TreeCount(); ++i )
            {
                const Int C = top_tree_rows.back()[i];
                tree_begin[i  ] = C;
                tree_end  [i  ] = C_next[C];
                tree_costs[i+1] = C_next[C] - C;
            }
            
            tree_costs.Accumulate();
            
            job_ptr = JobPointers<Int>( TreeCount(), tree_costs.data(), thread_count_ );
            
//            print(job_ptr.ToString());
//
//            print(tree_begin.ToString());
//
//            print(tree_end.ToString());
        }
        
        ~Forest_T() = default;
        
        Forest_T( const Forest_T & other ) noexcept
        :   top_tree_rows   ( other.top_tree_rows               )
        ,   tree_begin      ( other.tree_begin                  )
        ,   tree_end        ( other.tree_end                    )
        ,   job_ptr         ( other.job_ptr                     )
        {}

        Forest_T( Forest_T && other ) noexcept
        :   top_tree_rows   ( std::move(other.top_tree_rows)    )
        ,   tree_begin      ( std::move(other.tree_begin )      )
        ,   tree_end        ( std::move(other.tree_end )        )
        ,   job_ptr         ( std::move(other.job_ptr )         )
        {}

        Forest_T & operator=( const Forest_T & other ) noexcept
        {
            top_tree_rows = other.top_tree_rows;
            tree_begin    = other.tree_begin;
            tree_end      = other.tree_end;
            job_ptr       = other.job_ptr;
            
            return *this;
        }
        
        Forest_T & operator=( const Forest_T && other ) noexcept
        {
            top_tree_rows = std::move(other.top_tree_rows);
            tree_begin    = std::move(other.tree_begin);
            tree_end      = std::move(other.tree_end);
            job_ptr       = std::move(other.job_ptr);
            
            return *this;
        }
        
        const std::vector<Int> & Row( const Int level ) const
        {
            return top_tree_rows[level];
        }
        
        Int TopLevels() const
        {
            return static_cast<Int>(top_tree_rows.size()-1);
        }
        
        Int TreeCount() const
        {
            return static_cast<Int>(top_tree_rows.back().size());
        }
        
        Int TreeBegin( const Int i ) const
        {
            return tree_begin[i];
        }
        
        Int TreeEnd( const Int i ) const
        {
            return tree_end[i];
        }
        
        const JobPointers<Int> & JobPtr() const
        {
            return job_ptr;
        }
    };


public:

    const Forest_T & Forest() const
    {
        static std::string tag ( "Forest" );
        
        if( !this->InPersistentCacheQ( tag ) )
        {
            TOOLS_PTIC(ClassName()+"::" + tag );
            
            const Int max_depth = static_cast<Int>(settings.parallel_perc_depth);
            
            TopTreeRows_T top_tree_rows (max_depth+1);
            
            Tensor1<Int,Int> stack  ( 2*max_depth+2 );
            Tensor1<Int,Int> depths ( 2*max_depth+2 );
            
            Int stack_ptr     = null;
            stack [stack_ptr] = null;
            depths[stack_ptr] = null;
            
            while( stack_ptr >= null )
            {
                if( stack_ptr >= 2*max_depth+2 )
                {
                    eprint("Stack overflow.");
                }
                
                // We are at cluster C.
                const Int d = depths[stack_ptr  ];
                const Int C = stack [stack_ptr--]; //pop
                const Int L = C_left [C];
                const Int R = C_right[C];
                
                if( d > max_depth )
                {
                    eprint("max_depth overflow.");
                }
                
                top_tree_rows[d].push_back(C);
                
                // Using that children are either both interiors or both leaves.
                if( L >= null )
                {
                    if( d < max_depth  )
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
                }
            }
            
//            this->SetPersistentCache(tag,
//                std::make_any<Forest_T>( std::move(top_tree_rows), C_next, ThreadCount() )
//            );
            
            Forest_T forest ( std::move(top_tree_rows), C_next, ThreadCount() );

            this->SetPersistentCache(tag, std::move(forest) );
            
            TOOLS_PTOC( ClassName()+"::" + tag );
        }
        
        return this->template GetPersistentCache<Forest_T>(tag);
    }

    void PercolateDown_Parallel() const override
    {
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateDown_Parallel<1>();
                break;
            }
            case AMB_DIM:
            {
                percolateDown_Parallel<AMB_DIM>();
                break;
            }
            case AMB_DIM * AMB_DIM:
            {
                percolateDown_Parallel<AMB_DIM * AMB_DIM>();
                break;
            }
            case (AMB_DIM + 1) * AMB_DIM:
            {
                percolateDown_Parallel<(AMB_DIM + 1) * AMB_DIM>();
                break;
            }
            case FAR_DIM:
            {
                percolateDown_Parallel<FAR_DIM>();
                break;
            }
            default:
            {
                percolateDown_Parallel<0>();
            }
        }
        
    }


protected:

    template<Int BUFFER_DIM>
    void percolateDown_Parallel() const
    {
        if constexpr( BUFFER_DIM > 0 )
        {
            TOOLS_PTIC(ClassName()+"::PercolateDown_Parallel<"+ToString(BUFFER_DIM)+">");
        }
        else
        {
            TOOLS_PTIC(ClassName()+"::PercolateDown_Parallel<0> ("+ToString(buffer_dim)+")");
        }
        
        cref<Forest_T> forest = Forest();
        
        logprint("Breadth-first scan for the top levels of the tree.");
        for( Int level = 0; level < forest.TopLevels(); ++level )
        {
            for( Int C : forest.Row(level) )
            {
                PercolateDown_Step<BUFFER_DIM>(C);
            }
        }
        
        logprint("Parallel depth-first scan to for the subtrees.");
        ParallelDo(
            [this,&forest]( const Int tree )
            {
                // Exploiting that the subtrees are in depth-first order, so that the inverse ordering is a post-ordering.
                const Int begin = forest.TreeBegin(tree);
                const Int end   = forest.TreeEnd  (tree);
                
//                print( "{ begin, end } =  { " + ToString(begin) + ", " +ToString(end) + " }" );
                for( Int C = begin; C < end ; ++C )
                {
                    PercolateDown_Step<BUFFER_DIM>(C);
                }
            },
            forest.JobPtr()
        );
        
        if constexpr( BUFFER_DIM > 0 )
        {
            TOOLS_PTOC(ClassName()+"::PercolateDown_Parallel<"+ToString(BUFFER_DIM)+">");
        }
        else
        {
            TOOLS_PTOC(ClassName()+"::PercolateDown_Parallel<0> ("+ToString(buffer_dim)+")");
        }

    } // percolateDown_Parallel

    template<Int BUFFER_DIM>
    void PercolateDown_Step( const Int C ) const
    {
        const Int L = C_left [C];
        
        // Using that children are either both interiors or both leaves.
        if( L >= 0 )
        {
            const Int R = C_right[C];
            
            combine_buffers<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                Scalar::One<Real>, &C_out[buffer_dim * C],
                Scalar::One<Real>, &C_out[buffer_dim * L],
                buffer_dim
            );
            
            combine_buffers<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                Scalar::One<Real>, &C_out[buffer_dim * C],
                Scalar::One<Real>, &C_out[buffer_dim * R],
                buffer_dim
            );
        }
        
    } // PercolateDown_Step

public:

    void PercolateUp_Parallel() const override
    {
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateUp_Parallel<1>();
                break;
            }
            case AMB_DIM:
            {
                percolateUp_Parallel<AMB_DIM>();
                break;
            }
            case AMB_DIM * AMB_DIM:
            {
                percolateUp_Parallel<AMB_DIM * AMB_DIM>();
                break;
            }
            case (AMB_DIM + 1) * AMB_DIM:
            {
                percolateUp_Parallel<(AMB_DIM + 1) * AMB_DIM>();
                break;
            }
            case FAR_DIM:
            {
                percolateUp_Parallel<FAR_DIM>();
                break;
            }
            default:
            {
                percolateUp_Parallel<0>();
            }
        }
    }


protected:

    template<Int BUFFER_DIM>
    void percolateUp_Parallel() const
    {
        if constexpr( BUFFER_DIM > 0 )
        {
            TOOLS_PTIC(ClassName()+"::PercolateUp_Parallel<"+ToString(BUFFER_DIM)+">");
        }
        else
        {
            TOOLS_PTIC(ClassName()+"::PercolateUp_Parallel<0> ("+ToString(buffer_dim)+")");
        }
        
        cref<Forest_T> forest = Forest();
        
        logprint("Parallel reverse depth-first scan to for subtrees.");
        ParallelDo(
            [this,&forest]( const Int tree )
            {
                // Exploiting that the subtrees are in depth-first order, so that the inverse ordering is a post-ordering.
                const Int begin = forest.TreeBegin(tree);
                const Int end   = forest.TreeEnd  (tree);
                
//                print( "{ begin, end } =  { " + ToString(begin) + ", " +ToString(end) + " }" );
                for( Int C = end; C --> begin; )
                {
                    PercolateUp_Step<BUFFER_DIM>(C);
                }
            },
            forest.JobPtr()
        );
        
        
        logprint("Breadth-first scan for the top levels of the tree.");
        for( Int level = forest.TopLevels(); level --> 0 ; )
        {
            for( Int C : forest.Row(level) )
            {
                PercolateUp_Step<BUFFER_DIM>(C);
            }
        }

        if constexpr( BUFFER_DIM > 0 )
        {
            TOOLS_PTOC(ClassName()+"::PercolateUp_Parallel<"+ToString(BUFFER_DIM)+">");
        }
        else
        {
            TOOLS_PTOC(ClassName()+"::PercolateUp_Parallel<0> ("+ToString(buffer_dim)+")");
        }
        
    } // percolateUp_Parallel



    template<Int BUFFER_DIM>
    void PercolateUp_Step( const Int C ) const
    {
        const Int L = C_left [C];
        
        if( L >= 0 )
        {
            const Int R = C_right[C];
            
            combine_buffers3<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                Scalar::One<Real>, &C_in[buffer_dim * L],
                Scalar::One<Real>, &C_in[buffer_dim * R],
                                   &C_in[buffer_dim * C],
                buffer_dim
            );
        }
    } // PercolateUp_Step
