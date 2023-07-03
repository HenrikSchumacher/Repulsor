protected:
        
    void  ComputeClusterData() const
    {
        if( SplitThreshold() > 1 )
        {
            computeClusterData<false>();
        }
        else
        {
            computeClusterData<true>();
        }
    }

    template<bool leaves_are_singletons>
    void computeClusterData() const
    {
        ptic(ClassName()+"::ComputeClusterData<"+ToString(leaves_are_singletons)+">");
        
        const Forest_T & restrict forest = Forest();
        
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
                    ComputeClusterData_Step<leaves_are_singletons>(C);
                }
            },
            forest.JobPtr()
        );
        
        
        logprint("Breadth-first scan to for the top levels of the tree.");
        for( Int level = forest.TopLevels(); level --> 0 ; )
        {
            for( Int C : forest.Row(level) )
            {
                ComputeClusterData_Step<leaves_are_singletons>(C);
            }
        }

        ptoc(ClassName()+"::ComputeClusterData<"+ToString(leaves_are_singletons)+">");
        
    } // computeClusterData

    template<bool leaves_are_singletons>
    void ComputeClusterData_Step( const Int C ) const
    {
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        
        if( L >= 0 )
        {
            mut<Real> C_C = C_far.data(C);
            ptr<Real> C_L = C_far.data(L);
            ptr<Real> C_R = C_far.data(R);
            
            Real L_weight = C_L[0];
            Real R_weight = C_R[0];
            
            const Real C_mass = L_weight + R_weight;
            C_C[0] = C_mass;
            const Real C_invmass = Scalar::Inv<Real>(C_mass);
            
            L_weight *= C_invmass;
            R_weight *= C_invmass;
            
            for( Int k = 1; k < FAR_DIM; ++k )
            {
                C_C[k] = L_weight * C_L[k] + R_weight * C_R[k] ;
            }
        }
        else
        {
            //C points to leaf node; compute from primitive(s).
            
            if constexpr ( leaves_are_singletons )
            {
                copy_buffer<FAR_DIM,Sequential>( P_far.data(C_begin[C]), C_far.data(C) );
            }
            else
            {
                mut<Real> C_C = C_far.data(C);
                
                const Int begin = C_begin[C];
                const Int end   = C_end  [C];
                
                {
                    const Int i = begin;
                    ptr<Real> P = P_far.data(i);
                    
                    const Real a = P[0];
                    
                    C_C[0] = a;
                    
                    for( Int j = 1; j < FAR_DIM; ++j )
                    {
                        C_C[j] = a * P[j];
                    }
                }
                
                for( Int i = begin+1; i < end; ++i )
                {
                    ptr<Real> P = P_far.data(i);
                    
                    const Real a = P[0];
                    
                    C_C[0] += a;
                    
                    for( Int j = 1; j < FAR_DIM; ++j )
                    {
                        C_C[j] += a * P[j];
                    }
                }
                
                const Real C_invmass = Scalar::Inv<Real>(C_C[0]);
                
                for( Int j = 1; j < FAR_DIM; ++j )
                {
                    C_C[j] *= C_invmass;
                }
            }
        }
    } // ComputeClusterData_Step
