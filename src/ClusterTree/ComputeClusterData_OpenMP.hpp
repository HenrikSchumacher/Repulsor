protected:
        
    void ComputeClusterData() const
    {
        ptic(className()+"::ComputeClusterData (OpenMP)");
        
        // using the already serialized cluster tree
        #pragma omp parallel num_threads( ThreadCount() )
        {
            #pragma omp single nowait
            {
                computeClusterData( 0, ThreadCount() );
            }
        }
        ptoc(className()+"::ComputeClusterData (OpenMP)");
    }; //ComputeClusterData


    void computeClusterData( const Int C, const Int free_thread_count ) const // helper function for ComputeClusterData
    {
        const Int far_dim = FarDim();
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        mut<Real> C_C = C_far.data(C);
        
        if( L >= 0 && R >= 0 )
        {
            //C points to interior node.
            #pragma omp task final(free_thread_count<1)  shared( L )
            {
                computeClusterData( L, free_thread_count/2 );
            }
            #pragma omp task final(free_thread_count<1)  shared( R )
            {
                computeClusterData( R, free_thread_count-free_thread_count/2 );
            }
            #pragma omp taskwait

            ptr<Real> C_L = C_far.data(L);
            ptr<Real> C_R = C_far.data(R);
            
            Real L_weight = C_L[0];
            Real R_weight = C_R[0];
            
            const Real C_mass = L_weight + R_weight;
            C_C[0] = C_mass;
            const Real C_invmass = Scalar::Inv<Real>(C_mass);
            
            L_weight *= C_invmass;
            R_weight *= C_invmass;
            
            for( Int k = 1; k < far_dim; ++k )
            {
                C_C[k] = L_weight * C_L[k] + R_weight * C_R[k] ;
            }
        }
        else
        {
            //C points to leaf node.
            //compute from primitives
            const Int begin = C_begin[C];
            const Int end   = C_end  [C];
            
            zerofy_buffer( C_C, far_dim );
            
            for( Int i = begin; i < end; ++i )
            {
                ptr<Real> P = P_far.data(i);
                
                const Real a = P[0];
                
                C_C[0] += a;
                
                for( Int j = 1; j < far_dim; ++j )
                {
                    C_C[j] += a * P[j];
                }
            }
            
            const Real C_invmass = Scalar::Inv<Real>(C_C[0]);
    
            for( Int j = 1; j < far_dim; ++j )
            {
                C_C[j] *= C_invmass;
            }
        }
    }; //computeClusterData
