//#################################################################################
//##        Distance queries
//#################################################################################
    
    template<typename T, typename I>
    void FindNearestToPoints(
        cptr<T> query_points,
        Int     query_point_count,
        mptr<I> nearest_primitives,
        mptr<T> nearest_dist,
        mptr<T> nearest_points
    ) const
    {
        TOOLS_PTIC(ClassName()+"::FindNearestToPoints");
        using GJK_T    = GJK<AMB_DIM,GJK_Real,Int>;
        using Point_T  = Point<AMB_DIM,GJK_Real,Int,SReal,T>;
        using Vector_T = Tiny::Vector<AMB_DIM,GJK_Real,Int>;
        
        const Int thread_count = query_point_count >= Int(128) ? ThreadCount() : Int(1);
        
        ParallelDo(
            [=,this]( Int thread )
            {
                Int stack [128];
//                Tensor1<Int,Int> stack ( max_depth + 1 );
                
                GJK_T G;
                
                mptr<SReal> C_ser = C_serialized.data();
                mptr<SReal> P_ser = P_serialized.data();
                
                std::shared_ptr<BoundingVolume_T> C_box  = C_proto[thread];
                std::shared_ptr<Primitive_T     > P_prim = P_proto[thread];
                
                SReal x_buffer [AMB_DIM + 1] = {};
                Point_T x;
                x.SetPointer(&x_buffer[0]);
                
                const Int i_begin
                   = JobPointer( query_point_count, thread_count, thread     );
                const Int i_end
                   = JobPointer( query_point_count, thread_count, thread + 1 );
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    Vector_T P_witness;
                    Vector_T x_witness;
                    
                    x.FromCoordinates( query_points, i );

                    Int  best_idx = 0;
                    GJK_Real best_r2 = Scalar::Max<Real>;
                    Vector_T best_pt;
                    
                    Int stack_ptr = 0;
                    stack[0] = 0;
                    
                    while( (0 <= stack_ptr) && (stack_ptr < 126) )
                    {
                        const Int C = stack[stack_ptr];
                        --stack_ptr;
                        
                        const Int L = C_left [C];
                        const Int R = C_right[C];
                        
                        // Using that children are either both interiors or both leaves.
                        if( (L >= 0) /*|| (R >= 0)*/ )
                        {
                            C_box->SetPointer( C_ser, L );
                            const GJK_Real r2_to_L = G.SquaredDistance( *C_box, x );
                            
                            C_box->SetPointer( C_ser, R );
                            const GJK_Real r2_to_R = G.SquaredDistance( *C_box, x );
                            
                            if( r2_to_L < r2_to_R )
                            {
                                // Put the nearest cluster on top of the stack.
                                
                                if( r2_to_R < best_r2 )
                                {
                                    ++stack_ptr;
                                    stack[stack_ptr] = R;
                                }
                                
                                if( r2_to_L < best_r2 )
                                {
                                    ++stack_ptr;
                                    stack[stack_ptr] = L;
                                }
                            }
                            else
                            {
                                // Put the nearest cluster on top of the stack.
                                
                                if( r2_to_L < best_r2 )
                                {
                                    ++stack_ptr;
                                    stack[stack_ptr] = L;
                                }
                                
                                if( r2_to_R < best_r2 )
                                {
                                    ++stack_ptr;
                                    stack[stack_ptr] = R;
                                }
                            }
                        }
                        else // (L < null) && (R < null)
                        {
                            const Int P_begin = C_begin[C];
                            const Int P_end   = C_end  [C];
                            
                            for( Int P = P_begin; P < P_end; ++ P )
                            {
                                P_prim->SetPointer( P_ser, P );
                                
                                GJK_Real r2 = G.Witnesses(
                                    *P_prim, P_witness.data(),
                                    x      , x_witness.data()
                                );
                                
                                if( r2 < best_r2 )
                                {
                                    best_idx = P;
                                    best_r2  = r2;
                                    best_pt  = P_witness;
                                }
                            }
                        }
                        
                    } // while( (0 <= stack_ptr) && (stack_ptr < max_depth) )
                    
                    if( stack_ptr >= 126 )
                    {
                        eprint( ClassName()+"::FindNearestToPoints: Stack overflow.");
                    }
                    
                    nearest_primitives[i] = static_cast<I>(P_ordering[best_idx]);
                    nearest_dist[i] = static_cast<T>(Sqrt(best_r2));
                    best_pt.Write( &nearest_points[AMB_DIM * i] );
                    
                } // for( Int i = 0; i < query_point_count; ++i )
                
            },
            thread_count
        );
        
        TOOLS_PTOC(ClassName()+"::FindNearestToPoints");
    }


    virtual void FindNearestToPoints(
        cptr<Real> query_points,
        Int        query_point_count,
        mptr<Int>  nearest_primitives,
        mptr<Real> nearest_dist,
        mptr<Real> nearest_points
    ) const override
    {
        FindNearestToPoints<Real,Int>(
            query_points,query_point_count,nearest_primitives,nearest_dist,nearest_points
        );
    }
