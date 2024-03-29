namespace Repulsor
{
    template <int AMB_DIM, typename SReal, typename Int>
    void Compute_AABB_CollisionTimeInterval(
       cptr<SReal> p_0, cptr<SReal> p_1,
       cptr<SReal> q_0, cptr<SReal> q_1,
       mref<SReal> t_first,    // returning per reference to avoid std::pair
       mref<SReal> t_last      // returning per reference to avoid std::pair
    )
    {
        // Assuming p_0, p_1, q_0, q_1 are serizalizations in the format of some AABB<AMB_DIM,GJK_Real,Int,SReal>
        
        // Returns interval {t_first,t_last} of those t in [0,1] such that the AABBs defined by
        //     P(t) = (1-t) * p_0 + t * p_1
        // and
        //     Q(t) = (1-t) * q_0 + t * q_1
        // intersect.
        
        static constexpr SReal zero = Scalar::Zero<SReal>;
        static constexpr SReal one  = Scalar::One<SReal>;
        
        SReal t_0 = zero;
        SReal t_1 = one;
        
        for( Int k = 0; k < AMB_DIM; ++k )
        {
            // Compute the intervals of the four AABB in the k-th coordinate direction.
            
            const Int i = 1+k;
            const Int j = 1+k+AMB_DIM;
            
            const SReal P_min_0 = p_0[i] - p_0[j];
            const SReal P_max_0 = p_0[i] + p_0[j];
            
            const SReal P_min_1 = p_1[i] - p_1[j];
            const SReal P_max_1 = p_1[i] + p_1[j];
            
            const SReal Q_min_0 = q_0[i] - q_0[j];
            const SReal Q_max_0 = q_0[i] + q_0[j];
            
            const SReal Q_min_1 = q_1[i] - q_1[j];
            const SReal Q_max_1 = q_1[i] + q_1[j];
            
            bool initially_intersecting = (P_min_0 <= Q_max_0) && (Q_min_0 <= P_max_0);
            //            bool finally_intersecting   = (P_min_1 <= Q_max_1) && (Q_min_1 <= P_max_1);
            
            const SReal m_0 = Q_max_0 - Q_max_1 - P_min_0 + P_min_1;
            const SReal m_1 = P_max_0 - P_max_1 - Q_min_0 + Q_min_1;
            
            SReal A = ( m_0 != zero ) ? Min(one,(Q_max_0 - P_min_0)/m_0) : one;
            SReal B = ( m_1 != zero ) ? Min(one,(P_max_0 - Q_min_0)/m_1) : one;
            
            if( A > B )
            {
                std::swap(A,B);
            }
            
            SReal a;
            SReal b;
            
            //            GJK_DUMP(initially_intersecting);
            //            GJK_DUMP(A);
            //            GJK_DUMP(B);
            
            if( initially_intersecting )
            {
                a = zero;
                b = (A > zero) ? A : ( (B <= zero) ? one : B );
            }
            else
            {
                // initially_intersecting == false
                
                if( /*A < zero &&*/ B < zero )
                {
                    a = one;
                    b = one;
                }
                else
                {
                    if( A >= zero/* && B >= zero*/ )
                    {
                        a = A;
                        b = B;
                    }
                    else
                    {
                        a = B;
                        b = one;
                    }
                }
            }
            
            t_0 = Max(t_0,a);
            t_1 = Min(t_1,b);
        }
        
        if( t_1 < t_0 )
        {
            // empty intersection interval
            t_first = one;
            t_last  = one;
        }
        else
        {
            t_first = t_0;
            t_last  = t_1;
        }
    }

} // namespace Repulsor
