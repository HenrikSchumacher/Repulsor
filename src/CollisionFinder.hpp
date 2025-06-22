#pragma once

namespace Repulsor
{
    
    template <int AMB_DIM, typename Real, typename Int, typename SReal>
    class CollisionFinder
    {
        static_assert(FloatQ<Real>,"");
        static_assert(FloatQ<SReal>,"");
        
        static_assert(IntQ<Int>,"");
        
    public:
        
        using GJK_T = GJK<AMB_DIM+1,GJK_Real,Int>;
        
        using       Primitive_T =       PolytopeBase<AMB_DIM,GJK_Real,Int,SReal>;
        using MovingPrimitive_T = MovingPolytopeBase<AMB_DIM,GJK_Real,Int,SReal>;

        CollisionFinder(
            cref<MovingPrimitive_T> P_,
            cref<MovingPrimitive_T> Q_,
            const SReal TOL
        )
        :   P{ P_.Clone() }
        ,   Q{ Q_.Clone() }
        ,   eps(TOL)
        {}
        
        // Default constructor
        CollisionFinder() = default;
        
        // Destructor
        virtual ~CollisionFinder() = default;
        
        // Copy constructor
        CollisionFinder( const CollisionFinder & other )
        :   P { other.P->Clone() }
        ,   Q { other.Q->Clone() }
        ,   eps(other.eps)
        {}
        
        // Copy assignment
        CollisionFinder & operator=( CollisionFinder other )
        {
            swap(*this,other);
            return *this;
        }
        
        // Move constructor
        CollisionFinder( CollisionFinder && other )
        {
            swap(*this,other);
            return *this;
        }
        
        // Move assignment no needed.
        
        
        // Swap function
        friend void swap( CollisionFinder & X, CollisionFinder & Y ) noexcept
        {
            // for details see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function
            using std::swap;
            
            swap( X.P   , Y.P   );
            swap( X.Q   , Y.Q   );
            swap( X.eps , Y.eps );
        }
      
    protected:

        mutable std::shared_ptr<MovingPrimitive_T> P;
        mutable std::shared_ptr<MovingPrimitive_T> Q;
        
        mutable GJK_T G;
        
        const SReal eps = static_cast<SReal>(0.0625);
        
        static constexpr Int max_iter = 126;
        
        mutable SReal b_stack[max_iter+2] = {};

    public:
 
        SReal RelativeTolerance() const
        {
            return eps;
        }
        
//        void SetRelativeTolerance( const SReal eps_)
//        {
//            eps = eps_;
//        }
        
        SReal FindMaximumSafeStepSize(
            cptr<SReal> p, cptr<SReal> u,
            cptr<SReal> q, cptr<SReal> v,
            const SReal tinit,
            const bool reuse_direction = true
        ) const
        {
            SReal a = Scalar::Zero<SReal>;
            SReal b = tinit;

            P->ReadCoordinatesSerialized(p);
            P->ReadVelocitiesSerialized(u);
            P->SetFirstTime(a);
            P->SetSecondTime(b);
            P->SetTimeScale(Scalar::One<SReal>);

            Q->ReadCoordinatesSerialized(q);
            Q->ReadVelocitiesSerialized(v);
            Q->SetFirstTime(a);
            Q->SetSecondTime(b);
            Q->SetTimeScale(Scalar::One<SReal>);

//            wprint(ClassName()+"::FindMaximumSafeStepSize: This might not work correctly for tinit !=1?!");
            
            SReal T = Sqrt(
                G.InteriorPoints_SquaredDistance(*P,*Q)
                + P->SquaredRadius()
                + Q->SquaredRadius()
//                - static_cast<Real>(0.75) * tinit * tinit
            )/tinit;

            P->SetTimeScale(T);
            Q->SetTimeScale(T);
            
            Int iter = 0;
            Int stack_ptr = 0;
            b_stack[0] = iter;
            

//            while( (b-a > eps * a) && (iter < max_iter) )
            while( (b-a > eps * b) && (iter < max_iter) )
            {
                ++iter;
                
                bool intersecting = G.IntersectingQ( *P, *Q, 1, reuse_direction && (iter>0) );
                
                if( !intersecting )
                {
                    if( stack_ptr > 0 )
                    {
                        a = b;
                        P->SetFirstTime(a);
                        Q->SetFirstTime(a);
                        
                        // pop
                        b = b_stack[stack_ptr--];
                        P->SetSecondTime(b);
                        Q->SetSecondTime(b);
                    }
                    else
                    {
                        return b;
                    }
                }
                else
                {
                    // push
                    b_stack[++stack_ptr] = b;
                    b = Scalar::Half<SReal> * (a + b);
                    
                    P->SetSecondTime(b);
                    Q->SetSecondTime(b);
                }
            }

            if( stack_ptr >= max_iter )
            {
                eprint(ClassName()+"::FindMaximumSafeStepSize: Stack overflow.");
                return static_cast<SReal>(a);
            }
            
            if( iter >= max_iter )
            {
//                wprint(ClassName()+"::FindMaximumSafeStepSize: max_iter = "+ToString(max_iter)+" reached." );
                wprint(ClassName()+"::FindMaximumSafeStepSize: iter >= max_iter");
                return static_cast<SReal>(a);
            }
            
            if( a<= Scalar::Zero<SReal> )
            {
                return static_cast<SReal>(b);
            }
            else
            {
                return static_cast<SReal>(a);
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return std::string("CollisionFinder") + "<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+">";
        }
    };
    
    
} // namespace Repulsor
