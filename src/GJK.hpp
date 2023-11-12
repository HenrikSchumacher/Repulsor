#pragma once

#include "Primitives/PrimitiveBase.hpp"
#include "Primitives/PrimitiveSerialized.hpp"

//#include "Primitives/ConvexHull.hpp"
//#include "Primitives/Ellipsoid.hpp"
#include "Primitives/MovingPolytopeBase.hpp"
#include "Primitives/MovingPolytopeExt.hpp"
#include "Primitives/MovingPolytope.hpp"

#include "Primitives/Point.hpp"
#include "Primitives/PolytopeBase.hpp"
#include "Primitives/PolytopeExt.hpp"
#include "Primitives/Polytope.hpp"

#include "Primitives/SpaceTimePrism.hpp"

#include "Primitives/SpaceTimePrism.hpp"

#include "BoundingVolumes/BoundingVolumeBase.hpp"
#include "BoundingVolumes/AABB.hpp"
#include "BoundingVolumes/AABB_LongestAxisSplit.hpp"
#include "BoundingVolumes/AABB_MedianSplit.hpp"
#include "BoundingVolumes/AABB_PreorderedSplit.hpp"

//#include "BoundingVolumes/OBB.hpp"
//#include "BoundingVolumes/OBB_MedianSplit.hpp"
//#include "BoundingVolumes/OBB_PreorderedSplit.hpp"

namespace Repulsor
{
    // some template magic
    template<typename Int>
    constexpr Int face_count(Int amb_dim)
    {
        Int n = 2;
        for( Int k = 0; k < amb_dim; ++k )
        {
            n *= 2;
        }
        return n;
    }
    
    enum class GJK_Reason
    {
        NoReason,
        FullSimplex,
        InSimplex,
        MaxIteration,
        SmallProgress,
        SmallResidual,
        CollisionTolerance,
        Collision,
        Separated
    };
    
    template<int AMB_DIM, typename Real_, typename Int_>
    class alignas(ObjectAlignment) GJK
    {
        ASSERT_FLOAT(Real_);
        ASSERT_INT  (Int_ );

    public:
        
        using Int  = Int_;
        using Real = Real_;
        
        using Vector_T = Tiny::Vector<AMB_DIM  ,Real,Int>;
        using Lambda_T = Tiny::Vector<AMB_DIM+1,Real,Int>;
        
        using IntVec_T = Tiny::Vector<AMB_DIM+1,Int,Int>;
        
        using PrimitiveBase_T = PrimitiveBase<AMB_DIM,Real,Int>;
        using Op = Tensors::Op;
        using Flag = Scalar::Flag;
        
        static constexpr Real eps = cSqrt(std::numeric_limits<Real>::epsilon());
        static constexpr Real eps_squared = eps * eps;
        static constexpr Int max_iter = 100;
        
    protected:
        
        static constexpr Int FACE_COUNT = face_count(AMB_DIM);
        static constexpr Real zero      = Scalar::Zero<Real>;
        static constexpr Real one       = Scalar::One <Real>;
        
        std::array<Vector_T,AMB_DIM+1> coords = {{}};  //  position of the corners of the simplex; only the first simplex_size rows are defined.
        std::array<Vector_T,AMB_DIM+1> P_supp = {{}};  //  support points of the simplex in primitive P
        std::array<Vector_T,AMB_DIM+1> Q_supp = {{}};  //  support points of the simplex in primitive Q

        Tiny::Matrix<AMB_DIM+1,AMB_DIM+1,Real,Int> dots;  // simplex_size x simplex_size matrix of dots products of  the vectors coords[0],..., coords[simplex_size-1];
        Tiny::Matrix<AMB_DIM  ,AMB_DIM  ,Real,Int> G   ;  // Gram matrix of size = (simplex_size-1) x (simplex_size-1) of frame spanned by the vectors coords[i] - coords[simplex_size-1];
        
        Vector_T v;  // current direction vector (quite hot)
        
        Lambda_T Lambda;        // For the right hand sides of the linear equations to solve in DistanceSubalgorithm.
        Lambda_T best_lambda;
        
        std::array<Int     ,FACE_COUNT> facet_sizes;
        std::array<IntVec_T,FACE_COUNT> facet_vertices;
        std::array<IntVec_T,FACE_COUNT> facet_faces;
        
        std::array<bool    ,FACE_COUNT> visited = {};
        
        
        Real vv            = Scalar::Max<Real>;
        Real oldvv         = Scalar::Max<Real>;
        Real vw            = zero;
        Real TOL_squared   = zero;
        Real theta_squared = one;
        
        Int simplex_size = 0; //  index of last added point
        Int closest_facet;
        Int sub_calls = 0;
        GJK_Reason reason = GJK_Reason::NoReason;
        bool separatedQ = false;

    public:
        
        GJK()
        {
            // Initializing facet_sizes, facet_vertices, and facet_faces.
            // TODO: In principle, the compiler should be able to populate those at compile time. Dunno how to work this black magic.

            for( Int facet = 0; facet < FACE_COUNT; ++ facet )
            {
                Int i = 0;
                
                mref<IntVec_T> vertices = facet_vertices[facet];
                mref<IntVec_T> faces    = facet_faces   [facet];

                for( Int vertex = 0; vertex < AMB_DIM+1; ++vertex )
                {
                    if( get_bit(facet,vertex) )
                    {
                        vertices[i] = vertex;
                        faces[i] = set_bit_to_zero(facet,vertex);
                        ++i;
                    }
                }
                
                facet_sizes[facet] = i;
                
                for( Int j = i; j < AMB_DIM+1; ++j )
                {
                    vertices[i] = -1;
                    faces[i] = -1;
                }
            }
            
            visited[0] = true;
        }

        GJK( const GJK & other ) : GJK() {};
        
        GJK( GJK && other ) : GJK() {};
        
        ~GJK() = default;
        
        static constexpr Int AmbDim()
        {
            return AMB_DIM;
        }
        
        Int SimplexSize() const
        {
            return simplex_size;
        }
        
        Real LeastSquaredDistance() const
        {
            return vv;
        }
        
        bool SeparatedQ() const
        {
            return separatedQ;
        }
        
        template<typename ExtReal>
        void WriteClosestPoint( mptr<ExtReal> vec ) const
        {
            v.Write(vec);
        }

        Int SubCallCount() const
        {
            return sub_calls;
        }

    protected:
        
        
        static Int set_bit_to_zero( const Int n, const Int k )
        {
            return (n ^ (1 << k));
        }
        
        void Compute_dots()
        {
            // update matrix of dot products

            for( Int i = 0; i < simplex_size+1; ++i )
            {
                dots[i][simplex_size] = Dot( coords[i], coords[simplex_size] );
            }
        }
        
        int Compute_Gram()
        {
            // Computes Gram matrix G[i][i] = <coords[i] - coords[0], coords[j] - coords[0]>.
            // However, we do that in a convoluted way to save a few flops.
            
            for( Int i = 0; i < simplex_size; ++i )
            {
                const Real R1 = dots[i][simplex_size];
                const Real R2 = dots[simplex_size][simplex_size] - R1;

                Lambda[i] = R2;

                G[i][i] = dots[i][i] + R2 - R1;

                // This is to guarantee that the current facet has full rank.
                // If G is not of full rank, then most recently added point (w) was already contained in the old simplex.
                // (If another point were a problem, then the code below would have aborted already the previous GJK iteration.
                if( G[i][i] <= 0 )
                {
                    return 1;
                }

                for( Int j = i+1; j < simplex_size; ++j )
                {
                    G[i][j] = dots[i][j] - dots[j][simplex_size] + R2;
                }
            }
            return 0;
        }
        
        int Push()
        {
            // Pushes most recent difference of support points (i.e. w = P_supp[simplex_size] - Q_supp[simplex_size] )
            // into the simplex coords.

            LinearCombine<Flag::Plus,Flag::Minus>(
                one, P_supp[simplex_size], -one, Q_supp[simplex_size], coords[simplex_size]
            );
            
            Compute_dots();
            
            int stat = Compute_Gram();
            
            ++simplex_size;
            
            return stat;
        }

        Int PrepareDistanceSubalgorithm()
        {
            // Compute starting facet.
            const Int facet = (static_cast<Int>(1) << simplex_size) - static_cast<Int>(1) ;
            
            closest_facet = facet;
            oldvv = vv;
            
            vv = Scalar::Max<Real>;
            
            // Mark subsimplices of `facet` that do not contain the last vertex of the simplex as visited before starting.
            // Facet itself is marked as unvisited; if `facet` is a reasonable starting facet, it _must_ be visited.
            visited[facet] = false;
            
            const Int top_index = simplex_size - 1;
            // Subsimplices of `facet` have _lower_ index than `facet`.

            for( Int i = 0; i < facet; ++i )
            {
                visited[i] = !get_bit(i,top_index);
            }
            
            return facet;
        }
        
        void HandlePoints( const PrimitiveBase_T & P, const PrimitiveBase_T & Q )
        {
            // In the case that both primitives are points, we have to take care that witnesses are computed correctly.
            simplex_size = 1;
            best_lambda[1] = 1;
            
            P.InteriorPoint( &P_supp[0][0] );
            Q.InteriorPoint( &Q_supp[0][0] );

            
            LinearCombine<Flag::Plus,Flag::Minus>( one, P_supp[0], -one, Q_supp[0], v );
            
            vv = v.SquaredNorm();
            
            separatedQ = vv > zero;
        }
        
        // ################################################################
        // ###################   DistanceSubalgorithm   ###################
        // ################################################################
        
        template<Int FACET_SIZE>
        void DistanceSubalgorithm( const Int facet, const Int facet_size )
        {
            if( FACET_SIZE < facet_size )
            {
                if constexpr ( (1 <= FACET_SIZE) && (FACET_SIZE+1 <= AMB_DIM) )
                {
                    DistanceSubalgorithm<FACET_SIZE+1>( facet, facet_size );
                }
                return;
            }
            
            ++sub_calls;
            
            cref<IntVec_T> vertices = facet_vertices[facet];
            cref<IntVec_T> faces    = facet_faces   [facet];
            
            visited[facet] = true;
            
            bool interiorQ = true;
            
            // The local contiguous version of Lambda for facets of size > 1.
            Lambda_T lambda = {};
            
            // Compute lambdas.
            if constexpr ( FACET_SIZE  == 1 )
            {
                const Int i_0 = vertices[0];
                
                if( dots[i_0][i_0] < vv )
                {
                    closest_facet     = facet;
                    vv             = dots[i_0][i_0];
                    best_lambda[0]    = one;

                    v = coords[i_0];
                }
                return;
            }
            else if constexpr ( FACET_SIZE == 2 )
            {
                // Setting up linear system for the barycenter coordinates lambda.
                
                // Find first vertex in facet.
                const Int i_0 = vertices[0];

                lambda[0] = Lambda[i_0] / G[i_0][i_0];
                lambda[1] = one - lambda[0];

                interiorQ = (lambda[0] > eps) && (lambda[1] > eps);
            }
            else if constexpr ( FACET_SIZE == 3 )
            {
                // Setting up linear system for the barycenter coordinates lambda.
                
                // Find first two vertices in facet.
                const Int i_0 = vertices[0];
                const Int i_1 = vertices[1];

                // Using Cramer's rule to solve the linear system.
                const Real inv_det = one / ( G[i_0][i_0] * G[i_1][i_1] - G[i_0][i_1] * G[i_0][i_1] );

                lambda[0] = ( G[i_1][i_1] * Lambda[i_0] - G[i_0][i_1] * Lambda[i_1] ) * inv_det;
                lambda[1] = ( G[i_0][i_0] * Lambda[i_1] - G[i_0][i_1] * Lambda[i_0] ) * inv_det;
                lambda[2] = one - lambda[0] - lambda[1];

                // Check the barycentric coordinates for positivity ( and compute the 0-th coordinate).
                interiorQ = (lambda[0] > eps) && (lambda[1] > eps) && (lambda[2] > eps);
            }
            else if constexpr ( (3 < FACET_SIZE) && (FACET_SIZE <= AMB_DIM) )
            {
                constexpr Int n = FACET_SIZE - 1;
                
                Tiny::SelfAdjointMatrix<n,Real,Int> A;
                
                Tiny::Vector<n,Real,Int> mu;
                
                for( Int i = 0; i < n; ++i )
                {
                    const Int i_i = vertices[i];
                    
                    mu[i] = Lambda[i_i];
                    
                    A[i][i] = G[i_i][i_i];
                    
                    for( Int j = i+1; j < n; ++j )
                    {
                        const Int j_j = vertices[j];
                        
                        A[i][j] = G[i_i][j_j];
                    }
                }

                // Linear solve.
                
                A.Cholesky();
                
                A.CholeskySolve(mu);
                
                // Check the barycentric coordinates for positivity ( and compute the 0-th coordinate).
                
                Real mu_n = one;
                
                for( Int k = 0; k < n; ++k )
                {
                    mu_n -= lambda[k];
                    interiorQ = interiorQ && ( mu[k] > eps );
                }
                
                mu.Write( lambda.data() );
                lambda[n] = mu_n;
                
                interiorQ = interiorQ && (mu_n > eps );
            }
            
            // If we arrive here, then FACET_SIZE > 1.
            if( interiorQ )
            {
                // If the nearest point on `facet` lies in the interior, there is no point in searching its faces;
                // we simply mark them as visited.
                
                for( Int j = 0; j < FACET_SIZE; ++j )
                {
                    visited[ faces[j] ] = true;
                }
                
                // Computing closest point z on facet.
                
                Vector_T z = {};
                
                for( Int j = 0; j < facet_size; ++j )
                {
                    axpy( lambda[j], coords[ vertices[j] ], z );
                }
                
                const Real zz = z.SquaredNorm();
                
                if( zz < vv )
                {
                    closest_facet = facet;
                    vv         = zz;
                    v             = z;
                    best_lambda   = lambda;
                }
            }
            else
            {
                // We have to visit all faces of `facet` and lock for the minimizer there.
                for( Int j = 0; j < FACET_SIZE; ++j )
                {
                    const Int face = faces[j];
                    
                    if( !visited[face] )
                    {
                        if constexpr ( FACET_SIZE - 1 >= 1 )
                        {
                            if( (FACET_SIZE - 1) != facet_sizes[face] )
                            {
                                eprint("!!!");
                            }
                            
                            DistanceSubalgorithm<FACET_SIZE - 1>( face, FACET_SIZE - 1 );
                            
                        }
                    }
                }
            }
        }
        
    public:
    
        // ################################################################
        // ##########################  Compute  ###########################
        // ################################################################
        
        
        // If collision_only is set to true, then we abort if at any time, we found a point x in P and a point y in  Q such that:
        //
        //   - if TOL_squared_ > 0: theta_squared_ * |x-y|^2 <= TOL_squared_
        //
        //   - if TOL_squared_ <=0:
        //
        //          - if P.SquaredRadius() and Q.SquaredRadius() are positive:
        //
        //                 theta_squared_ * |x-y|^2 <= eps_squared * min( P.SquaredRadius(), Q.SquaredRadius() )
        //
        //          - if any of P.SquaredRadius(), Q.SquaredRadius() is nonpositive:
        //
        //                 theta_squared_ * |x-y|^2 <= eps_squared * max( P.SquaredRadius(), Q.SquaredRadius() )
        //
        // For intersections of thickened objects use
        //
        //    theta_squared = 1 and
        //
        //    TOL_squared_  = (r_P + r_Q)^2, where r_P and r_Q are the thickening radii of P and Q.
        //
        // For multipole acceptance criteria use
        //
        //    theta_squared = chi * chi, where chi >= 0 is the separation parameter and
        //
        //    TOL_squared_  = max(r_P^2,r_Q^2), where r_P and r_Q are the radii of P and Q and
        
        
        //Compute(P, Q, true, reuse_direction_, zero, theta_squared_ );
        void Compute(
            const PrimitiveBase_T & P,
            const PrimitiveBase_T & Q,
            const bool collision_only  = false,
            const bool reuse_direction = false,
            const Real TOL_squared_    = zero,
            const Real theta_squared_  = one
        )
        {
            separatedQ = false;
            
            Int iter = 0;

            int in_simplex;

            theta_squared = theta_squared_;

            if( TOL_squared_ > zero )
            {
                TOL_squared = TOL_squared_;
            }
            else
            {
                // Using a tolerance based on the radii of the primitives.
                
                TOL_squared = eps_squared * Min( P.SquaredRadius(), Q.SquaredRadius() );
                if( TOL_squared <= zero )
                {
                    // One of the primitives must be a point. Use the radius of the other one as tolerance.
                    
                    TOL_squared = eps_squared * Max( P.SquaredRadius(), Q.SquaredRadius() );
                    
                    if( TOL_squared <= zero )
                    {
                        // Both primitives are points.
                        
                        HandlePoints(P,Q);
                        
                        return;
                    }
                }
            }

            sub_calls = 0;
            reason = GJK_Reason::NoReason;
            simplex_size = 0;
            
            if( ! reuse_direction )
            {
                P.InteriorPoint(&v[0]);
                Q.InteriorPoint(&Q_supp[0][0]);
                
                v -= Q_supp[0];
            }
            
            oldvv = vv = v.SquaredNorm();

            // Unrolling the first iteration to avoid a call to DistanceSubalgorithm.
            
            // We use w = p-q, but do not define it explicitly.
            Real a = P.MinSupportVector( &v[0], &P_supp[0][0] );
            Real b = Q.MaxSupportVector( &v[0], &Q_supp[0][0] );
            vw  = a-b;
            
            Push();
            
            closest_facet = 1;
            vv = dots[0][0];
            best_lambda[0] = one;

            v = coords[0];
        
            while( true )
            {
                if( theta_squared * vv < TOL_squared )
                {
                    reason = GJK_Reason::CollisionTolerance;
                    break;
                }
                
                if( simplex_size >= AMB_DIM + 1 )
                {
                    reason = GJK_Reason::FullSimplex;
                    break;
                }
                
                if( iter >= max_iter)
                {
                    reason = GJK_Reason::MaxIteration;
                    break;
                }
                
                ++iter;
                
                // We use w = p-q, but do not define it explicitly.
                a = P.MinSupportVector( &v[0], &P_supp[simplex_size][0] );
                b = Q.MaxSupportVector( &v[0], &Q_supp[simplex_size][0] );
                vw = a-b;
            
                if( collision_only && (vw > zero) && (theta_squared * vw * vw > TOL_squared) )
                {
                    separatedQ = true;
                    reason = GJK_Reason::Separated;
                    break;
                }
                
                
                if( Abs(vv - vw) <= eps * vv )
                {
                    reason = GJK_Reason::SmallResidual;
                    break;
                }
                
                in_simplex = Push();
                
                if( in_simplex  )
                {
                    // The vertex added most recently to the simplex coincides with one of the previous vertices.
                    // So we stop here and use the old best_lambda. All we have to do is to reduce the simplex_size by 1.
                    --simplex_size;
                    reason = GJK_Reason::InSimplex;
                    break;
                }
                
                Int initial_facet = PrepareDistanceSubalgorithm();
                
                // This computes closest_facet, v, and vv of the current simplex.
                // If the simplex is degenerate, DistanceSubalgorithm terminates early and returns 1. Otherwise it returns 0.
                DistanceSubalgorithm<1>( initial_facet, facet_sizes[initial_facet] );
                
                simplex_size = facet_sizes[closest_facet];
                
                cref<IntVec_T> vertices = facet_vertices[closest_facet];
                
                // Deleting superfluous vertices in simplex and writing everything to the beginning of the array.
                
                for( Int i = 0; i < simplex_size; ++i )
                {
                    const Int i_i = vertices[i];
                    
                    coords[i] = coords[i_i];
                    Q_supp[i] = Q_supp[i_i];

                    for( Int j = i; j < simplex_size; ++j )
                    {
                        dots[i][j] = dots[i_i][vertices[j]];
                    }
                }

                if( Abs(oldvv - vv) <= eps * vv )
                {
                    reason = GJK_Reason::SmallProgress;
                    break;
                }
                
            } // while( true )

            separatedQ = separatedQ || ( TOL_squared < theta_squared * vv );
            
            if( iter >= max_iter)
            {
//                wprint(ClassName()+"::Compute: Stopped because iter = " + ToString(iter) + " >= " + ToString(max_iter) + " = max_iter iterations reached.");
            }
            
        } // Compute

        // ################################################################
        // #######################   IntersectingQ   ######################
        // ################################################################
        
        bool IntersectingQ(
            const PrimitiveBase_T & P,
            const PrimitiveBase_T & Q,
            const Real theta_squared_ = one,
            const bool reuse_direction_ = false
        )
        {
            // If both P.SquaredRadius() and Q.SquaredRadius() are positive, this routines checks
            // whether there are points x in P and y in Q such that
            //
            //     theta_squared_ * |x-y|^2 <= eps_squared * Min( P.SquaredRadius(), Q.SquaredRadius() );
            //
            // Otherwise it checks whether
            //
            //     theta_squared_ * |x-y|^2 <= eps_squared * Min( P.SquaredRadius(), Q.SquaredRadius() );
            //
            
            Compute(P, Q, true, reuse_direction_, zero, theta_squared_ );
            
            return !separatedQ;
        }
        
//        bool IntersectingQ(
//            const PrimitiveBase_T & P,
//            const PrimitiveBase_T & Q,
//            const Real theta_squared_
//        )
//        {
//            return IntersectingQ(P,Q,theta_squared_,false);
//        }
        
//        bool IntersectingQ(
//            const PrimitiveBase_T & P,
//            const PrimitiveBase_T & Q
//        )
//        {
//            return IntersectingQ(P,Q,one,false);
//        }
        
        // Faster overloads for AABBs.
        template<typename SReal>
        bool IntersectingQ(
            const AABB<AMB_DIM,Real,Int,SReal> & P,
            const AABB<AMB_DIM,Real,Int,SReal> & Q,
            const Real theta_squared_ = one
        )
        {
            return AABB_SquaredDistance(P,Q) <= zero;
        }
        
        // ################################################################
        // #####################   SquaredDistance   #####################
        // ################################################################
        
        Real SquaredDistance(
            const PrimitiveBase_T & P,
            const PrimitiveBase_T & Q,
            const bool reuse_direction_ = false
        )
        {
            Compute(P, Q, false, reuse_direction_, zero );
            
            return vv;
        }
        
        
        // Faster overload for AABBs.
        template<typename SReal>
        Real SquaredDistance(
            const AABB<AMB_DIM,Real,Int,SReal> & P,
            const AABB<AMB_DIM,Real,Int,SReal> & Q
        )
        {
            return AABB_SquaredDistance(P,Q);
        }
        
        // ##########################################################################
        // ##############################   Witnesses   #############################
        // ##########################################################################
        
        Real Witnesses(
            const PrimitiveBase_T & P, Real * restrict const x,
            const PrimitiveBase_T & Q, Real * restrict const y,
            const bool reuse_direction_ = false
        )
        {
            // x and y are the return parameters.
            // On return x lies in P while y lies in Q.
            // These points realize the minimal distance between any two points in these primitives.
            // Scalar return values is this minimal distance _squared_(!).
            
            Compute(P, Q, false, reuse_direction_, zero );

            
            // TODO: Unify?
            y  = Q_supp[0];
            y *= best_lambda[0];
            
            switch( simplex_size )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    axpy( best_lambda[1], Q_supp[1], y );
                    break;
                }
                case 3:
                {
                    for( Int j = 1; j < 3; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                    break;
                }
                case 4:
                {
                    for( Int j = 1; j < 4; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                    break;
                }
                case 5:
                {
                    for( Int j = 1; j < 5; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                    break;
                }
                default:
                {
                    for( Int j = 1; j < simplex_size; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                }
            }
            
            x  = y;
            x += v;
            
            return vv;
        } // Witnesses
        
        // ################################################################
        // ####################   Offset_IntersectingQ   ##################
        // ################################################################
        
        bool Offset_IntersectingQ(
            const PrimitiveBase_T & P, const Real P_offset,
            const PrimitiveBase_T & Q, const Real Q_offset,
            const bool reuse_direction_ = false
        )
        {
            const Real min_dist = P_offset + Q_offset;
            
            Compute(P, Q, true, reuse_direction_, min_dist * min_dist );
            
            return !separatedQ;
        }
        
        // ################################################################
        // #################   Offset_SquaredDistance   ##################
        // ################################################################
        
        Real Offset_SquaredDistance(
            const PrimitiveBase_T & P, const Real P_offset,
            const PrimitiveBase_T & Q, const Real Q_offset,
            const bool reuse_direction_ = false
        )
        {
            const Real min_dist = P_offset + Q_offset;
            
            Compute(P, Q, false, reuse_direction_, min_dist * min_dist );
            
            const Real dist = Ramp( Sqrt(vv) - P_offset - Q_offset );
            
            return dist * dist;
        }
        
        // ##########################################################################
        // ###########################   Offset_Witnesses   #########################
        // ##########################################################################
        
        Real Offset_Witnesses(
            const PrimitiveBase_T & P, const Real P_offset, Real * restrict const x,
            const PrimitiveBase_T & Q, const Real Q_offset, Real * restrict const y,
            const bool reuse_direction_ = false
        )
        {
            // x and y are the return variables.
            // On return,
            // x lies in the offset of P with offset P_offset and
            // y lies in the offset of Q with offset Q_offset.
            // The realize the minimal distance between two points in these offsets.
            // Scalar return values is this minimal distance _squared_(!).
            
            const Real min_dist = P_offset + Q_offset;
            
            Compute(P, Q, false, reuse_direction_, min_dist * min_dist );
            
            const Real dist0 = Sqrt(vv);
                  Real dist = dist0 - P_offset - Q_offset;
                  Real x_scale;
                  Real y_scale;
            
            if( dist > zero )
            {
                // If no intersection was detected, we have to set the witnesses onto the bounday of the thickened primitives.
                
                // We want y = y_0 + y_scale * v, where y_0 is constructed from Q_supp by barycentric coordinates
                y_scale =  Q_offset / dist0;
                // We want x = y + x_scale * v.
                x_scale = (dist0 - P_offset - Q_offset) / dist0;
            }
            else
            {
                // If an intersection was deteced, we just return the witnesses.
                dist = zero;
                // We want y = y_0 + y_scale * v, where y_0 is constructed from Q_supp by barycentric coordinates
                y_scale = zero;
                // We want x = y + x_scale * v.
                x_scale = one;
            }

            // Compute y = best_lambda * Q_supp;
            
            LinearCombine( y_scale, v, best_lambda[0], Q_supp[0], y );
            
            switch( simplex_size )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    axpy( best_lambda[1], Q_supp[1], y );
                    break;
                }
                case 3:
                {
                    for( Int j = 1; j < 3; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                    break;
                }
                case 4:
                {
                    for( Int j = 1; j < 4; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                    break;
                }
                case 5:
                {
                    for( Int j = 1; j < 5; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                    break;
                }
                default:
                {
                    for( Int j = 1; j < simplex_size; ++j )
                    {
                        axpy( best_lambda[j], Q_supp[j], y );
                    }
                }
            }
            

            LinearCombine<Flag::Plus,Flag::Generic>( one, y, x_scale, v, x );
            
            return dist * dist;
            
        } // Offset_Witnesses

        // ########################################################################
        // ##################   InteriorPoints_SquaredDistance   ##################
        // ########################################################################
        
        Real InteriorPoints_SquaredDistance(
            const PrimitiveBase_T & P,
            const PrimitiveBase_T & Q
        )
        {
            P.InteriorPoint( &coords[0][0] );
            Q.InteriorPoint( &coords[1][0] );
            
            Real r2 = 0;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                const Real diff = coords[1][k] - coords[0][k];
                r2 += diff * diff;
            }
            return r2;
            
        }
        
        
        // ########################################################################
        // ####################   MultipoleAcceptanceCriterion   ##################
        // ########################################################################
        
        bool MultipoleAcceptanceCriterion(
            const PrimitiveBase_T & P,
            const PrimitiveBase_T & Q,
            const Real theta_squared_,
            const bool reuse_direction_ = false
        )
        {
            Compute(P, Q, true, reuse_direction_, Max( P.SquaredRadius(), Q.SquaredRadius() ), theta_squared_ );
            
            return separatedQ;
        }
        
        // Faster overload for AABBs.
        template<typename SReal>
        bool MultipoleAcceptanceCriterion(
            const AABB<AMB_DIM,Real,Int,SReal> & P,
            const AABB<AMB_DIM,Real,Int,SReal> & Q,
            const Real theta_squared_
        )
        {
            return Max( P.SquaredRadius(), Q.SquaredRadius() ) < theta_squared_ * AABB_SquaredDistance( P, Q );
        }
        
        std::string ClassName() const
        {
            return std::string("GJK")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Real>+">";
        }
        
    }; // GJK

} // namespace Repulsor



