#pragma once

#include "../submodules/Tensors/Tensors.hpp"

#include "Primitives/PrimitiveBase.hpp"
#include "Primitives/PrimitiveSerialized.hpp"

#include "Primitives/ConvexHull.hpp"
#include "Primitives/Ellipsoid.hpp"
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
    using namespace Tensors;
    using namespace Tools;
    
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
    
    // TODO: This works well for AMB_DIM = 2 and AMB_DIM = 3. But AMB_DIM >= 4 seems to be problematic.
    
    template<int AMB_DIM, typename Real_, typename Int_>
    class alignas(ObjectAlignment) GJK
    {
        static_assert(FloatQ<Real_>,"");
        static_assert(IntQ<Int_>,"");
        
    public:
        
        using Int  = Int_;
        using Real = Real_;
        
        using PrimitiveBase_T = PrimitiveBase<AMB_DIM,Real,Int>;
        using Op = Tensors::Op;
        
        static constexpr Real eps = cSqrt(std::numeric_limits<Real>::epsilon());
        static constexpr Real eps_squared = eps * eps;
        static constexpr Int max_iter = 100;
        
    protected:
        
        static constexpr Int FACE_COUNT = face_count(AMB_DIM);
        static constexpr Real zero      = Scalar::Zero<Real>;
        static constexpr Real one       = Scalar::One <Real>;
        
        Real coords   [AMB_DIM+1][AMB_DIM  ] = {};  //  position of the corners of the simplex; only the first simplex_size rows are defined.
        Real P_supp   [AMB_DIM+1][AMB_DIM  ] = {};  //  support points of the simplex in primitive P
        Real Q_supp   [AMB_DIM+1][AMB_DIM  ] = {};  //  support points of the simplex in primitive Q

        Real dots     [AMB_DIM+1][AMB_DIM+1] = {};  // simplex_size x simplex_size matrix of dots products of  the vectors coords[0],..., coords[simplex_size-1];
        Real Gram     [AMB_DIM  ][AMB_DIM  ] = {};  // Gram matrix of size = (simplex_size-1) x (simplex_size-1) of frame spanned by the vectors coords[i] - coords[simplex_size-1];
        
        Real g        [AMB_DIM  ][AMB_DIM  ] = {}; // The local contiguous gram matrix for facets of size > 3.
        Real v                   [AMB_DIM  ] = {}; // current direction vector (quite hot)
        Real p                   [AMB_DIM  ] = {}; // just a buffer for storing current support point of p.
        Real Lambda              [AMB_DIM+1] = {}; // For the right hand sides of the linear equations to solve in DistanceSubalgorithm.
        Real best_lambda         [AMB_DIM+1] = {};
        Real facet_closest_point [AMB_DIM  ] = {};
        
        Int facet_sizes    [FACE_COUNT] = {};
        Int facet_vertices [FACE_COUNT][AMB_DIM+1] = {};
        Int facet_faces    [FACE_COUNT][AMB_DIM+1] = {};
        bool visited       [FACE_COUNT] = {};
        
        Real dotvv         = std::numeric_limits<Real>::max();
        Real olddotvv      = std::numeric_limits<Real>::max();
        Real dotvw         = zero;
        Real TOL_squared   = zero;
        Real theta_squared = one;
        
        Int simplex_size = 0; //  index of last added point
        Int closest_facet;
        Int sub_calls = 0;
        GJK_Reason reason = GJK_Reason::NoReason;
        bool separatedQ = false;

    public:
        
        // No default constructor
        GJK()
        {
            // Initializing facet_sizes, facet_vertices, and facet_faces.
            // TODO: In principle, the compiler should be able to populate those at compile time. Dunno how to work this black magic.

            for( Int facet = 0; facet < FACE_COUNT; ++ facet )
            {
                Int i = 0;
                
                mptr<Int> vertices = &facet_vertices[facet][0];
                mptr<Int> faces    = &facet_faces   [facet][0];

                for( Int vertex = 0; vertex < AMB_DIM+1; ++vertex )
                {
                    if( bit(facet,vertex) )
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

        ~GJK() = default;
        // Copy constructor
        GJK( const GJK & other ) = default;
        // Copy assignment operator
        GJK & operator=( const GJK & other ) = default;
        // Move constructor
        GJK( GJK && other ) = default;
        // Move assignment operator
        GJK & operator=( GJK && other ) = default;
        
        constexpr Int AmbDim() const
        {
            return AMB_DIM;
        }
        
        Int Size() const
        {
            return simplex_size;
        }
        
        Real LeastSquaredDistance() const
        {
            return dotvv;
        }
        
        bool SeparatedQ() const
        {
            return separatedQ;
        }
        
        template<typename ExtReal>
        void WriteClosestPoint( mptr<ExtReal> vec ) const
        {
            copy_buffer<AMB_DIM>( &v[0], vec );
        }

        Int SubCallCount() const
        {
            return sub_calls;
        }

    protected:
        
        static bool bit( const Int n, const Int k )
        {
            return ( (n & (1 << k)) >> k );
        }
        
        static  Int set_bit_to_zero( const Int n, const Int k )
        {
            return (n ^ (1 << k));
        }
        
        int Compute_Gram()
        {
            // Computes Gram matrix Gram[i][i] = <coords[i] - coords[0], coords[j] - coords[0]>.
            // However, we do that in a convoluted way to save a few flops.
            
            for( Int i = 0; i < simplex_size; ++i )
            {
                const Real R1 = dots[i][simplex_size];
                const Real R2 = dots[simplex_size][simplex_size] - R1;

                Lambda[i] = R2;

                Gram[i][i] = dots[i][i] + R2 - R1;

                // This is to guarantee that the current facet has full rank.
                // If g is not of full rank, then most recently added point (w) was already contained in the old simplex.
                // (If another point were a problem, then the code below would have aborted already the previous GJK iteration.
                if( Gram[i][i] <= 0 )
                {
                    return 1;
                }

                for( Int j = i+1; j < simplex_size; ++j )
                {
                    Gram[i][j] = dots[i][j] - dots[j][simplex_size] + R2;
                }
            }

            return 0;
        }
        
        
        template<Int SIMPLEX_SIZE>
        int Push()
        {
            if constexpr ( SIMPLEX_SIZE < AMB_DIM )
            {
                if( SIMPLEX_SIZE < simplex_size )
                {
                    return Push<SIMPLEX_SIZE+1>();
                }
            }
            
            // Pushes most recent difference of support points (i.e. w = P_supp[simplex_size] - Q_supp[simplex_size] )
            // into the simplex coords.

            for( Int k = 0; k < AMB_DIM; ++k)
            {
                coords[SIMPLEX_SIZE][k] = P_supp[SIMPLEX_SIZE][k] - Q_supp[SIMPLEX_SIZE][k];
            }

            // update matrix of dot products

            for( Int i = 0; i < SIMPLEX_SIZE+1; ++i )
            {
                dots[i][SIMPLEX_SIZE] = coords[i][0] * coords[SIMPLEX_SIZE][0];

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    dots[i][SIMPLEX_SIZE] += coords[i][k] * coords[SIMPLEX_SIZE][k];
                }
            }
            
            for( Int i = 0; i < SIMPLEX_SIZE; ++i )
            {
                const Real R1 = dots[i][SIMPLEX_SIZE];
                const Real R2 = dots[SIMPLEX_SIZE][SIMPLEX_SIZE] - R1;

                Lambda[i] = R2;

                Gram[i][i] = dots[i][i] + R2 - R1;

                // This is to guarantee that the current facet has full rank.
                // If g is not of full rank, then most recently added point (w) was already contained in the old simplex.
                // (If another point were a problem, then the code below would have aborted already the previous GJK iteration.
                if( Gram[i][i] <= 0 )
                {
                    ++simplex_size;
                    return 1;
                }

                for( Int j = i+1; j < SIMPLEX_SIZE; ++j )
                {
                    Gram[i][j] = dots[i][j] - dots[j][SIMPLEX_SIZE] + R2;
                }
            }
            
            ++simplex_size;
            
            return 0;
        }

        
        template<Int SIMPLEX_SIZE = 1>
        void Deflate()
        {
            if constexpr ( SIMPLEX_SIZE < AMB_DIM )
            {
                if( SIMPLEX_SIZE < simplex_size )
                {
                    return Deflate<SIMPLEX_SIZE+1>();
                }
            }
            
            cptr<Int> vertices = facet_vertices[closest_facet];
            
            // Deleting superfluous vertices in simplex and writing everything to the beginning of the array.
            
            for( Int i = 0; i < SIMPLEX_SIZE; ++i )
            {
                const Int i_i = vertices[i];
                
                copy_buffer<AMB_DIM>( &coords[i_i][0], &coords[i][0] );
                copy_buffer<AMB_DIM>( &Q_supp[i_i][0], &Q_supp[i][0] );

                for( Int j = i; j < SIMPLEX_SIZE; ++j )
                {
                    dots[i][j] = dots[i_i][vertices[j]];
                }
            }
        }
        Int PrepareDistanceSubalgorithm()
        {
            // Compute starting facet.
            Int facet = (static_cast<Int>(1) << simplex_size) - static_cast<Int>(1) ;
            
            closest_facet = facet;
            olddotvv = dotvv;
            
            dotvv = std::numeric_limits<Real>::max();
            
            // Mark subsimplices of `facet` that do not contain the last vertex of the simplex as visited before starting.
            // Facet itself is marked as unvisited; if `facet` is a reasonable starting facet, it _must_ be visited.
            visited[facet] = false;
            
            const Int top_index = simplex_size - 1;
            // Subsimplices of `facet` have _lower_ index than `facet`.

            for( Int i = 0; i < facet; ++i )
            {
                visited[i] = !bit(i,top_index);
            }
            return facet;
        }
        
        void HandlePoints(
            cref<PrimitiveBase_T> P,
            cref<PrimitiveBase_T> Q
        )
        {
            // In the case that both primitives are points, we have to take care that witnesses are computed correctly.
            simplex_size = 1;
            best_lambda[1] = 1;
            
            P.InteriorPoint(&P_supp[0][0]);
            Q.InteriorPoint(&Q_supp[0][0]);

            dotvv = 0;
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                v[k] = P_supp[0][k] - Q_supp[0][k];
                dotvv += v[k] * v[k];
            }
            
            separatedQ = dotvv > zero;
        }
        
/// ################################################################
/// ##                    DistanceSubalgorithm                    ##
/// ################################################################
        
        template<Int FACET_SIZE>
        void DistanceSubalgorithm( const Int facet, const Int facet_size )
        {
            if constexpr ( FACET_SIZE <= AMB_DIM )
            {
                if( FACET_SIZE < facet_size )
                {
                    DistanceSubalgorithm<FACET_SIZE+1>( facet, facet_size );
                    return;
                }
            }
            
            ++sub_calls;
            
//            const Int facet_size = facet_sizes[facet];
            cptr<Int> vertices = &facet_vertices[facet][0];
            cptr<Int> faces    = &facet_faces   [facet][0];
            
            Real lambda [AMB_DIM+1] = {}; // The local contiguous version of Lambda for facets of size > 3.
            
            visited[facet] = true;
            
            bool interior = true;
            
            constexpr Int LAST = FACET_SIZE - 1;

            
            // Compute lambdas.
            
            if constexpr ( FACET_SIZE == 1 )
            {
                
                const Int i_last = vertices[LAST];
                
                if( dots[i_last][i_last] < dotvv )
                {
                    closest_facet     = facet;
                    dotvv             = dots[i_last][i_last];
                    best_lambda[LAST] = one;

                    copy_buffer<AMB_DIM>( &coords[i_last][0], &v[0]);
                }
                
                return;
            }
            else if constexpr ( FACET_SIZE == 2 )
            {
                // Setting up linear system for the barycenter coordinates lambda.

                // Find first vertex in facet.
                const Int i_0 = vertices[0];

                lambda[0] = Lambda[i_0] / Gram[i_0][i_0];
                lambda[1] = one - lambda[0];

                interior = (lambda[0] > eps) && (lambda[1] > eps);
            }
            else if constexpr ( FACET_SIZE == 3 )
            {
                // Setting up linear system for the barycenter coordinates lambda.

                // Find first two vertices in facet.
                const Int i_0 = vertices[0];
                const Int i_1 = vertices[1];

                // Using Cramer's rule to solve the linear system.
                const Real inv_det = one / ( Gram[i_0][i_0] * Gram[i_1][i_1] - Gram[i_0][i_1] * Gram[i_0][i_1] );

                lambda[0] = ( Gram[i_1][i_1] * Lambda[i_0] - Gram[i_0][i_1] * Lambda[i_1] ) * inv_det;
                lambda[1] = ( Gram[i_0][i_0] * Lambda[i_1] - Gram[i_0][i_1] * Lambda[i_0] ) * inv_det;
                lambda[2] = one - lambda[0] - lambda[1];

                // Check the barycentric coordinates for positivity ( and compute the 0-th coordinate).
                interior = (lambda[0] > eps) && (lambda[1] > eps) && (lambda[2] > eps);
            }
            else
            {
                // Setting up linear system for the barycenter coordinates lambda.
                
                for( Int i = 0; i < LAST; ++i )
                {
                    const Int i_i = vertices[i];
                    
                    lambda[i] = Lambda[i_i];
                    
                    g[i][i] = Gram[i_i][i_i];
                    
                    for( Int j = i+1; j < LAST; ++j )
                    {
                        Int j_j = vertices[j];
                        g[j][i] = g[i][j] = Gram[i_i][j_j];
                    }
                }
                
                // Cholesky decomposition
                for( Int k = 0; k < LAST; ++k )
                {
                    const Real a = g[k][k] = Sqrt(g[k][k]);
                    const Real ainv = one/a;

                    for( Int j = k+1; j < LAST; ++j )
                    {
                        g[k][j] *= ainv;
                    }

                    for( Int i = k+1; i < LAST; ++i )
                    {
                        for( Int j = i; j < LAST; ++j )
                        {
                            g[i][j] -= g[k][i] * g[k][j];
                        }
                    }
                }

                // Lower triangular back substitution
                for( Int i = 0; i < LAST; ++i )
                {
                    for( Int j = 0; j < i; ++j )
                    {
                        lambda[i] -= g[j][i] * lambda[j];
                    }
                    lambda[i] /= g[i][i];
                }

                // Upper triangular back substitution
                for( Int i = LAST; i --> 0; )
                {
                    for( Int j = i+1; j < LAST; ++j )
                    {
                        lambda[i] -= g[i][j] * lambda[j];
                    }
                    lambda[i] /= g[i][i];
                }
                
                // Check the barycentric coordinates for positivity ( and compute the 0-th coordinate).

                lambda[LAST] = one;
                
                for( Int k = 0; k < LAST; ++k)
                {
                    lambda[LAST] -= lambda[k];
                    interior = interior && ( lambda[k] > eps );
                }
                interior = interior && (lambda[LAST] > eps );
            }
            
            // If we arrive here, then facet_size > 1.
            
            if( interior )
            {
                // If the nearest point on `facet` lies in the interior, there is no point in searching its faces;
                // we simply mark them as visited.
                
                for( Int j = 0; j < FACET_SIZE; ++j )
                {
                    visited[ faces[j] ] = true;
                }
                
                Real facet_squared_dist = zero;
                
                // Computing facet_closest_point.
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    Real x = lambda[0] * coords[ vertices[0] ][ k ];
                    
                    for( Int j = 1; j < FACET_SIZE; ++j )
                    {
                        x += lambda[j] * coords[ vertices[j] ][ k ];
                    }
                    
                    facet_squared_dist += x * x;
                    facet_closest_point[k] = x;
                }
                
                if( facet_squared_dist < dotvv )
                {
                    closest_facet = facet;
                    dotvv = facet_squared_dist;
                    copy_buffer<AMB_DIM>  ( &facet_closest_point[0], &v[0]           );
                    copy_buffer<AMB_DIM+1>( &lambda[0],              &best_lambda[0] );
                }
            }
            else
            {
                // Try to visit all faces of `facet`.
                for( Int j = 0; j < FACET_SIZE; ++j )
                {
                    const Int face = faces[j];
                    
                    if( !visited[face] )
                    {
                        // We may ignore the error code of DistanceSubalgorithm because the faces of the simplex should be nondegenerate if we arrive here.
                        if constexpr ( FACET_SIZE > 1 )
                        {
                            DistanceSubalgorithm<FACET_SIZE-1>( face, facet_sizes[face] );
                        }
                    }
                }
            }
        }
        
    public:
    
/// ################################################################
/// ##                          Compute                           ##
/// ################################################################
        
        
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
            cref<PrimitiveBase_T> P,
            cref<PrimitiveBase_T> Q,
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
                
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    v[k] -= Q_supp[0][k];
                }
            }
            
            olddotvv = dotvv = dot_buffers<AMB_DIM>(v,v);

            // Unrolling the first iteration to avoid a call to DistanceSubalgorithm.
            
            // We use w = p-q, but do not define it explicitly.
            Real a = P.MinSupportVector( &v[0], &P_supp[0][0] );
            Real b = Q.MaxSupportVector( &v[0], &Q_supp[0][0] );
            dotvw  = a-b;
            
            Push<0>();
            
            closest_facet = 1;
            dotvv = dots[0][0];
            best_lambda[0] = one;

            copy_buffer<AMB_DIM>( &coords[0][0], v );
        
            while( true )
            {
                if( theta_squared * dotvv < TOL_squared )
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
                dotvw = a-b;
            
                if( collision_only && (dotvw > zero) && (theta_squared * dotvw * dotvw > TOL_squared) )
                {
                    separatedQ = true;
                    reason = GJK_Reason::Separated;
                    break;
                }
                
                
                if( Abs(dotvv - dotvw) <= eps * dotvv )
                {
                    reason = GJK_Reason::SmallResidual;
                    break;
                }
                
                in_simplex = Push<1>();
                
                if( in_simplex  )
                {
                    // The vertex added most recently to the simplex coincides with one of the previous vertices.
                    // So we stop here and use the old best_lambda. All we have to do is to reduce the simplex_size by 1.
                    --simplex_size;
                    reason = GJK_Reason::InSimplex;
                    break;
                }
                
                Int initial_facet = PrepareDistanceSubalgorithm();
                
                // This computes closest_facet, v, and dotvv of the current simplex.
                // If the simplex is degenerate, DistanceSubalgorithm terminates early and returns 1. Otherwise it returns 0.
                DistanceSubalgorithm<2>( initial_facet, facet_sizes[initial_facet] );
                
                simplex_size = facet_sizes[closest_facet];
                
                Deflate();
                
                if( Abs(olddotvv - dotvv) <= eps * dotvv )
                {
                    reason = GJK_Reason::SmallProgress;
                    break;
                }
                
            } // while( true )
            
            separatedQ = separatedQ || ( TOL_squared < theta_squared * dotvv );
            
            if( iter >= max_iter)
            {
//                wprint(ClassName()+"::Compute: Stopped because iter = " + ToString(iter) + " >= " + ToString(max_iter) + " = max_iter iterations reached.");
            }
            
            
        } // Compute

/// ################################################################
/// ##                        IntersectingQ                       ##
/// ################################################################
        
        bool IntersectingQ(
            cref<PrimitiveBase_T> P,
            cref<PrimitiveBase_T> Q,
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
        
        // Faster overloads for AABBs.
        template<typename SReal>
        bool IntersectingQ(
            cref<AABB<AMB_DIM,Real,Int,SReal>> P,
            cref<AABB<AMB_DIM,Real,Int,SReal>> Q,
            const Real theta_squared_ = one
        )
        {
            (void)theta_squared_;
            
            return AABB_SquaredDistance(P,Q) <= zero;
        }
        
/// ################################################################
/// ##                       SquaredDistance                      ##
/// ################################################################
        
        Real SquaredDistance(
            cref<PrimitiveBase_T> P,
            cref<PrimitiveBase_T> Q,
            const bool reuse_direction_ = false
        )
        {
            Compute(P, Q, false, reuse_direction_, zero );
            
            return dotvv;
        }
        
        
        // Faster overload for AABBs.
        template<typename SReal>
        Real SquaredDistance(
            cref<AABB<AMB_DIM,Real,Int,SReal>> P,
            cref<AABB<AMB_DIM,Real,Int,SReal>> Q
        )
        {
            return AABB_SquaredDistance(P,Q);
        }
        
        // Faster overload for AABBs and Points.
        template<typename SReal>
        Real SquaredDistance(
            cref<AABB <AMB_DIM,Real,Int,SReal>> P,
            cref<Point<AMB_DIM,Real,Int,SReal>> Q
        )
        {
            return AABB_Point_SquaredDistance(P,Q);
        }
        
/// ##########################################################################
/// ##                               Witnesses                              ##
/// ##########################################################################
        
        Real Witnesses(
            cref<PrimitiveBase_T> P, mptr<Real> x,
            cref<PrimitiveBase_T> Q, mptr<Real> y,
            const bool reuse_direction_ = false
        )
        {
            // x and y are the return parameters.
            // On return x lies in P while y lies in Q.
            // These points realize the minimal distance between any two points in these primitives.
            // Scalar return values is this minimal distance _squared_(!).
            
            Compute(P, Q, false, reuse_direction_, zero );

            for( Int k = 0; k < AMB_DIM; ++k )
            {
                y[k] = best_lambda[0] * Q_supp[0][k];
            }
            
            
            switch( simplex_size )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 1 )
                    {
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            y[k] += best_lambda[1] * Q_supp[1][k];
                        }
                    }
                    break;
                }
                case 3:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 2 )
                    {
                        for( Int j = 1; j < 3; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                    break;
                }
                case 4:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 3 )
                    {
                        for( Int j = 1; j < 4; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                    break;
                }
                case 5:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 4 )
                    {
                        for( Int j = 1; j < 5; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                    break;
                }
                default:
                {
                    // Only needed in high-dimensional settings.
                    if constexpr ( AMB_DIM >= 5 )
                    {
                        for( Int j = 1; j < simplex_size; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                }
            }
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                x[k] = y[k] + v[k];
            }
            
            return dotvv;
        } // Witnesses
        
/// ################################################################
/// ##                     Offset_IntersectingQ                   ##
/// ################################################################
        
        bool Offset_IntersectingQ(
            cref<PrimitiveBase_T> P, const Real P_offset,
            cref<PrimitiveBase_T> Q, const Real Q_offset,
            const bool reuse_direction_ = false
        )
        {
            const Real min_dist = P_offset + Q_offset;
            
            Compute(P, Q, true, reuse_direction_, min_dist * min_dist );
            
            return !separatedQ;
        }
        
/// ################################################################
/// ##                   Offset_SquaredDistance                   ##
/// ################################################################
        
        Real Offset_SquaredDistance(
            cref<PrimitiveBase_T> P, const Real P_offset,
            cref<PrimitiveBase_T> Q, const Real Q_offset,
            const bool reuse_direction_ = false
        )
        {
            const Real min_dist = P_offset + Q_offset;
            
            Compute(P, Q, false, reuse_direction_, min_dist * min_dist );
            
            const Real dist = Ramp( Sqrt(dotvv) - P_offset - Q_offset );
            
            return dist * dist;
        }
        
/// ##########################################################################
/// ##                           Offset_Witnesses                           ##
/// ##########################################################################
        
        Real Offset_Witnesses(
            cref<PrimitiveBase_T> P, const Real P_offset, mptr<Real> x,
            cref<PrimitiveBase_T> Q, const Real Q_offset, mptr<Real> y,
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
            
            const Real dist0 = Sqrt(dotvv);
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
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                y[k] = y_scale * v[k] + best_lambda[0] * Q_supp[0][k];
            }
            
            switch( simplex_size )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 1 )
                    {
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            y[k] += best_lambda[1] * Q_supp[1][k];
                        }
                    }
                    break;
                }
                case 3:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 2 )
                    {
                        for( Int j = 1; j < 3; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                    break;
                }
                case 4:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 3 )
                    {
                        for( Int j = 1; j < 4; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                    break;
                }
                case 5:
                {
                    // This if constexpr is to silence some compiler warnings.
                    if constexpr ( AMB_DIM >= 4 )
                    {
                        for( Int j = 1; j < 5; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                    break;
                }
                default:
                {
                    // Only needed in high-dimensional settings.
                    if constexpr ( AMB_DIM >= 4 )
                    {
                        for( Int j = 1; j < simplex_size; ++j )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                y[k] += best_lambda[j] * Q_supp[j][k];
                            }
                        }
                    }
                }
            }
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                x[k] = y[k] + x_scale * v[k];
            }
            
            return dist * dist;
            
        } // Offset_Witnesses

/// ########################################################################
/// ##                   InteriorPoints_SquaredDistance                   ##
/// ########################################################################
        
        Real InteriorPoints_SquaredDistance(
            cref<PrimitiveBase_T> P,
            cref<PrimitiveBase_T> Q
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
        
        
/// ########################################################################
/// ##                     MultipoleAcceptanceCriterion                   ##
/// ########################################################################
        
        bool MultipoleAcceptanceCriterion(
            cref<PrimitiveBase_T> P,
            cref<PrimitiveBase_T> Q,
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
            cref<AABB<AMB_DIM,Real,Int,SReal>> P,
            cref<AABB<AMB_DIM,Real,Int,SReal>> Q,
            const Real theta_squared_
        )
        {
            return Max( P.SquaredRadius(), Q.SquaredRadius() ) < theta_squared_ * AABB_SquaredDistance( P, Q );
        }
        
        std::string ClassName() const
        {
            return std::string("GJK")+"<"+ToString(AMB_DIM) + "," + TypeName<Real> + "," + TypeName<Real>+">";
        }
        
    }; // GJK

} // namespace Repulsor
