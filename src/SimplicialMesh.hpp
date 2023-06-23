#pragma once

#include "SimplicialMesh/SimplicialMeshBase.hpp"

namespace Repulsor
{
//    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
//    class TangentPoint;
    
    template<typename Real_, typename Int_>
    class SimplicialRemesherBase;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_>
    class SimplicialRemesher;
    
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class SimplicialMesh : public SimplicialMeshBase<Real_,Int_,SReal_,ExtReal_>
    {
    private:
        
        using Base_T = SimplicialMeshBase<Real_,Int_,SReal_,ExtReal_>;
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        using LInt    = typename Base_T::LInt;
        
        using TangentVector_T      = typename Base_T::TangentVector_T;
        using CotangentVector_T    = typename Base_T::CotangentVector_T;
        
        using SparseMatrix_T       = typename Base_T::SparseMatrix_T;
        using SparseBinaryMatrix_T = typename Base_T::SparseBinaryMatrix_T;
        
        using Primitive_T = Polytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using MovingPrimitive_T = MovingPolytope<DOM_DIM+1,AMB_DIM,GJK_Real,Int,SReal,Real,Int>;
        using BoundingVolume_T = AABB<AMB_DIM,GJK_Real,Int,SReal>;
        
        
        using ClusterTree_T           = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        using BlockClusterTree_T      = BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal,true>;
        using ObstacleBlockClusterTree_T
                                      = BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal,false>;
        using CollisionTree_T         = CollisionTree<AMB_DIM,Real,Int,SReal,ExtReal,true>;
        using ObstacleCollisionTree_T = CollisionTree<AMB_DIM,Real,Int,SReal,ExtReal,false>;

        using Remesher_T         = SimplicialRemesher<DOM_DIM,AMB_DIM,Real,Int>;
        
        using RemesherBase_T     = typename Base_T::Remesher_T;
        
        using Obstacle_T         = Base_T;
        
        static constexpr Int  FAR_DIM = 1 + AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        static constexpr Int NEAR_DIM = 1 + (DOM_DIM+1) * AMB_DIM + (AMB_DIM * (AMB_DIM + 1)) / 2;
        
        static constexpr Int  SIZE      = DOM_DIM + 1;
        static constexpr Int  HULL_SIZE = AMB_DIM * SIZE;
        static constexpr Real nth       = Scalar::Inv<Real>( static_cast<Real>(SIZE));
        
        SimplicialMesh() = default;

        SimplicialMesh(
            const Tensor2<Real,Int> & V_coords_,
            // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Tensor2<Int,Int> & simplices_,
            // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T thread_count_ = 1
        )
        :   SimplicialMesh(
                  V_coords_.data(),
                  V_coords_.Dimension(0),
                  false,
                  simplices_.data(),
                  simplices_.Dimension(0),
                  false,
                  static_cast<Int>(thread_count_)
            )
        {
            ptic(className()+"()");
            if( V_coords_.Dimension(1) != AMB_DIM )
            {
                eprint(className()+" : V_coords.Dimension(1) != AMB_DIM");
                ptoc(className()+"()");
                return;
            }
            if( simplices_.Dimension(1) != DOM_DIM+1 )
            {
                eprint(className()+" : simplices_.Dimension(1) != DOM_DIM+1");
                ptoc(className()+"()");
                return;
            }
            ptoc(className()+"()");
        }
        
#ifdef LTEMPLATE_H
        
        SimplicialMesh(
            const mma::TensorRef<ExtReal> & V_coords_,   // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const mma::TensorRef<long long> & simplices_,   // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const long long thread_count_ = 1
        )
        :   SimplicialMesh(
                V_coords_.data(),
                static_cast<Int>(V_coords_.dimensions()[0]),
                simplices_.data(),
                static_cast<Int>(simplices_.dimensions()[0]),
                static_cast<Int>(thread_count_)
            )
        {
            ptic(className()+"() (from MTensor)");
            if( V_coords_.dimensions()[1] != AMB_DIM )
            {
                eprint(className()+" : V_coords_.dimensions()[1] != AMB_DIM");
                ptoc(className()+"() (from MTensor)");
                return;
            }
            if( simplices_.dimensions()[1] != DOM_DIM+1 )
            {
                eprint(className()+" : simplices_.dimensions()[1] != DOM_DIM+1");
                ptoc(className()+"() (from MTensor)");
                return;
            }
            ptoc(className()+"() (from MTensor)");
        }
        
#endif

        template<typename ExtReal_2, typename ExtInt>
        SimplicialMesh(
            const ExtReal_2 * V_coords_, // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Size_T vertex_count_,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T simplex_count_,
            const Size_T thread_count_ = 1
        )
        :   SimplicialMesh( V_coords_, vertex_count_, false, simplices_, simplex_count_, false, thread_count_)
        {}
        
        template<typename ExtReal_2, typename ExtInt>
        SimplicialMesh(
            const ExtReal_2 * vertex_coords_, // vertex coordinates; assumed to be of size vertex_count_ x AMB_DIM
            const Size_T vertex_count_,
            const bool vertex_coords_transpose,
            const ExtInt * simplices_, // simplices; assumed to be of size simplex_count_ x (DOM_DIM+1)
            const Size_T simplex_count_,
            const bool simplices_transpose,
            const Size_T thread_count_ = 1
        )
        :   Base_T    ( static_cast<Int>(thread_count_) )
        ,   V_coords  ( ToTensor2<Real,Int>(
                            vertex_coords_,
                            static_cast<Int>(vertex_count_),
                            static_cast<Int>(AMB_DIM),
                            vertex_coords_transpose
                        )
                      )
        ,   V_coords_frozen ( V_coords )
        ,   simplices ( ToTensor2<Int,Int>(
                            simplices_,
                            static_cast<Int>(simplex_count_),
                            static_cast<Int>(DOM_DIM+1),
                            simplices_transpose
                        )
                      )
        ,   details   ( static_cast<Int>(thread_count_) )
        {
            ptic(className()+" (pointer)");
            ptoc(className()+" (pointer)");
        }

        
        virtual ~SimplicialMesh() override = default;
        
    public:
        
        using Base_T::block_cluster_tree_settings;
        using Base_T::cluster_tree_settings;
        using Base_T::adaptivity_settings;
        using Base_T::ThreadCount;
        
        static constexpr Real StandardSimplexVolume()
        {
            return Scalar::Inv<Real>( Factorial(static_cast<Real>(DOM_DIM)) );
        }
        
    protected:

        
        mutable SReal max_update_step_size = 0;
        
        // Current vertex coordinates; will be updated by SemiStaticUpdate.
        mutable Tensor2<Real,Int> V_coords;
        // Frozen vertex coordinates to be used for the clustering, but not for the energy computation.
        // Won't be updated by SemiStaticUpdate so that energies that depend on the `cluster_tree` will look differentiable.
        Tensor2<Real,Int> V_coords_frozen;
        
        Tensor2<Int,Int>  simplices;
                
        mutable Primitive_T P_proto;
        
        mutable Tensor1<Int,Int> simplex_row_pointers;
        mutable Tensor1<Int,Int> simplex_column_indices;
        
        SimplicialMeshDetails<DOM_DIM,AMB_DIM,Real_,Int_> details;

        
    protected:
        
        void ComputeNearFarDataOps(
                  Tensor2<Real,Int> & P_coords,
                  Tensor3<Real,Int> & P_hull_coords,
                  Tensor2<Real,Int> & P_near,
                  Tensor2<Real,Int> & P_far,
            SparseMatrix_T & DiffOp,
            SparseMatrix_T & AvOp
        ) const
        {
            ptic(ClassName()+"::ComputeNearFarDataOps");
            
            ParallelDo(
                [&,this]( const Int thread )
                {
                    
                    mut<LInt> Av_outer = AvOp.Outer().data();
                    mut< Int> Av_inner = AvOp.Inner().data();
                    mut<Real> Av_value = AvOp.Values().data();
        
                    mut<LInt> Diff_outer = DiffOp.Outer().data();
                    mut< Int> Diff_inner = DiffOp.Inner().data();
                    mut<Real> Diff_value = DiffOp.Values().data();
        
                    Tiny::Vector<AMB_DIM,Real,Int> center;
                    
                    Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
                    Tiny::Matrix<AMB_DIM,SIZE,   Real,Int> Df;
                    
                    Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                    Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
                    Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
                    
                    Tiny::SelfAdjointMatrix    <DOM_DIM,Real,Int> g;
                    
                    Tiny::Vector<SIZE,Int,Int> simplex;
                    Tiny::Vector<SIZE,Int,Int> s_simplex;   // Sorted simplex.
        
                
                    Real std_simplex_volume = 1;
                    for( Int l = 2; l < SIZE; ++l )
                    {
                        std_simplex_volume *= l;
                    }
                    std_simplex_volume = Scalar::Inv<Real>(std_simplex_volume);
                    
//                    dump(std_simplex_volume);
                    
                    const Int simplex_count = simplices.Dimension(0);
                    const Int i_begin = JobPointer<Int>(simplex_count, ThreadCount(), thread  );
                    const Int i_end   = JobPointer<Int>(simplex_count, ThreadCount(), thread+1);
        
                    Sorter<SIZE,Int> sort;
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        mut<Real> near   = P_near.data(i);
                        mut<Real> far    = P_far.data(i);
                        
                          simplex.Read( simplices.data(i) );
                        s_simplex.Read( simplex.data() );
                        
//                        dump(simplex);
                      
                        // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                        sort( &s_simplex[0] );

//                        dump(s_simplex);
                        
                        Av_outer[i+1] = (i+1) * SIZE;
                        
                        s_simplex.Write( &Av_inner[SIZE * i] );
                        fill_buffer<SIZE>( &Av_value[SIZE * i], nth );
                        
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            const Int row = AMB_DIM * i + k;
                            
                            const Int rp = row * SIZE;
                            
                            Diff_outer[row+1] = rp + SIZE;
                            
                            s_simplex.Write( &Diff_inner[rp] );
                        }
           
                        for( Int l = 0; l < SIZE; ++l )
                        {
                            copy_buffer<AMB_DIM>( V_coords.data(simplex[l]), hull[l] );
                        }
                        
//                        dump(hull);

                        hull.Write( P_hull_coords.data(i) );
                        hull.Write( &near[1]              );
                        
                        
                        copy_buffer<AMB_DIM>(hull[0],center.data());
                        for( Int l = 1; l < SIZE; ++l )
                        {
                            add_to_buffer<AMB_DIM>( hull[l], center.data() );
                        }
                        center *= nth;
                        
                        center.Write( &far[1] );
                        center.Write( P_coords.data(i) );
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                df[k][l] =
                                    V_coords[s_simplex[l+1]][k]
                                    -
                                    V_coords[s_simplex[0  ]][k];
                            }
                        }
                        
//                        dump(df);
                        
                        df.Transpose( dfdagger );
                        
//                        dump(dfdagger);
                        
                        // g = df^T * df.
                        // At the moment dfdagger is just the transpose of df.
                        Tiny::gemm<Op::Id,Op::Id,DOM_DIM,DOM_DIM,AMB_DIM,
                            Scalar::Flag::Plus,Scalar::Flag::Zero
                        >(
                            Scalar::One<Real>,  dfdagger.data(), AMB_DIM,
                                                df.data(),       DOM_DIM,
                            Scalar::Zero<Real>, g.data(),        DOM_DIM
                        );
                        
//                        dump(g);
                        
                        // Factorize g in place.
                        g.Cholesky();
                        
//                        dump(g);
                        
                        Real a = StandardSimplexVolume();
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            a *= g[l][l];
                        }
                        
                        far[0] = near[0] = a;
                        
//                        dump(near[0]);
                        
                        //  dfdagger = g^{-1} * df^T
                        g.CholeskySolve( dfdagger );
                        
//                        dump(dfdagger);

                        // P = id - df * dfdagger
                        Tiny::gemm<Op::Id,Op::Id,AMB_DIM,AMB_DIM,DOM_DIM,
                            Scalar::Flag::Minus,Scalar::Flag::Zero
                        >(
                            -Scalar::One<Real>, df.data(),       DOM_DIM,
                                                dfdagger.data(), AMB_DIM,
                            Scalar::Zero<Real>, P.data(),        AMB_DIM
                        );
                        
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            P[k][k] += Scalar::One<Real>;
                        }
                        
//                        dump(P);
                        
                        {
                            // TODO: variable pos introduces a bad dependency in loop.
                            
                            Int pos = 0;
                            
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                for( Int l = k; l < AMB_DIM; ++l )
                                {
                                    ++pos;
                                    near[AMB_DIM * SIZE + pos] = far[AMB_DIM + pos] = P[k][l];
                                }
                            }
                        }
//
//                        for( Int k = 0; k < NearDim(); ++k )
//                        {
//                            valprint("near["+ToString(k)+"]",near[k]);
//                        }
//
//                        for( Int k = 0; k < FAR_DIM; ++k )
//                        {
//                            valprint("far["+ToString(k)+"]",far[k]);
//                        }
                        

                        // Create derivative operator  (AMB_DIM x SIZE matrix).
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            Df[k][0] = 0;
                            
                            for( Int l = 0; l < DOM_DIM; ++l )
                            {
                                Df[k][0  ] -= dfdagger[l][k];
                                Df[k][l+1]  = dfdagger[l][k];
                            }
                        }
                        
//                        dump(Df)
                        
                        Df.Write( &Diff_value[ HULL_SIZE * i ] );

                    }
                },
                ThreadCount()
            );
        
            ptoc(ClassName()+"::ComputeNearFarDataOps");
        };
        
        void ComputeNearFarData(
            Tensor2<Real,Int> & P_near,
            Tensor2<Real,Int> & P_far
        ) const
        {
            ptic(ClassName()+"::ComputeNearFarData");
            
            ParallelDo(
                [&]( const Int thread )
                {
                    Tiny::Vector<AMB_DIM,Real,Int> center;
                    
                    Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
                    
                    Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                    Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
                    Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
                    
                    Tiny::SelfAdjointMatrix<DOM_DIM,Real,Int> g;
                    
                    Tiny::Vector<SIZE,Int,Int> simplex;
                    Tiny::Vector<SIZE,Int,Int> s_simplex;   // Sorted simplex.
                    
                    const Int simplex_count = simplices.Dimension(0);
                    const Int i_begin = JobPointer<Int>(simplex_count, ThreadCount(), thread  );
                    const Int i_end   = JobPointer<Int>(simplex_count, ThreadCount(), thread+1);
        
                    Sorter<SIZE,Int> sort;
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        mut<Real> near   = P_near.data(i);
                        mut<Real> far    = P_far.data(i);
                        
                          simplex.Read( simplices.data(i) );
                        s_simplex.Read( simplex.data() );
                      
                        // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                        sort( &s_simplex[0] );

                        for( Int l = 0; l < SIZE; ++l )
                        {
                            copy_buffer<AMB_DIM>( V_coords.data(simplex[l]), hull[l] );
                        }

                        hull.Write( &near[1]              );
                        
                        
                        copy_buffer<AMB_DIM>(hull[0],center.data());
                        for( Int l = 1; l < SIZE; ++l )
                        {
                            add_to_buffer<AMB_DIM>( hull[l], center.data() );
                        }
                        center *= nth;
                        
                        center.Write( &far[1] );
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                df[k][l] =
                                    V_coords[s_simplex[l+1]][k]
                                    -
                                    V_coords[s_simplex[0  ]][k];
                            }
                        }
                        
                        df.Transpose( dfdagger );
                        
                        // g = df^T * df.
                        // At the moment dfdagger is just the transpose of df.
                        Tiny::gemm<Op::Id,Op::Id,DOM_DIM,DOM_DIM,AMB_DIM,
                            Scalar::Flag::Plus,Scalar::Flag::Zero
                        >(
                            Scalar::One<Real>,  dfdagger.data(), AMB_DIM,
                                                df.data(),       DOM_DIM,
                            Scalar::Zero<Real>, g.data(),        DOM_DIM
                        );
                        
                        // Factorize g in place.
                        g.Cholesky();
                        
                        Real a = StandardSimplexVolume();
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            a *= g[l][l];
                        }
                        
                        far[0] = near[0] = a;
                        
                        //  dfdagger = g^{-1} * df^T
                        g.CholeskySolve( dfdagger );

                        // P = id - df * dfdagger
                        Tiny::gemm<Op::Id,Op::Id,AMB_DIM,AMB_DIM,DOM_DIM,
                            Scalar::Flag::Minus,Scalar::Flag::Zero
                        >(
                            -Scalar::One<Real>, df.data(),       DOM_DIM,
                                                dfdagger.data(), AMB_DIM,
                            Scalar::Zero<Real>, P.data(),        AMB_DIM
                        );
                        
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            P[k][k] += Scalar::One<Real>;
                        }
                        
                        {
                            // TODO: variable pos introduces a bad dependency in loop.
                            
                            Int pos = 0;
                            
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                for( Int l = k; l < AMB_DIM; ++l )
                                {
                                    ++pos;
                                    near[AMB_DIM * SIZE + pos] = far[AMB_DIM + pos] = P[k][l];
                                }
                            }
                        }

                    }
                },
                ThreadCount()
            );
        
            ptoc(ClassName()+"::ComputeNearFarData");
    };
    
    protected:
        
        virtual SparseMatrix_T H1Metric( const Real c_1, const Real c_0 ) const override
        {
            Tensor2<Int, Int> ilist ( SimplexCount(), SIZE * SIZE );
            Tensor2<Int, Int> jlist ( SimplexCount(), SIZE * SIZE );
            Tensor2<Real,Int> alist ( SimplexCount(), SIZE * SIZE );
            
            ParallelDo(
                [this,&ilist,&jlist,&alist,c_1,c_0]( const Int thread )
                {
                    const Int n = SimplexCount();
                    const Int k_begin = JobPointer<Int>( n, ThreadCount(), thread     );
                    const Int k_end   = JobPointer<Int>( n, ThreadCount(), thread + 1 );
                    
                    Tiny::Vector<SIZE,Int,Int> simplex;
                    Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                    Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dftransp;
                    Tiny::Matrix<DOM_DIM,DOM_DIM,Real,Int> g_inv;
                    Tiny::Matrix<DOM_DIM,DOM_DIM,Real,Int> id;
                    
                    for( Int i = 0; i < DOM_DIM; ++i )
                    {
                        for( Int j = 0; j < DOM_DIM; ++j )
                        {
                            id[i][j] = Int(i==j);
                        }
                    }
                    
                    Tiny::Matrix<SIZE,SIZE,Real,Int> mass;
                    
                    const Real m_factor = c_0 / StandardSimplexVolume() / Factorial(static_cast<Real>(DOM_DIM + 2));
                    
                    for( Int i = 0; i < SIZE; ++i )
                    {
                        for( Int j = 0; j < SIZE; ++j )
                        {
                            mass[i][j] = m_factor * (Scalar::One<Real> + Delta<Real>(i,j) );
                        }
                    }
                    
                    Tiny::SelfAdjointMatrix<DOM_DIM,Real,Int> g;
                    
                    Tiny::Matrix<SIZE,SIZE, Int,Int> idx;
                    Tiny::Matrix<SIZE,SIZE, Int,Int> jdx;
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        simplex.Read( simplices.data(k) );
                        
                        Tiny::Matrix<SIZE,SIZE,Real,Int> val ( Scalar::Zero<Real> );
                        
                        for( Int i = 0; i < AMB_DIM; ++i )
                        {
                            for( Int j = 0; j < DOM_DIM; ++j )
                            {
                                df(i,j) = V_coords(simplex[j+1],i) - V_coords(simplex[0],i);
                            }
                        }
                        
                        df.Transpose( dftransp );

                        // g = df^T * df.
                        Dot<0>( dftransp, df, g_inv );
                        
                        g.Read( g_inv.data() );
                        
                        g.Cholesky();
                        
                        // Compute simplex area
                        Real a = StandardSimplexVolume();
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            a *= g[l][l];
                        }
                        
                        g_inv.Read( id.data() );
                        
                        g.CholeskySolve( g_inv );
                        
                        const Real l_factor = a * c_1;
                        
                        for( Int i = 0; i < DOM_DIM; ++i )
                        {
                            for( Int j = 0; j < DOM_DIM; ++j )
                            {
                                const Real a_ij = l_factor * g_inv[i][j];
                                
                                val[0  ][0  ] += a_ij;
                                val[i+1][0  ] -= a_ij;
                                val[0  ][j+1] -= a_ij;
                                val[i+1][j+1]  = a_ij;
                            }
                        }
                        
                        combine_buffers<
                            SIZE * SIZE, Scalar::Flag::Generic, Scalar::Flag::Plus
                        >(
                            a,                 mass.data(),
                            Scalar::One<Real>, val.data()
                        );
                        
                        val.Write( alist.data(k) );
                        
                        for( Int i = 0; i < SIZE; ++i )
                        {
                            for( Int j = 0; j < SIZE; ++j )
                            {
                                idx[i][j] = simplex[i];
                                jdx[i][j] = simplex[j];
                            }
                        }
                        
                        idx.Write( ilist.data(k) );
                        jdx.Write( jlist.data(k) );
                    }
                },
                ThreadCount()
            );
            
            SparseMatrix_T A (
                LInt(SimplexCount()) * LInt(SIZE * SIZE),
                ilist.data(),  jlist.data(),  alist.data(),
                VertexCount(), VertexCount(), ThreadCount(),
                true, false
            );
            
            return A;
        }
        
    public:
        
        virtual const SparseMatrix_T & StiffnessMatrix() const override
        {
            std::string tag ("StiffnessMatrix");
            if( !this->InCacheQ(tag))
            {
                ptic(ClassName()+"::"+tag);
                
                this->SetCache( tag,
                   std::any(std::move(H1Metric(1,0)))
                );
                
                ptoc(ClassName()+"::"+tag);
            }
            
            return std::any_cast<SparseMatrix_T &>( this->GetCache(tag) );
        }
        
        virtual const SparseMatrix_T & MassMatrix() const override
        {
            std::string tag ("MassMatrix");
            if( !this->InCacheQ(tag))
            {
                ptic(ClassName()+"::"+tag);
                
                this->SetCache( tag,
                    std::any(std::move(H1Metric(0,1)))
                );
                
                ptoc(ClassName()+"::"+tag);
            }
            
            return std::any_cast<SparseMatrix_T &>( this->GetCache(tag) );
        }
        
        virtual const Tensor1<Int,Int> & NestedDissectionOrdering() const override
        {
            std::string tag ("NestedDissectionOrdering");
            
            if( !this->InPersistentCacheQ(tag))
            {
                ptic(ClassName()+"::"+tag);
                
                auto & tree = GetClusterTree();
                
                // Single-threaded is actually faster here.
                // TODO: Maybe this is just because of the mutex?
                constexpr Int local_thread_count = 1;
                
                const Int m = VertexCount();
                const Int n = SimplexCount();
                
                Tensor1<Int,Int> perm_0 = iota<Int>( m );
                Tensor1<Int,Int> perm_1 ( m );
                
                const SparseMatrix_T & Adj = tree.LowOrderPostProcessor();
                
                const Tensor1<Int,Int> & C_begin = tree.ClusterBegin();
                const Tensor1<Int,Int> & C_end   = tree.ClusterEnd();
                const Tensor1<Int,Int> & C_left  = tree.ClusterLeft();
                const Tensor1<Int,Int> & C_right = tree.ClusterRight();
                
                Tensor1<Int,Int> V_begin ( tree.ClusterCount() );
                Tensor1<Int,Int> V_end   ( tree.ClusterCount() );
                
                Tensor2<Real,Int> v ( n, 2 );
                Tensor2<Real,Int> w ( m, 2 );
                
                constexpr Int16 s_zero = 0;
                constexpr Int16 s_one  = 1;
                constexpr Int16 s_two  = 2;
                
                Tensor1<Int16,Int> type ( m );
                
                V_begin[0] = 0;
                V_end  [0] = m;
                
                std::vector<Int> queue_0;
                std::vector<Int> queue_1;
                
                std::mutex mutex;
                
                queue_0.reserve( tree.ClusterCount() );
                queue_1.reserve( tree.ClusterCount() );
                
                // Push root cluster.
                queue_0.push_back( Scalar::Zero<Int> );
                
                Int level = 0;
                
                while( queue_0.size() > 0 )
                {
                    debug_print("Starting level " + ToString(level) + "." );
                    
                    queue_1.resize(0);
                    
                    debug_print("queue = "+ToString(queue));
                    
                    perm_1.Read( perm_0.data() );
                    
                    v.SetZero();
                    
                    debug_print("Compute indicators.");
                    ParallelDo(
                        [&]( const Int k )
                        {
                            const Int C = queue_0[k];
                            const Int L = C_left [C];
                            const Int R = C_right[C];
                            
                            if( L >= 0 )
                            {
                                
                                const Int L_begin = C_begin[L];
                                const Int L_end   = C_end  [L];
                                
                                const Int R_begin = C_begin[R];
                                const Int R_end   = C_end  [R];
                                
                                for( Int i = L_begin; i < L_end; ++i )
                                {
                                    v[i][0] = Scalar::One<Real>;
                                }
                                
                                for( Int i = R_begin; i < R_end; ++i )
                                {
                                    v[i][1] = Scalar::One<Real>;
                                }
                            }
                        },
                        static_cast<Int>(queue_0.size()), local_thread_count
                    );
                        

                    Adj.template Dot_<2>(
                        Scalar::One <Real>, v.data(), 2,
                        Scalar::Zero<Real>, w.data(), 2,
                        2
                    );

                    debug_print("Compute types.");
                    ParallelDo(
                        [&]( const Int i )
                        {
                            const Int j = perm_0[i];
                            
                            type[i] = - s_one
                                +
                                ( (w[j][0] > Scalar::Zero<Real>) ? s_one : s_zero )
                                +
                                ( (w[j][1] > Scalar::Zero<Real>) ? s_two : s_zero );
                        },
                        m, ThreadCount()
                    );
                    
                    debug_print("Modify permutation.");
                    ParallelDo(
                        [&]( const Int k )
                        {
                            const Int C = queue_0[k];
                            const Int L = C_left [C];
                            const Int R = C_right[C];

                            if( L >= 0 )
                            {
                                const Int i_begin = V_begin[C];
                                const Int i_end   = V_end  [C];
                                
                                Tiny::Vector<3,Int,Int> ctr ( Scalar::Zero<Int> );
                                
                                // Compute counters;
                                for( Int i = i_begin; i < i_end; ++i )
                                {
                                    const Int16 t = type[i];
                                    
                                    if( t >= s_zero )
                                    {
                                        ++ctr[t];
                                    }
                                }
                                
                                Tiny::Vector<4,Int,Int> pos;
                                
                                pos[0] = i_begin;
                                pos[1] = pos[0] + ctr[0];
                                pos[2] = pos[1] + ctr[1];
                                pos[3] = pos[2] + ctr[2];
                                
                                if( (ctr[0] > 0) && (ctr[1] > 0) && ( pos[3] == i_end ) )
                                {
                                    V_begin(L) = pos[0];
                                    V_end  (L) = pos[1];
                                    V_begin(R) = pos[1];
                                    V_end  (R) = pos[2];
                                    
                                    {
                                        const std::lock_guard<std::mutex> lock( mutex );
                                        queue_1.push_back( L );
                                        queue_1.push_back( R );
                                    }
                                    
                                    // Modify permutation;
                                    for( Int i = i_begin; i < i_end; ++i )
                                    {
                                        const Int16 t = type[i];
                                        
                                        if( t >= s_zero )
                                        {
                                            const Int p = pos[t]++;
                                            
                                            perm_1[p] = perm_0[i];
                                        }
                                    }
                                }
                            }
                        },
                        static_cast<Int>(queue_0.size()), local_thread_count
                    );
                
                    debug_print("Swapping.");
                    
                    std::swap( queue_0, queue_1 );
                    
                    swap( perm_0, perm_1 );
                    
                    debug_print("Done with level " + ToString(level) + "." );
                    
                    ++level;
                }
                
                this->SetPersistentCache( tag,
                   std::any( std::move(perm_0) )
                );
                
                ptoc(ClassName()+"::"+tag);
            }
            
            return std::any_cast<Tensor1<Int,Int> &>( this->GetPersistentCache(tag) );
        }
        
    public:
        
        virtual Int DomDim() const override
        {
            return DOM_DIM;
        }
        
        virtual Int AmbDim() const override
        {
            return AMB_DIM;
        }

        virtual const Tensor2<Real,Int> & VertexCoordinates() const override
        {
            return V_coords;
        }
        
        virtual const Tensor2<Int,Int> & Simplices() const override
        {
            return simplices;
        }

        virtual constexpr Int FarDim() const override
        {
            return FAR_DIM;
        }
        
        virtual constexpr Int NearDim() const override
        {
            return NEAR_DIM;
        }
        
        virtual Int VertexCount() const override
        {
            return V_coords.Dimension(0);
        }
        
        virtual Int SimplexCount() const override
        {
            return simplices.Dimension(0);
        }
    
        virtual Int DofCount() const override
        {
            return VertexCount() * AMB_DIM;
        }
        
        virtual const Real * Dofs() const override
        {
            return V_coords.data();
        }
        
        virtual void SemiStaticUpdate( ptr<ExtReal> V_coords_, const bool transp_ = false ) const override
        {
            ptic(className()+"::SemiStaticUpdate");
            
            // We read the new coordinates onto V_coords, but not into C_coords_frozen!
            if( transp_ )
            {
                V_coords.ReadTransposed(V_coords_);
            }
            else
            {
                V_coords.Read(V_coords_);
            }
            
            this->ClearCache();
            
            Tensor2<Real,Int> P_near( SimplexCount(), NEAR_DIM );
            Tensor2<Real,Int> P_far ( SimplexCount(), FAR_DIM  );

//            details.ComputeNearFarData( new_V_coords, simplices, P_near, P_far );
            ComputeNearFarData( P_near, P_far );
            
            GetClusterTree().SemiStaticUpdate( P_near, P_far );
            
            ptoc(className()+"::SemiStaticUpdate");
        }
       
        
        virtual void LoadUpdateVectors(
            ptr<ExtReal> vecs,
            const ExtReal max_time,
            const bool transp_
        ) const override
        {
            ptic(className()+"::LoadUpdateVectors");
            
            max_update_step_size = static_cast<SReal>(max_time);
            
            Tensor2<Real,Int> V_updates ( VertexCount(), AMB_DIM );

            if( vecs == nullptr )
            {
                V_updates.Fill(static_cast<Real>(0));
            }
            else
            {
                if( transp_ )
                {
                    V_updates.ReadTransposed(vecs);
                }
                else
                {
                    V_updates.Read(vecs);
                }
            }
            // ATTENTION: We reorder already here outside of the cluster tree to save a copy operation of a big Tensor2!
            
            MovingPrimitive_T P_moving;
            
            Tensor2<SReal,Int> P_velocities_serialized ( SimplexCount(), P_moving.VelocitySize(), 0 );
            mut<SReal> P_v_ser = P_velocities_serialized.data();
            
            const Tensor1<Int,Int> & P_ordering = GetClusterTree().PrimitiveOrdering();

            JobPointers<Int> job_ptr ( SimplexCount(), GetClusterTree().ThreadCount() );
            
            ParallelDo(
                [&]( const Int thread )
                {
                    MovingPrimitive_T P_mov;
                    
                    const Int i_begin = job_ptr[thread  ];
                    const Int i_end   = job_ptr[thread+1];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        const Int j = P_ordering[i];
                        P_mov.FromVelocitiesIndexList( V_updates.data(), simplices.data(), j );
                        P_mov.WriteVelocitiesSerialized( P_v_ser, i );
                    }
                },
                job_ptr.ThreadCount()
            );
            
            // P_moving and P_velocities_serialized will be swapped against potentiall empty containers.
            GetClusterTree().TakeUpdateVectors(
                P_moving, P_velocities_serialized, static_cast<SReal>(max_time)
            );
            
            ptoc(className()+"::LoadUpdateVectors");
        }
        
        virtual ExtReal MaximumSafeStepSize(
            ptr<ExtReal> vecs,
            const ExtReal max_time,
            const ExtReal TOL = scalar_cast<ExtReal>(0.0625),
            const bool transp_ = false
        ) override
        {
            ptic(className()+"::MaximumSafeStepSize");
            
            LoadUpdateVectors( vecs, max_time, transp_ );
            
            ExtReal t = max_time;

            if( this->InCacheQ("Obstacle") )
            {
                GetObstacle().LoadUpdateVectors( static_cast<ExtReal *>(nullptr), max_time );
                
                t = GetObstacleCollisionTree().MaximumSafeStepSize(t,TOL);
                
                logprint("GetObstacleCollisionTree().MaximumSafeStepSize(t) = "+ToString(t));
            }
                
            t = GetCollisionTree().MaximumSafeStepSize(t,TOL);
            
            logprint("GetCollisionTree().MaximumSafeStepSize(t)         = "+ToString(t));
            
            ptoc(className()+"::MaximumSafeStepSize");
            
            return t;
        }
        
        virtual const ClusterTree_T & GetClusterTree() const override
        {
            static std::string tag ( "ClusterTree" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                ptic(className()+"::GetClusterTree");
                if( (V_coords.Dimension(0) > 0) && (simplices.Dimension(0) > 0) )
                {
                    ptic("Allocations");
                    auto P_coords      = Tensor2<Real,Int> ( SimplexCount(), AMB_DIM, static_cast<Real>(0) );
                    auto P_hull_coords = Tensor3<Real,Int> ( SimplexCount(), SIZE, AMB_DIM );
                    auto P_near        = Tensor2<Real,Int> ( SimplexCount(), NEAR_DIM );
                    auto P_far         = Tensor2<Real,Int> ( SimplexCount(), FAR_DIM );

                    Tensor2<SReal,Int> P_serialized ( SimplexCount(), P_proto.Size() );

                    auto DiffOp = SparseMatrix_T(
                        SimplexCount() * AMB_DIM,
                        VertexCount(),
                        SimplexCount() * HULL_SIZE,
                        ThreadCount()
                    );

                    DiffOp.Outer()[SimplexCount() * AMB_DIM] = SimplexCount() * HULL_SIZE;

                    auto AvOp = SparseMatrix_T(
                        SimplexCount(),
                        VertexCount(),
                        SimplexCount() * SIZE,
                        ThreadCount()
                    );

                    AvOp.Outer()[SimplexCount()] = SimplexCount() * SIZE;

                    ptoc("Allocations");

                    // What remains is to compute P_coords, P_hull_coords, P_near and P_far and the nonzero values of DiffOp.
//                    details.ComputeNearFarDataOps( V_coords, simplices, P_coords, P_hull_coords, P_near, P_far, DiffOp, AvOp );
                    
                    ComputeNearFarDataOps( P_coords, P_hull_coords, P_near, P_far, DiffOp, AvOp );


                    const JobPointers<Int> job_ptr ( SimplexCount(), ThreadCount() );

                    ptic("Creating primitives");
                    ParallelDo(
                        [&]( const Int thread )
                        {
                            Primitive_T P;

                            const Int i_begin = job_ptr[thread];
                            const Int i_end   = job_ptr[thread+1];

                            for( Int i = i_begin; i < i_end; ++i )
                            {
                                P.SetPointer( P_serialized.data(), i );
                                // Beware, we use the frozen coordinates for clustering!
                                P.FromIndexList( V_coords_frozen.data(), simplices.data(), i );
                            }
                        },
                        job_ptr.ThreadCount()
                    );

                    ptoc("Creating primitives");

                    ptic("Initializing cluster prototypes");

                    std::shared_ptr<BoundingVolume_T> C_proto;

                    switch( cluster_tree_settings.bounding_volume_type )
                    {
                        case BoundingVolumeType::AABB_MedianSplit:
                        {
                            logprint("Using AABB_MedianSplit as bounding volume.");
                            C_proto = std::shared_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_MedianSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                            break;
                        }
                        case BoundingVolumeType::AABB_PreorderedSplit:
                        {
                            logprint("Using AABB_PreorderedSplit as bounding volume.");
                            C_proto = std::shared_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_PreorderedSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                            break;
                        }
                        default:
                        {
                            logprint("Using AABB_LongestAxisSplit as bounding volume.");
                            C_proto = std::shared_ptr<BoundingVolume_T>(
                                static_cast<BoundingVolume_T*>(new AABB_LongestAxisSplit<AMB_DIM,GJK_Real,Int,SReal>)
                            );
                        }
                    }

                    ptoc("Initializing cluster prototypes");

                    if( cluster_tree_settings.thread_count <= 0 )
                    {
                        cluster_tree_settings.thread_count = ThreadCount();
                    }

                    this->SetPersistentCache( tag,
                        std::make_any<ClusterTree_T>(
                            P_proto, P_serialized, *C_proto, Tensor1<Int,Int>(0),
                            P_near, P_far,
                            DiffOp, AvOp,
                            cluster_tree_settings
                        )
                    );
                }
                else
                {
                    eprint(ClassName()+"(V_coords.Dimension(0) <= 0) || (simplices.Dimension(0) <= 0)");

                    this->SetPersistentCache( tag, std::make_any<ClusterTree_T>() );
                }

                ptoc(className()+"::GetClusterTree");
            }

            return std::any_cast<ClusterTree_T &>( this->GetPersistentCache(tag) );
        }

        virtual const BlockClusterTree_T & GetBlockClusterTree() const override
        {
            static std::string tag ( "BlockClusterTree" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                ptic(className()+"::GetBlockClusterTree");
                
                block_cluster_tree_settings.near_field_separation_parameter = adaptivity_settings.theta;
                block_cluster_tree_settings.near_field_intersection_parameter  = adaptivity_settings.intersection_theta;
                
                block_cluster_tree_settings.max_refinement = adaptivity_settings.max_refinement;
                
                this->SetPersistentCache( tag,
                    std::make_any<BlockClusterTree_T>(
                        GetClusterTree(),
                        GetClusterTree(),
                        block_cluster_tree_settings
                    )
                );
                
                ptoc(className()+"::GetBlockClusterTree");
            }
            
            return std::any_cast<BlockClusterTree_T &>( this->GetPersistentCache(tag) );
        }
        
        virtual const CollisionTree_T & GetCollisionTree() const override
        {
            static std::string tag ( "CollisionTree" );
            
            if( !this->InCacheQ( tag ) )
            {
                ptic(className()+"::GetCollisionTree");

                this->SetCache( tag,
                    std::make_any<CollisionTree_T>( GetClusterTree(), GetClusterTree() )
                );

                ptoc(className()+"::GetCollisionTree");
            }
            
            return std::any_cast<CollisionTree_T &>( this->GetCache(tag) );
        }
        
        const SparseBinaryMatrix_T & DerivativeAssembler() const override
        {
            
            static std::string tag ( "DerivativeAssembler" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                ptic(className()+"::DerivativeAssembler");
            
                auto A = SparseBinaryMatrix_T(
                    SimplexCount() * SIZE,
                    VertexCount(),
                    SimplexCount() * SIZE,
                    ThreadCount()
                );
                
                A.Outer().iota();
                A.Inner().Read(simplices.data());
                
                this->SetPersistentCache( tag, std::any( std::move(A.Transpose()) ) );
                
                ptoc(className()+"::DerivativeAssembler");
            }
            
            return std::any_cast<SparseBinaryMatrix_T &>( this->GetPersistentCache(tag) );
            
        } // DerivativeAssembler
        
        void Assemble_ClusterTree_Derivatives( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_Derivatives");
            
            Tensor3<Real,Int> buffer ( SimplexCount(), SIZE, AMB_DIM, static_cast<Real>(0) );
            
            GetClusterTree().CollectDerivatives();
            
            
            details.DNearToHulls( V_coords, simplices, GetClusterTree().PrimitiveDNearFieldData(), buffer, false );
            
            details.DFarToHulls ( V_coords, simplices, GetClusterTree().PrimitiveDFarFieldData(), buffer, true );

            DerivativeAssembler().template Dot<AMB_DIM>(
                static_cast<Real>(weight),   buffer.data(),
                static_cast<ExtReal>(addTo), output,
                AMB_DIM
            );

            ptoc(className()+"::Assemble_ClusterTree_Derivatives");
        }
        
        void Assemble_ClusterTree_Density( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_Density");
            
            GetClusterTree().CollectDensity( output, weight, addTo );
            
            ptoc(className()+"::Assemble_ClusterTree_Density");
        }
        
        void Assemble_ClusterTree_SimplexEnergies( mut<ExtReal> output, const ExtReal weight, bool addTo = false ) const override
        {
            ptic(className()+"::Assemble_ClusterTree_SimplexEnergies");
            
            GetClusterTree().CollectPrimitiveEnergies( output, weight, addTo );
            
            ptoc(className()+"::Assemble_ClusterTree_SimplexEnergies");
        }

//##############################################################################################
//      Obstacle
//##############################################################################################

        virtual void LoadObstacle( std::unique_ptr<Obstacle_T> obstacle_ ) override
        {
            static std::string tag ( "Obstacle" );
            
            std::shared_ptr<Obstacle_T> obstacle;
            
            // Input obstacle is moved.
            if( obstacle_->AmbDim() != AmbDim() )
            {
                eprint(className()+"::LoadObstacle: Attempted to load obstacle of ambient dimension "+ToString(obstacle_->AmbDim())+" into mesh of ambient dimension "+ToString(AmbDim())+". Setting obstacle to nullptr."
                );
                
                // We have to construct an empy SimplicialMesh because Obstacle_T is abstract.
                obstacle = std::make_shared<SimplicialMesh>();
            
            }
            else
            {
                obstacle = std::move(obstacle_);
            }

            this->SetPersistentCache(tag, std::any( obstacle ) );

        }
        
        const Obstacle_T & GetObstacle() const override
        {
            static std::string tag ("Obstacle");

            if( !this->InPersistentCacheQ(tag) )
            {
                wprint( ClassName()+"::GetObstacle: Obstacle not initialized.");
                
                // We have to construct an empty SimplicialMesh because Obstacle_T is abstract.
                this->SetPersistentCache(tag, std::any(std::make_shared<SimplicialMesh>()) );
            }
            
            return *std::any_cast<std::shared_ptr<Obstacle_T>>( this->GetPersistentCache(tag) );
        }
        
        virtual const ClusterTree_T & GetObstacleClusterTree() const override
        {
            return *dynamic_cast<const ClusterTree_T *>( & (GetObstacle().GetClusterTree()) );
        }
        
        virtual const ObstacleBlockClusterTree_T & GetObstacleBlockClusterTree() const override
        {
            static std::string tag ( "ObstacleBlockClusterTree" );
            
            if( !this->InPersistentCacheQ( tag ) )
            {
                ptic(className()+"::GetObstacleBlockClusterTree");
                
                this->SetPersistentCache( tag,
                    std::make_any<ObstacleBlockClusterTree_T>(
                        GetClusterTree(), GetObstacleClusterTree(), block_cluster_tree_settings
                    )
                );
                
                ptoc(className()+"::GetObstacleBlockClusterTree");

            }
            
            return std::any_cast<ObstacleBlockClusterTree_T &>( this->GetPersistentCache(tag) );
        }
        
        virtual const ObstacleCollisionTree_T & GetObstacleCollisionTree() const override
        {
            static std::string tag ( "ObstacleCollisionTree" );
            
            if( !this->InPersistentCacheQ(tag) )
            {
                ptic(className()+"::GetObstacleCollisionTree");
                
                this->SetPersistentCache( tag,
                    std::make_any<ObstacleCollisionTree_T>( GetClusterTree(), GetObstacleClusterTree() )
                );

                ptoc(className()+"::GetObstacleCollisionTree");
            }
            
            return std::any_cast<ObstacleCollisionTree_T &>( this->GetPersistentCache(tag) );
        }

        
//##############################################################################################
//      IO
//##############################################################################################
     
    public:
        
        virtual void WriteToFile( const std::string & file_name ) const override
        {
            ptic(ClassName()+"::WriteToFile");
            
            print("Writing mesh to file "+file_name+".");
            
            std::ofstream s (file_name);
            
            valprint("std::numeric_limits<Real>::digits",std::numeric_limits<Real>::digits10);
            
            s << std::setprecision( std::numeric_limits<Real>::digits );
            
            const Int vertex_count = VertexCount();
            const Int simplex_count = SimplexCount();
            
            ptr<Real> V = V_coords.data();
            ptr<Int>  S = simplices.data();
            
            s << "domain_dimension" << "\t" << DOM_DIM << "\n";
            s << "ambient_dimension" << "\t" << AMB_DIM << "\n";
            s << "vertex_count" << "\t" << vertex_count << "\n";
            s << "simplex_count" << "\t" << simplex_count << "\n";
            
            for( Int i = 0; i < vertex_count; ++i )
            {
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    s << V[AMB_DIM * i + k] << "\t";
                }
                s <<"\n";
            }
            
            for( Int i = 0; i < simplex_count; ++i )
            {
                for( Int k = 0; k < SIZE; ++k )
                {
                    s << S[SIZE * i + k] << "\t";
                }
                s <<"\n";
            }
            
            ptoc(ClassName()+"::WriteToFile");
        }
        
//##############################################################################################
//      Remesher
//##############################################################################################
        
    public:
        
        virtual std::unique_ptr<RemesherBase_T> CreateRemesher() override
        {
            ptic(ClassName()+"::CreateRemesher");
            
            Real * null = nullptr;
            
            Remesher_T * R = new Remesher_T(
                VertexCoordinates().data(), VertexCoordinates().Dimension(0),
                Simplices().data(),         Simplices().Dimension(0),
                null,                       0,
                ThreadCount()
            );
            
            ptoc(ClassName()+"::CreateRemesher");
            
            return std::unique_ptr<RemesherBase_T>(R);
        }
        
//##############################################################################################
//      Standard interface
//##############################################################################################
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    
    private:
  
        static std::string className()
        {
            return "SimplicialMesh<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
    };
} // namespace Repulsor


#include "SimplicialMesh/SimplicialMesh_Factory.hpp"
