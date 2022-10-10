#pragma once

#define CLASS BlockClusterTree
#define BASE  BlockClusterTreeBase<Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);
        
    public:
        
        using BlockClusterTreeBase_T = BASE;
        
        using Setting_T = typename BASE::Setting_T;

        using Inter_T           = typename BASE::Inter_T;
        using VeryNear_T        = typename BASE::VeryNear_T;
        using Near_T            = typename BASE::Near_T;
        using Far_T             = typename BASE::Far_T;
        
        using ClusterTree_T     = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;
        
        using Primitive_T       = typename ClusterTree_T::Primitive_T;
        using BoundingVolume_T  = typename ClusterTree_T::BoundingVolume_T;
        
        using GJK_T             = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
        using BASE::AmbDim;
        using BASE::ThreadCount;
        using BASE::SeparationParameter;
        using BASE::AdaptivityParameter;
        using BASE::Settings;
        
    public:

        // In order to prevent GetS() and GetT() shooting a segfault, we have to initialize S and T here. This is the only case in which CLASS owns these raw pointers.
        
        virtual ~CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_,
            Setting_T settings_ = Setting_T()
        )
        :   S(S_)
        ,   T(T_)
        ,   settings(settings_)
        ,   near_theta2( pow(settings_.near_field_separation_parameter,2) )
        ,    far_theta2( pow(settings_. far_field_separation_parameter,2) )
        ,   thread_count( std::min(S_.ThreadCount(), T_.ThreadCount()) )
        ,   is_symmetric( std::addressof(S_) == std::addressof(T_) )
        {
            ptic(className());

            if( std::min( S.PrimitiveCount(), T.PrimitiveCount()) > 0 )
            {
                ComputeBlocks();
            }

            ptoc(className());
        } // Constructor
        
    protected:

        static constexpr Int zero = static_cast<Int>(0);
        
        // Not very elegant to use raw pointers here, but maybe acceptable due to constness.
        const ClusterTree_T & S; // "left"  BVH (output side of matrix-vector multiplication)
        const ClusterTree_T & T; // "right" BVH (input  side of matrix-vector multiplication)
        
        const Setting_T settings;
        
        mutable Inter_T    inter;
        mutable VeryNear_T verynear;
        mutable Near_T     near;
        mutable Far_T      far;
        
        const SReal near_theta2 = static_cast<SReal>(10);
        const SReal  far_theta2 = static_cast<SReal>(0.25);
        
        const Int thread_count = static_cast<Int>(1);
        const bool is_symmetric;
        
        mutable bool blocks_initialized = false;
        
        // Containers for parallel aggregation of {i,j}-index pairs
        mutable std::vector<std::vector<Int>> inter_i;
        mutable std::vector<std::vector<Int>> inter_j;
        mutable std::vector<std::vector<Int>> verynear_i;
        mutable std::vector<std::vector<Int>> verynear_j;
        mutable std::vector<std::vector<Int>> near_i;
        mutable std::vector<std::vector<Int>> near_j;
        mutable std::vector<std::vector<Int>> far_i;
        mutable std::vector<std::vector<Int>> far_j;
        
    public:

        Int AmbDim() const override
        {
            return AMB_DIM;
        }
        
        Int ThreadCount() const override
        {
            return thread_count;
        }
        
        Real SeparationParameter() const override
        {
            return sqrt(far_theta2);
        }
                        
        Real AdaptivityParameter() const override
        {
            return sqrt(near_theta2);
        }
        
        bool IsSymmetric() const override
        {
            return is_symmetric;
        }

        Int PrimitiveIntersectionCount() const override
        {
            return inter.NonzeroCount();
        }
        
        Int VeryNearFieldInteractionCount() const override
        {
            return verynear.NonzeroCount();
        }
        
        Int NearFieldInteractionCount() const override
        {
            return near.NonzeroCount();
        }
        
        Int FarFieldInteractionCount() const override
        {
            return far.NonzeroCount();
        }

        
        virtual const Inter_T & PrimitiveIntersectionMatrix() const override
        {
            return inter;
        }
        
        const VeryNear_T & VeryNear() const override
        {
            return verynear;
        }
        
        const Near_T & Near() const override
        {
            return near;
        }
        
        const Far_T & Far() const override
        {
            return far;
        }

        
        const ClusterTree_T & GetS() const override
        {
            return S;
        }
        
        const ClusterTree_T & GetT() const override
        {
            return T;
        }

        const Setting_T & Settings() const override
        {
            return settings;
        }
        
        virtual std::string Stats() const override
        {
            std::stringstream s;
            
            s
            << "\n==== "+className()+" Stats ====" << "\n\n"
            << " AmbDim()                     = " <<  AmbDim() << "\n"
            << " SeparationParameter()        = " <<  SeparationParameter() << "\n"
            << " AdaptivityParameter()        = " <<  AdaptivityParameter() << "\n"
            << " ThreadCount()                = " <<  ThreadCount() << "\n"
            << "\n"
            << "Tree S \n"
            << S.Stats() << "\n"
            << "Tree T \n"
            << T.Stats() << "\n"
            << " PrimitiveIntersectionCount()    = " <<  PrimitiveIntersectionCount() << "\n"
            << " VeryNearFieldInteractionCount() = " <<  VeryNearFieldInteractionCount() << "\n"
            << " NearFieldInteractionCount()     = " <<  NearFieldInteractionCount() << "\n"
            << " FarFieldInteractionCount()      = " <<  FarFieldInteractionCount() << "\n"

            << "\n---- bool data ----" << "\n"

            << " IsSymmetric()               = " <<  IsSymmetric() << "\n"

            << "\n==== "+className()+" Stats ====\n" << std::endl;

            return s.str();
        }
        

        
    
//#################################################################################################
//      Initialization
//#################################################################################################

    private:
        
        void ComputeBlocks()
        {
            if( blocks_initialized )
            {
                return;
            }
            
            ptic(className()+"::ComputeBlocks");

            inter_i    = std::vector<std::vector<Int>>(thread_count);
            inter_j    = std::vector<std::vector<Int>>(thread_count);
            verynear_i = std::vector<std::vector<Int>>(thread_count);
            verynear_j = std::vector<std::vector<Int>>(thread_count);
            near_i     = std::vector<std::vector<Int>>(thread_count);
            near_j     = std::vector<std::vector<Int>>(thread_count);
            far_i      = std::vector<std::vector<Int>>(thread_count);
            far_j      = std::vector<std::vector<Int>>(thread_count);
            
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                // TODO: Find a better initial guess for the number of admissable and inadmissable blocks.
                const Int expected = static_cast<Int>(10) * ( S.PrimitiveCount() + T.PrimitiveCount() );

                   inter_i[thread] = std::vector<Int>();
                   inter_j[thread] = std::vector<Int>();
                verynear_i[thread] = std::vector<Int>();
                verynear_j[thread] = std::vector<Int>();
                    near_i[thread] = std::vector<Int>();
                    near_j[thread] = std::vector<Int>();
                     far_i[thread] = std::vector<Int>();
                     far_j[thread] = std::vector<Int>();
                
                near_i[thread].reserve(expected);
                near_j[thread].reserve(expected);
                 far_i[thread].reserve(expected);
                 far_j[thread].reserve(expected);
            }
            
//            switch( settings.block_split_method )
//            {
//                case BlockSplitMethod::Parallel:
//                {
//                    Split_Parallel( static_cast<Int>(4) * thread_count * thread_count );
//                    break;
//                }
//                case BlockSplitMethod::Sequential:
//                {
//                    Split_Sequential_DFS( zero, zero );
//                    break;
//                }
//                default:
//                {
//                    Split_Parallel( static_cast<Int>(4) * thread_count * thread_count );
//                }
//            }
            
            // TODO: Reactivate switch here.
            Split_Sequential_DFS( zero, zero );

            ptic(className()+"  Far field interaction data");
            
            far = Far_T( far_i, far_j, S.ClusterCount(), T.ClusterCount(),
                    thread_count, false, is_symmetric );
            
            DUMP(far.Stats());
            
            ptoc(className()+"  Far field interaction data");

            ptic(className()+"  Near field interaction data");
            
            near = Near_T( near_i, near_j, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            DUMP(near.Stats());
            
            ptoc(className()+"  Near field interaction data");
            
            ptic(className()+"  Very near field interaction data");
            
            verynear = VeryNear_T( verynear_i, verynear_j, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            DUMP(verynear.Stats());
            
            ptoc(className()+"  Very near field interaction data");
            
            ptic(className()+"  Primitive intersection data");
            
            inter = Inter_T( inter_i, inter_j, S.PrimitiveCount(), T.PrimitiveCount(),
                thread_count, false, is_symmetric );
            
            DUMP(inter.Stats());
            
            ptoc(className()+"  Primitive intersection data");
            
            blocks_initialized = true;
            
            // Free memory that is no longer used.
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                   inter_i[thread] = std::vector<Int>();
                   inter_j[thread] = std::vector<Int>();
                verynear_i[thread] = std::vector<Int>();
                verynear_j[thread] = std::vector<Int>();
                    near_i[thread] = std::vector<Int>();
                    near_j[thread] = std::vector<Int>();
                     far_i[thread] = std::vector<Int>();
                     far_j[thread] = std::vector<Int>();
            }
            
            ptoc(className()+"::ComputeBlocks");
                
        }; // ComputeBlocks

        

//#################################################################################################
//      Initialization of matrices
//#################################################################################################

    public:

        void Split_Sequential_DFS( const Int i0, const Int j0 )
        {
    //        ptic(className()+"::Split_Sequential_DFS");
            
            const Int thread = omp_get_thread_num();
            
            const Int max_depth = 128;
            Int i_stack[128] = {};
            Int j_stack[128] = {};
//            Tensor1<Int,Int> i_stack (max_depth);
//            Tensor1<Int,Int> j_stack (max_depth);
            Int stack_ptr = 0;
            i_stack[0] = i0;
            j_stack[0] = j0;

            const SparseBinaryMatrixCSR<Int> & A = S.PrimitiveAdjacencyMatrix();
            const Real intersection_theta2 = AdaptivityParameter() * AdaptivityParameter();
            
            std::vector<Int> &    inter_idx =    inter_i[thread];
            std::vector<Int> &    inter_jdx =    inter_j[thread];
            std::vector<Int> & verynear_idx = verynear_i[thread];
            std::vector<Int> & verynear_jdx = verynear_j[thread];
            std::vector<Int> &     near_idx =     near_i[thread];
            std::vector<Int> &     near_jdx =     near_j[thread];
            std::vector<Int> &      far_idx =      far_i[thread];
            std::vector<Int> &      far_jdx =      far_j[thread];

            
            std::shared_ptr<BoundingVolume_T> S_C_proto_ptr = S.ClusterPrototype().Clone();
            std::shared_ptr<BoundingVolume_T> T_C_proto_ptr = T.ClusterPrototype().Clone();
            std::shared_ptr<Primitive_T>      S_P_proto_ptr = S.PrimitivePrototype().Clone();
            std::shared_ptr<Primitive_T>      T_P_proto_ptr = T.PrimitivePrototype().Clone();
            
            BoundingVolume_T & S_C_proto = *S_C_proto_ptr;
            BoundingVolume_T & T_C_proto = *T_C_proto_ptr;
            Primitive_T      & S_P_proto = *S_P_proto_ptr;
            Primitive_T      & T_P_proto = *T_P_proto_ptr;
            GJK_T G;
            
            const Int   * restrict const S_C_left       = S.ClusterLeft().data();
            const Int   * restrict const S_C_right      = S.ClusterRight().data();
            const Int   * restrict const S_C_begin      = S.ClusterBegin().data();
            const Int   * restrict const S_C_end        = S.ClusterEnd().data();
            
            const Int   * restrict const T_C_left       = T.ClusterLeft().data();
            const Int   * restrict const T_C_right      = T.ClusterRight().data();
            const Int   * restrict const T_C_begin      = T.ClusterBegin().data();
            const Int   * restrict const T_C_end        = T.ClusterEnd().data();
            
                  SReal * restrict const S_C_serialized = S.ClusterSerializedData().data();
                  SReal * restrict const T_C_serialized = T.ClusterSerializedData().data();
        
                  SReal * restrict const S_P_serialized = S.PrimitiveSerializedData().data();
                  SReal * restrict const T_P_serialized = T.PrimitiveSerializedData().data();
            
            while( (zero <= stack_ptr) && (stack_ptr < max_depth) )
            {
                const Int i = i_stack[stack_ptr];
                const Int j = j_stack[stack_ptr];
                stack_ptr--;
                
                S_C_proto.SetPointer(S_C_serialized,i);
                T_C_proto.SetPointer(T_C_serialized,j);
                
                const bool separatedQ = (!( IsSymmetric() && (i == j) )) && G.MultipoleAcceptanceCriterion( S_C_proto, T_C_proto, far_theta2 );
                
                
                if( !separatedQ )
                {
                    const Int left_i = S_C_left[i];
                    const Int left_j = T_C_left[j];

                    // Warning: This assumes that either both children are defined or empty.
                    if( left_i >= zero || left_j >= zero )
                    {

                        const Int right_i = S_C_right[i];
                        const Int right_j = T_C_right[j];
                        
                        const SReal score_i = (left_i>=zero) * S_C_proto.SquaredRadius();
                        const SReal score_j = (left_j>=zero) * T_C_proto.SquaredRadius();

                        if( score_i == score_j && score_i > static_cast<SReal>(0) /*&& score_j > static_cast<SReal>(0)*/ )
                        {
                            // tie breaker: split both clusters

                            if ( is_symmetric )
                            {
                                if( i == j )
                                {
                                    //  Creating 3 blockcluster children, since there is one block is just the mirror of another one.
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = left_i;
                                    j_stack[stack_ptr] = right_j;
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = right_i;
                                    j_stack[stack_ptr] = right_j;
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = left_i;
                                    j_stack[stack_ptr] = left_j;
                                }
                                else
                                {
                                    // This is a very seldom case; still required to preserve symmetry.
                                    // This happens only if i and j represent_diffent clusters with same radii.
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = right_i;
                                    j_stack[stack_ptr] = right_j;
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = left_i;
                                    j_stack[stack_ptr] = right_j;
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = right_i;
                                    j_stack[stack_ptr] = left_j;
                                    
                                    ++stack_ptr;
                                    i_stack[stack_ptr] = left_i;
                                    j_stack[stack_ptr] = left_j;
                                }
                            }
                            else
                            {
                                // Split both clusters
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = right_i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = right_i;
                                j_stack[stack_ptr] = left_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = left_j;
                            }
                        }
                        else
                        {
                            // split only larger cluster
                            if( score_i > score_j )
                            {
                                ++stack_ptr;
                                i_stack[stack_ptr] = right_i;
                                j_stack[stack_ptr] = j;
                                
                                //split cluster i
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = j;
                            }
                            else //score_i < score_j
                            {
                                //split cluster j
                                ++stack_ptr;
                                i_stack[stack_ptr] = i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = i;
                                j_stack[stack_ptr] = left_j;
                            }
                        }
                    }
                    else // left_i < zero && left_j < zero
                    {
                        // We know that i and j are leaf clusters and that they belong either to the near field, the very near field or contain intersecting primitives.
                        
                        // We have to go through all the primitive pairs to classify them.
                        
                        // For simplicity, we assume that split_threshold == 1, so that leaf clusters have only a single primitive contained in them:
                        
                        const Int ii = S_C_begin[i];
                        const Int jj = T_C_begin[j];

                        
                        if( is_symmetric )
                        {
                            if( ii != jj )
                            {
                                Int ii_;
                                Int jj_;
                                
                                if( ii < jj )
                                {
                                    ii_ = ii;
                                    jj_ = jj;
                                }
                                else
                                {
                                    ii_ = jj;
                                    jj_ = ii;
                                }
                                
                                const bool neighbor_found = A.FindNonzeroPosition(ii,jj) >= 0;
                                
                                S_P_proto.SetPointer(S_P_serialized,ii);
                                T_P_proto.SetPointer(T_P_serialized,jj);

                                const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
                                
                                if( admissable )
                                {
                                    near_idx.push_back(ii_);
                                    near_jdx.push_back(jj_);
                                }
                                else
                                {
                                    const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                                    
                                    if( intersecting )
                                    {
                                        inter_idx.push_back(ii_);
                                        inter_jdx.push_back(jj_);
                                    }
                                    else
                                    {
                                        verynear_idx.push_back(ii_);
                                        verynear_jdx.push_back(jj_);
                                    }
                                }
                            }
                        }
                        else // !is_symmetric
                        {
                            
                            S_P_proto.SetPointer(S_P_serialized,ii);
                            T_P_proto.SetPointer(T_P_serialized,jj);
                            
                            const bool admissable = G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
                            
                            if( admissable )
                            {
                                near_idx.push_back(ii);
                                near_jdx.push_back(jj);
                            }
                            else
                            {
                                const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                                
                                if( intersecting )
                                {
                                    inter_idx.push_back(ii);
                                    inter_jdx.push_back(jj);
                                }
                                else
                                {
                                    verynear_idx.push_back(ii);
                                    verynear_jdx.push_back(jj);
                                }
                            }
                        }
                    }
                }
                else // separatedQ
                {
                    //create far field leaf blockcluster
                    if( is_symmetric )
                    {
                        // For creating an upper triangle matrix we store the pair  { min(i,j), max(i,j) }.
                        if (i <= j)
                        {
                            far_idx.push_back(i);
                            far_jdx.push_back(j);
                        }
                        else
                        {
                            far_idx.push_back(j);
                            far_jdx.push_back(i);
                        }
                    }
                    else
                    {
                        // No symmetry exploited.
                        far_idx.push_back(i);
                        far_jdx.push_back(j);
                    }
                }
            }
            ptoc(className()+"::Split_Sequential_DFS");
            
        } // Split_Sequential_DFS


//        void Split_Parallel( const Int max_leaves )
//        {
//            ptic(className()+"::Split_Parallel");
//
//            std::deque<Int> i_queue;
//            std::deque<Int> j_queue;
//            i_queue.push_back(zero);
//            j_queue.push_back(zero);
//
//            std::vector<Int> &  sep_idx =  sep_i[0];
//            std::vector<Int> &  sep_jdx =  sep_j[0];
//            std::vector<Int> & mid_idx = mid_i[0];
//            std::vector<Int> & mid_jdx = mid_j[0];
//
//            BoundingVolume_T & P = *S_bv[0];
//            BoundingVolume_T & Q = *T_bv[0];
//            GJK_T            & G = *gjk [0];
//
//            const Int   * restrict const P_left       = S.ClusterLeft().data();
//            const Int   * restrict const P_right      = S.ClusterRight().data();
//
//            const Int   * restrict const Q_left       = T.ClusterLeft().data();
//            const Int   * restrict const Q_right      = T.ClusterRight().data();
//
//                  SReal * restrict const P_serialized = S.ClusterSerializedData().data();
//                  SReal * restrict const Q_serialized = T.ClusterSerializedData().data();
//
//            const Int   * restrict const S_lookup     = S.LeafClusterLookup().data();
//            const Int   * restrict const T_lookup     = T.LeafClusterLookup().data();
//
//            ptic(ClassName()+"::Split_Parallel: Singlethreaded BFS.");
//            while( !i_queue.empty() && ( static_cast<Int>(i_queue.size()) < max_leaves ) )
//            {
//                const Int i = i_queue.front();
//                const Int j = j_queue.front();
//                i_queue.pop_front();
//                j_queue.pop_front();
//
//                P.SetPointer( P_serialized, i );
//                Q.SetPointer( Q_serialized, j );
//
//                bool separatedQ = ( IsSymmetric() && (i == j) )
//                    ? false
//                    : G.MultipoleAcceptanceCriterion( P, Q, far_theta2 );
//
//
//                if( !separatedQ )
//                {
//                    const Int left_i = P_left[i];
//                    const Int left_j = Q_left[j];
//
//                    // Warning: This assumes that either both children are defined or empty.
//                    if( left_i >= zero || left_j >= zero )
//                    {
//
//                        const Int right_i = P_right[i];
//                        const Int right_j = Q_right[j];
//
//                        const SReal score_i = (left_i>=zero) * P.SquaredRadius();
//                        const SReal score_j = (left_j>=zero) * Q.SquaredRadius();
//
//                        if( score_i == score_j && score_i > static_cast<SReal>(0) /*&& score_j > static_cast<SReal>(0)*/ )
//                        {
//                            // tie breaker: split both clusters
//
//                            if( (exploit_symmetry) && (i == j) )
//                            {
//                                //  Creating 3 blockcluster children, since there is one block is just the mirror of another one.
//
//                                i_queue.push_back(left_i);
//                                j_queue.push_back(right_j);
//
//                                i_queue.push_back(right_i);
//                                j_queue.push_back(right_j);
//
//                                i_queue.push_back(left_i);
//                                j_queue.push_back(left_j);
//                            }
//                            else
//                            {
//                                // In case of exploit_symmetry !=0, this is a very seldom case; still requird to preserve symmetry.
//                                // This happens only if i and j represent_different clusters with same radii.
//
//                                i_queue.push_back(right_i);
//                                j_queue.push_back(right_j);
//
//                                i_queue.push_back(left_i);
//                                j_queue.push_back(right_j);
//
//                                i_queue.push_back(right_i);
//                                j_queue.push_back(left_j);
//
//                                i_queue.push_back(left_i);
//                                j_queue.push_back(left_j);
//                            }
//                        }
//                        else
//                        {
//                            // split only the larger cluster
//                            if (score_i > score_j)
//                            {
//                                i_queue.push_back(right_i);
//                                j_queue.push_back(j);
//
//                                //split cluster i
//                                i_queue.push_back(left_i);
//                                j_queue.push_back(j);
//                            }
//                            else //score_i < score_j
//                            {
//                                //split cluster j
//                                i_queue.push_back(i);
//                                j_queue.push_back(right_j);
//
//                                i_queue.push_back(i);
//                                j_queue.push_back(left_j);
//                            }
//                        }
//                    }
//                    else
//                    {
//                        // create nonsep leaf blockcluster
//
//                        // i and j must be leaves of each ClusterTree S and T, so we directly store their position in the list leaf_clusters. This is important for the sparse matrix generation.
//
//                        //            In know  this is a very deep branching. I optimized it a bit for the case exploit_symmetry != 0 and upper_triangular == 0, though. That seemed to work best in regard of the matrix-vector multiplication.
//
//                        if(exploit_symmetry)
//                        {
//                            if(!near_upper_triangular)
//                            {
//                                if(i != j)
//                                {
//                                    // Generate also the twin to get a full matrix.
//
//                                    // TODO: In the symmetric case two of these lookups are redundant.
//                                    mid_idx.push_back(S_lookup[i]);
//                                    mid_idx.push_back(S_lookup[j]);
//                                    mid_jdx.push_back(T_lookup[j]);
//                                    mid_jdx.push_back(T_lookup[i]);
//                                }
//                                else
//                                {
//                                    // This is a diagonal block; there is no twin to think about
//                                    mid_idx.push_back(T_lookup[i]);
//                                    mid_jdx.push_back(S_lookup[i]);
//                                }
//                            }
//                            else
//                            {
//                                // For creating an upper triangle matrix we store the pair  { min(i,j), max(i,j) }.
//                                if(i <= j)
//                                {
//                                    mid_idx.push_back(S_lookup[i]);
//                                    mid_jdx.push_back(T_lookup[j]);
//                                }
//                                else
//                                {
//                                    mid_idx.push_back(T_lookup[j]);
//                                    mid_jdx.push_back(S_lookup[i]);
//                                }
//                            }
//                        }
//                        else
//                        {
//                            // No symmetry exploited.
//                            mid_idx.push_back(S_lookup[i]);
//                            mid_jdx.push_back(T_lookup[j]);
//                        }
//                    }
//                }
//                else
//                {
//                    //create sep leaf blockcluster
//                    if(exploit_symmetry)
//                    {
//                        if(!far_upper_triangular)
//                        {
//                            // For creating an upper triangle matrix we store the pair  { min(i,j), max(i,j) }.
//                            if (i <= j)
//                            {
//                                sep_idx.push_back(i);
//                                sep_jdx.push_back(j);
//                            }
//                            else
//                            {
//                                sep_idx.push_back(j);
//                                sep_jdx.push_back(i);
//                            }
//                        }
//                        else
//                        {
//                            // Generate also the twin to get a full matrix
//                            sep_idx.push_back(i);
//                            sep_idx.push_back(j);
//                            sep_jdx.push_back(j);
//                            sep_jdx.push_back(i);
//                        }
//                    }
//                    else
//                    {
//                        // No symmetry exploited.
//                        sep_idx.push_back(i);
//                        sep_jdx.push_back(j);
//                    }
//                }
//            }
//            ptoc(ClassName()+"::Split_Parallel: Singlethreaded BFS.");
//
//            ptic(ClassName()+"::Split_Parallel: Parallel DFS.");
//            #pragma omp parallel for num_threads( thread_count ) schedule( dynamic )
//            for( Int k = 0; k < static_cast<Int>(i_queue.size()); ++k )
//            {
//                Split_Sequential_DFS(i_queue[k], j_queue[k]);
//            }
//            ptoc(ClassName()+"::Split_Parallel: Parallel DFS.");
//
//            ptoc(className()+"::Split_Parallel");
//
//        } // Split_Parallel

        
        std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS) + "<"+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
      
        
    };
    
} //namespace Repulsor

#undef BASE
#undef CLASS
