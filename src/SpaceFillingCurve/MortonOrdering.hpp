
public:


//#######################################################################################
//####                         MortonOrdering (hybrid sort)                         #####
//#######################################################################################





    template<typename Real, typename ExtInt>
    void MortonOrdering( cptr<Real> X, mptr<ExtInt> perm, const Int point_count )
    {
        ptic(ClassName()+"::MortonOrdering");

        Tensor1<Morton_T,Int> M = CoordToMorton( X, point_count );
        
        ptic("Ordering");
        Ordering( M.data(), perm, point_count );
        ptoc("Ordering");
        
        ptoc(ClassName()+"::MortonOrdering");
    }


    template<typename Real>
    Tensor1<Int,Int> MortonOrdering( cptr<Real> X, const Int point_count )
    {
        Tensor1<Int,Int> perm ( point_count );
        
        MortonOrdering( X, perm.data(), point_count );
        
        return perm;
    }

//
//    template<typename Real>
//    std::pair< Tensor1<Morton_T,Int>, Tensor1<Int,Int> > OrderedMortonCodes(
//        cptr<Real> X, const Int point_count
//    )
//    {
//        ptic(ClassName()+"::OrderedMortonCodes");
//
//        Tensor1<Morton_T,Int> code = CoordToMorton( X, point_count, thread_count );
//
//        Tensor1<Int,Int> perm = Sort( code );
//
//        ptoc(ClassName()+"::OrderedMortonCodes");
//        return std::tie( code, perm );
//    }
//
//
//    Tensor1<Int,Int> Sort( mref<Tensor1<Morton_T,Int>> code_0 )
//    {
//        ptic(ClassName()+"::Sort");
//        constexpr Int pass_count = bit_count / 16;
//
//        const Int point_count = code_0.Dimension(0);
//
//        Tensor1<Int,Int> perm ( point_count );
//
//        Tensor1<Int,Int> bin_ptr( bin_count + 1 );
//        bin_ptr[0] = 0;
//
//        // Using parallel counting sort for the most significant 16 bits to create 65536 buckets.
//        {
//            Tensor1<Morton_T,Int> code_1 ( point_count );
//
//            JobPointers<Int> job_ptr ( point_count, thread_count );
//
//            Tensor2<Int, Int> counters ( thread_count, bin_count, static_cast<Int>(0) );
//
//            const Int shift = 16 * (pass_count-1);
//
//            const UInt mask = static_cast<UInt>(bin_count-1) << shift;
//
//            ParallelDo(
//                [&]( const Int thread )
//                {
//                    cptr<Morton_T> c_0 = code_0.data();
//                    mptr<Int>      ctr = counters.data(thread);
//
//                    const Int i_begin = job_ptr[thread  ];
//                    const Int i_end   = job_ptr[thread+1];
//
//                    for( Int i = i_begin; i < i_end; ++ i )
//                    {
//                        const Morton_T code = c_0[i];
//
//                        const Int key = static_cast<Int>( ( code[0] & mask ) >> shift );
//
//                        ++ctr[key];
//                    }
//                },
//                thread_count
//            );
//
//            // Accumulate counters.
//            AccumulateAssemblyCounters_Parallel( counters );
//
//            for( Int i = 0; i < bin_count; ++i )
//            {
//                bin_ptr[i+1] = counters(thread_count-1,i);
//            }
//
//            ParallelDo(
//                [&]( const Int thread )
//                {
//                    cptr<Morton_T> c_0 = code_0.data();
//                    mptr<Morton_T> c_1 = code_1.data();
//                    mptr<Int   > p   = perm.data();
//                    mptr<Int   > ctr = counters.data(thread);
//
//                    const Int i_begin = job_ptr[thread  ];
//                    const Int i_end   = job_ptr[thread+1];
//
//                    for( Int i = i_end; i-->i_begin; )
//                    {
//                        const Morton_T code = c_0[i];
//                        const Int key = static_cast<Int>( ( code[0] & mask ) >> shift );
//
//                        const Int j = --ctr[key];
//
//                        c_1[j] = code;
//                        p  [j] = i;
//                    }
//                },
//                thread_count
//            );
//
//            std::swap(code_0,code_1);
//        }
//
//        // quicksort
//        {
//            JobPointers<Int> job_ptr ( bin_count, bin_ptr.data(), thread_count, false );
//
//            // Sort the buckets in parallel with quicksort.
//            ParallelDo(
//                [&]( const Int thread )
//                {
//                    TwoArrayQuickSort<Morton_T, Int, Int> Q;
//
//                    mptr<Morton_T> c = code_0.data();
//                    mptr<Int   > p = perm.data();
//
//                    const Int bin_begin = job_ptr[thread  ];
//                    const Int bin_end   = job_ptr[thread+1];
//
//                    for( Int bin = bin_begin; bin < bin_end; ++bin )
//                    {
//                        const Int i_begin = bin_ptr[bin  ];
//                        const Int i_end   = bin_ptr[bin+1];
//
//                        Q( &c[i_begin], &p[i_begin], i_end - i_begin );
//                    }
//                },
//                thread_count
//            );
//        }
//
//        ptoc(ClassName()+"::Sort");
//        return perm;
//    }

