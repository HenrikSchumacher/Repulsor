public:

    static constexpr Int SIZE = 1 + (DOM_DIM+2) * AMB_DIM;

    using Level_T  = int32_t;
    using Child_T  = int32_t;
    using Column_T = uint64_t;
//    using Column_T = int32_t;
 
protected:

    static constexpr SReal zero = Scalar::Zero<SReal>;
    static constexpr SReal one  = Scalar::One<SReal>;
    static constexpr SReal half = Scalar::Half<SReal>;
    static constexpr SReal two  = Scalar::Two<SReal>;
    static constexpr SReal nth  = Scalar::Inv<SReal>(DOM_DIM+1);

    SReal root_serialized    [SIZE] = {};
    SReal simplex_serialized [SIZE] = {};
    
    SReal corners [DOM_DIM+1][DOM_DIM+1] = {};
    SReal center  [DOM_DIM+1] = {};
    
    Polytope<DOM_DIM+1,AMB_DIM,Real,Int,SReal,SReal,Int> P;
    
    Int      simplex_id      = -1;
    Level_T  level           =  0;
    Column_T column          =  0;
    Child_T  child_id        =  0;
    Child_T  former_child_id = -1;
    
    SReal scale  = one;
    SReal weight = one;
    
    bool current_simplex_computed = false;
    
public:
    
    SimplexHierarchy()
    {
        P.SetPointer(&simplex_serialized[0]);
    }

    virtual ~SimplexHierarchy() = default;
    
    bool Check()
    {
        if( !(0 <= simplex_id && simplex_id <= std::numeric_limits<Int>::max()) )
        {
            eprint("Overflow of simplex_id detected.");
            return false;
        };
        
        if( !(0 <= column && column <= std::numeric_limits<Column_T>::max()) )
        {
            eprint("Overflow of column detected.");
            return false;
        };
        
        if( !(0 <= child_id && child_id <= ChildCount()) )
        {
            eprint("Overflow of child_id detected.");
            return false;
        };
        
        if( !(0 <= former_child_id && former_child_id <= ChildCount()) )
        {
            eprint("Overflow of former_child_id detected.");
            return false;
        };
        
        return true;
    }


    constexpr Int VertexCount() const
    {
        return DOM_DIM+1;
    }

    constexpr Int DomDim() const
    {
        return DOM_DIM;
    }
    
    constexpr Int AmbDim() const
    {
        return AMB_DIM;
    }
    

    void RequireSimplex( const SReal * input_serialized, const Int k )
    {
//        if( simplex_id != k )
//        {
        
            current_simplex_computed = false;
            simplex_id = k;
            level  = 0;
            column = 0;
            child_id = 0;
            former_child_id = -1;
            
            scale  = one;
            weight = one;
            
            P.SetPointer(&root_serialized[0]);
            P.Read( input_serialized,k);
            P.SetPointer(&simplex_serialized[0]);
            P.Read( input_serialized,k);
        
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                center[i] = nth;
            }
        
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                for( Int j = 0; j < DOM_DIM+1; ++j )
                {
                    corners[i][j] = static_cast<SReal>(i == j);
                }
            }

//        }
    }
    
    void RequireSubsimplex()
    {
        if( !current_simplex_computed )
        {
            // root_serialized[0] is the squared radius!
            simplex_serialized[0] = scale * scale * root_serialized[0];
            
            cptr<SReal> root    =    &root_serialized[1 + AMB_DIM];
            mptr<SReal> c       = &simplex_serialized[1          ];
            mptr<SReal> current = &simplex_serialized[1 + AMB_DIM];

            for( Int k = 0; k < AMB_DIM; ++k )
            {
                c[k] = center[0] * root[AMB_DIM * 0 + k];
            }
        
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                for( Int i = 1; i < DOM_DIM+1; ++i )
                {
                    c[k] += center[i] * root[ AMB_DIM * i + k];
                }
            }
        
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                for( Int i = 0; i < DOM_DIM+1; ++i )
                {
                    current[ AMB_DIM * i + k] = corners[i][0] * root[ AMB_DIM * 0 + k];
                    
                    for( Int j = 1; j < DOM_DIM+1; ++j )
                    {
                        current[ AMB_DIM * i + k] += corners[i][j] * root[ AMB_DIM * j + k];
                    }
                }
            }
            
            P.SetPointer( &simplex_serialized[0] );
            
            current_simplex_computed = true;
        }
    }
    


    void RequireSimplexFromCoordinates( const SReal * coords, const Int k )
    {
    //        if( simplex_id != k )
    //        {
            simplex_id = k;
            level  = 0;
            column = 0;
            
            scale  = Scalar::One<SReal>;
            weight = Scalar::One<SReal>;

            P.SetPointer(&root_serialized[0]);

            P.FromCoordinates(coords,k);

            P.SetPointer(&simplex_serialized[0]);

            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                center[i] = Scalar::Inv<SReal>(DOM_DIM+1);
            
                for( Int j = 0; j < DOM_DIM+1; ++j )
                {
                    corners[i][j] = static_cast<SReal>(i == j);
                }
            }
            
            current_simplex_computed = false;
            child_id = 0;
            former_child_id = -1;
    //        }
    }

    const Polytope<DOM_DIM+1,AMB_DIM,Real,Int,SReal,SReal,Int> & SimplexPrototype()
    {
        RequireSubsimplex();
        return P;
    }
    
    SReal Scale() const
    {
        return scale;
    }
    
    SReal Weight() const
    {
        return weight;
    }
    
    Level_T Level() const
    {
        return level;
    }
    
    Column_T Column() const
    {
        return column;
    }

    Int SimplexID() const
    {
        return simplex_id;
    }

    Child_T ChildID() const
    {
        return child_id;
    }
    Child_T FormerChildID() const
    {
        return former_child_id;
    }

    const SReal * Center() const
    {
        return &center[0];
    }

    const SReal * BarycentricCoordinates() const
    {
        return &corners[0][0];
    }

    const SReal * SimplexSerialized() const
    {
        return simplex_serialized;
    }

    std::string CornerString() const
    {
        std::stringstream s;
        
        s << "{ ";
        
        s << "{ " << corners[0][0];

        for( Int j = 1; j < DOM_DIM+1; ++j )
        {
            s << ", " << corners[0][j];
        }
        s << " }";
        

        for( Int i = 1; i < DOM_DIM+1; ++i )
        {
            s << ", " << "{ " << corners[i][0];
            
            for( Int j = 1; j < DOM_DIM+1; ++j )
            {
                s << ", " << corners[i][j];
            }
            s << " }";
        }
        
        s << " }";
        return s.str();
    }

    std::string CenterString() const
    {
        std::stringstream s;

        s << "{ " << center[0];
        
        for( Int j = 1; j < DOM_DIM+1; ++j )
        {
            s << ", " << center[j];
        }
        s << " }";
        return s.str();
    }

    std::string EmbeddedSimplexString() const
    {
        std::stringstream s;
        
        s << "{ { ";

        const SReal * X = simplex_serialized+1+AMB_DIM;

        s << X[0];

        for( Int j = 1; j < AMB_DIM; ++j )
        {
            s << ", " << X[(DOM_DIM+1)*0+j];
        }

        for( Int i = 1; i < DOM_DIM+1; ++i )
        {
            s << " }" << ", " << "{ ";

            s << X[(DOM_DIM+1)*i+0];

            for( Int j = 1; j < AMB_DIM; ++j )
            {
                s << ", " << X[(DOM_DIM+1)*i+j];
            }
        }

        s << " } }";
        
        return s.str();
    }

    
#ifdef LTEMPLATE_H

    mma::TensorRef<mreal> BarycentricCoordinates_TensorRef()
    {
        auto A = mma::makeMatrix<mreal>( DOM_DIM+1, DOM_DIM+1 );
        
        mreal * a = A.data();
        
        for( Int i = 0; i < DOM_DIM+1; ++i )
        {
            for( Int j = 0; j < DOM_DIM+1; ++j )
            {
                a[ DOM_DIM+1 * i + j ] = static_cast<mreal>(corners[i][j]);
            }
        }
        
        return A;
    }

    mma::TensorRef<mreal> VertexCoordinates_TensorRef()
    {
        auto A = mma::makeMatrix<mreal>( DOM_DIM+1, AMB_DIM );
        
        mreal * a = A.data();
        
        RequireSubsimplex();
        
        mptr<SReal> q = &simplex_serialized[0] + 1 + AMB_DIM;
        
        for( Int i = 0; i < DOM_DIM+1 * AMB_DIM; ++i )
        {
            a[i] = static_cast<mreal>(q[i]);
        }
        
        return A;
    }

    mma::TensorRef<mreal> RootSerialized_TensorRef()
    {
        auto A = mma::makeVector<mreal>( P.Size() );

        P.SetPointer(&root_serialized[0]);
        
        P.Write(A.data());
        
        P.SetPointer(&simplex_serialized[0]);
        
        return A;
    }
    
    mma::TensorRef<mreal> SimplexSerialized_TensorRef()
    {
        auto A = mma::makeVector<mreal>( P.Size() );
    
        RequireSubsimplex();
        
        P.SetPointer(&simplex_serialized[0]);
        
        P.Write(A.data());
        
        return A;
    }
    
#endif
    
public:
    
    std::string ClassName() const
    {
        return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
    }
