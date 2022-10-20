//##############################################################################################
//      Factory-like code
//##############################################################################################

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal, typename ExtInt>
    std::unique_ptr<BASE> CONCAT(Make_,CLASS)(
         const ExtReal * const vertex_coords_, // vertex coordinates; size = vertex_count_ x amb_dim
         const long long vertex_count_,
         const int amb_dim,
         const bool vertex_coords_transpose, // whether to transpose input; set to false for row-major
         const ExtInt  * const simplices_, // simplices; size = simplex_count_ x simplex_size
         const long long simplex_count_,
         const int simplex_size,
         const bool simplices_transpose, // whether to transpose input; set to false for row-major
         const long long thread_count_ = 1
     )
    {
        const int dom_dim = simplex_size-1;
        ptic("Make_"+TO_STD_STRING(CLASS)+" ("+ToString(dom_dim)+","+ToString(amb_dim)+")");
        
        BASE * r = nullptr;
        
        switch(dom_dim)
        {
    #if ENABLE_CURVES == 1
            case 1:
            {
                switch(amb_dim)
                {
    #if ENABLE_2D == 1
                    case 2:
                    {
                        r = new CLASS<1,2,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
    #if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<1,3,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
    #if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<1,4,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
                    default:
                    {
                        eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported for domain dimension "+ToString(dom_dim)+". Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
    #endif
    #if ENABLE_SURFACES == 1
            case 2:
            {
                switch(amb_dim)
                {
    #if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<2,3,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
    #if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<2,4,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
                    default:
                    {
                        eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported for domain dimension "+ToString(dom_dim)+". Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
    #endif
    #if ENABLE_POINTCLOUDS == 1
            case 0:
            {
                switch(amb_dim)
                {
    #if ENABLE_1D == 1
                    case 1:
                    {
                        r = new CLASS<0,1,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
    #if ENABLE_2D == 1
                    case 2:
                    {
                        r = new CLASS<0,2,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
    #if ENABLE_3D == 1
                    case 3:
                    {
                        r = new CLASS<0,3,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
    #if ENABLE_4D == 1
                    case 4:
                    {
                        r = new CLASS<0,4,Real,Int,SReal,ExtReal>(
                            vertex_coords_, vertex_count_,  vertex_coords_transpose,
                            simplices_,     simplex_count_, simplices_transpose,
                            thread_count_
                        );
                        break;
                    }
    #endif
                    default:
                    {
                        eprint("Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported for domain dimension "+ToString(dom_dim)+". Returning nullptr.");
                        r = nullptr;
                        break;
                    }
                }
                break;
            }
    #endif
            default:
            {
                eprint("Make_"+TO_STD_STRING(CLASS)+": domain dimension "+ToString(dom_dim)+" not supported. Returning nullptr.");
                r = nullptr;
                break;
            }
                
        }
        
        ptoc("Make_"+TO_STD_STRING(CLASS)+" ("+ToString(dom_dim)+","+ToString(amb_dim)+")");
        
        return std::unique_ptr<BASE>(r);
    }

    template<typename Real, typename Int, typename SReal, typename ExtReal, typename ExtInt>
    std::unique_ptr<BASE> CONCAT(Make_,CLASS)(
        const ExtReal * const V_coords_, // vertex coordinates; size = vertex_count_ x amb_dim
        const long long vertex_count_,
        const int amb_dim,
        const ExtInt  * const simplices_, // simplices; size = simplex_count_ x simplex_size
        const long long simplex_count_,
        const int simplex_size,
        const long long thread_count_ = 1
    )
    {
        // For backward-compatibility.
        // Assuming that V_coords_ and simplices_ are both in row-major format.
        return CONCAT(Make_,CLASS)<Real,Int,SReal,ExtReal, ExtInt>(
            V_coords_,  vertex_count_,  amb_dim,      false,
            simplices_, simplex_count_, simplex_size, false,
            thread_count_
        );
    }


    template<typename Real, typename Int, typename SReal, typename ExtReal, typename ExtInt>
    std::unique_ptr<BASE> CONCAT3(Make_,CLASS,_FromFile)( const std::string & file_name, const Int thread_count )
    {
        ptic("Make_"+TO_STD_STRING(CLASS)+"FromFile");
        
        print("Reading mesh from file "+file_name+".");
        
        std::ifstream s (file_name);
        
        std::string str;
        int amb_dim;
        int dom_dim;
        Int vertex_count;
        Int simplex_count;
        s >> str;
        
        //        print(str);
        
        s >> dom_dim;
        
        valprint("dom_dim",dom_dim);
        
        s >> str;
        
        //        print(str);
        
        s >> amb_dim;
        
        valprint("amb_dim",amb_dim);
        
        
        s >> str;
        
        //        print(str);
        
        s >> vertex_count;
        
        valprint("vertex_count",vertex_count);
        
        s >> str;
        
        //        print(str);
        
        s >> simplex_count;
        
        valprint("simplex_count",simplex_count);
        
        const int simplex_size = dom_dim+1;
        
        valprint("simplex_size",simplex_size);
        
        Tensor2<ExtReal,Int> coords    (vertex_count, amb_dim);
        Tensor2<ExtInt,Int>  simplices (simplex_count,simplex_size);
        
        ExtReal * restrict const V = coords.data();
        ExtInt  * restrict const S = simplices.data();
        
        
        for( Int i = 0; i < vertex_count; ++i )
        {
            for( Int k = 0; k < amb_dim; ++k )
            {
                s >> V[amb_dim * i + k];
            }
        }
        
        for( Int i = 0; i < simplex_count; ++i )
        {
            for( Int k = 0; k < simplex_size; ++k )
            {
                s >> S[simplex_size * i + k];
            }
        }
        
        
        auto M = CONCAT(Make_,CLASS)<Real,Int,SReal,ExtReal,ExtInt>(
            V, vertex_count,  amb_dim,      false,
            S, simplex_count, simplex_size, false,
            thread_count
        );

        ptoc("Make_"+TO_STD_STRING(CLASS)+"FromFile");
        
        return M;
    }
    
} // namespace Repulsor

