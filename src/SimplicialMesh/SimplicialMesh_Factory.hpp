#pragma once


namespace Repulsor
{
    
    
    // Some generic template for all the factories to come.
    template < class T, int ...N> class SimplicialMesh_Factory {};
    
    template<
        int MinDomDim_, int MaxDomDim_,
        int MinAmbDim_, int MaxAmbDim_,
        typename Real_, typename Int_, typename SReal_, typename ExtReal_
    >
    class SimplicialMesh_Factory< SimplicialMeshBase<Real_,Int_,SReal_,ExtReal_>, MinDomDim_, MaxDomDim_, MinAmbDim_, MaxAmbDim_ >
    {
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        
        using Base_T  = SimplicialMeshBase<Real,Int,SReal,ExtReal>;
        
        static constexpr Int MinAmbDim = std::max( Int(1), Int(MinAmbDim_) );
        static constexpr Int MaxAmbDim = std::max( Int(1), Int(MaxAmbDim_) );
        
        static constexpr Int MinDomDim = std::max( Int(0), Int(MinDomDim_) );
        static constexpr Int MaxDomDim = std::min( Int(MaxAmbDim-1), Int(MaxDomDim_) );
        
        SimplicialMesh_Factory() = default;
        
        ~SimplicialMesh_Factory() = default;

        template<typename ExtInt>
        std::unique_ptr<Base_T> Make(
            const ExtReal * const vertex_coords, // vertex coordinates; size = vertex_count_ x amb_dim
            const Int             vertex_count,
            const Int             amb_dim,
            const bool            vertex_coords_transpose, // whether to transpose input; set to false for row-major
            const ExtInt  * const simplices, // simplices; size = simplex_count_ x simplex_size
            const Int             simplex_count,
            const Int             simplex_size,
            const bool            simplices_transpose, // whether to transpose input; set to false for row-major
            const Int             thread_count = 1
        )
        {
            const Int dom_dim = simplex_size-1;
            
            if( (dom_dim < MinDomDim) || (dom_dim > MaxDomDim) )
            {
                eprint(ClassName()+"::Make: dom_dim is out of range.");
                return Error( dom_dim, amb_dim );
            }
            
            if( (amb_dim < MinAmbDim) || (amb_dim > MaxAmbDim) )
            {
                eprint(ClassName()+"::Make: amb_dim is out of range.");
                return Error( dom_dim, amb_dim );
            }
            
            if( dom_dim >= amb_dim )
            {
                return Error( dom_dim, amb_dim );
            }
            
            return make_1<MaxAmbDim>(
                vertex_coords, vertex_count, amb_dim, vertex_coords_transpose,
                simplices, simplex_count, dom_dim, simplices_transpose,
                thread_count
            );
        }
        
    private:
        
        template<Int AmbDim, typename ExtInt>
        std::unique_ptr<Base_T> make_1(
            const ExtReal * const v,
            const Int             v_cnt,
            const Int             amb_dim,
            const bool            v_transp,
            const ExtInt  * const s,
            const Int             s_cnt,
            const Int             dom_dim,
            const bool            s_transp,
            const Int             thread_count = 1
        )
        {
            if( amb_dim == AmbDim )
            {
                return make_2<std::min(MaxDomDim,AmbDim-1),AmbDim>(
                    v, v_cnt, amb_dim, v_transp, s, s_cnt, dom_dim, s_transp, thread_count
                );
            }
            else if constexpr ( AmbDim > MinAmbDim )
            {
                return make_1<AmbDim-1>(
                    v, v_cnt, amb_dim, v_transp, s, s_cnt, dom_dim, s_transp, thread_count
                );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
        
        template<Int DomDim, Int AmbDim, typename ExtInt>
        std::unique_ptr<Base_T> make_2(
            const ExtReal * const v,
            const Int             v_cnt,
            const Int             amb_dim,
            const bool            v_transp,
            const ExtInt  * const s,
            const Int             s_cnt,
            const Int             dom_dim,
            const bool            s_transp,
            const Int             thread_count = 1
        )
        {
            if( dom_dim == DomDim )
            {
                return std::unique_ptr<Base_T>(
                    new SimplicialMesh<DomDim,AmbDim,Real,Int,SReal,ExtReal>(
                        v, v_cnt, v_transp, s, s_cnt, s_transp, thread_count
                    )
                );
            }
            else if constexpr ( DomDim > MinDomDim )
            {
                return make_2<DomDim-1,AmbDim>(
                    v, v_cnt, amb_dim, v_transp, s, s_cnt, dom_dim, s_transp, thread_count
                );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
            
        std::unique_ptr<Base_T> Error( const Int dom_dim, const Int amb_dim )
        {
            eprint(ClassName()+" cannot create SimplicialMesh with domain dimension "+ToString(dom_dim)+" and  ambient dimension "+ToString(amb_dim)+".");
            
            return std::unique_ptr<Base_T>( nullptr );
        }
        
    public:
        
        template<typename Real, typename Int, typename SReal, typename ExtReal>
        std::unique_ptr<Base_T> Make_FromFile( const std::string & file_name, const Int thread_count )
        {
            ptic(ClassName()+"Make_FromFile");
            
            print("Reading mesh from file "+file_name+".");
            
            std::ifstream s (file_name);
            
            if( !s.good() )
            {
                eprint("Make_FromFile: File "+file_name+" could not be opened.");
                
                return nullptr;
            }
            
            std::string str;
            Int amb_dim;
            Int dom_dim;
            Int vertex_count;
            Int simplex_count;
            s >> str;
            s >> dom_dim;
            valprint("dom_dim",dom_dim);
            s >> str;
            s >> amb_dim;
            valprint("amb_dim",amb_dim);
            s >> str;
            s >> vertex_count;
            valprint("vertex_count",vertex_count);
            s >> str;
            s >> simplex_count;
            valprint("simplex_count",simplex_count);
            
            const Int simplex_size = dom_dim+1;
            
            valprint("simplex_size",simplex_size);
            
            Tensor2<ExtReal,Int> coords    (vertex_count, amb_dim);
            Tensor2<Int,Int>     simplices (simplex_count,simplex_size);
            
            
            mut<ExtReal> V = coords.data();
            mut<Int>     S = simplices.data();
            
            
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
            
            
            std::unique_ptr<Base_T> M = Make(
                V, vertex_count,  amb_dim,      false,
                S, simplex_count, simplex_size, false,
                thread_count
            );

            ptoc(ClassName()+"Make_FromFile");
            
            return M;
        }
        
    public:
        
        std::string ClassName()
        {
            return std::string("SimplicialMesh_Factory")+"<"
            + ToString(MinDomDim)+ ","
            + ToString(MaxDomDim)+ ","
            + ToString(MinAmbDim)+ ","
            + ToString(MaxAmbDim)+ ","
            + TypeName<Real>     + ","
            + TypeName<Int>      + ","
            + TypeName<SReal>    + ","
            + TypeName<ExtReal>  + ","
            + ">";
        }
        
    }; // SimplicialMesh_Factory
}
