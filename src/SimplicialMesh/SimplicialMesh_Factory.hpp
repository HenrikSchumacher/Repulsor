#pragma once

namespace Repulsor
{
    
    // Some generic template for all the factories to come.
    template < class T, int ...N> class SimplicialMesh_Factory {};
    
    template<
        int MinDomDim_, int MaxDomDim_,
        int MinAmbDim_, int MaxAmbDim_,
        typename Real_, typename Int_, typename LInt_, typename SReal_, typename ExtReal_
    >
    class SimplicialMesh_Factory< SimplicialMeshBase<Real_,Int_,LInt_,SReal_,ExtReal_>, MinDomDim_, MaxDomDim_, MinAmbDim_, MaxAmbDim_ >
    {
        static_assert( MinDomDim_ <= MaxDomDim_, "MinDomDim_ <= MaxDomDim_ required." );
        static_assert( MinAmbDim_ <= MaxAmbDim_, "MinAmbDim_ <= MaxAmbDim_ required." );
        
        static_assert(IntQ< Int_>,"");
        static_assert(IntQ<LInt_>,"");
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<SReal_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using LInt    = LInt_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        
        using MeshBase_T  = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
        
        static constexpr Int MinAmbDim = Max( Int(1), Int(MinAmbDim_) );
        static constexpr Int MaxAmbDim = Max( Int(1), Int(MaxAmbDim_) );
        
        static constexpr Int MinDomDim = Max( Int(0), Int(MinDomDim_) );
        static constexpr Int MaxDomDim = Min( Int(MaxAmbDim-1), Int(MaxDomDim_) );
        
        SimplicialMesh_Factory() = default;

        template<typename ExtInt>
        [[nodiscard]] std::unique_ptr<MeshBase_T> Make(
            const ExtReal * const vertex_coords, // vertex coordinates; size = vertex_count_ x amb_dim
            const Int             vertex_count,
            const Int             amb_dim,
            const bool            vertex_coords_ColMajorQ, // whether to transpose input; set to false for row-major
            const ExtInt  * const simplices, // simplices; size = simplex_count_ x simplex_size
            const Int             simplex_count,
            const Int             simplex_size,
            const bool            simplices_ColMajorQ, // whether to transpose input; set to false for row-major
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
                vertex_coords, vertex_count, amb_dim, vertex_coords_ColMajorQ,
                simplices, simplex_count, dom_dim, simplices_ColMajorQ,
                thread_count
            );
        }
        
    private:
        
        template<Int AmbDim, typename ExtInt>
        std::unique_ptr<MeshBase_T> make_1(
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
                return make_2<Min(MaxDomDim,AmbDim-1),AmbDim>(
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
        std::unique_ptr<MeshBase_T> make_2(
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
                return std::unique_ptr<MeshBase_T>(
                    new SimplicialMesh<DomDim,AmbDim,Real,Int,LInt,SReal,ExtReal>(
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
            
        std::unique_ptr<MeshBase_T> Error( const Int dom_dim, const Int amb_dim )
        {
            eprint(ClassName()+" cannot create SimplicialMesh with domain dimension "+ToString(dom_dim)+" and  ambient dimension "+ToString(amb_dim)+".");
            
            return std::unique_ptr<MeshBase_T>( nullptr );
        }
        
    public:
        
        [[nodiscard]] std::unique_ptr<MeshBase_T> Make_FromFile( 
            const std::filesystem::path & file, const Int thread_count
        )
        {
            TOOLS_PTIC(ClassName()+"::Make_FromFile");
            
            print("Reading mesh from file " + file.string() + ".");
            
            std::ifstream s ( file );
            
            if( !s.good() )
            {
                eprint("Make_FromFile: File " + file.string() + " could not be opened.");
                
                TOOLS_PTOC(ClassName()+"::Make_FromFile");
                
                return nullptr;
            }
            
            std::string str;
            Int amb_dim;
            Int dom_dim;
            Int vertex_count;
            Int simplex_count;
            s >> str;
            s >> dom_dim;
            TOOLS_DUMP(dom_dim);
            s >> str;
            s >> amb_dim;
            TOOLS_DUMP(amb_dim);
            s >> str;
            s >> vertex_count;
            TOOLS_DUMP(vertex_count);
            s >> str;
            s >> simplex_count;
            TOOLS_DUMP(simplex_count);
            
            const Int simplex_size = dom_dim+1;
            
            TOOLS_DUMP(simplex_size);
            
            Tensor2<ExtReal,Int> coords    (vertex_count, amb_dim);
            Tensor2<Int,Int>     simplices (simplex_count,simplex_size);
            
            
            mptr<ExtReal> V = coords.data();
            mptr<Int>     S = simplices.data();
            
            
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
            
            
            std::unique_ptr<MeshBase_T> M = Make(
                V, vertex_count,  amb_dim,      false,
                S, simplex_count, simplex_size, false,
                thread_count
            );

            TOOLS_PTOC(ClassName()+"::Make_FromFile");
            
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
