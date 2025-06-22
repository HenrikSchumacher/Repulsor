#pragma once

namespace Repulsor
{
    
    // Some generic template for all the factories to come.
    template < class T, int ...N> class SimplicialRemesher_Factory {};
    
    template<
        int MinDomDim_, int MaxDomDim_,
        int MinAmbDim_, int MaxAmbDim_,
        typename Real_, typename Int_, typename ExtReal_, typename ExtInt_
    >
    class SimplicialRemesher_Factory< SimplicialRemesherBase<Real_,Int_,ExtReal_,ExtInt_>, MinDomDim_, MaxDomDim_, MinAmbDim_, MaxAmbDim_ >
    {
        static_assert( MinDomDim_ <= MaxDomDim_, "MinDomDim_ <= MaxDomDim_ required." );
        static_assert( MinAmbDim_ <= MaxAmbDim_, "MinAmbDim_ <= MaxAmbDim_ required." );
        
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
        static_assert(IntQ<Int_>,"");
        static_assert(IntQ<ExtInt_>,"");
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using ExtReal = ExtReal_;
        using ExtInt  = ExtInt_;
        
        using RemesherBase_T  = SimplicialRemesherBase<Real,Int,ExtReal,ExtInt>;
        
        static constexpr Int MinAmbDim = Max( Int(1), Int(MinAmbDim_) );
        static constexpr Int MaxAmbDim = Max( Int(1), Int(MaxAmbDim_) );
        
        static constexpr Int MinDomDim = Max( Int(0), Int(MinDomDim_) );
        static constexpr Int MaxDomDim = Min( Int(MaxAmbDim-1), Int(MaxDomDim_) );
        
        SimplicialRemesher_Factory() = default;

        [[nodiscard]] std::unique_ptr<RemesherBase_T> Make(
            cptr<ExtReal> vertex_coords, // vertex coordinates; size = vertex_count_ x amb_dim
            const Int     vertex_count,
            const Int     amb_dim,
            const bool    vertex_coords_ColMajorQ, // whether to transpose input; set to false for row-major
            cptr<ExtInt>  simplices, // simplices; size = simplex_count_ x simplex_size
            const Int     simplex_count,
            const Int     simplex_size,
            const bool    simplices_ColMajorQ, // whether to transpose input; set to false for row-major
            cptr<ExtReal> vertex_data,
            const Int     data_dim,
            const bool    vertex_data_ColMajorQ,
            const Int     thread_count = 1
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
                vertex_coords, vertex_count,  amb_dim,      vertex_coords_ColMajorQ,
                simplices,     simplex_count, dom_dim,      simplices_ColMajorQ,
                vertex_data,                  data_dim,     vertex_data_ColMajorQ,
                thread_count
            );
        }
        
    private:
        
        template<Int AmbDim>
        std::unique_ptr<RemesherBase_T> make_1(
            cptr<ExtReal> vertex_coords,
            const Int     vertex_count,
            const Int     amb_dim,
            const bool    vertex_coords_ColMajorQ,
            cptr<ExtInt>  simplices,
            const Int     simplex_count,
            const Int     dom_dim,
            const bool    simplices_ColMajorQ,
            cptr<ExtReal> vertex_data,
            const Int     data_dim,
            const bool    vertex_data_ColMajorQ,
            const Int     thread_count = 1
        )
        {
            if( amb_dim == AmbDim )
            {
                return make_2<Min(MaxDomDim,AmbDim-1),AmbDim>(
                    vertex_coords, vertex_count,  amb_dim,  vertex_coords_ColMajorQ,
                    simplices,     simplex_count, dom_dim,  simplices_ColMajorQ,
                    vertex_data,                  data_dim, vertex_data_ColMajorQ,
                    thread_count
                );
            }
            else if constexpr ( AmbDim > MinAmbDim )
            {
                return make_1<AmbDim-1>(
                    vertex_coords, vertex_count,  amb_dim,  vertex_coords_ColMajorQ,
                    simplices,     simplex_count, dom_dim,  simplices_ColMajorQ,
                    vertex_data,                  data_dim, vertex_data_ColMajorQ,
                    thread_count
                );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
        
        template<Int DomDim, Int AmbDim>
        std::unique_ptr<RemesherBase_T> make_2(
            cptr<ExtReal> vertex_coords,
            const Int     vertex_count,
            const Int     amb_dim,
            const bool    vertex_coords_ColMajorQ,
            cptr<ExtInt>  simplices,
            const Int     simplex_count,
            const Int     dom_dim,
            const bool    simplices_ColMajorQ,
            cptr<ExtReal> vertex_data,
            const Int     data_dim,
            const bool    vertex_data_ColMajorQ,
            const Int     thread_count = 1
        )
        {
            if( dom_dim == DomDim )
            {
                return std::unique_ptr<RemesherBase_T>(
                    new SimplicialRemesher<DomDim,AmbDim,Real,Int,ExtReal,ExtInt>(
                        vertex_coords, vertex_count,            vertex_coords_ColMajorQ,
                        simplices,     simplex_count,           simplices_ColMajorQ,
                        vertex_data,                  data_dim, vertex_data_ColMajorQ,
                        thread_count
                    )
                );
            }
            else if constexpr ( DomDim > MinDomDim )
            {
                return make_2<DomDim-1,AmbDim>(
                    vertex_coords, vertex_count,  amb_dim,  vertex_coords_ColMajorQ,
                    simplices,     simplex_count, dom_dim,  simplices_ColMajorQ,
                    vertex_data,                  data_dim, vertex_data_ColMajorQ,
                    thread_count
                );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
            
        std::unique_ptr<RemesherBase_T> Error( const Int dom_dim, const Int amb_dim )
        {
            eprint(ClassName()+" cannot create SimplicialRemesher with domain dimension "+ToString(dom_dim)+" and  ambient dimension "+ToString(amb_dim)+".");
            
            return std::unique_ptr<RemesherBase_T>( nullptr );
        }
        
        
    public:
        
        static std::string ClassName()
        {
            return std::string("SimplicialRemesher_Factory")+"<"
            + ToString(MinDomDim)+ ","
            + ToString(MaxDomDim)+ ","
            + ToString(MinAmbDim)+ ","
            + ToString(MaxAmbDim)+ ","
            + TypeName<Real>     + ","
            + TypeName<Int>      + ","
            + ">";
        }
        
    }; // SimplicialRemesher_Factory
}

