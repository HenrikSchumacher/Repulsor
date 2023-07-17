#pragma once

namespace Repulsor
{
    
    // Some generic template for all the factories to come.
    template < class T, int ...N> class SimplicialRemesher_Factory {};
    
    template<
        int MinDomDim_, int MaxDomDim_,
        int MinAmbDim_, int MaxAmbDim_,
        typename Real_, typename Int_
    >
    class SimplicialRemesher_Factory< SimplicialRemesherBase<Real_,Int_>, MinDomDim_, MaxDomDim_, MinAmbDim_, MaxAmbDim_ >
    {
        static_assert( MinDomDim_ <= MaxDomDim_, "MinDomDim_ <= MaxDomDim_ required." );
        static_assert( MinAmbDim_ <= MaxAmbDim_, "MinAmbDim_ <= MaxAmbDim_ required." );
        
        ASSERT_FLOAT(Real_);
        ASSERT_INT(Int_);
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        
        using RemesherBase_T  = SimplicialRemesherBase<Real,Int>;
        
        static constexpr Int MinAmbDim = std::max( Int(1), Int(MinAmbDim_) );
        static constexpr Int MaxAmbDim = std::max( Int(1), Int(MaxAmbDim_) );
        
        static constexpr Int MinDomDim = std::max( Int(0), Int(MinDomDim_) );
        static constexpr Int MaxDomDim = std::min( Int(MaxAmbDim-1), Int(MaxDomDim_) );
        
        SimplicialRemesher_Factory() = default;
        
        ~SimplicialRemesher_Factory() = default;

        template<typename ExtReal, typename ExtInt>
        std::unique_ptr<RemesherBase_T> Make(
            ptr<ExtReal> vertex_coords, // vertex coordinates; size = vertex_count_ x amb_dim
            const Int    vertex_count,
            const Int    amb_dim,
            const bool   vertex_coords_ColMajorQ, // whether to transpose input; set to false for row-major
            ptr<ExtInt>  simplices, // simplices; size = simplex_count_ x simplex_size
            const Int    simplex_count,
            const Int    simplex_size,
            const bool   simplices_ColMajorQ, // whether to transpose input; set to false for row-major
            const Int    thread_count = 1
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
        
        template<Int AmbDim, typename ExtReal, typename ExtInt>
        std::unique_ptr<RemesherBase_T> make_1(
            ptr<ExtReal> v,
            const Int    v_cnt,
            const Int    amb_dim,
            const bool   vertex_coords_ColMajorQ,
            ptr<ExtInt>  s,
            const Int    s_cnt,
            const Int    dom_dim,
            const bool   simplices_ColMajorQ,
            const Int    thread_count = 1
        )
        {
            if( amb_dim == AmbDim )
            {
                return make_2<std::min(MaxDomDim,AmbDim-1),AmbDim>(
                    v, v_cnt, amb_dim, vertex_coords_ColMajorQ,
                    s, s_cnt, dom_dim, simplices_ColMajorQ,
                    thread_count
                );
            }
            else if constexpr ( AmbDim > MinAmbDim )
            {
                return make_1<AmbDim-1>(
                    v, v_cnt, amb_dim, vertex_coords_ColMajorQ,
                    s, s_cnt, dom_dim, simplices_ColMajorQ,
                    thread_count
                );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
        
        template<Int DomDim, Int AmbDim, typename ExtInt>
        std::unique_ptr<RemesherBase_T> make_2(
            ptr<ExtReal> v,
            const Int    v_cnt,
            const Int    amb_dim,
            const bool   vertex_coords_ColMajorQ,
            ptr<ExtInt>  s,
            const Int    s_cnt,
            const Int    dom_dim,
            const bool   simplices_ColMajorQ,
            const Int    thread_count = 1
        )
        {
            if( dom_dim == DomDim )
            {
                return std::unique_ptr<RemesherBase_T>(
                    new SimplicialRemesher<DomDim,AmbDim,Real,Int>(
                        v, v_cnt, vertex_coords_ColMajorQ,
                        s, s_cnt, simplices_ColMajorQ,
                        thread_count
                    )
                );
            }
            else if constexpr ( DomDim > MinDomDim )
            {
                return make_2<DomDim-1,AmbDim>(
                    v, v_cnt, amb_dim, vertex_coords_ColMajorQ,
                    s, s_cnt, dom_dim, simplices_ColMajorQ,
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
        
        std::string ClassName()
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

