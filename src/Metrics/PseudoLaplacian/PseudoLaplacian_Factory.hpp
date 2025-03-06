namespace Repulsor
{
    template < class T, int ...N> class TOOLS_CONCAT2(CLASS,_Factory) {};
    
    template<
        int MIN_DOM_DIM_, int MAX_DOM_DIM_,
        int MIN_AMB_DIM_, int MAX_AMB_DIM_,
        typename Real, typename Int, typename LInt, typename SReal, typename ExtReal
    >
    class TOOLS_CONCAT2(CLASS,_Factory)<BESH, MIN_DOM_DIM_, MAX_DOM_DIM_, MIN_AMB_DIM_, MAX_AMB_DIM_ >
    {
    public:

        static constexpr Int MIN_AMB_DIM = Max( Int(1), Int(MIN_AMB_DIM_) );
        static constexpr Int MAX_AMB_DIM = Max( Int(1), Int(MAX_AMB_DIM_) );
        
        static constexpr Int MIN_DOM_DIM = Max( Int(0), Int(MIN_DOM_DIM_) );
        static constexpr Int MAX_DOM_DIM = Min( Int(MAX_AMB_DIM-1), Int(MAX_DOM_DIM_) );
        
        TOOLS_CONCAT2(CLASS,_Factory)() = default;
        
        ~TOOLS_CONCAT2(CLASS,_Factory)() = default;

        [[nodiscard]] std::unique_ptr<ROOT> Make( const Int dom_dim, const Int amb_dim, const Real s )
        {
            if( (dom_dim < MIN_DOM_DIM) || (dom_dim > MAX_DOM_DIM) )
            {
                eprint(ClassName()+"::Make: dom_dim is out of range.");
                return Error( dom_dim, amb_dim );
            }
            
            if( (amb_dim < MIN_AMB_DIM) || (amb_dim > MAX_AMB_DIM) )
            {
                eprint(ClassName()+"::Make: amb_dim is out of range.");
                return Error( dom_dim, amb_dim );
            }
            
            if( dom_dim >= amb_dim )
            {
                return Error( dom_dim, amb_dim );
            }
            
            return make_1<MAX_AMB_DIM>( dom_dim, amb_dim, s );
        }
        
    private:
        
        template<Int AMB_DIM>
        std::unique_ptr<ROOT> make_1( const Int dom_dim, const Int amb_dim, const Real s )
        {
            if( amb_dim == AMB_DIM )
            {
                return make_2<Min(MAX_DOM_DIM,AMB_DIM-1),AMB_DIM>( dom_dim, amb_dim, s );
            }
            else if constexpr ( AMB_DIM > MIN_AMB_DIM )
            {
                return make_1<AMB_DIM-1>( dom_dim, amb_dim, s );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
        
        template<Int DOM_DIM, Int AMB_DIM>
        std::unique_ptr<ROOT> make_2( const Int dom_dim, const Int amb_dim, const Real s )
        {
            if( dom_dim == DOM_DIM )
            {
                return make_3<DOM_DIM,AMB_DIM>( s );
            }
            else if constexpr ( DOM_DIM > MIN_DOM_DIM )
            {
                return make_2<DOM_DIM-1,AMB_DIM>( dom_dim, amb_dim, s );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim, amb_dim );
            }
        }
        
        template<Int DOM_DIM, Int AMB_DIM>
        std::unique_ptr<ROOT> make_3( const Real s )
        {
            if( (1 < s) && (s < 2) )
            {
                return std::unique_ptr<ROOT>(
                    new CLASS<
                        SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>,
                        true
                    >( s ) );
            }
            else if( (0 < s) && (s < 1) )
            {
                return std::unique_ptr<ROOT>(
                    new CLASS<
                        SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>,
                        false
                    >( s ) );
            }
            else
            {
                eprint(ClassName()+": parameter s = "+ToString(s)+" is out of bounds.");
                return std::unique_ptr<ROOT>( nullptr );
            }
        }
            
        std::unique_ptr<ROOT> Error( const Int dom_dim, const Int amb_dim )
        {
            eprint(ClassName()+" cannot create "+TOOLS_TO_STD_STRING(CLASS)+" with domain dimension "+ToString(dom_dim)+" and  ambient dimension "+ToString(amb_dim)+".");
            
            return std::unique_ptr<ROOT>( nullptr );
        }
        
    public:
        
        std::string ClassName()
        {
            return TOOLS_TO_STD_STRING(CLASS)+"_Factory<"
            + ToString(MIN_DOM_DIM)+ ","
            + ToString(MAX_DOM_DIM)+ ","
            + ToString(MIN_AMB_DIM)+ ","
            + ToString(MAX_AMB_DIM)+ ","
            + TypeName<Real>       + ","
            + TypeName<Int>        + ","
            + TypeName<LInt>       + ","
            + TypeName<SReal>      + ","
            + TypeName<ExtReal>    + ","
            + ">";
        }
        
    }; // TOOLS_CONCAT2(CLASS,_Factory)
}


