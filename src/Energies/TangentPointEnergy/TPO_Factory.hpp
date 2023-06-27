namespace Repulsor
{
    template < class T, int ...N> class CONCAT(CLASS,_Factory) {};
    
    template<
        int MIN_DOM_DIM_S_, int MAX_DOM_DIM_S_,
        int MIN_DOM_DIM_T_, int MAX_DOM_DIM_T_,
        int MIN_AMB_DIM_,  int MAX_AMB_DIM_,
        typename Real, typename Int, typename LInt, typename SReal, typename ExtReal
    >
    class CONCAT(CLASS,_Factory)<BESH, MIN_DOM_DIM_S_, MAX_DOM_DIM_S_, MIN_DOM_DIM_T_, MAX_DOM_DIM_T_, MIN_AMB_DIM_, MAX_AMB_DIM_ >
    {
    public:
        static constexpr Int MIN_AMB_DIM = std::max( Int(1), Int(MIN_AMB_DIM_) );
        static constexpr Int MAX_AMB_DIM = std::max( Int(1), Int(MAX_AMB_DIM_) );
        
        static constexpr Int MIN_DOM_DIM_S = std::max( Int(0), Int(MIN_DOM_DIM_S_) );
        static constexpr Int MAX_DOM_DIM_S = std::min( Int(MAX_AMB_DIM-1), Int(MAX_DOM_DIM_S_) );

        static constexpr Int MIN_DOM_DIM_T = std::max( Int(0), Int(MIN_DOM_DIM_T_) );
        static constexpr Int MAX_DOM_DIM_T = std::min( Int(MAX_AMB_DIM-1), Int(MAX_DOM_DIM_T_) );

        
        CONCAT(CLASS,_Factory)() = default;
        
        ~CONCAT(CLASS,_Factory)() = default;

        std::unique_ptr<ROOT> Make( const Int dom_dim_S, const Int dom_dim_T, const Int amb_dim, const Real q, const Real p )
        {
            if( (dom_dim_S<MIN_DOM_DIM_S) || (dom_dim_S > MAX_DOM_DIM_S) )
            {
                eprint(ClassName()+"::Make: dom_dim_S is out of range.");
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
            
            if( (dom_dim_T<MIN_DOM_DIM_T) || (dom_dim_T > MAX_DOM_DIM_T) )
            {
                eprint(ClassName()+"::Make: dom_dim_T is out of range.");
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
            
            if( (amb_dim < MIN_AMB_DIM) || (amb_dim > MAX_AMB_DIM) )
            {
                eprint(ClassName()+"::Make: amb_dim is out of range.");
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
            
            
            if( dom_dim_S >= amb_dim || dom_dim_T >= amb_dim )
            {
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
            
            return make_1<MAX_AMB_DIM>( dom_dim_S, dom_dim_T, amb_dim, q, p );
        }
        
    private:
        
        template<Int AMB_DIM>
        std::unique_ptr<ROOT> make_1( const Int dom_dim_S, const Int dom_dim_T, const Int amb_dim, const Real q, const Real p )
        {
            if( amb_dim == AMB_DIM )
            {
                return make_2<std::min(MAX_DOM_DIM_T,AMB_DIM-1),AMB_DIM>( dom_dim_S, dom_dim_T, amb_dim, q, p );
            }
            else if constexpr ( AMB_DIM > MIN_AMB_DIM )
            {
                return make_1<AMB_DIM-1>( dom_dim_S, dom_dim_T, amb_dim, q, p );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
        }
        
        template<Int DOM_DIM_T, Int AMB_DIM>
        std::unique_ptr<ROOT> make_2( const Int dom_dim_S, const Int dom_dim_T, const Int amb_dim, const Real q, const Real p )
        {
            if( dom_dim_T == DOM_DIM_T )
            {
                return make_3<std::min(MAX_DOM_DIM_S,AMB_DIM-1),DOM_DIM_T,AMB_DIM>( dom_dim_S, dom_dim_T, amb_dim, q, p );
            }
            else if constexpr ( DOM_DIM_T > MIN_DOM_DIM_T )
            {
                return make_2<DOM_DIM_T-1,AMB_DIM>( dom_dim_S, dom_dim_T, amb_dim, q, p );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
        }
        
        template<Int DOM_DIM_S, Int DOM_DIM_T, Int AMB_DIM>
        std::unique_ptr<ROOT> make_3( const Int dom_dim_S, const Int dom_dim_T, const Int amb_dim, const Real q, const Real p )
        {
            if( dom_dim_S == DOM_DIM_S )
            {
                return std::unique_ptr<BASE>( new CLASS<MESH,DOM_DIM_T>( q, p ) );
            }
            else if constexpr ( DOM_DIM_S > MIN_DOM_DIM_S )
            {
                return make_3<DOM_DIM_S-1,DOM_DIM_T,AMB_DIM>( dom_dim_S, dom_dim_T, amb_dim, q, p );
            }
            else
            {
                // Actually, this should not be reachable.
                return Error( dom_dim_S, dom_dim_T, amb_dim );
            }
        }
            
        std::unique_ptr<ROOT> Error( const Int dom_dim_S, const Int dom_dim_T, const Int amb_dim )
        {
            eprint(ClassName()+" cannot create "+TO_STD_STRING(CLASS)+" with domain dimension "+ToString(dom_dim_S)+",obstacle's domain dimension "+ToString(dom_dim_T)+", and ambient dimension "+ToString(amb_dim)+".");
            
            return std::unique_ptr<ROOT>( nullptr );
        }
        
    public:
        
        std::string ClassName()
        {
            return TO_STD_STRING(CLASS)+"_Factory<"
            + ToString(MIN_DOM_DIM_S)    + ","
            + ToString(MAX_DOM_DIM_S)    + ","
            + ToString(MIN_DOM_DIM_T)    + ","
            + ToString(MAX_DOM_DIM_T)    + ","
            + ToString(MIN_AMB_DIM)      + ","
            + ToString(MAX_AMB_DIM)      + ","
            + TypeName<Real>    + ","
            + TypeName<Int>     + ","
            + TypeName<LInt>    + ","
            + TypeName<SReal>   + ","
            + TypeName<ExtReal> + ","
            + ">";
        }
        
    }; // CONCAT(CLASS,_Factory)
}


