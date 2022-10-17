#pragma once

namespace Repulsor
{
    template<typename S, typename T, typename Int>
    class alignas(OBJECT_ALIGNMENT) PairAggregator
    {
        Int current_size = 0;
        Int capacity     = 0;
        Tensor1<S, Int> buffer_0;
        Tensor1<T, Int> buffer_1;
        
    public:
        
        PairAggregator( const Int n )
        :   buffer_0 (n)
        ,   buffer_1 (n)
        ,   capacity (n)
        ,   current_size(0)
        {}
        
        void Push( const S a, const T b )
        {
            if( current_size < capacity )
            {
                buffer_0[current_size] = a;
                buffer_1[current_size] = b;
                ++current_size;
            }
            else
            {
                const Int new_capacity = 2 * capacity;
                Tensor1<S, Int> new_buffer_0 (new_capacity);
                Tensor1<T, Int> new_buffer_1 (new_capacity);
                
                copy_buffer( buffer_0.data(), new_buffer_0.data(), capacity );
                copy_buffer( buffer_1.data(), new_buffer_1.data(), capacity );
                
                swap( buffer_0, new_buffer_0 );
                swap( buffer_1, new_buffer_1 );
                
                buffer_0[0] = a;
                buffer_1[0] = b;
                
                current_size = 1;
                capacity     = new_capacity;
            }
        }
        
        template<>
        Tensor1<S,Int> & Get<0>()
        {
            return buffer_0;
        }
        
        template<>
        const Tensor1<S,Int> & Get<0>() const
        {
            return buffer_0;
        }
        
        template<>
        Tensor1<T,Int> & Get<1>()
        {
            return buffer_1;
        }
        
        template<>
        const Tensor1<T,Int> & Get<1>() const
        {
            return buffer_1;
        }
    };
