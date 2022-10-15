#pragma once

namespace Repulsor
{
    
    template<typename Int>
    struct Cluster // slim POD container to hold only the data relevant for the construction phase in the tree, before it is serialized
    {
    public:
        
        Cluster() = default;

        Cluster( Int thread_, Int ID_, Int begin_, Int end_, Int depth_ )
        :   thread(thread_)
        ,   ID(ID_)
        ,   begin(begin_)
        ,   end(end_)
        ,   depth(depth_)
        ,   max_depth(depth_)
        ,   descendant_count(0)
        ,   descendant_leaf_count(0)
        ,   left(nullptr)
        ,   right(nullptr)
        {}
        
        ~Cluster() = default;
        
        Int thread = -1;               // thread that created this cluster
        Int ID = -1;                   // ID of cluster within this thread
        Int begin = 0;                 // position of first primitive in cluster relative to array ordering
        Int end = 0;                   // position behind last primitive in cluster relative to array ordering
        Int depth = 0;                 // depth within the tree -- not absolutely necessary but nice to have for plotting images
        Int max_depth = 0;             // maximal depth of all descendants of this cluster
        Int descendant_count = 0;      // number of descendents of cluster, _this cluster included_
        Int descendant_leaf_count = 0; // number of leaf descendents of cluster
        Cluster<Int> *left  = nullptr; // left child
        Cluster<Int> *right = nullptr; // right child
        
    }; //Cluster
        
} // namespace Repulsor
