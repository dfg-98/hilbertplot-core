#ifndef PARALLEL_ALGORITHM_H
#define PARALLEL_ALGORITHM_H

#include "threads_utility.h"
#include <algorithm>
#include <cassert>

#include <iostream>

class parallel_algorithm
{
    public:
        parallel_algorithm() {}
};


template <typename Iter>
void swap_range_reverse(Iter first_beg, Iter first_end, Iter last_beg, Iter last_end)
{
   uint lenght = first_end - first_beg;
   assert (lenght == (last_end - last_beg));
   for(uint i = 0; i < lenght; ++i)
   {
       std::iter_swap(first_beg++, --last_end);
   }
}

template <typename Iter>
void reverse_parallel(Iter first, Iter last)
{
    auto lenght = std::distance(first, last);
    if(!lenght)
        return;

    thread_pool pool;
    const unsigned long mid_lenght = lenght / 2;
    //const unsigned long block_size = mid_lenght / (pool.get_thread_count () * 5) + 1;
    const unsigned long block_size = 10000;
    if(mid_lenght < block_size)
    {
        return std::reverse(first, last);
    }

    const unsigned long num_blocks = (mid_lenght + block_size - 1) / block_size;


    Iter start = first;
    Iter end = last;
    for(unsigned long i = 0; i < num_blocks-1; ++i)
    {
        std::function<void()> task = std::bind(swap_range_reverse<Iter>, start, start+block_size, end-block_size, end);
        pool.push_task (task);
        start += block_size;
        end -= block_size;
    }
    swap_range_reverse(start, first + mid_lenght, last-mid_lenght, end);
    while (pool.isWorking ())
    {
        pool.run_task ();
    }
}

template <typename Iter, typename Func>
void for_each_parallel(Iter first, Iter last, Func f)
{
    unsigned long const length=std::distance(first,last);
    if(!length)
        return;
    unsigned long const min_per_thread=10000;

    if(length<(2*min_per_thread))
    {
        std::for_each(first,last,f);
    }
    else
    {
        Iter const mid_point=first+length/2;
        std::future<void> first_half= std::async(&for_each_parallel<Iter,Func>, first,mid_point,f);
        for_each_parallel(mid_point,last,f);
        first_half.get();
    }
}
#endif // PARALLEL_ALGORITHM_H
