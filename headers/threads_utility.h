#ifndef SCOPED_THREAD_H
#define SCOPED_THREAD_H

#include <thread>
#include <atomic>
#include <mutex>
#include <list>
#include <functional>
#include <vector>
#include <future>

class scoped_thread
{
    public:
        explicit scoped_thread(std::thread t_):
            t(std::move(t_))
        {
            if(!t.joinable ())throw std::logic_error("No thread");
        }
        ~scoped_thread()
        {
            t.join ();
        }
        scoped_thread(scoped_thread&)= delete;
        scoped_thread& operator=(scoped_thread const &)= delete;

    private:
        std::thread t;
};

class function_wrapper
{
        struct impl_base
        {
                virtual void call()=0;
                virtual ~impl_base() {}
        };

        std::unique_ptr<impl_base> impl;
        template<typename F>
        struct impl_type: impl_base
        {
                F f;
                impl_type(F&& f_): f(std::move(f_)) {}
                void call() { f(); }
        };
    public:
        template<typename F>
        function_wrapper(F&& f):
            impl(new impl_type<F>(std::move(f)))
        {}
        void operator()() { impl->call(); }
        function_wrapper() = default;
        function_wrapper(function_wrapper&& other):
            impl(std::move(other.impl))
        {}
        function_wrapper& operator=(function_wrapper&& other)
        {
            impl=std::move(other.impl);
            return *this;
        }
        function_wrapper(const function_wrapper&)=delete;
        function_wrapper(function_wrapper&)=delete;
        function_wrapper& operator=(const function_wrapper&)=delete;
};


class thread_pool
{
    public:
        thread_pool(int nr_threads = 0):
           done(false), unfinish_tasks(0)
        {
            if(nr_threads <= 0)
            {
                int hardware_concurrency = std::thread::hardware_concurrency ();
                thread_count = hardware_concurrency > 1?  hardware_concurrency - 1: 1;
            }
            else
                thread_count = nr_threads;
            for (unsigned int i = 0; i < thread_count; ++i)
            {
                threads.push_back (std::thread(&thread_pool::worker_thread, this));
            }
        }
        virtual ~thread_pool()
        {
            done = true;
            for(auto & th: threads)
            {
                if(th.joinable ())
                {
                    th.join ();
                }
            }
        }
        //template<typename F, typename ...Args>
        //std::future<typename std::result_of<F()>::type> push_task(F func)
        void push_task(std::function<void()> func)
        {
//            typedef typename  std::result_of<F()>::type result_t;

//            std::packaged_task<result_t()> packed_task(std::move(func));
//            std::future<result_t> future = packed_task.get_future();
            std::unique_lock<std::mutex> lck(m_mutex);
            work_queue.push_back (func);
            ++unfinish_tasks;
            //return future;
        }
        int get_work_queue_size()
        {
            std::unique_lock<std::mutex> lck(m_mutex);
            return work_queue.size ();
        }
        uint get_thread_count()
        {
            return thread_count;
        }
        void run_task()
        {
            m_mutex.lock ();
            if(work_queue.empty ())
            {
                m_mutex.unlock ();
                std::this_thread::yield ();
            }
            else
            {
                //function_wrapper task (std::move(work_queue.front ()));
                auto task = work_queue.front ();
                work_queue.pop_front ();
                m_mutex.unlock ();
                task ();

                --unfinish_tasks;

            }
        }
        bool isWorking()
        {
            return unfinish_tasks != 0;
        }
        // Global instance
        static thread_pool& instance()
        {
            static thread_pool instance;
            return instance;
        }
    private:
        std::atomic_bool done;
        std::atomic_int unfinish_tasks;
        unsigned int thread_count;
        std::mutex m_mutex;
        //std::list<function_wrapper> work_queue;
        std::list<std::function<void()>> work_queue;
        std::vector<std::thread> threads;
        void worker_thread()
        {
            while (!done)
            {
                run_task();
            }
        }
};

#endif // SCOPED_THREAD_H




