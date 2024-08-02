#include "bench_template.hpp"
#include "../hourglass.cpp"

template <typename t_itr, typename... Args>
inline Hourglass init_hourglass(const t_itr begin, const t_itr end, const double bpk, Args... args)
{
    auto&& t = std::forward_as_tuple(args...);
    auto max_range_size = std::get<0>(t);
    auto corr_degree = std::get<1>(t);
    start_timer(build_time);
    Hourglass f(begin, end, max_range_size, corr_degree * 0.1, bpk);
    stop_timer(build_time);
    return f;
}

template <typename value_type>
inline bool query_hourglass(Hourglass &f, const value_type left, const value_type right)
{
    return f.Lookup_Adaptive(left, right);
}


inline size_t size_hourglass(Hourglass &f)
{
    return f.size() >> 3;
}

int main(int argc, char const *argv[])
{
    argparse::ArgumentParser parser("bench-hourglass");
    init_parser(parser);
    try
    {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    auto [ keys, queries, arg , range_size, corr_degree ] = read_parser_arguments(parser);
    experiment(pass_fun(init_hourglass),pass_ref(query_hourglass),
               pass_ref(size_hourglass), arg, keys, queries, range_size, corr_degree);
    print_test();

    return 0;
}