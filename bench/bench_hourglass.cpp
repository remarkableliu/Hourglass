#include "./bench_template.hpp"
#include "../hourglass.cpp"


inline std::vector<std::pair<uint64_t, uint64_t>> SampleQueries(std::vector<std::pair<uint64_t, uint64_t>> queries, double sample_rate) {
    std::vector<std::pair<uint64_t, uint64_t>> sample_queries;

    std::mt19937 generator(std::random_device{}());
    std::shuffle(queries.begin(), queries.end(), generator);
    auto sample_size = queries.size() * sample_rate;
    
    for (size_t i = 0; i < sample_size; i++)
    {
        sample_queries.push_back(queries[i]);
    }
    return sample_queries;
}

template <typename t_itr>
inline int EsitimateCorrDegree(const t_itr begin, const t_itr end, std::vector<std::pair<uint64_t, uint64_t>> queries, size_t range_size) {
    auto interval = 1 << range_size;
    uint64_t corr = 0;
    for (auto pair : queries)
    {
        auto l = pair.first;
        auto r = pair.second;

        uint64_t distance = 0;
        auto low = std::lower_bound(begin, end, l);
        if (low == end) {
            distance = l - *(end - 1);
        } else {
            distance = min(l - *(low - 1), *low - r);
        }
        if (distance <= interval)
            corr++;
    }
    if (corr == 0)
        return 0;
    auto D = 1 - std::log2(queries.size() / (double) corr * interval) / 30;
    return std::round(D * 10);
} 

template <typename t_itr, typename... Args>
inline Hourglass init_hourglass(const t_itr begin, const t_itr end, const double bpk, Args... args)
{
    auto&& t = std::forward_as_tuple(args...);
    auto queries_temp = std::get<0>(t);
    auto max_range_size = std::get<1>(t);

    auto queries = std::vector<std::pair<uint64_t, uint64_t>>(queries_temp.size());
    std::transform(queries_temp.begin(), queries_temp.end(), queries.begin(), [](auto x) {
        auto [left, right, result] = x;
        return std::make_pair(left, right);
    });

    auto sample_queries = SampleQueries( queries, 0.2);
    start_timer(modelling_time);
    auto D = EsitimateCorrDegree(begin, end, sample_queries, max_range_size);
    stop_timer(modelling_time);
    start_timer(build_time);
    Hourglass f(begin, end, max_range_size, D * 0.1, bpk);
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
               pass_ref(size_hourglass), arg, keys, queries, queries, range_size);
    print_test();

    return 0;
}