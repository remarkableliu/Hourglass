#include "hourglass.cpp"

int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> distr_value(0, UINT64_MAX - 1);

    std::set<uint64_t> key_set;
    std::vector<uint64_t> keys;
    std::vector<std::pair<uint64_t, uint64_t>> queries;
    size_t n_keys = 15'938'355; // 2^24 * 0.95
    size_t n_queries = 1'000'000; 
    size_t range_size = 32;

    // generate keys
    while (key_set.size() < n_keys) {
        auto key = distr_value(gen);
        key_set.insert(key);
    }
    keys = {key_set.begin(), key_set.end()};

    // generate queries
    std::set<uint64_t> middle_points;
    while (middle_points.size() < n_queries) {
        auto middle = distr_value(gen);
        middle_points.insert(middle);
    }

    for (auto m : middle_points)
    {
        queries.push_back(std::pair<uint64_t, uint64_t>(m, m + range_size - 1));
    }
    
    // build Hourglass
    Hourglass h(keys.begin(), keys.end(), 5, 0, 12);

    // query on Hourglass
    size_t fp = 0;
    for (size_t i = 0; i < n_queries; i++)
    {
        if (h.Lookup_Adaptive(queries[i].first, queries[i].second))
            fp++;
    }
    std::cout << "FPR: " << fp / (double)n_queries << std::endl;
    std::cout << "BPK: " << h.size() / (double)n_keys << std::endl;

    return 0;
}