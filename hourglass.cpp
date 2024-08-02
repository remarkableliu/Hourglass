#include <iostream>
#include <cstdint>
#include <vector>
#include <set>
#include "hashutil.h"
#include "permencoding.h"
#include "rank_util.cpp"
#include <cassert>
#include <fstream>

#define MASK(i) ((1ULL << i) - 1)
#define MASK128(i) (((__uint128_t)1 << i) - 1)
#define POW(i) (1ULL << i)

uint64_t getNextPow2(uint64_t n)
{
    n--;
	n |= n >> 1;
	n |= n >> 2;
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16;
	n |= n >> 32;
	n++;
	return n;
}

uint64_t getTrailingZero(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n = n >> 1;
        i++;
    }
    return i - 1;
}

using namespace std;

const size_t kMaxCuckooCount = 500;

class Hourglass
{   
    public:
    size_t num_items_ = 0;
    static const size_t bucket_size_ = 4;
    size_t num_layers_;
    size_t threshold_;
    size_t num_buckets_;
    size_t fp_len_;
    size_t kBitsPerBucket_;
    size_t keys_tag_len_;   // log(max_suffix_size) + 1
    size_t kBytesKeys_; // get from preprocessing
    size_t kBytesTags_;

    double corr_degree_;
    size_t pruned_layers_;
    bool correlated_ = false;

    char *c_packed_tags_;
    rank_multi rm;
    char *c_keys_;

    // adaptive lookup table
    uint64_t** value_table_;
    uint64_t current_value_;
#ifdef ADAPT
    static const size_t selector_len_ = 1;
    uint64_t* hash_selector_multi_[selector_len_];
    rank_multi rm_selector_;
#endif
    size_t current_selector_;
    
    typedef struct {
        size_t index;
        uint32_t tag;
        set<uint16_t> keys;
        size_t selector;
        bool used;
    } VictimCache;

    VictimCache victim_;

    PermEncoding perm_;

    ~Hourglass() {
        
    }

    template <typename t_itr>
    Hourglass(const t_itr begin, const t_itr end, const size_t max_range_size, const float corr_degree, const size_t bpk) {
        if (begin == end)
            return;
        if (!std::is_sorted(begin, end))
            throw std::runtime_error("error, the input is not sorted");
        if (max_range_size > 16)
            throw std::runtime_error("error, range size is too big to store(uint16_t).");
        if (corr_degree > 1 || corr_degree < 0)
            throw std::runtime_error("error, illegal correlated degree.");
#ifdef ADAPT
        std::cout << "[+] adaptivity used." << std::endl;
#endif

        auto num_keys = std::distance(begin, end);

        corr_degree_ = corr_degree;
        auto query_scale = floor(30 * (1 - corr_degree_ + 0.01));
        auto suffix_layer = max_range_size;

        if (corr_degree_ == 0) {
            correlated_ = false;
            num_layers_ = max_range_size;
            build_uncorr(begin, end, num_keys, bpk);
            return;
        }

        uint64_t prev_prefix = *begin >> suffix_layer;
        size_t num_unique_prefix = 1;
        for (auto it = begin + 1; it != end; it++) {
            uint64_t cur_prefix = *it >> suffix_layer; 
            if (prev_prefix != cur_prefix) {
                prev_prefix = cur_prefix;
                num_unique_prefix++;
            } 
        }
        num_layers_ = corr_aware_model(suffix_layer, query_scale, bpk, num_keys, num_unique_prefix);
        pruned_layers_ = suffix_layer - num_layers_;
        correlated_ = num_layers_ != 0;
        if (correlated_) {
            threshold_ = POW(num_layers_) / num_layers_ + 1;
            build_corr(begin, end, num_keys, bpk);
        }
        else {
            num_layers_ = max_range_size;
            build_uncorr(begin, end, num_keys, bpk);
        }
        
    }

    template <typename t_itr>
    void build_uncorr(const t_itr begin, const t_itr end, const size_t num_keys, const size_t bpk) {
        num_items_ = num_keys;
        uint64_t *prefixes = new uint64_t[num_keys]();
        uint64_t prev_prefix = *begin >> num_layers_;
        size_t num_unique_prefix = 0;
        prefixes[num_unique_prefix] = prev_prefix;
        num_unique_prefix++;
        for (auto it = begin + 1; it != end; it++) {
            uint64_t cur_prefix = *it >> num_layers_; 
            if (prev_prefix != cur_prefix) {
                prev_prefix = cur_prefix;
                prefixes[num_unique_prefix] = prev_prefix;
                num_unique_prefix++;
            } 
        }
        
        size_t max_num_items = num_unique_prefix++;
        size_t assoc = 4;
        size_t num_buckets = getNextPow2(std::max<uint64_t>(1, max_num_items / assoc));
        double frac = (double)max_num_items / num_buckets / assoc;
        if (frac > 0.95) {
            num_buckets <<= 1;
        }
        num_buckets_ = num_buckets;
        
        fp_len_ = bpk + 1;
        if (fp_len_ > 32) {
            // throw std::runtime_error("error, bucket length exceeds 128.");
            fp_len_ = 32;
        }

        if (fp_len_ <= 5)
            fp_len_ = 6;
            // throw std::runtime_error("error, fp is too short.");

        std::cout << "[+] prefix length: " << fp_len_ << ", suffix length: 0." << std::endl;
        
        kBitsPerBucket_ = (fp_len_ - 1) * 4;
        kBytesTags_ = (num_buckets_ * bucket_size_ * (fp_len_ - 1) + 127) >> 3;
        c_packed_tags_ = new char[kBytesTags_]();

        
        value_table_ = new uint64_t*[num_buckets_]();
        for (size_t i = 0; i < num_buckets_; i++) {
            value_table_[i] = new uint64_t[bucket_size_]();
        }

#ifdef ADAPT
        size_t num_words = (num_buckets_ * bucket_size_ + 63) / 64;
        for (size_t i = 0; i < selector_len_; i++)
        {
            hash_selector_multi_[i] = new uint64_t[num_words]();
        }
#endif

        victim_.used = false;

        for (size_t i = 0; i < max_num_items; i++)
        {
            if (!Insert(prefixes[i]))
                throw std::runtime_error("error, failed to insert(fp is too short).");
        }

        Compress();

        delete[] prefixes;
    }

    template <typename t_itr>
    void build_corr(const t_itr begin, const t_itr end, const size_t num_keys, const size_t bpk) {

        uint64_t key;
        
        uint64_t kBitsKeys = 0;
        size_t *num_ones_len, *num_suffix_table;
        num_ones_len = new size_t[num_layers_ + 1]();
        num_suffix_table = new size_t[num_keys]();

        uint64_t prev_prefix = *begin >> (num_layers_ + pruned_layers_);
        uint16_t prev_suffix = *begin >> pruned_layers_ & MASK(num_layers_);
        size_t num_unique_prefix = 1;
        size_t max_num_suffix = 1;
        size_t num_suffix = 1;
        for (auto it = begin + 1; it != end; it++) {
            key = *it;
            uint64_t cur_prefix = key >> (num_layers_ + pruned_layers_);
            uint16_t cur_suffix = key >> pruned_layers_ & MASK(num_layers_); 

            if (prev_prefix != cur_prefix) {
                num_suffix_table[num_unique_prefix - 1] = num_suffix;
                num_unique_prefix++;
                max_num_suffix = max_num_suffix > num_suffix ? max_num_suffix : num_suffix;

                for (size_t i = 0; i < ceil(log2(num_suffix + 1)); i++)
                {
                    num_ones_len[i] += num_suffix >> i & 0x1;
                }

                kBitsKeys += num_suffix < threshold_ ? num_suffix * num_layers_ : threshold_ * num_layers_;
                
                prev_prefix = cur_prefix;
                prev_suffix = cur_suffix;
                num_suffix = 1;
            } else if (cur_suffix != prev_suffix) {
                num_suffix++;
                prev_suffix = cur_suffix;
            }
        }
        num_suffix_table[num_unique_prefix - 1] = num_suffix;
        max_num_suffix = max_num_suffix > num_suffix ? max_num_suffix : num_suffix;
        for (size_t i = 0; i < ceil(log2(num_suffix + 1)); i++)
        {
            num_ones_len[i] += num_suffix >> i & 0x1;
        }
        kBitsKeys += num_suffix < threshold_ ? num_suffix * num_layers_ : threshold_ * num_layers_;
        
          
        keys_tag_len_ = ceil(log2(max_num_suffix + 1));
        kBytesKeys_ = (kBitsKeys + 63) >> 3;

        uint64_t *prefixes;
        uint16_t **suffixes;
        prefixes = new uint64_t[num_unique_prefix]();
        suffixes = new uint16_t*[num_unique_prefix]();
        for (size_t i = 0; i < num_unique_prefix; i++)
        {
            suffixes[i] = new uint16_t[max_num_suffix]();
        }
        
        prev_prefix = *begin >> (num_layers_ + pruned_layers_);
        prev_suffix = *begin >> pruned_layers_ & MASK(num_layers_);
        size_t p_i = 0;
        size_t s_i = 0;
        prefixes[p_i] = prev_prefix;
        suffixes[p_i][s_i] = prev_suffix + 1;
        
        for (auto it = begin + 1; it != end; it++) {
            key = *it;
            uint64_t cur_prefix = key >> (num_layers_ + pruned_layers_);
            uint16_t cur_suffix = key >> pruned_layers_ & MASK(num_layers_);    

            if (prev_prefix != cur_prefix) {
                s_i = 0;
                p_i++;
                prefixes[p_i] = cur_prefix;
                suffixes[p_i][s_i] = cur_suffix + 1;
                prev_prefix = cur_prefix;
                prev_suffix = cur_suffix;
            } else if (cur_suffix != prev_suffix) {
                s_i++;
                suffixes[p_i][s_i] = cur_suffix + 1;
                prev_suffix = cur_suffix;
            }
        }

        size_t max_num_items = num_unique_prefix;
        size_t assoc = 4;
        size_t num_buckets = getNextPow2(std::max<uint64_t>(1, max_num_items / assoc));
        double frac = (double)max_num_items / num_buckets / assoc;
        if (frac > 0.95) {
            num_buckets <<= 1;
        }
        num_buckets_ = num_buckets;

        size_t num_entries = num_buckets_ * bucket_size_;
        size_t kBitsLen = 0;
        for (size_t i = 0; i < keys_tag_len_; i++)
        {
            float proportion = num_ones_len[i] / (float)num_entries;
            proportion = proportion < 0.7 ? proportion : 1 - proportion;
            if (proportion > 0.3) 
                kBitsLen += 1.25 * num_entries;
            else
                kBitsLen += (34.72 * pow(proportion, 3) - 21.85 * pow(proportion, 2) + 7.23 * proportion + 0.02) * num_entries;
                // kBitsLen += 34718296.5 * pow(proportion, 3) - 21846772.9 * pow(proportion, 2) + 7228624.33 * proportion + 19916.5641;
        }
        if (bpk * num_keys < kBitsKeys + kBitsLen)
            fp_len_ = 0;
        else {
            uint64_t kBitsPrefix = bpk * num_keys - kBitsKeys - kBitsLen;
            fp_len_ = kBitsPrefix / (float)num_entries + 1.7;
        }
        if (fp_len_ > 32) {
            // throw std::runtime_error("error, bucket length exceeds 128.");
            fp_len_ = 32;
        } else if (fp_len_ <= 5) {
            //  throw std::runtime_error("error, fp is too short.");
            fp_len_ = 6;
        }

        std::cout << "[+] prefix length: " << fp_len_ << ", suffix length: " << num_layers_ << std::endl;
     
        // shorten 127 to 7
        kBitsPerBucket_ = (fp_len_ - 1) * 4;
        kBytesTags_ = (num_buckets_ * bucket_size_ * (fp_len_ - 1) + 127) >> 3;
        c_packed_tags_ = new char[kBytesTags_]();

        
        value_table_ = new uint64_t*[num_buckets_]();
        uint64_t** keys_table_ = new uint64_t*[num_buckets_]();
        for (size_t i = 0; i < num_buckets_; i++) {
            value_table_[i] = new uint64_t[bucket_size_]();
            keys_table_[i] = new uint64_t[bucket_size_]();
        }

#ifdef ADAPT
        size_t num_words = (num_buckets_ * bucket_size_ + 63) / 64;
        for (size_t i = 0; i < selector_len_; i++)
        {
            hash_selector_multi_[i] = new uint64_t[num_words]();
        }
#endif

        victim_.used = false;


        for (size_t i = 0; i < max_num_items; i++)
        {
            // index = 0 left for empty keys
            if(!Insert(prefixes[i], i + 1, keys_table_, suffixes[i], num_suffix_table[i]))
                throw std::runtime_error("error, failed to insert(fp is too short).");
        }

        Compress(suffixes, keys_table_, num_suffix_table);

        delete[] prefixes;
        for (size_t i = 0; i < num_unique_prefix; i++)
        {
            delete[] suffixes[i];
        }
        delete[] suffixes;
        delete[] num_ones_len;
        delete[] num_suffix_table;

    }

    inline size_t corr_aware_model(int l_s, int l_q, int b, int n, int n_0) {
        int l;
        l = l_s - l_q - 2 + 0.95 * n / n_0 * (b - 1) - log2(0.95 * n / n_0);
        l = l < 0 ? 0 : l;
        l = l > l_s ? l_s : l;
        if (l == 0 || (double)n / n_0 <= (1 << l) / l) {
            if ((double)n / n_0 * (b - l) < 4) 
                l = b > (double)n_0 / n * 4 ? b - (double)n_0 / n * 4 : 0;
            return l;
        }
        auto gap = 2000;
        for (size_t i = 0; i <= l_s; i++)
        {
            auto derivative = 0.95 * (1 << i) + 2 * i - l_s - l_q - 2 + 0.95 * b * n / n_0 + 0.53;
            if (abs(derivative) < gap) {
                l = i;
                gap = abs(derivative);
            }
        }
        return l;
    }


    inline bool LookupImpl_Interval(size_t i, size_t j, uint16_t bound_l, uint16_t bound_r) {
        size_t key_idx = rm.rank(i * bucket_size_ + j) * num_layers_;
        size_t key_num = rm.getBit(i * bucket_size_ + j);
        if (key_num < threshold_) {
            return BinarySearch(key_idx, bound_l, bound_r, key_num);
        } else {
            return BitmapSearch(key_idx, bound_l, bound_r);
        }
    }

    inline bool LookupImpl(size_t i, size_t j, uint16_t bound, bool left) {
        size_t key_idx = rm.rank(i * bucket_size_ + j) * num_layers_;
        size_t key_num = rm.getBit(i * bucket_size_ + j);

        if (key_num < threshold_) {
            if (left) {
                uint16_t max = ReadKeys(key_idx, key_num - 1);
                return bound <= max;
            } else {
                uint16_t min = ReadKeys(key_idx, 0);
                return bound >= min;
            }
        } else {
            if (left)
                return BitmapSearch(key_idx, bound, MASK(num_layers_));
            else
                return BitmapSearch(key_idx, 0, bound);
        }
        return false;
    }

    inline bool BinarySearch(size_t base_pos, uint16_t target_l, uint16_t target_r, size_t key_num) {
        int l = 0, r = key_num - 1;
        while(l <= r) {
            uint16_t key_l = ReadKeys(base_pos, l);
            uint16_t key_r = ReadKeys(base_pos, r);
            if ((key_l > target_r) || key_r < target_l)
                return false;
            if (((key_l <= target_l) && (key_r <= target_r)) || ((key_r >= target_r) && (key_l >= target_l)))
                return true;
            size_t m = l + (r - l) / 2;
            uint16_t key_m = ReadKeys(base_pos, m);
            if (key_m < target_l) {
                l = m + 1;
            } else if (key_m > target_r) {
                r = m - 1;
            } else
                return true;
        }
        return false;
    }

    inline bool BitmapSearch(size_t base_pos, uint16_t l, uint16_t r) {
        size_t begin_pos = base_pos + l;
        size_t count = 0;
        const char *p = c_keys_ + (begin_pos >> 3);
        size_t start = begin_pos % 8;

        while(count < r - l + 1) {
            uint64_t temp;
            if (count + 64 - start > r - l + 1) {
                temp = *((uint64_t *)p) & (MASK(r - l + 1 - count) << start);
                count += r - l + 1 - count;
            } else {
                temp = *((uint64_t *)p) & ~MASK(start);
                count += 64 - start;
            }
            if (__builtin_popcountll(temp) != 0)
                return true;
            start = 0;
            p += 8;
        }       
        return false;
    }

    inline void CleanState() {
        current_selector_ = 0;
        current_value_ = 0;
    }

    inline bool Insert(uint64_t prefix, uint64_t keys_index, uint64_t** keys_table_, const uint16_t* keys, const size_t num_keys)
    {
        CleanState();

        size_t i;
        uint32_t tag;

        if (victim_.used) {
            return false;
        }

        GenerateIndexTagHash(prefix, &i, &tag, &current_value_);
        
        return InsertImpl(i, keys_index, keys_table_, keys, num_keys);
    }

    inline bool Insert(uint64_t prefix)
    {
        CleanState();

        size_t i;
        uint32_t tag;

        if (victim_.used) {
            return false;
        }

        GenerateIndexTagHash(prefix, &i, &tag, &current_value_);
        
        return InsertImpl(i);
    }


    inline bool InsertImpl(const size_t i, uint64_t keys_index, uint64_t** keys_table_, const uint16_t* keys, const size_t num_keys)
    {
        size_t curindex = i;
        uint64_t o_keys_index = 0;

        for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
            bool kickout = count > 0;
            if (InsertTagToBucket(curindex, keys_index, kickout, o_keys_index, keys_table_)) {
                num_items_ += num_keys;
                return true;
            }
            if (kickout) {
                keys_index = o_keys_index;
            }
            curindex = AltIndex(curindex, TagHash(current_value_));
        }

        victim_.index = curindex;
        victim_.tag = TagHash(current_value_, current_selector_);
        for (size_t j = 0; j < num_keys; j++)
        {
            victim_.keys.emplace(*(keys + j) - 1);
        }
        victim_.selector = current_selector_;
        victim_.used = true;
        return true;
    }

    inline bool InsertImpl(const size_t i)
    {
        size_t curindex = i;

        for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
            bool kickout = count > 0;
            if (InsertTagToBucket(curindex, kickout)) {
                return true;
            }
            curindex = AltIndex(curindex, TagHash(current_value_));
        }

        victim_.index = curindex;
        victim_.tag = TagHash(current_value_, current_selector_);
        victim_.selector = current_selector_;
        victim_.used = true;
        return true;
    }

    inline bool InsertTagToBucket(const size_t i, uint64_t keys_index, 
                        const bool kickout, uint64_t &o_keys_index, uint64_t** keys_table_) {
        for (size_t j = 0; j < bucket_size_; j++) {
#ifndef ADAPT
            if (value_table_[i][j] == 0) {
                keys_table_[i][j] = keys_index;
                value_table_[i][j] = current_value_;
                return true;
            }
#endif

#ifdef ADAPT
            size_t r_selector = ReadSelector(i, j);
            uint32_t r_tag = value_table_[i][j] == 0 ? 0 : TagHash(value_table_[i][j], r_selector);
            uint32_t tag = TagHash(current_value_, current_selector_);
            if (r_tag == 0) {
                keys_table_[i][j] = keys_index;
                value_table_[i][j] = current_value_;
                WriteSelector(i, j, current_selector_);
                return true;
            }
            if ((r_selector == current_selector_) && r_tag == tag) {
                //adaptive: update selector and tag, unchange keys
                ++r_selector;
                WriteSelector(i, j, r_selector);
                // change original one too
                ++current_selector_;
            } else if (r_selector != current_selector_) {
                for (size_t k = 0; k <= current_selector_; k++)
                {
                    if (TagHash(current_value_, k) == TagHash(value_table_[i][j], k)) {
                        if (r_selector < current_selector_) {
                            ++r_selector;
                            WriteSelector(i, j, r_selector);
                        } else {
                            ++current_selector_;
                        }
                        break;
                    }
                }
            } 
#endif
        }
        if (kickout) {
            size_t r = rand() % bucket_size_;

            // swap keys
            o_keys_index = keys_table_[i][r];
            keys_table_[i][r] = keys_index;

            // swap value
            uint64_t value = current_value_;
            current_value_ = value_table_[i][r];
            value_table_[i][r] = value;

#ifdef ADAPT
            // swap selector
            size_t selector = current_selector_;
            current_selector_ = ReadSelector(i, r);
            WriteSelector(i, r, selector);
#endif
        }
        return false;
    }

    inline bool InsertTagToBucket(const size_t i, const bool kickout) {
        for (size_t j = 0; j < bucket_size_; j++) {
#ifdef ADAPT
            size_t r_selector = ReadSelector(i, j);
            uint32_t r_tag = value_table_[i][j] == 0 ? 0 : TagHash(value_table_[i][j], r_selector);
            uint32_t tag = TagHash(current_value_, current_selector_);
            if (r_tag == 0) {
                value_table_[i][j] = current_value_;
                WriteSelector(i, j, current_selector_);
                return true;
            }
            if ((r_selector == current_selector_) && r_tag == tag) {
                //adaptive: update selector and tag, unchange keys
                ++r_selector;
                WriteSelector(i, j, r_selector);
                // change original one too
                ++current_selector_;
            } else if (r_selector != current_selector_) {
                for (size_t k = 0; k <= current_selector_; k++)
                {
                    if (TagHash(current_value_, k) == TagHash(value_table_[i][j], k)) {
                        if (r_selector < current_selector_) {
                            ++r_selector;
                            WriteSelector(i, j, r_selector);
                        } else {
                            ++current_selector_;
                        }
                        break;
                    }
                }
            } 
#else
            if (value_table_[i][j] == 0) {
                value_table_[i][j] = current_value_;
                return true;
            }
#endif
        }
        if (kickout) {
            size_t r = rand() % bucket_size_;

            // swap value
            uint64_t value = current_value_;
            current_value_ = value_table_[i][r];
            value_table_[i][r] = value;

#ifdef ADAPT
            // swap selector
            size_t selector = current_selector_;
            current_selector_ = ReadSelector(i, r);
            WriteSelector(i, r, selector);
#endif
        }
        return false;
    }

#ifdef ADAPT
    inline void WriteSelector(size_t i, size_t j, uint64_t selector) {
        size_t pos = i * bucket_size_ + j;
        for (size_t k = 0; k < selector_len_; k++)
        {
            hash_selector_multi_[k][pos / 64] &= ~(1ULL << pos % 64);
            hash_selector_multi_[k][pos / 64] |= (selector >> k & 0x1) << pos % 64;
        }
    }

    inline size_t ReadSelector(size_t i, size_t j) {
        size_t pos = i * bucket_size_ + j;
        size_t sel = 0;
        for (size_t k = 0; k < selector_len_; k++)
        {
            sel += (hash_selector_multi_[k][pos / 64] >> pos % 64 & 0x1) << k;
        }
        return sel;
    }
#endif

    inline void WriteKeys(size_t begin_pos, uint16_t *keys, size_t num_keys, bool raw) {
        if (raw) {
            return WriteKeysImpl1(begin_pos, keys, num_keys);
        } else {
            return WriteKeysImpl2(begin_pos, keys, num_keys);
        }
    }

    inline uint16_t ReadKeys(size_t base_pos, size_t i){
        uint16_t res = 0;
        size_t begin_pos = base_pos + i * num_layers_;
        size_t start = begin_pos % 8;
        const char *p = c_keys_ + (begin_pos >> 3);
        return *((uint64_t*)p) >> start & MASK(num_layers_);
    }

    inline void WriteKeysImpl1(size_t begin_pos, uint16_t *keys, size_t num_keys) {
        for (size_t i = 0; i < num_keys; i++)
        {
            size_t start = begin_pos % 8;
            const char *p = c_keys_ + (begin_pos >> 3);
            // *((uint64_t *)p) &= ~MASK(start + num_layers_) + MASK(start);
            *((uint64_t *)p) |= keys[i] - 1 << start;
            begin_pos += num_layers_;
        }
        
        // for (auto key = keys->begin(); key != keys->end(); key++) {
        //     size_t start = begin_pos % 8;
        //     const char *p = c_keys_ + (begin_pos >> 3);
        //     // *((uint64_t *)p) &= ~MASK(start + num_layers_) + MASK(start);
        //     *((uint64_t *)p) |= *key << start;
        //     begin_pos += num_layers_;
        // }
    }

    inline void WriteKeysImpl2(size_t begin_pos, uint16_t *keys, size_t num_keys) {
        size_t idx = 0, set_pos = 0;
        for (size_t i = 0; i < num_keys; i++)
        {
            idx = (begin_pos + *(keys + i) - 1) >> 3;
            set_pos = (begin_pos + keys[i] - 1) & MASK(3);
            c_keys_[idx] |= 1U << set_pos;
        }
        
        // for (auto it = keys->begin(); it != keys->end(); it++)
        // {
        //     idx = (begin_pos + *it) >> 3;
        //     set_pos = (begin_pos + *it) & MASK(3);
        //     c_keys_[idx] |= 1U << set_pos;
        // }
    }


    inline void GenerateIndexTagHash(uint64_t item, size_t* index,
                                   uint32_t* tag) const {
        // const uint64_t hash = hasher_(item);
        // *index = IndexHash(hash >> 32);
        // *tag = TagHash(hash);

        uint32_t pc = 0, pb = 0;
        HashUtil::BobHash(&item, 8, &pc, &pb);
        *tag = TagHash(pc);
        *index = IndexHash(pb);

        // uint64_t hash = HashUtil::MurmurHash3(&item, 8, 13);
        // *tag = TagHash(hash);
        // *index = IndexHash(hash >> 32);
    }

    inline void GenerateIndexTagHash(uint64_t item, size_t* index,
                                   uint32_t* tag, uint64_t* cur_hash) const {
        // const uint64_t hash = hasher_(item);
        // *cur_hash = hash;
        // *index = IndexHash(hash >> 32);
        // *tag = TagHash(hash);

        uint32_t pc = 0, pb = 0;
        HashUtil::BobHash(&item, 8, &pc, &pb);
        *cur_hash = pc + (((uint64_t)pb)<<32);
        *tag = TagHash(pc);
        *index = IndexHash(pb);

        // *cur_hash = HashUtil::MurmurHash3(&item, 8, 13);
        // *tag = TagHash(*cur_hash);
        // *index = IndexHash(*cur_hash >> 32);
    }

    inline size_t AltIndex(const size_t index, const uint32_t tag) const {
        // NOTE(binfan): originally we use:
        // index ^ HashUtil::BobHash((const void*) (&tag), 4)) & table_->INDEXMASK;
        // now doing a quick-n-dirty way:
        // 0x5bd1e995 is the hash constant from MurmurHash2
        return IndexHash((uint32_t)(index ^ (tag * 0x5bd1e995)));

        // return (index ^ HashUtil::BobHash((const void*) (&tag), 4)) % num_buckets_;
    }

    inline size_t IndexHash(uint32_t hv) const {
        // table_->num_buckets is always a power of two, so modulo can be replaced
        // with
        // bitwise-and:
        return hv & (num_buckets_ - 1);
        // return hv % num_buckets_; 
    }

    // inline uint32_t TagHash(uint32_t hv) const {
    //     uint32_t tag;
    //     tag = hv & ((1ULL << fp_len_) - 1);
    //     tag += (tag == 0);
    //     return tag;
    // }

    inline uint32_t TagHash(uint32_t hv) const {
        uint32_t tag;
#ifdef ADAPT
        tag = hv & ((1ULL << fp_len_) - 1);
        if (tag == 0) {
            tag |= 1ULL << 4;
        }
        return tag;
#endif
        tag = hv & ((1ULL << fp_len_) - 1);
        tag += (tag == 0);
        return tag;
    }

    inline uint32_t TagHash(uint32_t hv, size_t selector) const {
        uint32_t tag;
#ifdef ADAPT
        uint32_t low_bits = hv & 0xf;
        selector &= MASK(selector_len_);
        tag = (hv >> selector * (fp_len_ - 4) + 4) & MASK(fp_len_ - 4);
        if (tag == 0 && low_bits == 0)
            tag++;
        tag = tag << 4 | low_bits;
        // tag += (tag == 0);
        return tag;
#endif

        tag = (hv >> selector * fp_len_) & ((1ULL << fp_len_) - 1);
        tag += (tag == 0);
        return tag;
    }

    inline bool SortPair(uint32_t &a, uint32_t &b) {
        if ((a & 0x0f) > (b & 0x0f)) {
            std::swap(a, b);
            return true;
        }
        return false;
    }

    inline void SortTags(uint32_t *tags, uint32_t *orders) {
        // SortPair(tags[0], tags[2]);
        // SortPair(tags[1], tags[3]);
        // SortPair(tags[0], tags[1]);
        // SortPair(tags[2], tags[3]);
        // SortPair(tags[1], tags[2]);

        if (SortPair(tags[0], tags[2]))
            std::swap(orders[0], orders[2]);
        if (SortPair(tags[1], tags[3]))
            std::swap(orders[1], orders[3]);
        if (SortPair(tags[0], tags[1]))
            std::swap(orders[0], orders[1]);
        if (SortPair(tags[2], tags[3]))
            std::swap(orders[2], orders[3]);
        if (SortPair(tags[1], tags[2]))
            std::swap(orders[1], orders[2]);
    }

    inline void SortPayloads(const size_t i, uint64_t** keys_table_, uint32_t *orders, uint64_t *num_ones) {
#ifdef ADAPT
        /* adjust the orders of hash selectors*/
        size_t pos = i * bucket_size_;
        for (size_t j = 0; j < selector_len_; j++)
        {
            uint8_t o_selectors = 0;
            uint64_t n_selectors = 0;
            o_selectors = hash_selector_multi_[j][pos / 64] >> pos % 64 & 0xf;
            num_ones[j] += (o_selectors & 0x1) + (o_selectors >> 1 & 0x1) + (o_selectors >> 2 & 0x1) + (o_selectors >> 3 & 0x1); 
            n_selectors |= (o_selectors >> orders[0] & 0x1);
            n_selectors |= (o_selectors >> orders[1] & 0x1) << 1;
            n_selectors |= (o_selectors >> orders[2] & 0x1) << 2;
            n_selectors |= (o_selectors >> orders[3] & 0x1) << 3;
            hash_selector_multi_[j][pos / 64] &= ~(0xfULL << pos % 64);
            hash_selector_multi_[j][pos / 64] |= n_selectors << pos % 64;
        }
#endif

        /* the same effect as above*/
        // size_t ss[4] = {};
        // ss[0] = ReadSelector(i, 0);
        // ss[1] = ReadSelector(i, 1);
        // ss[2] = ReadSelector(i, 2);
        // ss[3] = ReadSelector(i, 3);
        // WriteSelector(i, 0, ss[orders[0]]);
        // WriteSelector(i, 1, ss[orders[1]]);
        // WriteSelector(i, 2, ss[orders[2]]);
        // WriteSelector(i, 3, ss[orders[3]]);

        /* adjust the orders of value_table_ and keys_table_*/
        uint64_t values[4] = { value_table_[i][0], value_table_[i][1], value_table_[i][2], value_table_[i][3] };
        uint64_t suffix_index[4] = { keys_table_[i][0], keys_table_[i][1], keys_table_[i][2], keys_table_[i][3] };
        if (orders[0] != 0) {
            value_table_[i][0] = values[orders[0]];
            keys_table_[i][0] = suffix_index[orders[0]];
        }
        if (orders[1] != 1) {
            value_table_[i][1] = values[orders[1]];
            keys_table_[i][1] = suffix_index[orders[1]];
        }
        if (orders[2] != 2) {
            value_table_[i][2] = values[orders[2]];
            keys_table_[i][2] = suffix_index[orders[2]];
        }
        if (orders[3] != 3) {
            value_table_[i][3] = values[orders[3]];
            keys_table_[i][3] = suffix_index[orders[3]];
        }
    }

    inline void SortPayloads(const size_t i, uint32_t *orders, uint64_t *num_ones) {
#ifdef ADAPT
        /* adjust the orders of hash selectors*/
        size_t pos = i * bucket_size_;
        for (size_t j = 0; j < selector_len_; j++)
        {
            uint8_t o_selectors = 0;
            uint64_t n_selectors = 0;
            o_selectors = hash_selector_multi_[j][pos / 64] >> pos % 64 & 0xf;
            num_ones[j] += (o_selectors & 0x1) + (o_selectors >> 1 & 0x1) + (o_selectors >> 2 & 0x1) + (o_selectors >> 3 & 0x1); 
            n_selectors |= (o_selectors >> orders[0] & 0x1);
            n_selectors |= (o_selectors >> orders[1] & 0x1) << 1;
            n_selectors |= (o_selectors >> orders[2] & 0x1) << 2;
            n_selectors |= (o_selectors >> orders[3] & 0x1) << 3;
            hash_selector_multi_[j][pos / 64] &= ~(0xfULL << pos % 64);
            hash_selector_multi_[j][pos / 64] |= n_selectors << pos % 64;
        }
#endif

        uint64_t values[4] = { value_table_[i][0], value_table_[i][1], value_table_[i][2], value_table_[i][3] };
        if (orders[0] != 0) {
            value_table_[i][0] = values[orders[0]];
        }
        if (orders[1] != 1) {
            value_table_[i][1] = values[orders[1]];
        }
        if (orders[2] != 2) {
            value_table_[i][2] = values[orders[2]];
        }
        if (orders[3] != 3) {
            value_table_[i][3] = values[orders[3]];
        }
    }

    inline void WriteBucket(const size_t i, uint32_t tags[4]) {
        /* put in direct bits for each tag*/
        uint8_t lowbits[4];
        uint32_t highbits[4];

        lowbits[0] = tags[0] & 0x0f;
        lowbits[1] = tags[1] & 0x0f;
        lowbits[2] = tags[2] & 0x0f;
        lowbits[3] = tags[3] & 0x0f;

        highbits[0] = tags[0] & 0xfffffff0;
        highbits[1] = tags[1] & 0xfffffff0;
        highbits[2] = tags[2] & 0xfffffff0;
        highbits[3] = tags[3] & 0xfffffff0;

        // note that :  tags[j] = lowbits[j] | highbits[j]

        uint16_t codeword = perm_.encode(lowbits);

        /* write out the bucketbits to its place*/
        const char *p = c_packed_tags_ + ((kBitsPerBucket_ * i) >> 3);
        const size_t dir_bits = kBitsPerBucket_ / 4 - 3;

        if (kBitsPerBucket_ % 8 == 0 || (i & 0x0001) == 0) {
            if (kBitsPerBucket_ <= 64) {
                *((uint64_t *)p) &= MASK(64 - kBitsPerBucket_) << kBitsPerBucket_;
                *((uint64_t *)p) |= codeword | ((uint64_t)highbits[0] << 8) | ((uint64_t)highbits[1] << 8 + dir_bits) |
                            ((uint64_t)highbits[2] << 8 + dir_bits * 2) | ((uint64_t)highbits[3] << 8 + dir_bits * 3);
            } else {
                __uint128_t temp = codeword | ((__uint128_t)highbits[0] << 8) | ((__uint128_t)highbits[1] << 8 + dir_bits) |
                            ((__uint128_t)highbits[2] << 8 + dir_bits * 2) | ((__uint128_t)highbits[3] << 8 + dir_bits * 3);
                *((uint64_t *)p) |= (uint64_t)temp;
                *((uint64_t *)(p + 8)) &= ~MASK(kBitsPerBucket_ - 64);
                *((uint64_t *)(p + 8)) |= temp >> 64;
            }
        } else {
            if (kBitsPerBucket_ <= 64) {
                *((uint64_t *)p) &= (MASK(60 - kBitsPerBucket_) << kBitsPerBucket_ + 4) + 0x000f;
                *((uint64_t *)p) |= ((codeword << 4) | ((uint64_t)highbits[0] << 12) | ((uint64_t)highbits[1] << 12 + dir_bits) |
                            ((uint64_t)highbits[2] << 12 + dir_bits * 2) | ((uint64_t)highbits[3] << 12 + dir_bits * 3));
            } else {
                __uint128_t temp = ((codeword << 4) | ((__uint128_t)highbits[0] << 12) | ((__uint128_t)highbits[1] << 12 + dir_bits) |
                            ((__uint128_t)highbits[2] << 12 + dir_bits * 2) | ((__uint128_t)highbits[3] << 12 + dir_bits * 3));
                *((uint64_t *)p) &= 0x000f;
                *((uint64_t *)p) |= (uint64_t)temp;
                *((uint64_t *)(p + 8)) &= ~MASK(kBitsPerBucket_ - 60);
                *((uint64_t *)(p + 8)) |= temp >> 64;
            }
        }
    }

    inline void WriteBucket(const size_t i, uint32_t tags[4], uint64_t** keys_table_, uint32_t *orders, uint64_t *num_ones, bool sort = true) {
        /* first sort the tags in increasing order is arg sort = true*/
        if (sort) {
            SortTags(tags, orders);

            /* adjust the orders of hash selectors*/
            SortPayloads(i, keys_table_, orders, num_ones);
        }

        WriteBucket(i, tags);
        
    }

    inline void WriteBucket(const size_t i, uint32_t tags[4], uint32_t *orders, uint64_t *num_ones, bool sort = true) {
        /* first sort the tags in increasing order is arg sort = true*/
        if (sort) {
            SortTags(tags, orders);

            /* adjust the orders of hash selectors*/
            SortPayloads(i, orders, num_ones);
        }

        WriteBucket(i, tags);
        
    }

    /* read and decode the bucket i, pass the 4 decoded tags to the 2nd arg
    * bucket bits = 12 codeword bits + dir bits of tag1 + dir bits of tag2 ...
    */
    inline void ReadBucket(const size_t i, uint32_t tags[4]) const {

        const char *p =  c_packed_tags_ + ((kBitsPerBucket_ * i) >> 3);
        const uint32_t kDirBitsMask = MASK(fp_len_ - 4) << 4;
        const size_t dir_bits = kBitsPerBucket_ / 4 - 3;
        uint16_t codeword;
        uint8_t lowbits[4];

        // __uint128_t bucketbits = *((__uint128_t *)p);
        __uint128_t bucketbits = *((uint64_t *)p) + ((__uint128_t)*((uint64_t *)(p + 8)) << 64);

        if (kBitsPerBucket_ % 8 == 0) {

            codeword = bucketbits & 0x0fff;
            tags[0] = ((bucketbits >> 8) & kDirBitsMask);
            tags[1] = ((bucketbits >> 8 + dir_bits) & kDirBitsMask);
            tags[2] = ((bucketbits >> 8 + dir_bits * 2) & kDirBitsMask);
            tags[3] = ((bucketbits >> 8 + dir_bits * 3) & kDirBitsMask);

        } else {
            
            codeword = bucketbits >> ((i & 1) << 2) & 0x0fff;
            tags[0] = (bucketbits >> (8 + ((i & 1) << 2))) & kDirBitsMask;
            tags[1] = (bucketbits >> (8 + dir_bits + ((i & 1) << 2))) & kDirBitsMask;
            tags[2] = (bucketbits >> (8 + dir_bits * 2 + ((i & 1) << 2))) & kDirBitsMask;
            tags[3] = (bucketbits >> (8 + dir_bits * 3 + ((i & 1) << 2))) & kDirBitsMask;
        }

        /* codeword is the lowest 12 bits in the bucket */
        uint16_t v = perm_.dec_table[codeword];
        lowbits[0] = (v & 0x000f);
        lowbits[2] = ((v >> 4) & 0x000f);
        lowbits[1] = ((v >> 8) & 0x000f);
        lowbits[3] = ((v >> 12) & 0x000f);

        tags[0] |= lowbits[0];
        tags[1] |= lowbits[1];
        tags[2] |= lowbits[2];
        tags[3] |= lowbits[3];
    }

    inline void Compress(uint16_t** key_sets, uint64_t** keys_table_, size_t* num_suffix_table) {
        uint64_t keys_pos = 0;

        // allocate memory for keysLen
        uint64_t* ranks[keys_tag_len_];
        uint64_t num_ones[keys_tag_len_] = {};
        uint64_t num_bits = num_buckets_ * bucket_size_ + 63;
        size_t rank_size = num_bits / 64;
        for (size_t i = 0; i < keys_tag_len_; i++)
        {
            ranks[i] = new uint64_t[rank_size];
        }

        // allocate memory for keys
        c_keys_ = new char[kBytesKeys_];

#ifdef ADAPT
        uint64_t num_ones_selector[selector_len_] = {};
        uint64_t num_bits_selector = num_buckets_ * bucket_size_ + 63;
#endif

        
        for (size_t i = 0; i < num_buckets_; i++)
        {
            uint32_t tags[4];
            uint32_t orders[4] = {0, 1, 2, 3};
#ifdef ADAPT
            tags[0] = TagHash(value_table_[i][0], ReadSelector(i, 0));
            tags[1] = TagHash(value_table_[i][1], ReadSelector(i, 1));
            tags[2] = TagHash(value_table_[i][2], ReadSelector(i, 2));
            tags[3] = TagHash(value_table_[i][3], ReadSelector(i, 3));
            WriteBucket(i, tags, keys_table_, orders, num_ones_selector, true);
#else
            tags[0] = TagHash(value_table_[i][0], 0);
            tags[1] = TagHash(value_table_[i][1], 0);
            tags[2] = TagHash(value_table_[i][2], 0);
            tags[3] = TagHash(value_table_[i][3], 0);
            WriteBucket(i, tags, keys_table_, orders, nullptr, true);
#endif
            
            for (size_t j = 0; j < bucket_size_; j++)
            {
                size_t num_keys = 0;
                uint16_t* keys;
                if (keys_table_[i][j] == 0) {
                    num_keys = 0;
                } else {
                    keys = key_sets[keys_table_[i][j] - 1];
                    num_keys = num_suffix_table[keys_table_[i][j] - 1];
                }


                if (num_keys == 0) {}  // skip empty keys 
                else if (num_keys < threshold_) {   // raw keys
                    WriteKeys(keys_pos, keys, num_keys, true);  // insert raw keys
                    keys_pos += num_keys * num_layers_;
                } else {    // bitmap
                    WriteKeys(keys_pos, keys, num_keys, false);  // insert bitmap keys
                    num_keys = threshold_;  // align for key position calculation
                    keys_pos += num_keys * num_layers_;
                }

                // store keys len
                for (size_t k = 0; k < keys_tag_len_; k++)
                {
                    num_ones[k] += num_keys >> k & 0x1;
                    ranks[k][(i * bucket_size_ + j) / 64] |= (num_keys >> k & 0x1) << (i * bucket_size_ + j) % 64;
                }
            } 
        }

        rm = rank_multi(ranks, num_bits, num_ones, keys_tag_len_);
#ifdef ADAPT
        rm_selector_ = rank_multi(hash_selector_multi_, num_bits_selector, num_ones_selector, selector_len_);
#endif
             
    }

    inline void Compress() {
#ifdef ADAPT
        uint64_t num_ones_selector[selector_len_]={};
        uint64_t num_bits_selector = num_buckets_ * bucket_size_ + 63;
        
        for (size_t i = 0; i < num_buckets_; i++)
        {
            uint32_t tags[4];
            uint32_t orders[4] = {0, 1, 2, 3};
            tags[0] = TagHash(value_table_[i][0], ReadSelector(i, 0));
            tags[1] = TagHash(value_table_[i][1], ReadSelector(i, 1));
            tags[2] = TagHash(value_table_[i][2], ReadSelector(i, 2));
            tags[3] = TagHash(value_table_[i][3], ReadSelector(i, 3));
            WriteBucket(i, tags, orders, num_ones_selector, true);
        }
        rm_selector_ = rank_multi(hash_selector_multi_, num_bits_selector, num_ones_selector, selector_len_); 
#else 
        for (size_t i = 0; i < num_buckets_; i++)
        {
            uint32_t tags[4];
            uint32_t orders[4] = {0, 1, 2, 3};
            tags[0] = TagHash(value_table_[i][0], 0);
            tags[1] = TagHash(value_table_[i][1], 0);
            tags[2] = TagHash(value_table_[i][2], 0);
            tags[3] = TagHash(value_table_[i][3], 0);
            WriteBucket(i, tags, orders, nullptr, true);
        }
#endif
    }

    bool Lookup_Adaptive(uint64_t left, uint64_t right) {
        CleanState();
        if (correlated_) {
            uint64_t prefix_l = left >> (num_layers_  + pruned_layers_);
            uint64_t prefix_r = right >> (num_layers_  + pruned_layers_);
            uint16_t bound_l = left >> pruned_layers_ & MASK(num_layers_);
            uint16_t bound_r = right >> pruned_layers_ & MASK(num_layers_);
            if (prefix_l == prefix_r)
                return LookupImpl_Adaptive_Corr(prefix_l, bound_l, bound_r);
            else
                return LookupImpl_Adaptive_Corr(prefix_l, prefix_r, bound_l, bound_r);
        }
        uint64_t prefix_l = left >> num_layers_;
        uint64_t prefix_r = right >> num_layers_;
        if (prefix_l == prefix_r)
            return LookupImpl_Adaptive_Uncorr(prefix_l);
        else
            return LookupImpl_Adaptive_Uncorr(prefix_l, prefix_r);
    }

    bool LookupImpl_Adaptive_Corr(uint64_t prefix, uint16_t bound_l, uint16_t bound_r) {
        size_t i1, i2, j;
        uint32_t tag;

        GenerateIndexTagHash(prefix, &i1, &tag, &current_value_);
        i2 = AltIndex(i1, tag);

        assert(i1 == AltIndex(i2, tag));

        if (victim_.used && (TagHash(current_value_, victim_.selector) == victim_.tag) &&
                (i1 == victim_.index || i2 == victim_.index)) {
            auto it = victim_.keys.upper_bound(bound_l);
            return (it != victim_.keys.end()) && (*it <= bound_r);
        }

        vector<size_t> j1, j2;
        FindTagIndexInBuckets(i1, &j1);
        FindTagIndexInBuckets(i2, &j2);

        for (auto j : j1) {
            if (LookupImpl_Interval(i1, j, bound_l, bound_r)) {
#ifdef ADAPT
                // check remote
                if (value_table_[i1][j] != current_value_) {
                    Adapt(i1, j, value_table_[i1][j]);
                }
#endif
                return true;
            }
        }
        for (auto j : j2) {
            if(LookupImpl_Interval(i2, j, bound_l, bound_r)) {
#ifdef ADAPT
                // check remote
                if (value_table_[i2][j] != current_value_) {
                    Adapt(i2, j, value_table_[i2][j]);
                }
#endif
                return true;
            }
        }
        return false;
    }

    bool LookupImpl_Adaptive_Uncorr(uint64_t prefix) {
        size_t i1, i2, j;
        uint32_t tag;

        GenerateIndexTagHash(prefix, &i1, &tag, &current_value_);
        i2 = AltIndex(i1, tag);

        assert(i1 == AltIndex(i2, tag));

        if (victim_.used && (TagHash(current_value_, victim_.selector) == victim_.tag) &&
                (i1 == victim_.index || i2 == victim_.index)) {
            return true;
        }

#ifdef ADAPT
        vector<size_t> j1, j2;
        bool f = false;
        f |= FindTagIndexInBuckets(i1, &j1);
        f |= FindTagIndexInBuckets(i2, &j2);
        if (f)  {
            Adapt(i1, j1);
            Adapt(i2, j2);
            return true;
        }
        return false;
#endif

#ifndef ADAPT
        return FindTagInBuckets(i1) || FindTagInBuckets(i2);
#endif
    }

#ifdef ADAPT
    inline void Adapt(size_t i, size_t j, uint64_t hash) {
        size_t selector = ReadSelector(i, j);
        ++selector;
        uint32_t tag = TagHash(hash, selector);
        WriteSelector(i, j, selector);
        uint32_t tags[4] = {};
        ReadBucket(i, tags);
        tags[j] = tag;
        WriteBucket(i, tags);

        // we exclude the updating time of hash selector in benchmark using code above, because
        // it only occurs at false positive, if you want to test end-to-end latency, use code below 
        // uint64_t num_ones[selector_len_];
        // rm_selector_.getNumOnes(num_ones);
        // uint64_t num_bits = num_buckets_ * bucket_size_ + 63;
        // size_t selector = rm_selector_.getBit(i * bucket_size_ + j);
        // ++selector;
        // if (selector & 0x1 == 1) {
        //     num_ones[0]++;
        // } else {
        //     for (size_t k = 0; k < selector_len_; k++)
        //     {
        //         if (selector >> k & 0x1) {
        //             num_ones[k]++;
        //             break;
        //         } else {
        //             num_ones[k]--;
        //         }
        //     }
        // }
        // uint32_t tag = TagHash(hash, selector);
        // WriteSelector(i, j, selector);
        // rm_selector_ = rank_multi(hash_selector_multi_, num_bits, num_ones, selector_len_);
        // uint32_t tags[4] = {};
        // ReadBucket(i, tags);
        // tags[j] = tag;
        // WriteBucket(i, tags);
    }

    inline void Adapt(size_t i, vector<size_t> js) {
        for (auto j : js) {
            if (value_table_[i][j] != current_value_) {
                    Adapt(i, j, value_table_[i][j]);
            }
        }
    }

#endif

    bool LookupImpl_Adaptive_Corr(uint64_t prefix_l, uint64_t prefix_r, uint16_t bound_l, uint16_t bound_r) {
        size_t i11, i12, i21, i22;
        uint32_t tag1, tag2;

        GenerateIndexTagHash(prefix_l, &i11, &tag1, &current_value_);
        i12 = AltIndex(i11, tag1);
        assert(i11 == AltIndex(i12, tag1));

        if (victim_.used && (TagHash(current_value_, victim_.selector) == victim_.tag) &&
                (i11 == victim_.index || i12 == victim_.index)) {
            uint16_t max = *(--victim_.keys.end());
            if (bound_l <= max)
                return true;
        }

        vector<size_t> j1, j2;
        FindTagIndexInBuckets(i11, &j1);
        FindTagIndexInBuckets(i12, &j2);

        for (auto j : j1) {
            if (LookupImpl(i11, j, bound_l, true)) {
#ifdef ADAPT
                // check remote
                if (value_table_[i11][j] != current_value_) {
                    Adapt(i11, j, value_table_[i11][j]);
                }
#endif
                return true;
            }
        }
        for (auto j : j2) {
            if(LookupImpl(i12, j, bound_l, true)) {
#ifdef ADAPT
                // check remote
                if (value_table_[i12][j] != current_value_) {
                    Adapt(i12, j, value_table_[i12][j]);
                }
#endif
                return true;
            }
        }

        GenerateIndexTagHash(prefix_r, &i21, &tag2, &current_value_);
        i22 = AltIndex(i21, tag2);
        assert(i21 == AltIndex(i22, tag2));

        if (victim_.used && (TagHash(current_value_, victim_.selector) == victim_.tag) &&
                (i21 == victim_.index || i22 == victim_.index)) {
            uint16_t min = *(victim_.keys.begin());
            if (bound_r >= min)
                return true;
        }

        j1.clear();
        j2.clear();
        FindTagIndexInBuckets(i21, &j1);
        FindTagIndexInBuckets(i22, &j2);

        for (auto j : j1) {
            if (LookupImpl(i21, j, bound_r, false)) {
#ifdef ADAPT
                // check remote
                if (value_table_[i21][j] != current_value_) {
                    Adapt(i21, j, value_table_[i21][j]);
                }
#endif
                return true;
            }
        }
        for (auto j : j2) {
            if(LookupImpl(i22, j, bound_r, false)) {
#ifdef ADAPT
                // check remote
                if (value_table_[i22][j] != current_value_) {
                    Adapt(i22, j, value_table_[i22][j]);
                }
#endif
                return true;
            }
        }

        return false;
    }

    bool LookupImpl_Adaptive_Uncorr(uint64_t prefix_l, uint64_t prefix_r) {
        size_t i11, i12, i21, i22;
        uint32_t tag1, tag2;

        GenerateIndexTagHash(prefix_l, &i11, &tag1, &current_value_);
        i12 = AltIndex(i11, tag1);
        assert(i11 == AltIndex(i12, tag1));

        if (victim_.used && (TagHash(current_value_, victim_.selector) == victim_.tag) &&
                (i11 == victim_.index || i12 == victim_.index)) {
            return true;
        }

#ifdef ADAPT
        vector<size_t> j1, j2;
        bool f = false;
        f |= FindTagIndexInBuckets(i11, &j1);
        f |= FindTagIndexInBuckets(i12, &j2);
        if (f)  {
            Adapt(i11, j1);
            Adapt(i12, j2);
            return true;
        }
#endif

#ifndef ADAPT
        if (FindTagInBuckets(i11) || FindTagInBuckets(i12))
            return true;
#endif
        

        GenerateIndexTagHash(prefix_r, &i21, &tag2, &current_value_);
        i22 = AltIndex(i21, tag2);
        assert(i21 == AltIndex(i22, tag2));

        if (victim_.used && (TagHash(current_value_, victim_.selector) == victim_.tag) &&
                (i21 == victim_.index || i22 == victim_.index)) {
            return true;
        }

#ifdef ADAPT
        j1.clear();
        j2.clear();
        f |= FindTagIndexInBuckets(i21, &j1);
        f |= FindTagIndexInBuckets(i22, &j2);
        if (f)  {
            Adapt(i21, j1);
            Adapt(i22, j2);
            return true;
        }
        else return false;
#endif

#ifndef ADAPT
        return FindTagInBuckets(i21) || FindTagInBuckets(i22);
#endif
        
    }

    inline bool FindTagIndexInBuckets(size_t i, vector<size_t>* j_) {
        bool found = false;
        uint32_t tags[4] = {};
        size_t selector;
        uint32_t rtag;
        ReadBucket(i, tags);
        for (size_t j = 0; j < bucket_size_; j++) {
#ifdef ADAPT
            selector = ReadSelector(i, j);
            // selector = rm_selector_.getBit(bucket_size_ * i + j);
            rtag = tags[j];
            if (rtag == TagHash(current_value_, selector)) {
                j_->emplace_back(j);
                found = true;
            }
#else
            rtag = tags[j];
            if (rtag == TagHash(current_value_, 0)) {
                j_->emplace_back(j);
                found = true;
            }
#endif
        }
        return found;
    }

    inline bool FindTagInBuckets(size_t i) {
        uint32_t tags[4] = {0,0,0,0};
        uint32_t rtag;
        ReadBucket(i, tags);
        for (size_t j = 0; j < bucket_size_; j++) {
            rtag = tags[j];
            if (rtag == TagHash(current_value_, 0)) {
                return true;
            }
        }
        return false;
    }

    /**
     * return Hourglass size in bits
     */
    inline uint64_t size() {
        double prefix_usage_ = kBytesTags_ * 8;
        double suffix_usage_ = kBytesKeys_ * 8;
        double length_usage_ = rm.bitCount();
        double selector_usage_ = 0;
#ifdef ADAPT
        selector_usage_ = rm_selector_.bitCount();
#endif
        return prefix_usage_ + suffix_usage_ + length_usage_ + selector_usage_ + sizeof(*this) * 8;
    }

};