#include "./sux/bits/EliasFano.hpp"
#include "./sux/bits/Rank9.hpp"

class rank_vector {
    private:
        sux::bits::Rank* ds;
        bool isOne;
        bool isEmpty;

    public:
        rank_vector(const uint64_t *const bits, const uint64_t num_bits, const uint64_t num_one) {
            if (num_one == 0) {
                isEmpty = true;
                return;
            }
            isEmpty = false;
            float ratio = (float)num_one / num_bits;
            if (ratio < 0.3) {
                ds = new sux::bits::EliasFano<>(bits, num_bits);
                isOne = true;
            } else if (ratio > 0.7) {
                uint64_t *reverse_bits = new uint64_t[num_bits / 64 + 1];
                for (size_t i = 0; i < num_bits / 64 + 1; i++)
                {
                    reverse_bits[i] = ~bits[i];
                }
                ds = new sux::bits::EliasFano<>(reverse_bits, num_bits);
                isOne = false;
            } else {
                ds = new sux::bits::Rank9<>(bits, num_bits);
                isOne = true;
            }
        }

        size_t getBit(const size_t pos) { 
            if (isEmpty) 
                return 0;
            else
                return isOne ? ds->getBit(pos) : ds->getBit(pos) ^ 0x1; 
        }

        uint64_t getNumOnes() { return isEmpty ? 0 : ds->getNumOnes(); }

        uint64_t rank(std::size_t pos) {
            if (isEmpty) 
                return 0;
            else
                return isOne ? ds->rank(pos) : ds->rankZero(pos);
        }

        std::size_t bitCount() const {
            if (isEmpty)
                return sizeof(*this) * 8;
            else {
                return ds->bitCount() - sizeof(ds) * 8 + sizeof(*this) * 8;
            }
        }
};

class rank_multi {
    private:
        rank_vector* ds;
        size_t rank_num;

    public:
        rank_multi() {}
        
        rank_multi(const uint64_t *const bits[], const uint64_t num_bits[], const uint64_t num_ones[], const size_t rank_num) {

            this->rank_num = rank_num;

            ds = (rank_vector*)malloc(sizeof(rank_vector) * rank_num);
            for (size_t i = 0; i < rank_num; i++)
                ds[i] = rank_vector(bits[i], num_bits[i], num_ones[i]);
        }

        rank_multi(const uint64_t *const bits[], const uint64_t num_bits, const uint64_t num_ones[], const size_t rank_num) {

            this->rank_num = rank_num;

            ds = (rank_vector*)malloc(sizeof(rank_vector) * rank_num);
            for (size_t i = 0; i < rank_num; i++)
                ds[i] = rank_vector(bits[i], num_bits, num_ones[i]);
        }

        size_t getBit(const size_t pos) { 
            uint64_t r = 0;
            for (size_t i = 0; i < rank_num; i++)
                r += ds[i].getBit(pos) << i;
            return r;
        }

        void getNumOnes(uint64_t *num_ones) { 
            for (size_t i = 0; i < rank_num; i++)
            {
                num_ones[i] = ds[i].getNumOnes();
            }
        }

        uint64_t rank(std::size_t pos) {
            uint64_t r = 0;
            for (size_t i = 0; i < rank_num; i++)
                r += ds[i].rank(pos) << i;
            return r;
        }
        
        std::size_t bitCount() const {
            size_t bc = 0;
            for (size_t i = 0; i < rank_num; i++) {
                bc += ds[i].bitCount();
            }
            return bc + sizeof(*this) * 8 - sizeof(ds) * 8;
        }
};