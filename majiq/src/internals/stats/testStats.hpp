#ifndef TEST_STATS_H
#define TEST_STATS_H
#include <algorithm>
#include <string>


namespace MajiqStats{

    struct tTRecord{
        int neg ;
        int pos ;
        double score ;

        bool operator<(const tTRecord& rhs) const{
                return (neg < rhs.neg || (neg == rhs.neg && pos < rhs.pos) ||
                        (neg == rhs.neg && pos == rhs.pos && score < rhs.score)) ;
        }

        bool operator==(const tTRecord& rhs) const{
                return (neg == rhs.neg && pos == rhs.pos && score == rhs.score) ;
        }

    } ;



    template <class T>
    inline void hash_combine(std::size_t & s, const T & v){
        std::hash<T> h;
        s^= h(v) + 0x9e3779b9 + (s<< 6) + (s>> 2);
    }


    struct hash_fn{

        size_t operator()(tTRecord const& s) const {
            size_t res = 0;
            hash_combine(res,s.neg);
            hash_combine(res,s.pos);
            hash_combine(res,s.score);
            return res;
        }
    } ;
    class TestStat{
        public:
            TestStat() {}
            ~TestStat() {}
            virtual double Calc_pval(std::vector<float>& data, std::vector<int>& labels) { return 0 ;}
    } ;
}

#endif
