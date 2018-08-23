#ifndef TEST_STATS_H
#define TEST_STATS_H
#include <algorithm>
#include <string>


namespace MajiqStats{

    class TestStat{
        public:
            TestStat() {}
            ~TestStat() {}
            virtual double Calc_pval(std::vector<float> data, std::vector<int> labels) { return 0 ;}
    } ;
}

#endif
