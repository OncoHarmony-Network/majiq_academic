#ifndef TEST_STATS_H
#define TEST_STATS_H
#include <algorithm>
#include <string>

using namespace std ;
namespace stats{
    class TestStat{
        public:
            TestStat() ;
            ~TestStat() ;
            virtual double Calc_pval(vector<float> data, vector<int> labels) { return 0 ;}
    } ;
}
#endif
//} ;