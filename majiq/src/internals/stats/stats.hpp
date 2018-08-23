#ifndef stats_H
#define stats_H
#include "TNOM.hpp"
#include "Wilcoxon.hpp"
#include "Info.hpp"
#include "testStats.hpp"
#include <map>

class HetStats{
    public:
        HetStats() {}
        ~HetStats() {}

        std::vector<MajiqStats::TestStat *> statistics ;

        bool initialize_statistics(std::vector<string> list_stats){

            for (const auto & st: list_stats){
                if (st == "TNOM" ) statistics.push_back((MajiqStats::TestStat*) new MajiqStats::TNOM()) ;
                if (st == "WILCOXON") statistics.push_back((MajiqStats::TestStat*) new MajiqStats::Wilcoxon()) ;
                if (st == "INFOSCORE") statistics.push_back((MajiqStats::TestStat*) new MajiqStats::InfoScore()) ;
                if (st == "ALL"){
                    statistics.clear() ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::TNOM()) ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::Wilcoxon()) ;
//                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::InfoScore()) ;
                    break ;
                }
            }
            return statistics.size()>0 ;
       }
} ;

#endif
//} ;