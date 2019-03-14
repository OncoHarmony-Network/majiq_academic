#ifndef stats_H
#define stats_H
#include "TNOM.hpp"
#include "Wilcoxon.hpp"
#include "Info.hpp"
#include "tTest.hpp"
#include "testStats.hpp"
#include <set>
#include <map>

class HetStats{

    private:

    public:
        HetStats() {}
        ~HetStats() {}
        std::vector<string> names ;
        std::vector<MajiqStats::TestStat *> statistics ;

        bool initialize_statistics(std::vector<string> list_stats){
            sort(list_stats.begin(), list_stats.end()) ;
            std::vector<string>::iterator it;
            it = std::unique (list_stats.begin(), list_stats.end());
            list_stats.resize( std::distance(list_stats.begin(),it) );

            for (it=list_stats.begin(); it!=list_stats.end(); ++it){
                string st = *it ;
                if (st == "TNOM")
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::TNOM()) ;
                if (st == "WILCOXON")
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::Wilcoxon()) ;
                if (st == "INFOSCORE")
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::InfoScore()) ;
                if (st == "TTEST")
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::tTest()) ;
                if (st == "ALL"){
                    statistics.clear() ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::InfoScore()) ;
                    names.push_back("INFOSCORE") ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::TNOM()) ;
                    names.push_back("TNOM") ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::tTest()) ;
                    names.push_back("TTEST") ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::Wilcoxon()) ;
                    names.push_back("WILCOXON") ;
                    break ;
                } else
                    names.push_back(st) ;
            }
            return statistics.size()>0 ;
        }

        int get_number_stats(){
            return statistics.size() ;
        }
} ;

#endif
//} ;