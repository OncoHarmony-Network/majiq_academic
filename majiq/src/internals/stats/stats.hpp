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
        std::set<string> names ;
        std::vector<MajiqStats::TestStat *> statistics ;

        bool initialize_statistics(std::vector<string> list_stats){

            for (const auto & st: list_stats){
                if (st == "TNOM" && names.count("TNOM") == 0)
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::TNOM()) ;
                if (st == "WILCOXON" && names.count("WILCOXON") == 0)
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::Wilcoxon()) ;
                if (st == "INFOSCORE" && names.count("INFOSCORE") == 0)
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::InfoScore()) ;
                if (st == "TTEST" && names.count("TTEST") == 0)
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::tTest()) ;
                if (st == "ALL"){
                    statistics.clear() ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::TNOM()) ;
                    names.insert("TNOM") ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::Wilcoxon()) ;
                    names.insert("WILCOXON") ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::InfoScore()) ;
                    names.insert("INFOSCORE") ;
                    statistics.push_back((MajiqStats::TestStat*) new MajiqStats::tTest()) ;
                    names.insert("TTEST") ;
                    break ;
                } else
                    names.insert(st) ;
            }
            return statistics.size()>0 ;
        }

        int get_number_stats(){
            return statistics.size() ;
        }
} ;

#endif
//} ;