#ifndef stats_H
#define stats_H
#include "TNOM.hpp"
#include "Wilcoxon.hpp"
#include "testStats.hpp"
#include <map>

class HetStats{
    public:
        HetStats() {};
        ~HetStats() {};

//        std::map<string, std::function> statistics ;
        std::vector<stats::TestStat *> statistics ;

        bool initialize_statistics(vector<string> list_stats){

            for (const auto & st: list_stats){
                if (st == "TNOM" ) statistics.push_back((stats::TestStat*) new stats::TNOM()) ;
                if (st == "WILCOXON") statistics.push_back((stats::TestStat*) new stats::Wilcoxon()) ;


//                    case "TNOM": statistics["TNOM"] = new stats::TNOM() ;
//                                 break ;
//                    case "WILCOXON" : statistics["WILCOXON"] = new stats::Wilcoxon() ;
//                                      break ;
//                    default:
//                        return false ;
                }
            return true ;
       }
} ;

#endif
//} ;