#include "TNOM.hpp"
#include "Wilcoxon.hpp"
#include <functional>

class HetStats{
    public:
        std::map<string, std::function> statistics ;

        bool initialize_statistics(vector<string> list_stats){

            for (const auto & st: list_stats){
                switch(st){
                    case "TNOM": statistics["TNOM"] = new stats::TNOM() ;
                                 break ;
                    case "WILCOXON" : statistics["WILCOXON"] = new stats::Wilcoxon() ;
                                      break ;
                    default:
                        return false ;
                }
            }
            return true ;
        }
} ;