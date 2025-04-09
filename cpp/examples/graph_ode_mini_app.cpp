#include <vector>
#include <iostream>
#include <random>
#include <nvtx3/nvToolsExt.h>


class Node
{
    Node(){
        //            S, E, I, R,   S, E, I, R,     S, E, I, R 
        population = {1000, 0, 0, 0, 1000, 0, 0, 0, 1000, 0, 0, 0};
    }



    public:
    std::vector<int> population;
};

void exchange_person(std::vector<int>& from, std::vector<int>& to){
    for(size_t i = 0; i < from.size(); i++){
        if (from[i] > 1){
            from[i] -= 1;
            to[i] += 1;
        }
    }
}

int main(){
    //Initialise nodes
    std::vector<Node> nodes(4000);

    for (size_t time = 0; time < 1; time ++){
        
        //Infect people (not relevant)
        for(size_t i = 0; i < nodes.size(); i++){
            for(size_t j = 0; j < nodes.size(); j++){
                if(nodes[i].population[j] > 1 && rand() < 0.5 && (i+1)%4 == 0){
                    nodes[i].population[j] -= 1;
                    nodes[i].population[j+1] += 1;
                }
            }
        }

        nvtxRangePushA("Exchange");
        //Exchange people (to be optimized)
        for(size_t i = 0; i < nodes.size(); i++){
            #pragma acc parallel loop
            for(size_t j = 0; j < nodes.size(); j++){
                exchange_person(nodes[i].population, nodes[j].population);
            }
        }
        nvtxRangePop();

    }
    // for(size_t j = 0; j < 100; j++){
    //     std::vector<double> a(1000000, 1.0);
    //     std::vector<double> b(1000000, 1.0);
    //     std::vector<double> c(1000000);
    //     #pragma acc 
    //     for (size_t i = 0; i < a.size(); i++){
    //         c[i] = a[i] + b[i];
    //     }

    //     double result = 0;
    //     for (size_t i = 0; i < c.size(); i++){
    //         result += c[i];
    //     }
    // }
//37 mit OpenACC, 38 ohne

}