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

void exchange_person(Node from, Node to){
    for(size_t i = 0; i < from.population.size(); i++){
        if (from.population[i] > 1){
            from.population[i] -= 1;
            to.population[i] += 1;
        }
    }
}

int main(){
    //Initialise nodes
    std::vector<Node> nodes(4000);

    for (size_t time = 0; time < 28; time ++){
        
        nvtxRangePushA("Exchange");
        //Infect people (not relevant)
        for(size_t i = 0; i < nodes.size(); i++){
            for(size_t j = 0; j < nodes[i].population.size()-1; j++){
                if(nodes[i].population[j] > 1 && rand() < 0.5 && (i+1)%4 == 0){
                    nodes[i].population[j] -= 1;
                    nodes[i].population[j+1] += 1;
                }
            }
        }
        nvtxRangePop();

        //Exchange people (to be optimized)
        for(size_t i = 0; i < nodes.size(); i++){
            for(size_t j = 0; j < nodes.size(); j++){
                exchange_person(nodes[i], nodes[j]);
            }
        }
    }


}