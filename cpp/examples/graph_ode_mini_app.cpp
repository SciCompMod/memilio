#include <vector>
#include <iostream>
#include <random>
#include <nvtx3/nvToolsExt.h>

size_t num_states = 5*1000;

class Node
{
    Node(){
        //            S, E, I, R,   S, E, I, R,     S, E, I, R 
        // population = {1000, 0, 0, 0, 1000, 0, 0, 0, 1000, 0, 0, 0};
        population = std::vector<int>(num_states);
        for (size_t i = 0; i < num_states; i+=5){
            population[i] = 1000;
            population[i+1] = 0;
            population[i+2] = 0;
            population[i+3] = 0;
            population[i+4] = 1000;
        }
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
    size_t num_nodes = 10000;
    std::vector<Node> nodes(num_nodes);

    for (size_t time = 0; time < 1; time ++){
        
        //Infect people (not relevant)
        for(size_t i = 0; i < nodes.size(); i++){
            for(size_t j = 0; j < nodes.size(); j++){
                if(nodes[i].population[j] > 1 && rand() < 0.5 && (i+1)%5 == 0){
                    nodes[i].population[j] -= 1;
                    nodes[i].population[j+1] += 1;
                }
            }
        }

        nvtxRangePushA("Exchange");
        //Exchange people (to be optimized)
        // #pragma acc parallel loop
        for(size_t i = 0; i < num_nodes; i++){
            #pragma acc parallel loop
            for(size_t j = 0; j < num_nodes; j+=2){
                exchange_person(nodes[(i+j) % num_nodes].population, nodes[(i+1+j) % num_nodes].population);
            }
            #pragma acc parallel loop
            for(size_t j = 1; j < num_nodes; j+=2){
                exchange_person(nodes[(i+j) % num_nodes].population, nodes[(i+1+j) % num_nodes].population);
            }
        }
        nvtxRangePop();

    }
}