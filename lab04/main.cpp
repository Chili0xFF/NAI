#include <iostream>
#include <functional>
#include <math.h>
#include <random>
#include <time.h>
#include <algorithm>
std::random_device rd;
std::mt19937 mt_generator(rd());
using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;
auto beal_f = [](double x, double y) {
    double firstPart = pow((1.5-x+(x*y)),2);
    double secondPart = pow(2.25-x+(x*pow(y,2)),2);
    double thirdPart = pow(2.625-x+x*pow(y,3),2);
    return firstPart+secondPart+thirdPart;
};
population_t populate(int pop_size, int chrom_size){
    srand(time(nullptr));
    population_t population;
    for(int i=0;i<pop_size;i++){
        chromosome_t chromosome;
        for(int j=0;j<chrom_size;j++){
            chromosome.push_back(rand()%2);
        }
        population.push_back(chromosome);
    }
    return population;
}
double translate(chromosome_t chromosome){
    using namespace std;
    double result=0;
    reverse(chromosome.begin(),chromosome.end());
    for(int i=0;i<chromosome.size();i++){
        result+=(chromosome.at(i)*pow(2,i));                            //DODAJ TUTAJ ABY MOGŁO DAWAĆ WARTOŚCI UJEMNE i ułamki, podzielic translate na x i y w połwoei
    }
    return result;
}
auto genetic_algorithm = [](
        auto initial_population, auto fitness, auto term_condition,
        auto selection, double p_crossover,
        auto crossover, double p_mutation,  auto mutation) {
    using namespace std;
    uniform_real_distribution<double> uniform(0.0,1.0);
    auto population = initial_population;
    vector<double> population_fit = fitness(population);
    while (!term_condition(population,population_fit)) {
        auto parents_indexes = selection(population_fit);
        decltype(population) new_population;
        for (int i = 0 ; i < parents_indexes.size(); i+=2) {
            decltype(initial_population) offspring = {population[i],population[i+1]};
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring);
            }
            for (auto chromosome : offspring) new_population.push_back(chromosome);
        }
        for (auto & chromosome : new_population) {
            chromosome = mutation(chromosome,p_mutation);
        }
        population = new_population;
        population_fit = fitness(population);
    }
    return population;
};

std::vector<double> fitness_function(population_t population){


    std::vector<double> top10;
    std::vector<double> translated;
    std::vector<double> results;
    for(chromosome_t chrom : population){
        double a = fmod(translate(chrom),4.5);
        translated.push_back(a);
    }

    for(int i=1;i<translated.size();i+=2){
        results.push_back(
                (beal_f(translated.at(i),translated.at(i-1)))
                );
        //std::cout<<results.back()<<std::endl;
    }

    for(int i=0;i<10;i++){
        int best=0;
        int bestId=0;
        for(int j=0;j<results.size();j++){
            if(best<results.at(j)) {
                best = results.at(j);
                bestId = j;
            }
        }
        top10.push_back(best);
        results.at(bestId)=0;
    }
    return top10;
}
std::vector<int> selection_empty(std::vector<double> fitnesses) {
    return {};
}
std::vector<chromosome_t > crossover_empty(std::vector<chromosome_t > parents) {
    return parents;
}
chromosome_t mutation_empty(chromosome_t parents, double p_mutation) {
    return parents;
}


int main() {

    auto himmel_f = [](double x, double y) { return pow(pow(x,2)+y-11,2) + pow(x+pow(y,2)-7,2); };
    auto tcamel_f = [](double x, double y) {
        double firstPart = 2*pow(x,2);
        double secondPart = 1.05*pow(x,4);
        double thirdPart = (pow(x,6))/6;
        double fourthPart = x*y;
        double fifthPart = pow(y,2);
        return (firstPart-secondPart+thirdPart+fourthPart+fifthPart);
    };

    using namespace std;
    population_t population = populate(100,100+(23049%10)*2);
//    for (chromosome_t chromosome: result) {
//        cout << "[";
//        for (int p: chromosome) {
//            cout << p;
//        }
//        cout << "]\n";
//    }
//    cout << endl;
//    auto result = genetic_algorithm(population,
//                                    fitness_function,
//                                    [](auto a, auto b){return true;},
//                                    selection_empty, 1.0,
//                                    crossover_empty,
//                                    0.01, mutation_empty);
    for(chromosome_t chrom : population){
        cout<<translate(chrom)<<endl;
    }
    auto result = fitness_function(population);
    cout<<"top10\n";
    for(double a:result){
        cout<<a<<endl;
    }

    return 0;
}
