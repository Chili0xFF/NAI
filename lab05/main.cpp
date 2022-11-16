#include <iostream>
#include <functional>
#include <math.h>
#include <random>
#include <time.h>
#include <algorithm>
#include <vector>
std::random_device rd;
std::mt19937 mt_generator(rd());
using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;
using result_t = std::vector<std::vector<double>>;
using fitness_f = std::function<double(const chromosome_t &)>;
using term_condition_f = std::function<bool(const population_t &, const std::vector<double> &)>;
using selection_f = std::function<int(const std::vector<double> &)>;
using crossover_f = std::function<std::vector<chromosome_t>(const std::vector<chromosome_t> &)>;
using mutation_f = std::function<chromosome_t(const chromosome_t, double)>;


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
    bool flagNegative = false;
    if(chromosome.at(0)==1)flagNegative=true;

    double twos=1;
    for(int i=2;i>=1;--i){
        result+=(chromosome.at(i)*twos);                            //DODAJ TUTAJ ABY MOGŁO DAWAĆ WARTOŚCI UJEMNE i ułamki, podzielic translate na x i y w połwoei
        twos *= 2.0;
    }
    twos = 2;
    for (int i = 3; i < chromosome.size(); ++i) {
        result += (chromosome.at(i)/twos);
        twos *= 2.0;
    }

    if(flagNegative)result*=-1;
    return result;
}
int PopulationIndex=0;

void writePopulation(population_t vector1, std::vector<double> vector2);

population_t genetic_algorithm(population_t initial_population,
                               fitness_f fitness,
                               term_condition_f term_condition,
                               selection_f selection,
                               double p_crossover, crossover_f crossover,
                               double p_mutation, mutation_f mutation) {
    using namespace std;
    uniform_real_distribution<double> uniform(0.0, 1.0);
    auto population = initial_population;
    vector<double> population_fit(population.size());
    transform(population.begin(), population.end(), population_fit.begin(),fitness);
    while (!term_condition(population, population_fit)) {
        writePopulation(population,population_fit);
        vector<int> parents_indexes(population.size());
        population_t new_population(population.size());
        // calculate fitness
        transform(population_fit.begin(), population_fit.end(),
                  parents_indexes.begin(),
                  [&](auto e) { return selection(population_fit); });
        // perform crossover operations
        for (int i = 0; i < parents_indexes.size() - 1; i += 2) {
            vector<chromosome_t> offspring = {population[parents_indexes[i]], population[parents_indexes[i + 1]]};
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring);
            }
            new_population[i] = offspring[0];
            new_population[i + 1] = offspring[1];
        }
        for (auto &chromosome : new_population) {
            chromosome = mutation(chromosome, p_mutation);
        }
        population = new_population;
        std::transform(population.begin(), population.end(), population_fit.begin(),
                       fitness);
    PopulationIndex++;
    }
    return population;
}

void writePopulation(population_t vector1, std::vector<double> vector2) {
    std::cout<<"Population"<<PopulationIndex<<": \n";1
    for(int i=0;i<vector1.size();i++){
        std::cout<<vector2.at(i)<<" ";
    }
    std::cout<<std::endl;
};


int selection_roulette(std::vector<double> fitnesses);

double fitness(const chromosome_t chromosome){
    using namespace std;
    vector<double> translated;
    chromosome_t chrom_a, chrom_b;
    for (int i = 0; i < (chromosome.size()/2)-1; ++i) {
        chrom_a.push_back(chromosome.at(i));
    }
    for (int i=chromosome.size()/2;i<chromosome.size()-1;i++){
        chrom_b.push_back(chromosome.at(i));
    }
    double a = translate(chrom_a);
    double b = translate(chrom_b);
    double fit = 1/(beal_f(a,b));
    return fit;
}

std::vector<chromosome_t > crossover_1p(std::vector<chromosome_t > parents) {
    chromosome_t child1, child2;
    for(int i=0;i<parents.at(0).size()/2;i++){
        child1.push_back(parents.at(0).at(i));
    }
    for(int i=parents.at(0).size()/2;i<parents.at(0).size();i++){
        child1.push_back(parents.at(1).at(i));
    }
    for(int i=0;i<parents.at(0).size()/2;i++){
        child2.push_back(parents.at(1).at(i));
    }
    for(int i=parents.at(0).size()/2;i<parents.at(0).size();i++){
        child2.push_back(parents.at(0).at(i));
    }
    return {child1,child2};
}
chromosome_t mutation_singular(chromosome_t chrom, double p_mutation) {
    std::uniform_real_distribution<double> uniform(0.0001, 1.0);
    chromosome_t afterMutations = chrom;
    for(int i=0;i<afterMutations.size();i++){
        if(p_mutation>=uniform(mt_generator)){
            afterMutations.at(i)=1-afterMutations.at(i);
        }
    }
    return afterMutations;
}
std::ostream &operator<<(std::ostream &o, const chromosome_t chromosome) {
    for (const int p : chromosome) {
        o << p;
    }
    return o;
}
std::ostream &operator<<(std::ostream &o,
                         std::pair<population_t, fitness_f> pop) {
    for (const auto p : pop.first) {
        o << "{" << p << " " << (pop.second(p)) << "} ";
    }
    return o;
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
    population_t population = populate(1000,118);
    population_t result = genetic_algorithm(
            population, fitness,
            [](auto a, auto b) {
                static int i = 0;
                i++;
                //cout << i << ": " << make_pair(a, fitness) << endl;
                return i > 10;
            },
            selection_roulette, 0.1, crossover_1p, 0.5, mutation_singular);

    //cout << make_pair(result, fitness);
    //cout << endl;
    return 0;
}

int selection_roulette(std::vector<double> fitnesses) {
    int resultIndex;
    double sum=0;
    for(int i=0;i<fitnesses.size();i++){
        sum+=fitnesses.at(i);
    }
    double tempsum = 0;
    std::uniform_real_distribution<double> uniform(0.0, sum);
    std::vector<double> probabilities;
    double random = uniform(mt_generator);
    for(int i=0;i<fitnesses.size();i++){
        tempsum+= fitnesses.at(i);
        if(tempsum>=random){
            resultIndex=i;
            break;
        }
    }
    return resultIndex;
}
/*
 * // zaczynamy od tego
#include <algorithm>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
std::random_device rd;
std::mt19937 mt_generator(rd());
using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;
using fitness_f = std::function<double(const chromosome_t &)>;
using term_condition_f =
        std::function<bool(const population_t &, const std::vector<double> &)>;
using selection_f = std::function<int(const std::vector<double> &)>;
using crossover_f =
        std::function<std::vector<chromosome_t>(const std::vector<chromosome_t> &)>;
using mutation_f = std::function<chromosome_t(const chromosome_t, double)>;
population_t genetic_algorithm(population_t initial_population,
                               fitness_f fitness,
                               term_condition_f term_condition,
                               selection_f selection, double p_crossover,
                               crossover_f crossover, double p_mutation,
                               mutation_f mutation) {
    using namespace std;
    uniform_real_distribution<double> uniform(0.0, 1.0);
    auto population = initial_population;
    vector<double> population_fit(population.size());
    transform(population.begin(), population.end(), population_fit.begin(),fitness);
    while (!term_condition(population, population_fit)) {
        vector<int> parents_indexes(population.size());
        population_t new_population(population.size());
        // calculate fitness
        transform(population_fit.begin(), population_fit.end(),
                  parents_indexes.begin(),
                  [&](auto e) { return selection(population_fit); });
        // perform crossover operations
        for (int i = 0; i < parents_indexes.size() - 1; i += 2) {
            vector<chromosome_t> offspring = {population[parents_indexes[i]], population[parents_indexes[i + 1]]};
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring);
            }
            new_population[i] = offspring[0];
            new_population[i + 1] = offspring[1];
        }
        for (auto &chromosome : new_population) {
            chromosome = mutation(chromosome, p_mutation);
        }
        population = new_population;
        std::transform(population.begin(), population.end(), population_fit.begin(),
                       fitness);
    }
    return population;
};
using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;
int selection_empty(std::vector<double> fitnesses) {
    return 0;
}
int selection_tournament_2(std::vector<double> fitnesses) {
    std::uniform_int_distribution<int> uniform(0, fitnesses.size()-1);
    int a = uniform(mt_generator);
    int b = uniform(mt_generator);
    return (fitnesses[a]>fitnesses[b])?a:b;
}
std::vector<chromosome_t> crossover_empty(std::vector<chromosome_t> parents) {
    return parents;
}
std::vector<chromosome_t> crossover_two_point(std::vector<chromosome_t> parents) {
    using namespace std;
    uniform_int_distribution<int> locus(0,parents.at(0).size()-1);
    int a = locus(mt_generator);
    int b = locus(mt_generator);
    if (a > b) swap(a,b);
    auto children = parents;
    for (int i = a; i < b; i++) {
        swap(children[0][i],children[1][i]);
    }
    return children;
}
chromosome_t mutation_empty(const chromosome_t parent, double p_mutation) {
    return parent;
}
chromosome_t mutation_one_point(const chromosome_t parent, double p_mutation) {
    using namespace std;
    uniform_real_distribution<double> uni(0.0,1.0);
    // mutation??
    if (uni(mt_generator) < p_mutation) {
         uniform_int_distribution<int> locus(0,parent.size()-1);
         chromosome_t child = parent;
         auto l = locus(mt_generator);
         child[l] = 1 - child[l];
        return child;
    } else
        return parent;
}

double fitness_function(const chromosome_t &chromosome) {
    return std::accumulate(chromosome.begin(), chromosome.end(), 0);
}
std::vector<chromosome_t> generate_initial_population(int n) {
    std::vector<chromosome_t> ret(n);
    std::uniform_int_distribution<int> uniform(0, 1);
    std::transform(ret.begin(), ret.end(), ret.begin(), [&](auto e) {
        chromosome_t c(8);
        for (int i = 0; i < c.size(); i++) c[i] = uniform(mt_generator);
        return c;
    });
    return ret;
}
std::ostream &operator<<(std::ostream &o, const chromosome_t chromosome) {
    for (const int p : chromosome) {
        o << p;
    }
    return o;
}
std::ostream &operator<<(std::ostream &o,
                         std::pair<population_t, fitness_f> pop) {
    for (const auto p : pop.first) {
        o << "{" << p << " " << (pop.second(p)) << "} ";
    }
    return o;
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
    population_t population = populate(10000,100+(23049%10)*2);
    auto result = genetic_algorithm(
            population, fitness,
            [](auto a, auto b) {
                static int i = 0;
                i++;
                cout << i << ": " << make_pair(a, fitness_function) << endl;
                return i > 10;
            },
            selection_roulette, 0.5, crossover_1p, 0.01, mutation_singular);

    cout make_pair(result, fitness_function);
    cout << endl;
    return 0;

}
 */