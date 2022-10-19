#include <iostream>
#include <functional>
#include <math.h>
#include <random>
#include <time.h>
using namespace std;
using domain_t = std::vector<double>;
std::random_device rd;
std::mt19937 mt_generator(rd());

void result(double time, double result){
    using namespace std;
    cout<<"\nTime: "<<time<<" Result: "<<result;
}

void climb(auto function, auto domain, int maxIterations=1000){

}

vector<double> neighbours(){
    vector<double> neightbours
}

void anneal(auto function, auto domain, int maxIterations=1000, int step=5){
    double V;
    std::uniform_real_distribution<double> dist(domain.at(0), domain.at(1));
    vector<double[]> sk;
    sk[0]=dist(mt_generator),dist(mt_generator);
    double rand1 = dist(mt_generator);
    double rand2 = dist(mt_generator);
    sk.at(0)=function(rand1, rand2);



    //czym jest K-.ilośc iteracji, czym jest uk->losujemy za każdym razem gdy potrzebujemy, czym jest tk

    for(int k=1;k<maxIterations;k+=step){
        double Tk = 1.0/step;           //Wyznaczamy temperature

        vector<double> NK = neighbours();

        double rand1 = dist(mt_generator);
        double rand2 = dist(mt_generator);
        if(function(rand1,rand2)>=sk){

        }
    }

}



int main() {
    auto beal_f = [](double x, double y) {
        double firstPart = pow((1.5-x+(x*y)),2);
        double secondPart = pow(2.25-x+(x*pow(y,2)),2);
        double thirdPart = pow(2.625-x+x*pow(y,3),2);
        return firstPart+secondPart+thirdPart;
    };
    auto himmel_f = [](double x, double y) { return pow(pow(x,2)+y-11,2) + pow(x+pow(y,2)-7,2); };
    auto tcamel_f = [](double x, double y) {
        double firstPart = 2*pow(x,2);
        double secondPart = 1.05*pow(x,4);
        double thirdPart = (pow(x,6))/6;
        double fourthPart = x*y;
        double fifthPart = pow(y,2);
        return (firstPart-secondPart+thirdPart+fourthPart+fifthPart);
    };
    vector<double> domain={-4.5,4.5};
    for(int i=0;i<20;i++){
        climb(beal_f,domain,10000);
    }

    return 0;
}
