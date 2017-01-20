#include <cmath>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

typedef double function(double);

std::pair<std::vector<double>, std::vector<double>> eulersMethod(function f, double x0, double t0, double h, int maxt);
std::pair<std::vector<double>, std::vector<double>> improvedEulersMethod(function f, double x0, double t0, double h, int maxt);
std::pair<std::vector<double>, std::vector<double>> rungeKuttaMethod(function f, double x0, double t0, double h, int maxt);

void problem3(bool debug = false);
void problem4(bool debug = false);
void problem5(bool debug = false);

int main(){
    problem3();
    problem4();
    problem5();

    return EXIT_SUCCESS;
}

/**
* Euler's Method for numerical approximation of a first order differential equation
*
* xnp1 = xn+f(xn)*deltaT
*
* @param f - x'=f(x)
* @param x0 - the intial condition for x
* @param t0 - the intial condition for t
* @param h - the step size of t
* @param maxt - the maximum value of t
* @returns <x[],t[]>
*/
std::pair<std::vector<double>, std::vector<double>> eulersMethod(function f, double x0, double t0, double h, int maxt){
    auto numSteps = maxt / h;
    std::vector<double> t(numSteps, t0);
    std::vector<double> x(numSteps, x0);
    for(auto k = 1; k < numSteps; k++){
        auto slope = f(x[k-1]);
        t[k] = t[k-1] + h;
        x[k] = x[k-1] + slope * h;
    }

    return std::make_pair(x,t);
}

/**
* Improved Euler's Method for numerical approximation of a first order differential equation
*
* x~ = xn + f(xn)*deltaT
* xnp1 = xn + .5 * (f(xn) + f(x~)) * deltaT
*
* @param f - x'=f(x)
* @param x0 - the intial condition for x
* @param t0 - the intial condition for t
* @param h - the step size of t
* @param maxt - the maximum value of t
* @returns <x[],t[]>
*/
std::pair<std::vector<double>, std::vector<double>> improvedEulersMethod(function f, double x0, double t0, double h, int maxt){
    auto numSteps = maxt / h;
    std::vector<double> t(numSteps, t0);
    std::vector<double> x(numSteps, x0);
    for(auto k = 1; k < numSteps; k++){
        auto slope = f(x[k-1]);
        t[k] = t[k-1] + h;
        auto xT = x[k-1] + slope * h;
        x[k] = x[k-1] + .5 * (slope + f(xT)) * h;
    }

    return std::make_pair(x,t);
}

/**
* Runge-Kutta Method for numerical approximation of a first order differential equation
*
* k1 = f(xn) * deltaT
* k2 = f(xn + .5 * k1) * deltaT
* k3 = f(xn +. 5 * k2) * deltaT
* k4 = f(xn + k3) * deltaT
* xnp1 = xn + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
*
* @param f - x'=f(x)
* @param x0 - the intial condition for x
* @param t0 - the intial condition for t
* @param h - the step size of t
* @param maxt - the maximum value of t
* @returns <x[],t[]>
*/
std::pair<std::vector<double>, std::vector<double>> rungeKuttaMethod(function f, double x0, double t0, double h, int maxt){
        auto numSteps = maxt / h;
    std::vector<double> t(numSteps, t0);
    std::vector<double> x(numSteps, x0);
    for(auto k = 1; k < numSteps; k++){
        auto slope = f(x[k-1]);
        t[k] = t[k-1] + h;
        auto k1 = slope * h;
        auto k2 = f(x[k-1] + .5 * k1) * h;
        auto k3 = f(x[k-1] + .5 * k2) * h;
        auto k4 = f(x[k-1] + k3) * h;
        x[k] = x[k-1] + (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4);
    }

    return std::make_pair(x,t);
}

void problem3(bool debug){
    std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "++++  Problem 3 (Euler's Method)  ++++" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;

    std::cout << "\nx' = -x" << std::endl;
    std::cout << "x  = exp(-t)" << std::endl;   
    auto f = [](double x){ return -1 * x; };
    auto x = std::exp(-1);

    std::cout << "x(1) = " << x << std::endl;

    std::cout << "\nh, E" << std::endl;
    for(auto n = 0; n <= 4; n++){
        auto h = std::pow(10, -1*n);
        std::pair<std::vector<double>, std::vector<double>> x_t = eulersMethod(f, 1, 0, h, 1);
        std::cout << h << ", " << std::abs(x - x_t.first.back()) << std::endl;
        if(debug){
            std::cout << "x, t" << std::endl;
            for(int i = 0; i < x_t.first.size(); i++){
                std::cout << x_t.first[i] << ", " << x_t.second[i] << std::endl;
            }
        }
    }
}

void problem4(bool debug){
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "++++  Problem 4 (Improved Euler's Method)  ++++" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    std::cout << "\nx' = -x" << std::endl;
    std::cout << "x  = exp(-t)" << std::endl;   
    auto f = [](double x){ return -1 * x; };
    auto x = std::exp(-1);

    std::cout << "x(1) = " << x << std::endl;

    std::cout << "\nh, E" << std::endl;
    for(auto n = 0; n <= 4; n++){
        auto h = std::pow(10, -1*n);
        std::pair<std::vector<double>, std::vector<double>> x_t = improvedEulersMethod(f, 1, 0, h, 1);
        std::cout << h << ", " << std::abs(x - x_t.first.back()) << std::endl;
        if(debug){
            std::cout << "x, t" << std::endl;
            for(int i = 0; i < x_t.first.size(); i++){
                std::cout << x_t.first[i] << ", " << x_t.second[i] << std::endl;
            }
        }
    }
}

void problem5(bool debug){
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "++++  Problem 5 (Runge-Kutta Method)  ++++" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    std::cout << "\nx' = -x" << std::endl;
    std::cout << "x  = exp(-t)" << std::endl;   
    auto f = [](double x){ return -1 * x; };
    auto x = std::exp(-1);

    std::cout << "x(1) = " << x << std::endl;

    std::cout << "\nh, E" << std::endl;
    for(auto n = 0; n <= 4; n++){
        auto h = std::pow(10, -1*n);
        std::pair<std::vector<double>, std::vector<double>> x_t = rungeKuttaMethod(f, 1, 0, h, 1);
        std::cout << h << ", " << std::abs(x - x_t.first.back()) << std::endl;
        if(debug){
            std::cout << "x, t" << std::endl;
            for(int i = 0; i < x_t.first.size(); i++){
                std::cout << x_t.first[i] << ", " << x_t.second[i] << std::endl;
            }
        }
    }
}