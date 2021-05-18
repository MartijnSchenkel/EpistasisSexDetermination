#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

#include <random>

int randomize();
extern std::mt19937 rng;

// random integer [0,n)

template <typename T>
T rn(const T n){
    static std::uniform_int_distribution<T> d{0,n-1};
    return d(rng);
}


// random uniform [0,1)
double ru();

// random standard normal
double rnorm(const double &mean, const double &stddev);

// random binary
bool r2(const double& p = 0.5);

// random Poisson
int rpois(const double&);

// random exponential
double rexp(const double&);


extern std::mt19937 rng;


template <typename T1 = double, typename T2 = size_t>
class ridx{
public:
    ridx(const std::vector<T1>& w) : w(w) {}
    T2 operator()() const {
        static std::discrete_distribution<T2> d(w.begin(),w.end());
        return d(rng);
    }
private:
    const std::vector<T1>& w;
};




#endif // RANDOMNUMBERS_H
