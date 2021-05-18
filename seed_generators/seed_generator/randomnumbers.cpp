#include "randomnumbers.h"

std::mt19937 rng;

int randomize() {
    static std::random_device rd{};
    auto seed = rd();
    rng.seed(seed);
    return seed;
}

// random double [0,1)
double ru()
{
    std::uniform_real_distribution<> d{};
    return d(rng);
}

// random standard normal
double rnorm(const double &mean, const double &stddev)
{
    std::normal_distribution<> d{mean, stddev};
    return d(rng);
}

// random bernoulli {0,1}
bool r2(const double& p)
{
    static std::bernoulli_distribution d{};
    using parm_t = decltype(d)::param_type;
    return d(rng, parm_t{p});
}

int rpois(const double &lambda)
{
    std::poisson_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d(rng, parm_t{lambda});
}

double rexp(const double &lambda)
{
    std::exponential_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d(rng, parm_t{ lambda });
}

