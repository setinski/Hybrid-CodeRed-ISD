#include "hspace.hpp"

std::vector<long double> convolve(const std::vector<long double>& a, const std::vector<long double>& b) {

    std::vector<long double> result(a.size() + b.size() - 1, 0.0L);

    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            if (i >= j && (i - j) < a.size()) {
                result[i] += a[i - j] * b[j];
            }
        }
    }
    return result;
}

std::vector<long double> HSpace::weight_dist_ball(size_t r) {

    // Check r
    if (r > _d) {
        throw std::invalid_argument("Radius exceeds dimension");
    }

    // Vector to store the weight distribution
    std::vector<long double> dist(r + 1);

    // Compute binomial coefficients C(d, i) for i = 0 to _r
    for (size_t i = 0; i <= r; ++i) {
        dist[i] = binop::util::count_combinations(_d, i);
    }

    return dist;
}

std::vector<long double> HSpace::weight_dist_fundamental_ball() {

    size_t half_range = (_d + 1) / 2;

    // Vector to store the weight distribution (initialized with the correct size)
    std::vector<long double> dist(half_range + (_d % 2 == 0 ? 1 : 0));

    // Compute binomial coefficients for i = 0 to half_range
    for (size_t i = 0; i < half_range; ++i) {
        dist[i] = binop::util::count_combinations(_d, i);
    }

    // Handle the middle term for even _d
    if (_d % 2 == 0) {
        dist[half_range] = binop::util::count_combinations(_d, half_range) / 2;
    }

    return dist;
}

std::vector<long double> weight_dist_fundamental_domain(const std::vector<size_t>& profile, size_t k1) {

    // Return an empty vector if the profile is empty
    if (profile.empty()) {
        return {};
    }

    // Store distributions for each profile value
    std::vector<std::vector<long double>> dist_array;
    for (size_t i = 0; i <= k1; ++i) {
        HSpace space(profile[i]);
        dist_array.push_back(space.weight_dist_fundamental_ball());
    }

    // Initialize result with the first distribution
    std::vector<long double> dist_rez = dist_array[0];

    // Convolve all distributions iteratively
    for (size_t i = 1; i < dist_array.size(); ++i) {
        dist_rez = convolve(dist_rez, dist_array[i]);
    }

    return dist_rez;
}
