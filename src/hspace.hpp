#ifndef HSPACE_H
#define HSPACE_H

#include <iostream>
#include <vector>
#include "binop.hpp"

/**
 * @brief Represents a Hamming space of dimension d.
 *
 */

class HSpace {

protected:
    size_t _d; ///< Hamming space dimension

public:

    /**
     * @brief Constructs an HSpace instance.
     *
     * @param d Dimension of the Hamming space.
     */

    explicit HSpace(size_t d) {
        _d = d;
    }

    /**
     * @brief Sets the Hamming space dimension.
     *
     * @param d Dimension of the Hamming space.
     */

    void set_d(size_t d) {
        _d = d;
    }

    /**
     * @brief Gets the Hamming space dimension.
     *
     * @return The Hamming space dimension.
     */

    size_t get_d() const {
        return _d;
    }

    /**
     * @brief Computes the weight distribution of all vectors within a Hamming ball of radius r.
     *
     * This function calculates the binomial coefficients C(d, i) for all i in [0, r],
     * representing the number of binary vectors with weight i in d-dimensional space.
     *
     * @param r Radius of the Hamming ball.
     * @return std::vector<long double> containing the weight distribution.
     */

    std::vector<long double> weight_dist_ball(size_t r);

    /**
     * @brief Computes the weight distribution of a fundamental ball of a
     * d-dimensional Hammin space.
     *
     * @return A vector containing the weight distribution of the fundamental ball.
     */

    std::vector<long double> weight_dist_fundamental_ball();

};

/**
 * @brief Performs convolution of two discrete distributions.
 *
 * This function takes two input vectors representing discrete distributions
 * and returns their convolution. The convolution result is obtained by
 * multiplying and summing corresponding elements, resulting in a new
 * distribution representing the combined effect of the inputs.
 *
 * @param a The first input vector representing a discrete distribution.
 * @param b The second input vector representing a discrete distribution.
 * @return A vector containing the convolved distribution as long doubles.
 */

std::vector<long double> convolve(const std::vector<long double>& a, const std::vector<long double>& b);

/**
* @brief Computes the weight distribution of a Hamming fundamental domain.
*
* This function calculates the number of binary vectors in a fundamental domain
* of d-dimensional Hamming space.
*
* @return A vector containing the weight distribution of a fundamental domain.
*/

std::vector<long double> weight_dist_fundamental_domain(const std::vector<size_t>& profile, size_t k1);

#endif // HSPACE_H