#ifndef ISD_HPP
#define ISD_HPP

#include <iostream>
#include <vector>
#include <map>
#include "binop.hpp"
#include "lcode.hpp"

#ifndef MAX_NUM_LEVELS
#define MAX_NUM_LEVELS 3
#endif

using binop::operator*;
using binop::operator+;

/**
 * @brief Represents the Information Set Decoding (ISD) algorithm.
 *
 * Provides functionalities to systematize a generator matrix, permute
 * target vectors, and perform error decoding.
 */

class ISD {
protected:
    LCode _C;                ///< Linear code
    binop::binvec _t;        ///< Target vector
    std::string _name;       ///< Algorithm's name

public:

    /**
     * @brief Constructs an ISD instance.
     *
     * Initializes the linear code, target vector, weight, and name.
     * Performs initial systematization and prints the ISD input.
     *
     * @param C Linear code (generator matrix).
     * @param t Target vector.
     * @param w Hamming weight limit.
     * @throws std::invalid_argument If w is negative or greater than the code length n.
     */

    ISD(const LCode& C, const binop::binvec& t, const std::string& name = "Prange")
        : _C(C), _t(t), _name(name) {}

    /**
     * @brief Sets the linear code.
     *
     * Updates the generator matrix associated with the code.
     *
     * @param C Linear code to set.
     */

    void set_C(const LCode& C) { _C = C; }

    /**
     * @brief Sets the target vector.
     *
     * Updates the target vector `t` used in decoding.
     *
     * @param t Target vector to set.
     */

    void set_t(const binop::binvec& t) { _t = t; }

    /**
     * @brief Sets the name of the algorithm.
     *
     * Updates the name of the algorithm.
     *
     * @param name Name of the algorithm to be set.
     */

    void set_name(const std::string& name) { _name = name; }

    /**
     * @brief Gets the linear code.
     *
     * @return A constant reference to the linear code.
     */

    const LCode& get_C() const { return _C; }

    /**
     * @brief Gets the target vector.
     *
     * @return A constant reference to the target vector.
     */

    const binop::binvec& get_t() const { return _t; }

    /**
     * @brief Gets the algorithm name.
     *
     * @return The algorithm name.
     */

    const std::string& get_name() const { return _name; }

    /**
     * @brief Prints the ISD input.
     *
     * Displays the linear code, target vector, and weight limit.
     */

    void print() const {
        std::cout << "\nLinear Code (C):\n";
        _C.print();

        std::cout << "\nTarget (t): ";
        binop::util::print(_t);
    }

    /**
     * @brief Semi-systematizes the generator matrix and target vector.
     *
     * This function permutes the rows and columns of the augmented matrix `[G | t]`,
     * where `G` is the generator matrix and `t` is the target vector, and then
     * transforms it into a reduced row-echelon form (RREF) below a specified row index `k1`.
     * It iteratively adjusts the permutation vector `pi` and retries the process until
     * the RREF condition is successfully achieved.
     *
     * @param pi A permutation vector of size `binop::n`, specifying the order in
     *           which rows and columns should be permuted. Throws an exception if
     *           the size of `pi` is not equal to `binop::n`.
     * @param k1 The row index below which RREF is computed. Default is 0.
     *
     * @throws std::invalid_argument If the size of `pi` is not `binop::n`.
     *
     * @details
     * - The function uses `std::next_permutation` to explore permutations of `pi`
     *   systematically and ensures all possible configurations are attempted.
     * - If all permutations are exhausted without success, the permutation vector
     *   is reset, and the process is retried.
     * - Upon successful systematization, the target vector `t` and generator matrix
     *   `G` are updated with the modified versions from the augmented matrix.
     */

    bool systematize(std::vector<size_t> pi, size_t k1 = 0);

    /**
     * @brief Finds the minimum-weight error vector using a specific decoding strategy.
     *
     * @param e Output binary vector representing the error with the minimum weight.
     * @param pi Permutation vector used in decoding.
     * @return True if the error vector is successfully found; otherwise, false.
     */

    virtual void find_min_error(binop::binvec& e, std::vector<size_t> pi,
        std::map<size_t, std::pair<double, long double>>& stats) = 0;

    /**
     * @brief Predicts the error weight distribution in ISD algorithms.
     *
     * This function estimates the distribution of error weights obtained by
     * an ISD algorithm based on the given profile.
     *
     * @param profile A vector representing the profile used for prediction.
     *                If not provided, a default empty profile is used.
     * @return A vector containing the normalized weight distribution.
     */

    virtual std::vector<long double> errors_dist_predict(std::vector<size_t> profile = {}) = 0;

    /**
     * @brief Virtual destructor.
     */

    virtual ~ISD() = default;

};


/**
 * @brief Represents the Prange decoding algorithm.
 *
 * Implements a basic ISD strategy to find errors using the Prange method.
 */

class Prange : public ISD {
public:

    /**
     * @brief Constructs a Prange instance using the ISD constructor.
     */

    using ISD::ISD;

    /**
     * @brief Finds the minimum-weight error vector using the Prange algorithm.
     *
     * @param e Output binary vector representing the error with the minimum weight.
     * @param pi Permutation vector used in decoding.
     * @return True if the error vector is successfully found; otherwise, false.
     */

    void find_min_error(binop::binvec& e, std::vector<size_t> pi,
        std::map<size_t, std::pair<double, long double>>& stats) override;

    /**
     * @brief Predicts the error weight distribution in Prange's algorithm.
     *
     * This function estimates the distribution of error weights obtained by
     * Prange's algorithm based on the given profile.
     *
     * @param profile A vector representing the profile used for prediction.
     *                If not provided, a default empty profile is used.
     * @return A vector containing the normalized weight distribution.
     */

    std::vector<long double> errors_dist_predict(std::vector<size_t> profile = {}) override;

    /**
    * @brief Virtual destructor.
    */

    virtual ~Prange() = default;

};


/**
 * @brief Represents the Lee-Brickell decoding algorithm (with potential preprocessing).
 *
 * Implements an ISD strategy that splits the error weight into two components
 * (w1 and w2) and performs decoding using using the Lee-Brickell approach.
 */

class LeeBrickell : public ISD {

protected:

    size_t _w2; ///< Secondary weight parameter for Lee-Brickell algorithm.
    bool _pre;  ///< Determines if preprocessing is applied.

public:

    /**
     * @brief Constructs a Lee-Brickell ISD instance.
     *
     * Initializes the Lee-Brickell ISD algorithm with a linear code, target vector,
     * total Hamming weight limit, and secondary weight limit (w2). Inherits the
     * base ISD class to initialize the linear code, target vector, and total weight.
     *
     * @param C Linear code (generator matrix).
     * @param t Target vector.
     * @param w Total Hamming weight limit.
     * @param w2 Secondary weight limit for partial search.
     * @throws std::invalid_argument If w2 < 0 or w2 exceeds w.
     */

    LeeBrickell(const LCode& C, const binop::binvec& t, size_t w2, bool pre = false)
        : ISD(C, t, "LB"), _w2(w2), _pre(pre) {
        if (w2 > binop::n) {
            throw std::invalid_argument("w2 cannot be smaller than 0 and cannot exceed n.");
        }

        if (pre) {
            set_name("LBB");
        }
    }

    /**
     * @brief Sets the secondary weight limit (w2).
     *
     * Ensures that the new value of w2 does not exceed the total weight limit (w)
     * and that it is not smaller than 0.
     *
     * @param w2 The new value for w2.
     * @throws std::invalid_argument If w2 < 0 or w2 exceeds w.
     */

    void set_w2(size_t w2) {
        if (w2 > binop::n) {
            throw std::invalid_argument("w2 cannot be smaller than 0 and cannot exceed n.");
        }
        _w2 = w2;
    }

    /**
     * @brief Gets the secondary weight limit (w2).
     *
     * @return The current value of w2.
     */

    size_t get_w2() const { return _w2; }

    /**
     * @brief Finds the minimum-weight error vector using the Lee-Brickell method.
     *
     * @param e Output binary vector representing the error with the minimum weight.
     * @param pi Permutation vector used in decoding.
     * @return True if the error vector is successfully found; otherwise, false.
     */

    void find_min_error(binop::binvec& e, std::vector<size_t> pi,
        std::map<size_t, std::pair<double, long double>>& stats) override;

    /**
     * @brief Predicts the error weight distribution in Lee-Brickell's algorithm.
     *
     * This function estimates the distribution of error weights obtained by
     * Lee-Brickell's algorithm based on the given profile.
     *
     * @param profile A vector representing the profile used for prediction.
     *                If not provided, a default empty profile is used.
     * @return A vector containing the normalized weight distribution.
     */

    std::vector<long double> errors_dist_predict(std::vector<size_t> profile = {}) override;

    /**
    * @brief Virtual destructor.
    */

    virtual ~LeeBrickell() = default;

};


/**
 * @brief Represents the Stern decoding algorithm (with potential preprocessing).
 *
 * Implements an ISD strategy that splits the error weight into two components
 * (w1 and w2) and performs decoding using using the Stern approach.
 */

class Stern : public ISD {

protected:
    size_t _w2; ///< Secondary weight parameter for the Stern algorithm.
    size_t _c;   ///< Collision parameter.
    bool _pre;  ///< Determines if preprocessing is applied.

public:
    /**
     * @brief Constructs a Stern ISD instance.
     *
     * Initializes the Stern ISD algorithm with a linear code, target vector,
     * total Hamming weight limit, secondary weight limit (w2), and collision parameter.
     * Inherits the base ISD class to initialize the linear code, target vector,
     * and total weight.
     *
     * @param C Linear code (generator matrix).
     * @param t Target vector.
     * @param w Total Hamming weight limit.
     * @param w2 Secondary weight limit for partial search.
     * @param c Collision parameter.
     * @throws std::invalid_argument If w2 < 0, w2 exceeds w, c < 0.
     */

    Stern(const LCode& C, const binop::binvec& t, size_t w2, size_t c, bool pre = false)
        : ISD(C, t, "Stern"), _w2(w2), _c(c), _pre(pre) {
        if (w2 > binop::n) {
            throw std::invalid_argument("w2 cannot be smaller than 0 and cannot exceed n.");
        }

        if (c > binop::n - _C.get_k()) {
            throw std::invalid_argument("c cannot be smaller than 0 and cannot exceed n - k.");
        }

        if (pre) {
            set_name("SternBabai");
        }
    }

    /**
     * @brief Sets the secondary weight limit (w2).
     *
     * Ensures that the new value of w2 does not exceed the total weight limit (w)
     * and that it is not smaller than 0.
     *
     * @param w2 The new value for w2.
     * @throws std::invalid_argument If w2 < 0 or w2 exceeds w.
     */

    void set_w2(size_t w2) {
        if (w2 > binop::n) {
            throw std::invalid_argument("w2 cannot be smaller than 0 and cannot exceed n.");
        }
        _w2 = w2;
    }

    /**
     * @brief Sets the collision parameter.
     *
     * Ensures that the new value of c does not exceed n - k
     * and that it is not smaller than 0.
     *
     * @param c The new value for c.
     * @throws std::invalid_argument If c < 0 or c exceeds n - k.
     */

    void set_c(size_t c) {
        if (c > binop::n - _C.get_k()) {
            throw std::invalid_argument("c cannot be smaller than 0 and cannot exceed w.");
        }
        _c = c;
    }

    /**
     * @brief Gets the secondary weight limit (w2).
     *
     * @return The current value of w2.
     */

    size_t get_w2() const { return _w2; }

    /**
     * @brief Gets the collision parameter (c).
     *
     * @return The current value of c.
     */

    size_t get_c() const { return _c; }

    /**
     * @brief Finds the minimum-weight error vector using the Stern method.
     *
     * @param e Output binary vector representing the error with the minimum weight.
     * @param pi Permutation vector used in decoding.
     * @return True if the error vector is successfully found; otherwise, false.
     */

    void find_min_error(binop::binvec& e, std::vector<size_t> pi,
        std::map<size_t, std::pair<double, long double>>& stats) override;

    /**
     * @brief Predicts the error weight distribution in Stern's algorithm.
     *
     * This function estimates the distribution of error weights obtained by
     * Stern's algorithm based on the given profile.
     *
     * @param profile A vector representing the profile used for prediction.
     *                If not provided, a default empty profile is used.
     * @return A vector containing the normalized weight distribution.
     */

    std::vector<long double> errors_dist_predict(std::vector<size_t> profile = {}) override;

    /**
    * @brief Virtual destructor.
    */

    virtual ~Stern() = default;

};

/**
 * @brief Represents the combination of Wagner's generalized birthday decoding with
 * the Shamir-Schroeppel approach (with potential preprocessing).
 *
 * Implements an ISD strategy that splits the error weight into power of 2 components
 * (w1, w2, ...) and performs decoding using using the Wagner + Shamir-Schroeppel approach.
 */

class Wagner : public ISD {

protected:
    size_t _w2; ///< Secondary weight parameter for the algorithm.
    size_t _c;  ///< Collision parameter.
    size_t _a;  ///< Number of levels of decoding algorithm.
    size_t _m;  ///< Number of lists at the lowest merging level. 
    bool _pre;  ///< Determines if preprocessing is applied.

public:
    /**
     * @brief Constructs a Wagner ISD instance.
     *
     * Initializes the Wagner ISD algorithm with a linear code, target vector,
     * total Hamming weight limit, secondary weight limit (w2), collision parameter
     * number of levels in the decoding algorithm.
     * Inherits the base ISD class to initialize the linear code, target vector,
     * and total weight.
     *
     * @param C Linear code (generator matrix).
     * @param t Target vector.
     * @param w Total Hamming weight limit.
     * @param w2 Secondary weight limit for partial search.
     * @param c Collision parameter.
     * @param a Number of levels in the decoding algorithm.
     * @throws std::invalid_argument If w2 < 0, w2 exceeds w, c < 0.
     */

    Wagner(const LCode& C, const binop::binvec& t, size_t w2, size_t c, size_t a, bool pre = false)
        : ISD(C, t, "Wagner"), _w2(w2), _c(c), _a(a), _pre(pre) {
        if (w2 > binop::n) {
            throw std::invalid_argument("w2 cannot be smaller than 0 and cannot exceed n.");
        }

        if (c > binop::n - _C.get_k()) {
            throw std::invalid_argument("c cannot be smaller than 0 and cannot exceed n - k.");
        }

        if (a > MAX_NUM_LEVELS) {
            throw std::invalid_argument("a is too high.");
        }
        set_m();

        if (pre) {
            set_name("WagnerBabai");
        }
    }

    /**
     * @brief Sets the secondary weight limit (w2).
     *
     * Ensures that the new value of w2 does not exceed the total weight limit (w)
     * and that it is not smaller than 0.
     *
     * @param w2 The new value for w2.
     * @throws std::invalid_argument If w2 < 0 or w2 exceeds w.
     */

    void set_w2(size_t w2) {
        if (w2 > binop::n) {
            throw std::invalid_argument("w2 cannot be smaller than 0 and cannot exceed n.");
        }
        _w2 = w2;
    }

    /**
     * @brief Sets the collision parameter.
     *
     * Ensures that the new value of c does not exceed n - k
     * and that it is not smaller than 0.
     *
     * @param c The new value for c.
     * @throws std::invalid_argument If c < 0 or c exceeds n - k.
     */

    void set_c(size_t c) {
        if (c > binop::n - _C.get_k()) {
            throw std::invalid_argument("c cannot be smaller than 0 and cannot exceed w.");
        }
        _c = c;
    }

    /**
     * @brief Sets the number of levels in the decoding algorithm.
     *
     * Ensures that the number of levels is not too high.
     *
     * @param a The new number of levels.
     * @throws std::invalid_argument If a exceeds MAX_NUM_LEVELS.
     */

    void set_a(size_t a) {
        if (a > MAX_NUM_LEVELS) {
            throw std::invalid_argument("a is too high.");
        }
        _a = a;
    }

    /**
     * @brief Sets the number of lists in the lowest level of the decoding algorithm.
     *
     */

    void set_m() {
        _m = static_cast<size_t>(pow(2, _a));
    }

    /**
     * @brief Gets the secondary weight limit (w2).
     *
     * @return The current value of w2.
     */

    size_t get_w2() const {
        return _w2;
    }

    /**
     * @brief Gets the collision parameter (c).
     *
     * @return The current value of c.
     */

    size_t get_c() const {
        return _c;
    }

    /**
     * @brief Gets the number levels (a).
     *
     * @return The current value of a.
     */

    size_t get_a() const {
        return _a;
    }

    /**
     * @brief Gets the number lists at the lowest level of the algorithm (m).
     *
     * @return The current value of m.
     */

    size_t get_m() const {
        return _m;
    }

    /**
     * @brief Finds the minimum-weight error vector using the Wagner method.
     *
     * @param e Output binary vector representing the error with the minimum weight.
     * @param pi Permutation vector used in decoding.
     * @return True if the error vector is successfully found; otherwise, false.
     */

    void find_min_error(binop::binvec& e, std::vector<size_t> pi,
        std::map<size_t, std::pair<double, long double>>& stats) override;

    /**
     * @brief Predicts the error weight distribution in Wagner's algorithm.
     *
     * This function estimates the distribution of error weights obtained by
     * Wagner's algorithm based on the given profile.
     *
     * @param profile A vector representing the profile used for prediction.
     *                If not provided, a default empty profile is used.
     * @return A vector containing the normalized weight distribution.
     */

    std::vector<long double> errors_dist_predict(std::vector<size_t> profile = {}) override;

    /**
    * @brief Virtual destructor.
    */

    virtual ~Wagner() = default;

};

/**
     * @brief Retrieves a list of positions based on the given parameter k2.
     *
     * This function determines and returns a vector of positions based on the value
     * of `k2` and the internal state of the `Stern` object. If `_pre` is enabled,
     * it extracts positions using bitwise operations on the error set `_C.get_E()`.
     * Otherwise, it generates a sequence of positions using the `_c` parameter.
     *
     * @param k2 The threshold index used to determine which positions to include.
     * @return A vector containing the selected positions.
     */

std::vector<size_t> get_positions(bool pre, const LCode& code, size_t coll, size_t k2);

/**
     * @brief Merges two lists of binary vectors based on collision detection. 
     *
     * Combines vectors from two lists (`L1` and `L2`) by checking for collisions
     * using specified positions. If a collision is found between a vector from `L1`
     * and `L2`, the vectors are combined (bitwise XOR) and added to the output list.
     *
     * @param L1 The first list of binary vectors (std::vector<binop::binvec>).
     * @param L2 The second list of binary vectors (std::vector<binop::binvec>).
     * @param positions Positions at which the merge is performed.
     * @return A vector containing the merged lists (with the collisions at the given positions).
     */

std::vector<binop::binvec> merge_join(const std::vector<binop::binvec>& L1, const std::vector<binop::binvec>& L2,
    const std::vector<size_t>& positions);

template <typename Callback>
void merge_lists(std::vector<std::vector<binop::binvec>>& L, size_t num_levels,
    const std::vector<size_t>& positions, Callback process_pair);

#endif // ISD_HPP