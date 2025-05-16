#ifndef BINOP_HPP
#define BINOP_HPP

#include <iostream>
#include <vector>
#include <bitset>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>

#ifndef NUM_BITS
#define NUM_BITS 1280
#endif

namespace binop {

    // Fixed size for binary vectors and matrices
    constexpr std::size_t n = NUM_BITS;

    using binvec = std::bitset<n>;         // Alias for binary vector
    using binmat = std::vector<binvec>;    // Alias for binary matrix

    /**
     * @brief Add two binary vectors (binvec) using XOR.
     *
     * @param a The first binary vector (binvec).
     * @param b The second binary vector (binvec).
     * @return The resulting binary vector after XOR (binvec).
     */

    binvec operator+(const binvec& a, const binvec& b);

    /**
     * @brief Add two binary vectors (binvec) inplace using XOR.
     *
     * @param a The first binary vector (binvec).
     * @param b The second binary vector (binvec).
     * @return The resulting binary vector after XOR (binvec).
     */

    binvec& operator+=(binvec& a, const binvec& b);

    /**
     * @brief Add two binary vectors (vector<bool>) using XOR.
     *
     * @param a The first binary vector (vector<bool>).
     * @param b The second binary vector (vector<bool>).
     * @return The resulting binary vector after XOR (vector<bool>).
     */

    std::vector<bool> operator+(const std::vector<bool>& a, const std::vector<bool>& b);
    
    /**
     * @brief Multiply a binary vector (std::vector<bool>) with a binary matrix (binmat).
     *
     * @param v The binary vector (std::vector<bool>).
     * @param M The binary matrix (binmat).
     * @return The resulting binary vector (std::vector<bool>).
     */

    binvec operator*(const std::vector<bool>& v, const binmat& M);

    /**
     * @brief Convert a binary vector (binvec) to a binary vector in std::vector<bool> format.
     *
     * @param b The input binary vector of type `binvec`.
     * @return The converted binary vector in `std::vector<bool>` format.
     */

    std::vector<bool> to_vec(const binvec& b);

    /**
     * @brief Check if two binary vectors (binvec) are orthopodal.
     *
     * This function determines whether two binary vectors are orthopodal
     * by computing their bitwise AND and checking if the result is zero.
     *
     * @param x The first binary vector (binvec).
     * @param y The second binary vector (binvec).
     * @return True if the vectors are orthopodal (bitwise AND is zero), false otherwise.
     */

    bool is_orthopodal(const binvec& x, const binvec& y);

    namespace util {

        /**
         * @brief Find the index of the first set bit in a binary vector.
         *
         * This function scans the binary vector `v` and returns the index of the first bit
         * that is set to 1. If no bits are set, it returns `v.size()`.
         *
         * @param v The binary vector to scan.
         * @return The index of the first set bit, or `v.size()` if all bits are zero.
         */

        std::size_t find_first(const binop::binvec& v);

        /**
         * @brief Find the index of the last set bit in a binary vector.
         *
         * This function scans the binary vector `v` and returns the index of the last bit
         * that is set to 1. If no bits are set, it returns 0.
         *
         * @param v The binary vector to scan.
         * @return The index of the last set bit, or 0 if all bits are zero.
         */

        std::size_t find_last(const binop::binvec& v);


        /**
         * @brief Generates a random binary vector.
         *
         * Populates a binary vector with random 0s and 1s using a
         * Mersenne Twister random number generator.
         *
         * @param x A reference to the binary vector to fill with random values.
         */

        void gen_rnd(binvec& x);

        /**
         * @brief Generates a random binary matrix of dimension d.
         *
         * Creates a binary matrix with `d` rows, where each row is a
         * randomly generated binary vector.
         *
         * @param d The number of rows in the binary matrix.
         * @param X A reference to the binary matrix to populate.
         */

        void gen_rnd(std::size_t d, binmat& X);

        /**
         * @brief Generates a random permutation of integers from 0 to n-1.
         *
         * Creates a random permutation vector of integers that can be
         * used to permute binary vectors or matrices.
         *
         * @param pi A reference to the permutation vector to populate.
         */

        void gen_rnd_perm(std::vector<std::size_t>& pi);

        /**
         * @brief Permutes a binary vector using the provided permutation.
         *
         * Rearranges the elements of the binary vector according to the
         * permutation vector. Throws an exception if the permutation
         * vector size is not equal to `n`.
         *
         * @param x A reference to the binary vector to permute.
         * @param pi The permutation vector of size `n`.
         * @throws std::invalid_argument if pi.size() != n.
         */

        void perm(binvec& x, const std::vector<std::size_t>& pi);

        /**
         * @brief Permutes a binary matrix by permuting each row.
         *
         * Applies the given permutation to each row of the binary matrix.
         * Throws an exception if the permutation vector size is not equal to `n`.
         *
         * @param X A reference to the binary matrix to permute.
         * @param pi The permutation vector of size `n`.
         * @throws std::invalid_argument if pi.size() != n.
         */

        void perm(binmat& X, const std::vector<std::size_t>& pi);

        /**
         * @brief Inversely permutes a binary vector using the provided permutation.
         *
         * Restores the original order of a permuted binary vector using
         * the inverse of the given permutation vector. Throws an exception
         * if the permutation vector size is not equal to `n`.
         *
         * @param x_pi A reference to the permuted binary vector to inverse.
         * @param pi The permutation vector of size `n`.
         * @throws std::invalid_argument if pi.size() != n.
         */

        void inv_perm(binvec& x_pi, const std::vector<std::size_t>& pi);

        /**
         * @brief Inversely permutes a binary matrix by inversely permuting each row.
         *
         * Restores the original order of a permuted binary matrix by applying
         * the inverse of the permutation vector to each row. Throws an exception
         * if the permutation vector size is not equal to `n`.
         *
         * @param X_pi A reference to the permuted binary matrix to inverse.
         * @param pi The permutation vector of size `n`.
         * @throws std::invalid_argument if pi.size() != n.
         */

        void inv_perm(binmat& X_pi, const std::vector<std::size_t>& pi);

        /**
         * @brief Compute the Reduced Row Echelon Form (RREF) of a binary matrix.
         *
         * Uses the Gauss-Jordan elimination algorithm to transform the
         * binary matrix into its Reduced Row Echelon Form (RREF). Operates
         * in-place and assumes the matrix has a fixed number of columns `n`.
         *
         * @param X A reference to the binary matrix to transform into RREF.
         */

        bool get_rref(binmat& X, std::size_t r);

        /**
         * @brief Compute the binomial coefficient C(n, k).
         *
         * @param n The size of the set.
         * @param k The size of the subset.
         * @return The binomial coefficient C(n, k).
         * @throws std::invalid_argument if k > n.
         */

        long double count_combinations(std::size_t n, std::size_t k);

        /**
         * @brief Checks for a collision between two binary vectors at specific positions.
         *
         * This function compares two binary vectors (`e1` and `e2`) at a set of given
         * positions. A collision is detected if the bits at all specified positions are
         * equal in both vectors.
         *
         * @param e1 The first binary vector (std::vector<bool>).
         * @param e2 The second binary vector (std::vector<bool>).
         * @param positions A vector of indices specifying the positions to compare.
         * @return True if all bits at the specified positions are equal, false otherwise.
         */

        bool find_collision(const binop::binvec& e1, const binop::binvec& e2, const std::vector<std::size_t>& positions);

        /**
         * @brief Prints a binary vector to the specified output stream.
         *
         * Converts the binary vector to a string and writes it to the
         * output stream. The default output stream is `std::cout`.
         *
         * @param x The binary vector to print.
         * @param os The output stream to print to (default is std::cout).
         */

        void print(const binvec& x, std::ostream& os = std::cout);

        /**
         * @brief Prints a binary matrix to the specified output stream.
         *
         * Prints each row of the binary matrix on a separate line in the
         * output stream. The default output stream is `std::cout`.
         *
         * @param X The binary matrix to print.
         * @param os The output stream to print to (default is std::cout).
         */

        void print(const binmat& X, std::ostream& os = std::cout);

        /**
         * @brief Class for enumerating binary vectors of lenght d
         * and the Hamming weight r.
         */

        class BinaryEnumerator {
            std::size_t _d;
            std::size_t _r;
            bool _init;
            std::vector<bool> curr_;
        public:

            /**
             * @brief Constructor to initialize enumerator.
             *
             * @param d Length of a binary vector.
             * @param r The Hamming weight of the vector.
             * @throws std::invalid_argument if r > d.
             */

            BinaryEnumerator(std::size_t d, std::size_t r)
                : _d(d), _r(r), _init(false) {
                if (r > d) {
                    throw std::invalid_argument("Hamming weight cannot exceed vector length.");
                }
            }

            /**
             * @brief Reset the enumerator with a new weight.
             *
             * @param r New weight r.
             * @throws std::invalid_argument if r > d.
             */

            void reset(std::size_t r) {
                if (r > _d) {
                    throw std::invalid_argument("Hamming weight cannot exceed vector length.");
                }
                _r = r;
                _init = false;
            }

            /**
             * @brief Generate the next binary vector with the specified Hamming weight.
             *
             * @param x A reference to store the next binary vector.
             * @param i Number of zeros in the offset.
             * @return `true` if a new vector was generated, `false` if all combinations are exhausted.
             */

            bool next(std::vector<bool>& x, std::size_t i = 0);

            /**
             * @brief Enumerate all the binary vectors of radius r in the d-dimensional Hamming space.
             *
             * @return A list of binary vectors of raidus r in the d-dimensional Hamming space.
             */

            std::vector<std::vector<bool>> enumerate();

        };

    } // namespace util

} // namespace binop

#endif // BINOP_HPP
