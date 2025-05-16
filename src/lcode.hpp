#ifndef LCODE_H
#define LCODE_H

#include "binop.hpp"
#include <iostream>
#include <stdexcept>
#include <optional>

/**
 * @brief Represents a linear code using its generator matrix.
 *
 * Provides functionalities to set and retrieve the generator matrix,
 * as well as to print details about the code.
 */

class LCode {

protected:
    binop::binmat _G; ///< Generator matrix defining the linear code.
    binop::binmat _E; ///< Corresponding epipodal matrix.
    binop::binmat _P; ///< Corresponding cumulative projector matrix.
    std::vector<size_t> _profile; ///< Weights of the epipodal vectors.
    size_t _k1; ///< The index of the last epipodal vector of weight 
    /// greater than 1 in sorted epipodal matrix.

public:
    /**
     * @brief Constructs an LCode instance.
     *
     * Initializes the linear code with the provided generator matrix
     * and sets the corresponding cumulative projector matrix and
     * epipodal matrix.
     *
     * @param G Generator matrix for the linear code.
     * @throws std::invalid_argument If the generator matrix is empty.
     */

    explicit LCode(const binop::binmat& G) : _G(G) {
        if (G.empty()) {
            throw std::invalid_argument("Generator matrix is empty.");
        }
        set_P();
        set_E();
    }

    /**
     * @brief Sets the generator matrix.
     *
     * @param G New generator matrix to set.
     */

    void set_G(const binop::binmat& G) {
        if (G.empty()) {
            throw std::invalid_argument("Generator matrix is empty.");
        }
        _G = G;
        set_P();
        set_E();
    }

    /**
     * @brief Sets the cumulative projector matrix.
     *
     */

    void set_P();

    /**
     * @brief Sets the epipodal matrix.
     *
     */

    void set_E();

    /**
     * @brief Sets the profile corresponding to the epipodal matrix.
     *
     */

    void set_profile() {
        _profile.clear();
        for (const auto& e : _E) _profile.push_back(e.count());
    }

    /**
     * @brief Sets the k1 corresponding to the epipodal matrix.
     *
     */

    void set_k1() {

        // Reset k1
        _k1 = 0;

        for (size_t prof : _profile) {
            if (prof > 1)
                ++_k1;
            else
                break;
        }

        if (_k1 > 0)
            --_k1;
    }

    /**
     * @brief Gets the generator matrix.
     *
     * @return A constant reference to the generator matrix.
     */

    const binop::binmat& get_G() const {
        return _G;
    }

    /**
     * @brief Gets the cumulative projector matrix.
     *
     * @return A constant reference to the cumulative projector matrix.
     */

    const binop::binmat& get_P() const {
        return _P;
    }

    /**
     * @brief Gets the epipodal matrix.
     *
     * @return A constant reference to the epipodal matrix.
     */

    const binop::binmat& get_E() const {
        return _E;
    }

    /**
     * @brief Gets the code dimension.
     *
     * @return The dimension of the given code.
     */

    size_t get_k() const {
        return _G.size();
    }

    /**
     * @brief Gets the profile corresponding
     * to the epipodal matrix, namely, the epipodal vectors of the
     * Hamming weight greater than 1.
     *
     * @return The Hamming weight of epipodal vectors of the
     * weight greater than 1.
     */

    std::vector<size_t> get_profile() const {
        return _profile;
    }

    /**
     * @brief Gets the number of epipodal vectors of the Hamming
     * weight greater than 1.
     *
     * @return The number of vectors of the Hamming weight
     * greater than 1.
     */

    size_t get_k1() const {
        return _k1;
    }

    /**
     * @brief Gets the index of the epipodal vector after which
     * weights of the epipodal vector k2 such that the weight of
     * the epipodal vectors between k1 and k2 sums up to the
     * collision constant c. It also updates collision constant c
     * to the newly obtained value.
     *
     * @param c The collision constant c.
     * @return The index of the epipodal vector k2 such that the
     * weight of the epipodal vectors between k1 and k2 sums up
     * to the collision constant c.
     */

    size_t get_k2(size_t& c) const {

        binop::binvec bin_positions;
        for (size_t k2 = _k1; k2 > 0; --k2) {
            bin_positions |= _E[k2];

            size_t c_new = bin_positions.count();
            if (c_new >= c + _k1 - k2 + 1) {
                c = c_new;
                return k2;
            }
        }

        return 0;
    }

    /**
     * @brief Prints details about the linear code.
     *
     * Displays the code length (n), code dimension (k), and
     * the generator matrix.
     */

    void print() const {
        std::cout << "Epipodal matrix: " << "\n";
        binop::util::print(_E);

        std::cout << "Generator matrix: " << "\n";
        binop::util::print(_G);
    }

    /**
     * @brief Reduces a binary vector using the epipodal matrix of the code.
     *
     * This function modifies the binary vector `x` by reducing it with respect to the
     * rows of the epipodal matrix `_E`, which corresponds to the generator matrix _G.
     *
     * The reduction is performed iteratively in reverse order, starting from the last row
     * of the matrix up to the specified index `ind`.
     *
     * For each row of the matrix `_E[i - 1]` up to index `ind`, the function checks
     * if the reduction condition is satisfied:
     * - If the weighted overlap between `x` and the current row of `_E`, adjusted
     *   by the tie-breaking condition, exceeds the weight of the current row, the vector
     *   `x` is updated by adding the corresponding row from the generator matrix _G.
     *
     * @param x The binary vector to be reduced. It is modified in-place during the reduction process.
     * @param ind The index up to which rows of the epipodal matrix `_E` are considered for reduction.
     *
     * @throws std::invalid_argument If `ind` is out of range.
     */

    void size_red(binop::binvec& x, std::size_t ind);

    /**
     * @brief Performs size reduction of the basis vectors i.e. vectors of the generator matrix.
     *
     * This function modifies the generator matrix (`_G`) of the code to obtain a size-reduced basis.
     *
     */

    void size_red_basis();

    /**
     * @brief Semi-systematizes the generator matrix of the code.
     *
     * This function modifies the generator matrix (`_G`) of the code to achieve a semi-systematic form. Rows in the
     * epipodal matrix `E` with a Hamming weight of 1 are moved after rows with a Hamming weight greater than 1 using
     * a bubble sort approach.
     *
     * @note This function assumes that the generator matrix (`_G`) and the epipodal matrix (`_E`) are valid and properly
     *       initialized. The semi-systematization process ensures that the generator matrix reflects the desired ordering
     *       based on the properties of the epipodal matrix.
     *
     */

    void semi_systematize(size_t m);

    /**
     * @brief Performs the Lenstra�Lenstra�Lov�sz (LLL) basis reduction on the generator matrix of the `LCode` instance.
     *
     * This function applies the LLL algorithm to reduce the basis of the generator matrix (`_G`) of the code.
     * The reduction ensures a more orthogonal basis representation, optimizing the generator matrix for certain
     * computational tasks. The function operates iteratively, modifying both the generator matrix (`_G`) and the
     * epipodal matrix (`_E`) to reflect the reduced basis.
     *
     *
     * @details The function performs the following steps:
     * - For each pair of consecutive rows in the generator matrix:
     *   - Compute the projection of the current row onto the basis.
     *   - Evaluate the reduction condition using the Hamming weights of the rows and a tie-breaking function.
     *   - Update and potentially swap rows in both the generator and epipodal matrices to ensure the basis
     *     satisfies the LLL reduction conditions.
     *
     * @note This function assumes that the generator matrix (`_G`) and the epipodal matrix (`_E`) are valid and
     *       properly initialized. The tie-breaking logic ensures deterministic results when reduction conditions
     *       are ambiguous.
     */

    void LLL(std::size_t begin, std::size_t end);

    void epi_sort();

    /**
     * @brief Sorts the columns of the internal binary matrix and updates related structures.
     *
     * This function sorts the columns of the matrix `_G` using `sort_cols`.
     *
     * @parameter k1 The number of rows according to which `_G` is sorted.
     */

    void sort_cols(size_t k1);

    /**
     * @brief Grants the ISD class access to private members.
     *
     * Allows the Information Set Decoding (ISD) class to access the
     * private members of LCode.
     */

    friend class ISD;

};

/**
     * @brief Sorts the columns of the given binary matrix based on the first k1 rows.
     *
     * The sorting follows the rule where in each row, the block of 1s is left-aligned
     * and starts where the previous row�s block ended. Additionally, it adjusts the
     * lower part of the matrix (rows k1 to n-k) to ensure having identity matrix.
     *
     * @param mat The binary matrix to be sorted in-place.
     * @param m The number of rows to process for sorting.
     * @return A vector representing the final column permutation.
     */

std::vector<size_t> sort_cols(binop::binmat& mat, size_t m);

/**
* @brief Project the vector x to the complement of the support for the
* first ind dimensions of the linear code (i.e. to its ind-dimensional
* subspace).
*
* This function modifies the given binary vector `x` by zeroing out bits
* that are covered by the bitwise OR of the first `ind` rows of the
* generator matrix contained in the provided `LCode` instance (`lc`).
*
* @param x The binary vector to be modified.
* @param lc The `LCode` instance containing the generator matrix used
*           to determine the epipodal part of the vector.
* @param ind The index up to which rows of the generator matrix are
*           considered. Must be in the range [0, lc.G_.size()).
*
* @throws std::invalid_argument If `ind` is out of range.
*/

void project(binop::binvec& x, const LCode& lc, std::size_t ind);

#endif // LCODE_H