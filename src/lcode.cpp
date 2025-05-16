#include "lcode.hpp"
#include "binop.hpp"
#include <stdexcept>
#include <iostream>

using binop::operator+;
using binop::operator+=;
using binop::operator*;

void LCode::set_P() {

    // Resize _P to accommodate _G.size() elements
    _P.resize(_G.size());

    // Initialize the first element of _P with all zeros
    _P[0].set();

    // Compute the remaining elements of _P
    for (std::size_t i = 0; i < _G.size() - 1; ++i)
    {
        _P[i + 1] = _P[i] & ~_G[i];
    }
}

void LCode::set_E() {

    // Resize _P to accommodate _G.size() elements
    _E.resize(_G.size());

    // Compute each entry of _E as the bitwise AND of _G and _P
    for (std::size_t i = 0; i < _E.size(); ++i)
    {
        _E[i] = _G[i] & _P[i];
    }

    // Update profile
    set_profile();

    // Update k1
    set_k1();
}

void LCode::size_red(binop::binvec& x, std::size_t ind) {
    if (ind > _G.size()) {
        throw std::invalid_argument("Index ind cannot be smaller than 0 or bigger than k.");
    }

    // If ind is equal 0, reduction is not performed
    if (ind == 0) {
        return;
    }

    // Perform size reduction
    for (std::size_t i = ind; i > 0; --i) {
        const binop::binvec& e = _E[i - 1];
        if (2 * (x & e).count() > e.count()) {
            x += _G[i - 1];
        }
    }
}

void LCode::size_red_basis() {

    // Perform size reduction: version 4 (from the paper)
    for (size_t i = 1; i < _G.size(); ++i) {
        size_red(_G[i], i - 1);
    }

    // Set the cumulative projector matrix and the epipodal matrix
    set_P();
    set_E();
}

void LCode::epi_sort() {

    std::vector<std::size_t> indices(_G.size());
    bool sorted = false;
    const size_t max_iter = 1000;
    size_t iter = 0;

    while (!sorted) {
        std::iota(indices.begin(), indices.end(), 0);

        // Stable sort indices based on profile sizes
        std::stable_sort(indices.begin(), indices.end(), [=](std::size_t i, std::size_t j) {
            return _profile[i] > _profile[j];
            });

        // Apply sorted order to _G, _E, and _profile
        binop::binmat sorted_G(_G.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            sorted_G[i] = _G[indices[i]];
        }

        // Move sorted vectors back
        _G = std::move(sorted_G);

        // Set the cumulative projector matrix
        set_P();
        set_E();

        // Check if it's indeed sorted properly
        sorted = true;
        for (size_t i = 0; i < _profile.size() - 1; ++i) {
            if (_profile[i + 1] > _profile[i])
                sorted = false;
        }
        ++iter;
    }

    if (iter == max_iter && !sorted)
        throw std::logic_error("Profile is not sorted.");
}

void LCode::semi_systematize(size_t m) {

    binop::util::get_rref(_G, m);

    // Set the cumulative projector matrix
    set_P();
    set_E();
}

void LCode::LLL(std::size_t begin, std::size_t end) {

    // Obtain the LLL reduced basis: version 2 (from Leo's code)
    // Loop invariant: the basis is LLL-reduced from beg to i.
    for (std::size_t i = begin; i < end - 1; ++i) {
        // Check size condition: if size is larger than a treshold, reduce it
        if (((_G[i + 1] ^ _G[i]) & _P[i]).count() < ((_G[i + 1]) & _P[i]).count()) {
            _G[i + 1] += _G[i];
        }

        // Check Lovasz condition: if size is larger than a treshold, swap i+i and i basis' vectors
        // and update projector and epipodal matrix
        if ((_G[i + 1] & _P[i]).count() < (_G[i] & _P[i]).count())
        {
            std::swap(_G[i + 1], _G[i]);

            // Update auxiliary data
            _E[i] = _G[i] & _P[i];
            _P[i + 1] = _P[i] & ~_G[i];

            _E[i + 1] = _G[i + 1] & _P[i + 1];
            _P[i + 2] = _P[i + 1] & ~_G[i + 1];

            if (i > 0)
            {
                --i;
                continue;
            }
        }
    }
}

std::vector<size_t> sort_cols(binop::binmat& mat, size_t m) {

    if (mat.empty() || m > mat.size()) {
        return {};
    }

    // Initialize indices with {0, 1, ..., n-1}
    std::vector<size_t> col_perm(binop::n);
    std::iota(col_perm.begin(), col_perm.end(), 0);

    // Position where the next block of 1s should start
    size_t start_col = 0;

    // Phase 1: Sort top m rows
    for (size_t row = 0; row < m; ++row) {
        size_t ones_count = 0;
        for (size_t j = start_col; j < binop::n; ++j) {
            if (mat[row][col_perm[j]]) {
                std::swap(col_perm[start_col + ones_count], col_perm[j]);
                ones_count++;
            }
        }
        start_col += ones_count;
    }

    // Phase 2: Sort lower part (rows k1 to n-k)
    for (size_t row = m; row < mat.size(); ++row) {
        size_t target_col = binop::n - mat.size() + row;

        if (mat[row][col_perm[target_col]] == 0) {
            for (size_t j = target_col + 1; j < binop::n; ++j) {
                if (mat[row][col_perm[j]] == 1) {
                    std::swap(col_perm[target_col], col_perm[j]);
                    break;
                }
            }
        }
    }

    // Apply the final column order to the matrix using swaps
    for (size_t i = 0; i < mat.size(); ++i) {
        binop::binvec temp;
        for (size_t j = 0; j < binop::n; ++j) {
            temp[j] = mat[i][col_perm[j]];
        }
        mat[i] = temp;
    }

    return col_perm;
}

void LCode::sort_cols(size_t k1) {

    // Sort columns according to k1
    ::sort_cols(_G, k1);

    // Update the cummulative projector and 
    // epipodal matrices
    set_P();
    set_E();
}

void project(binop::binvec& x, const LCode& lc, std::size_t ind) {
    if (ind >= lc.get_k()) {
        throw std::invalid_argument("Invalid index.");
    }

    binop::binvec y = lc.get_G()[0];
    for (std::size_t i = 1; i < ind; ++i) {
        y |= lc.get_G()[i];
    }
    x &= ~y;
}