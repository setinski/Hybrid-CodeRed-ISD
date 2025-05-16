#include "binop.hpp"

namespace binop {

    binvec operator+(const binvec& a, const binvec& b) {
        return a ^ b;
    }

    binvec& operator+=(binvec& a, const binvec& b) {
        a ^= b;
        return a;
    }

    std::vector<bool> operator+(const std::vector<bool>& a, const std::vector<bool>& b) {

        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must be of the same size for bitwise XOR");
        }

        std::vector<bool> result;
        result.reserve(a.size());

        for (std::size_t i = 0; i < a.size(); ++i) {
            result.push_back(a[i] ^ b[i]); // Avoid indexing overhead
        }
        return result;
    }


    binvec operator*(const std::vector<bool>& v, const binmat& M) {
        if (M.empty() || v.size() != M.size()) {
            throw std::invalid_argument("Vector and matrix dimensions are incompatible for multiplication.");
        }

        binvec result;
        result.reset();

        for (std::size_t i = 0; i < M.size(); ++i) {
            if (v[i]) {
                result += M[i];
            }
        }

        return result;
    }

    std::vector<bool> to_vec(const binvec& b) {
        std::vector<bool> vec(binop::n, 0);
        for (size_t i = 0; i < binop::n; ++i) {
            vec[i] = b[i];
        }
        return vec;
    }

    bool is_orthopodal(const binvec& x, const binvec& y) {
        if ((x & y).none()) {
            return true;
        }

        return false;
    }

    namespace util {

        std::size_t find_first(const binop::binvec& v) {
            for (size_t i = 0; i < n; ++i) {
                if (v[i]) {
                    return i;
                }
            }
            return v.size();
        }

        std::size_t find_last(const binop::binvec& v) {
            for (size_t i = n; i > 0; --i) {
                if (v[i]) {
                    return i;
                }
            }
            return 0;
        }

        void gen_rnd(binvec& x) {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            static std::uniform_int_distribution<int> dis(0, 1);

            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] = dis(gen);
            }
        }

        void gen_rnd(std::size_t d, binmat& X) {
            X.clear();
            binvec x;
            for (std::size_t i = 0; i < d; ++i) {
                gen_rnd(x);
                X.push_back(x);
            }
        }

        void gen_rnd_perm(std::vector<std::size_t>& pi) {
            static std::random_device rd;
            static std::mt19937 gen(rd());

            pi.resize(n);
            std::iota(pi.begin(), pi.end(), 0);
            std::shuffle(pi.begin(), pi.end(), gen);
        }

        void perm(binvec& x, const std::vector<std::size_t>& pi) {
            if (pi.size() != n) {
                throw std::invalid_argument("Permutation vector size must be n.");
            }

            binvec x_pi;
            for (std::size_t i = 0; i < n; ++i) {
                x_pi[pi[i]] = x[i];
            }
            x = std::move(x_pi);
        }

        void perm(binmat& X, const std::vector<std::size_t>& pi) {
            for (auto& row : X) {
                perm(row, pi);
            }
        }

        void inv_perm(binvec& x_pi, const std::vector<std::size_t>& pi) {
            if (pi.size() != n) {
                throw std::invalid_argument("Permutation vector size must equal n.");
            }

            binvec x;
            for (std::size_t i = 0; i < n; ++i) {
                x[i] = x_pi[pi[i]];
            }
            x_pi = std::move(x);
        }

        void inv_perm(binmat& X_pi, const std::vector<std::size_t>& pi) {
            for (auto& row : X_pi) {
                inv_perm(row, pi);
            }
        }

        bool get_rref(binmat& X, size_t r) {
            size_t row_num = X.size() - 1; //Due to specific needs of ISD
            size_t col_num = n;

            // Tracks the pivot position
            std::size_t pivot_pos = col_num - row_num + r;

            size_t row = r;
            while (pivot_pos < col_num) {

                // Find the first row with a non-zero entry in the pivot_pos column
                size_t i = row;
                while (i < row_num && X[i][pivot_pos] == 0) {
                    ++i;
                }

                // If there is no 1 in column corresponding to pivot_pos, the submatrix is not full rank
                if (i == row_num) {
                    return false;
                }

                // Swap row `i` with row `row` if i != row
                if (i != row) {
                    std::swap(X[i], X[row]);
                }

                // Eliminate all other entries in the lead column
                for (std::size_t j = 0; j < row_num + 1; ++j) {
                    if (j != row && X[j][pivot_pos] == 1) {
                        X[j] ^= X[row]; // Binary addition (XOR) to eliminate the entry
                    }
                }

                // Move to the next column
                ++pivot_pos;
                ++row;
            }

            // Check the rank of the last part of matrix
            size_t rank = 0;
            for (size_t i = r; i < row_num; ++i) {
                rank += X[i][col_num - row_num + i];
            }
            if (rank < row_num - r) {
                return false;
            }

            return true;
        }

        long double count_combinations(std::size_t n, std::size_t k) {
            if (k > n) {
                throw std::invalid_argument("k cannot be greater than n in count_combinations.");
            }

            long double result = 1.0L;
            for (std::size_t i = 1; i <= k; ++i) {
                result *= static_cast<long double>(n - i + 1) / static_cast<long double>(i);
            }
            return result;
        }

        bool find_collision(const binop::binvec& e1, const binop::binvec& e2, const std::vector<std::size_t>& positions) {
            // Verify positions
            if (!std::all_of(positions.begin(), positions.end(), [=](std::size_t pos) { return pos < n; })) {
                throw std::invalid_argument("Positions need to be in range [0, n].");
            }

            // Check for collisions at the given positions
            for (std::size_t pos : positions) {
                if (e1[pos] != e2[pos]) {
                    return false;
                }
            }

            // If collisions exist at all given positions
            return true;
        }

        void print(const binvec& x, std::ostream& os) {
            std::string bit_string = x.to_string();
            std::reverse(bit_string.begin(), bit_string.end());
            os << bit_string << '\n';
        }

        void print(const binmat& X, std::ostream& os) {
            for (const auto& row : X) {
                print(row, os);
            }
        }

        bool BinaryEnumerator::next(std::vector<bool>& x, std::size_t i) {
            if (!_init) {
                // Initialize the first combination
                curr_.assign(_d, false);
                std::fill(curr_.begin() + i, curr_.begin() + _r + i, true);
                _init = true;
            }
            else {
                // Generate the next permutation
                if (!std::prev_permutation(curr_.begin() + i, curr_.end())) {
                    _init = false;
                    return false; // All combinations are exhausted
                }
            }

            x = curr_;
            return true;
        }

        std::vector<std::vector<bool>> BinaryEnumerator::enumerate() {

            std::vector<std::vector<bool>> L;

            std::vector<bool> l(_d, false);
            for (size_t i = 0; i <= _r; ++i) {
                while (next(l)) {
                    L.push_back(l);
                }
                reset(i);
            }

            return L;
        }
    } // namespace util
} // namespace binop
