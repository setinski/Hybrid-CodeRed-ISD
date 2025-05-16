#include "hspace.hpp"
#include "isd.hpp"
#include <thread>
#include <mutex>
#include <unordered_map>

using binop::operator*;
using binop::operator+;

binop::binvec get_mask(const std::vector<size_t>& positions) {

    binop::binvec mask;
    for (auto pos : positions) {
        if (pos < mask.size()) {
            mask.set(pos);
        }
    }
    return mask;
}

std::pair<binop::binvec, binop::binvec> get_hash(const binop::binvec& e, const binop::binvec& mask) {

    binop::binvec key = e & mask;
    return { key, e };
}

std::vector<size_t> get_positions(bool pre, const LCode& code, size_t coll, size_t k2) {

    std::vector<size_t> positions;

    if (pre) {
        binop::binvec bin_positions;
        for (size_t i = k2; i <= code.get_k1(); ++i) {
            bin_positions |= code.get_E()[i];
        }

        for (size_t i = 0; i < binop::n; ++i) {
            if (bin_positions[i]) {
                positions.push_back(i);
            }
        }

        for (size_t j = 0; j <= code.get_k1() - k2; ++j)
            positions.pop_back();
    }
    else {
        positions.resize(coll);
        std::iota(positions.begin(), positions.end(), binop::n - code.get_k() - coll);
    }

    return positions;

}

std::vector<binop::binvec> merge_join(const std::vector<binop::binvec>& L1, const std::vector<binop::binvec>& L2,
    const std::vector<size_t>& positions) {

    std::vector<binop::binvec> L;

    // Create mask
    binop::binvec mask = get_mask(positions);

    // Obtain hash_map
    std::unordered_multimap<binop::binvec, binop::binvec> hash_map;
    for (const binop::binvec& elem : L1) {
        hash_map.insert(get_hash(elem, mask));
    }

    // Merge the two lists
    for (const auto& l2_elem : L2) {
        auto hash_el = get_hash(l2_elem, mask);
        auto bucket = hash_map.equal_range(hash_el.first);

        for (auto it = bucket.first; it != bucket.second; ++it) {
            L.push_back(hash_el.second + it->second);
        }
    }

    return L;
}

template <typename Callback>
void merge_lists(std::vector<std::vector<binop::binvec>>& L, size_t num_levels, 
    const std::vector<size_t>& positions, Callback process_pair) {
    
    // Merging list hiearchically
    size_t level = num_levels;

    while (level > 1) {
        // Set the merge positions
        size_t end_pos = (num_levels - level +1) * positions.size() / num_levels;
        std::vector<size_t> pos(positions.begin(), positions.begin() + end_pos);

        // Update list L with the new (merged) list
        size_t num_pairs = static_cast<size_t>(pow(2, level)) / 2;
        std::vector<std::vector<binop::binvec>> L_new;
        for (size_t j = 0; j < num_pairs; ++j) {
            L_new.push_back(merge_join(L[2 * j], L[2 * j + 1], pos));
        }
        L = std::move(L_new);

        // Move to the previous level
        --level;
    }

    std::vector<binop::binvec> L_final;
    L_final = merge_join(L[0], L[1], positions);
    
    // Process each element of the list
    for (binop::binvec& l : L_final) {
        process_pair(l);
    }
}

bool ISD::systematize(std::vector<std::size_t> pi, std::size_t k1) {

    // Ensure the permutation vector has the correct size
    if (pi.size() != binop::n) {
        throw std::invalid_argument("The length of the permutation vector pi must be equal to binop::n.");
    }

    // Create an augmented generator matrix by appending the target vector _t
    binop::binmat augmented_G = _C._G;
    augmented_G.push_back(_t);

    // Flag to indicate if a valid permutation is found
    bool success = false;

    // Counter for the number of attempts
    std::size_t attempts = 0, max_attempts = 10000;

    // Try different permutations until a successful one is found 
    // or attempts are exhausted
    while (!success && attempts < max_attempts) {

        // Generate the previous lexicographical permutation of pi
        std::prev_permutation(pi.begin(), pi.end());

        // Apply the permutation to the augmented generator matrix
        binop::util::perm(augmented_G, pi);

        // Attempt to put the matrix in row-reduced echelon form (RREF) 
        // with the required number of pivot elements
        success = binop::util::get_rref(augmented_G, k1);

        // Increment the attempt counter
        attempts++;
    }

    // If a valid permutation was found, update the class members
    if (success) {

        // Extract and update the target vector _t from the
        //  augmented matrix
        _t = augmented_G.back();
        augmented_G.pop_back();

        // Update the generator matrix with the modified version
        _C._G = std::move(augmented_G);

        // Set the corresponding cumulative projector matrix
        // and the epipodal matrix
        _C.set_P();
        _C.set_E();
        return true;
    }
    else {
        return false;
    }
}

void Prange::find_min_error(binop::binvec& e, std::vector<std::size_t> pi,
    std::map<std::size_t, std::pair<double, long double>>& stats) {

    // Prepare the generator matrix in systematic form
    if (!systematize(pi)) {
        throw std::runtime_error("Failed to systematize the generator matrix.");
    }

    // Calculate the candidate error vector
    e = _t;

    // Write error statistics
    stats[e.count()].first = 1;
}

std::vector<long double> Prange::errors_dist_predict(std::vector<size_t> profile) {

    // Determine the code dimension
    size_t k = _C.get_k();

    // Create an HSpace object representing the Hamming space of dimension (n - k)
    HSpace space(binop::n - k);

    // Compute the weight distribution of a Hamming ball of radius (n - k) in that space
    std::vector<long double> dist = space.weight_dist_ball(binop::n - k);

    // Compute the normalization coefficient (total number of cosets: 2^(n-k))
    long double coeff = std::pow(2.0L, static_cast<long double>(binop::n - k));
    std::vector<long double> dist_d(dist.size());

    // Normalize each entry by dividing by the total number of cosets
    for (size_t i = 0; i < dist.size(); ++i) {
        dist_d[i] = dist[i] / coeff;
    }

    return dist_d;
}

void LeeBrickell::find_min_error(binop::binvec& e, std::vector<std::size_t> pi,
    std::map<std::size_t, std::pair<double, long double>>& stats) {

    // Prepare the generator matrix in systematic form
    if (!systematize(pi)) {
        throw std::runtime_error("Failed to systematize the generator matrix.");
    }

    // Special case (reverts back to Prange)
    if (_w2 == 0) {
        e = _t;
        stats[e.count()].first = 1;

        return;
    }

    // Define parameters for enumerator and Hamming ball traversal
    size_t k = _C.get_k(), k1 = 0;

    // If Babai is required, preprocess
    if (_pre) {

        //// LLL basis reduction
        //_C.LLL(0, _C.get_k());
        
        // Sort based on the length of epipodal vectors
        _C.epi_sort();

        // Put the matrix in a nice form
        k1 = _C.get_k1();
        _C.sort_cols(k1);
    }
    std::cout << "\nk1: " << k1 << "\n";

    // Explore vectors in the Hamming ball
    std::vector<bool> x(k);
    binop::util::BinaryEnumerator enumerator(k, _w2);

    e.set();
    for (size_t i = 0; i <= _w2; ++i) {
        while (enumerator.next(x, k1)) {

            // Calculate the candidate error vector
            binop::binvec candidate = _t + (x * _C.get_G());

            // Size-reduce e
            if (_pre) {
                _C.size_red(candidate, k1);
            }

            // Writes down statistics
            stats[candidate.count()].first++;

            // Updating e
            if (candidate.count() < e.count()) e = candidate;
        }
        enumerator.reset(i);
    }
}

std::vector<long double> LeeBrickell::errors_dist_predict(std::vector<size_t> profile) {

    // Determine the code dimension
    size_t k = _C.get_k(), k1 = 0;

    // If the final vectors are preprocessed and reduced,
    // obtain corresponding k1
    if (_pre) {
        k1 = _C.get_k1();
    }

    // Compute the weight distribution of a Hamming ball of radius (n - k) or
    // fundamental domain distribution in that space (for the given profile)
    std::vector<long double> dist1;
    if (!_pre) {
        HSpace space1(binop::n - k);
        dist1 = space1.weight_dist_ball(binop::n - k);
    }
    else {
        dist1 = weight_dist_fundamental_domain(profile, k1);
    }

    // Create an HSpace object representing the Hamming space of dimension k - k1
    HSpace space2(k - k1);

    // Compute the weight distribution of a Hamming ball in that space
    std::vector<long double> dist2 = space2.weight_dist_ball(_w2);

    // Compute joined distribution
    std::vector<long double> dist = convolve(dist1, dist2);

    // Compute the normalization coefficient (total number of cosets: 2^(n-k))
    long double coeff = std::pow(2.0L, static_cast<long double>(binop::n - k));
    std::vector<long double> dist_d(dist.size());

    // Normalize each entry by dividing by the total number of cosets
    for (size_t i = 0; i < dist.size(); ++i) {
        dist_d[i] = dist[i] / coeff;
    }

    return dist_d;
}

void Stern::find_min_error(binop::binvec& e, std::vector<std::size_t> pi,
    std::map<std::size_t, std::pair<double, long double>>& stats) {

    // Prepare the generator matrix in systematic form
    if (!systematize(pi)) {
        throw std::runtime_error("Failed to systematize the generator matrix.");
    }

    // Special case (reverts back to Prange)
    if (_w2 == 0) {
        e = _t;
        stats[e.count()].first = 1;

        return;
    }

    // Define parameters for enumerator and Hamming ball traversal
    size_t k = _C.get_k(), k1 = 0, k2 = 0;
    std::vector<bool> l(k, false);
    std::vector<std::vector<binop::binvec>> L_t;

    if (_pre) {

        //// LLL basis reduction
        //_C.LLL(0, _C.get_k());

        // Sort based on the length of epipodal vectors
        _C.epi_sort();

        // Update values of k1, k2 and c
        size_t c_new = get_c();
        k1 = _C.get_k1();
        k2 = _C.get_k2(c_new);
        set_c(c_new);

        // Put the matrix in a nice form
        _C.sort_cols(k2);
    }
    std::cout << "\nk1: " << k1 << ", k2: " << k2 << ", c_new: " << get_c() << "\n";

    // Obtain the first list
    binop::util::BinaryEnumerator enumerator1((k - k2) / 2, _w2);
    std::vector<std::vector<bool>> L1 = enumerator1.enumerate();
    std::vector<binop::binvec> L1_t(L1.size());
    for (size_t i = 0; i < L1.size(); ++i) {
        std::copy(L1[i].begin(), L1[i].end(), l.begin() + k2);
        L1_t[i] = l * _C.get_G();
        l.assign(l.size(), false);
    }
    L_t.push_back(L1_t);

    // Obtain the second list
    binop::util::BinaryEnumerator enumerator2(k - k2 - (k - k2) / 2, _w2);
    std::vector<std::vector<bool>> L2 = enumerator2.enumerate();
    std::vector<binop::binvec> L2_t(L2.size());
    for (size_t j = 0; j < L2.size(); ++j) {
        std::copy(L2[j].begin(), L2[j].end(), l.begin() + k2 + (k - k2) / 2);
        L2_t[j] = _t + (l * _C.get_G());
        l.assign(l.size(), false);
    }
    L_t.push_back(L2_t);

    // Perform merge-join with on-the-fly processing
    e.set();
    std::vector<size_t> positions = get_positions(_pre, _C, _c, k2);

    merge_lists(L_t, 1, positions, [&](binop::binvec& candidate) {

        // Size-reduce e
        if (_pre) {
            _C.size_red(candidate, k2);
        }

        // Write statistics
        stats[candidate.count()].first++;

        // Check if a valid error vector is found
        if (candidate.count() < e.count()) e = candidate;

        });
}

std::vector<long double> Stern::errors_dist_predict(std::vector<size_t> profile) {

    // Determine the code dimension
    size_t k = _C.get_k(), k2 = 0;

    // If the final vectors are preprocessed and reduced,
    // obtain corresponding k1, k2 and c
    if (_pre) {
        size_t c_new = get_c();
        k2 = _C.get_k2(c_new);
    }

    // Compute the weight distribution of a Hamming ball in that space
    std::vector<long double> dist1;
    if (!_pre) {
        HSpace space1(binop::n - k - _c);
        dist1 = space1.weight_dist_ball(binop::n - k - _c);
    }
    else {
        dist1 = weight_dist_fundamental_domain(profile, k2);
    }

    // Create an HSpace object representing the Hamming space of dimension k - k2
    HSpace space2((k - k2) / 2);

    // Compute the weight distribution of a Hamming ball in that space
    std::vector<long double> dist2 = space2.weight_dist_ball(_w2);

    // Compute joined distribution
    std::vector<long double> dist = convolve(dist1, dist2);
    dist = convolve(dist, dist2);

    // Compute the normalization coefficient (total number of cosets: 2^(n-k))
    long double coeff = std::pow(2.0L, static_cast<long double>(binop::n - k));
    std::vector<long double> dist_d(dist.size());

    // Normalize each entry by dividing by the total number of cosets
    for (size_t i = 0; i < dist.size(); ++i) {
        dist_d[i] = static_cast<long double>(dist[i]) / coeff;
    }

    return dist_d;
}

void Wagner::find_min_error(binop::binvec& e, std::vector<std::size_t> pi,
    std::map<std::size_t, std::pair<double, long double>>& stats) {

    // Prepare the generator matrix in systematic form
    if (!systematize(pi)) {
        throw std::runtime_error("Failed to systematize the generator matrix.");
    }

    // Special case (reverts back to Prange)
    if (_w2 == 0) {
        e = _t;
        stats[e.count()].first = 1;

        return;
    }

    // Define parameters for enumerator and Hamming ball traversal
    size_t k = _C.get_k(), k1 = 0, k2 = 0;
    std::vector<bool> l(k, false);
    std::vector<std::vector<binop::binvec>> L_t;

    if (_pre) {

        //// LLL basis reduction
        //_C.LLL(0, _C.get_k());

        // Sort based on the length of epipodal vectors
        _C.epi_sort();

        // Update values of k1, k2 and c
        size_t c_new = get_c();
        k1 = _C.get_k1();
        k2 = _C.get_k2(c_new);
        set_c(c_new);

        // Put the matrix in a nice form
        _C.sort_cols(k2);
    }
    std::cout << "\nk1: " << k1 << ", k2: " << k2 << ", c_new: " << get_c() << ", m: " << _m << "\n";

    // Obtain the first 2^a - 1 list
    binop::util::BinaryEnumerator enumerator1((k - k2) / _m, _w2);
    std::vector<std::vector<bool>> L1 = enumerator1.enumerate();

    for (size_t i = 0; i < _m - 1; ++i) {
        std::vector<binop::binvec> L1_t(L1.size());
        for (size_t j = 0; j < L1.size(); ++j) {
            std::copy(L1[j].begin(), L1[j].end(), l.begin() + k2 + i * (k - k2) / _m);
            L1_t[j] = l * _C.get_G();
            l.assign(l.size(), false);
        }
        L_t.push_back(L1_t);
    }

    // Obtain the second list
    binop::util::BinaryEnumerator enumerator2(k - k2 - (_m - 1) * (k - k2) / _m, _w2);
    std::vector<std::vector<bool>> L2 = enumerator2.enumerate();

    std::vector<binop::binvec> L2_t(L2.size());
    for (size_t j = 0; j < L2.size(); ++j) {
        std::copy(L2[j].begin(), L2[j].end(), l.begin() + k2 + (_m - 1) * (k - k2) / _m);
        L2_t[j] = _t + (l * _C.get_G());
        l.assign(l.size(), false);
    }
    L_t.push_back(L2_t);

    // Perform merge-join with on-the-fly processing
    e.set();
    std::vector<size_t> positions = get_positions(_pre, _C, _c, k2);

    merge_lists(L_t, _a, positions, [&](binop::binvec& candidate) {

        // Size-reduce e
        if (_pre) {
            _C.size_red(candidate, k2);
        }

        // Write statistics
        stats[candidate.count()].first++;

        // Check if a valid error vector is found
        if (candidate.count() < e.count()) e = candidate;

        });
}

std::vector<long double> Wagner::errors_dist_predict(std::vector<size_t> profile) {

    // Determine the code dimension
    size_t k = _C.get_k(), k2 = 0;

    // If the final vectors are preprocessed and reduced,
    // obtain corresponding k1, k2 and c
    if (_pre) {
        size_t c_new = get_c();
        k2 = _C.get_k2(c_new);
    }

    // Compute the weight distribution of a Hamming ball in that space
    std::vector<long double> dist1;
    if (!_pre) {
        HSpace space1(binop::n - k - _c);
        dist1 = space1.weight_dist_ball(binop::n - k - _c);
    }
    else {
        dist1 = weight_dist_fundamental_domain(profile, k2);
    }

    // Create an HSpace object representing the Hamming space of dimension k - k2
    HSpace space2((k - k2) / _m);

    // Compute the weight distribution of a Hamming ball in that space
    std::vector<long double> dist2 = space2.weight_dist_ball(_w2);

    // Compute joined distribution
    std::vector<long double> dist = convolve(dist1, dist2);
    for (size_t i = 1; i < _m; ++i) {
        dist = convolve(dist, dist2);
    }
    
    // Compute the normalization coefficient (total number of cosets: 2^(n-k))
    long double coeff = std::pow(2.0L, static_cast<long double>(binop::n - k));
    std::vector<long double> dist_d(dist.size());

    // Normalize each entry by dividing by the total number of cosets
    for (size_t i = 0; i < dist.size(); ++i) {
        dist_d[i] = static_cast<long double>(dist[i]) / coeff;
    }

    return dist_d;
}
