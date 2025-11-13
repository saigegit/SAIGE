#include "Duplicate_vec.hpp"
#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_map>

// ---- Hash for arma::vec ----
std::size_t DenseVecHash::operator()(const arma::vec& v) const {
    std::size_t h = 0;
    for (arma::uword i = 0; i < v.n_elem; ++i) {
        std::size_t value_hash = std::hash<double>{}(v[i]);
        h ^= value_hash + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
}

// ---- Equality comparator with tolerance ----
bool DenseVecEqual::operator()(const arma::vec& a, const arma::vec& b) const {
    if (a.n_elem != b.n_elem) return false;
    for (arma::uword i = 0; i < a.n_elem; ++i) {
        if (std::fabs(a[i] - b[i]) > 1e-12) return false;
    }
    return true;
}

// ---- DenseVecDeduplicator class methods ----

// Check vector and add if unique
bool DenseVecDeduplicator::check_and_add(const arma::vec& v, int idx) {
    auto it = seen_map.find(v);
    if (it == seen_map.end()) {
        seen_map[v] = idx; // store first occurrence index
        return true;       // unique
    } else {
        return false;      // duplicate
    }
}

// Number of unique vectors seen
int DenseVecDeduplicator::size() const {
    return static_cast<int>(seen_map.size());
}

// Get index of first occurrence of a vector
int DenseVecDeduplicator::get_first_index(const arma::vec& v) const {
    auto it = seen_map.find(v);
    if (it != seen_map.end()) return it->second;
    return -1; // not found
}
