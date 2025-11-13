#ifndef DUPLICATE_VEC_HPP
#define DUPLICATE_VEC_HPP

#include <RcppArmadillo.h>
#include <unordered_map>
#include <cmath>

// ---- Hash for arma::vec ----
struct DenseVecHash {
    std::size_t operator()(const arma::vec& v) const;
};

// ---- Equality comparator with tolerance ----
struct DenseVecEqual {
    bool operator()(const arma::vec& a, const arma::vec& b) const;
};

// ---- Deduplicator class ----
class DenseVecDeduplicator {
public:
    DenseVecDeduplicator() = default;

    // Check vector and add if unique
    bool check_and_add(const arma::vec& v, int idx);

    // Number of unique vectors seen
    int size() const;

    // Get index of first occurrence of a vector
    int get_first_index(const arma::vec& v) const;

private:
    std::unordered_map<arma::vec, int, DenseVecHash, DenseVecEqual> seen_map;
};

#endif // DUPLICATE_VEC_HPP

