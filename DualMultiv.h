/**
 * \file DualMultiv.h
 * \brief Multi-variable dual number for first-order automatic differentiation.
 *
 * \details Tracks partial derivatives per variable ID via a sorted
 * std::map<int, double>. Supports all common arithmetic and math
 * functions, plus weighted-average merging and correlated-systematic
 * uncertainty computation via a static correlation matrix.
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   - Renamed MyDualMultiv -> DualMultiv; snake_case with get_/set_ prefix
 *   - Removed const from operator return-by-value (enables move semantics)
 *   - Direct computation in operators — no intermediate DualNumber temporaries
 *   - std::string fully qualified; fixed param name (symboltype -> opt)
 *   - operator/(double, dm) computes derivative analytically instead of
 *     constructing DualNumber objects per variable
 *   - Added noexcept qualifiers on trivial operations
 *   - Added operator<<, to_string(); replaced Print() and StrLatex()
 *   - Added reset_id_counter(), get_current_id()
 *   - Added static correlation matrix for correlated systematics
 *     (get_uncertainty() includes cross-terms: sigma^2 = Sum d_i^2 + 2*Sum_{i<j} d_i d_j rho_ij)
 *   - Added missing headers (climits, sstream)
 *   - Iterators replace C-style (int)list.size() casts
 */

#ifndef DualMultiv_H
#define DualMultiv_H

#include <climits>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "DualNumber.h"

/// \class DualMultiv
/// \brief Multi-variable dual number with automatic differentiation.
///
/// Stores a value and a map of partial derivatives keyed by variable ID.
/// Operators implement the multi-variable chain rule. A static correlation
/// matrix enables correlated-uncertainty propagation via get_uncertainty().
class DualMultiv
{
private:
    double value_;                         ///< Function value
    std::map<int, double> duals_;         ///< Map: variable ID -> partial derivative
    inline static int id_counter_ = 0;    ///< Auto-incrementing variable ID counter
    inline static std::map<std::pair<int,int>, double> correlations_; ///< Sparse correlation matrix

public:
    static const int ID_MAX = INT_MAX;     ///< Maximum allowed variable ID

    // --- auto-incrementing variable ID ---

    /// \brief Allocate a new unique variable ID.
    /// \return the new ID (1-based).
    static int new_id() {
        if (id_counter_ >= ID_MAX - 1) {
            std::cerr << "DualMultiv::new_id: counter exhausted at "
                      << id_counter_ << " / " << (ID_MAX - 1) << std::endl;
        } else {
            ++id_counter_;
        }
        return id_counter_;
    }

    /// \brief Reset the ID counter to 0.
    static void reset_id_counter() noexcept { id_counter_ = 0; }

    /// \brief Advance the ID counter to at least \p id.
    /// \param id the minimum counter value (useful when reading objects back
    ///           from file to prevent ID collisions).
    static void set_id_counter(int id) noexcept {
        if (id > id_counter_) id_counter_ = id;
    }

    /// \brief Get the current ID counter value.
    static int get_current_id() noexcept { return id_counter_; }

    // --- correlation matrix (global, static) ---

    /// \brief Set the correlation coefficient \f$ \rho_{ij} \f$ between
    ///        systematic sources i and j.
    ///
    /// The key is always normalized so that i < j.
    /// \param i first variable ID
    /// \param j second variable ID
    /// \param rho correlation coefficient in [-1, 1]
    static void set_correlation(int i, int j, double rho) {
        if (i > j) std::swap(i, j);
        correlations_[{i, j}] = rho;
    }

    /// \brief Get the correlation coefficient between sources i and j.
    /// \return \f$ \rho_{ij} \f$, or 0 if not set.
    static double get_correlation(int i, int j) noexcept {
        if (i > j) std::swap(i, j);
        auto it = correlations_.find({i, j});
        return (it != correlations_.end()) ? it->second : 0.0;
    }

    /// \brief Clear all stored correlations.
    static void clear_correlations() noexcept { correlations_.clear(); }

    // --- constructors ---

    /// \brief Default constructor — value 0, no partial derivatives.
    DualMultiv() noexcept : value_(0.0) {}

    /// \brief Construct a constant (no derivatives).
    explicit DualMultiv(double value) noexcept : value_(value) {}

    /// \brief Construct with a value and derivative, auto-assigning a new ID.
    DualMultiv(double value, double deriv) : value_(value) {
        duals_.insert({new_id(), deriv});
    }

    /// \brief Construct from an existing ID and a DualNumber.
    DualMultiv(int id, const DualNumber& dn) : value_(dn.get_value()) {
        duals_.insert({id, dn.get_derivative()});
    }

    /// \brief Construct from a DualNumber, auto-assigning a new ID.
    explicit DualMultiv(const DualNumber& dn) : value_(dn.get_value()) {
        duals_.insert({new_id(), dn.get_derivative()});
    }

    // --- accessors ---

    /// \brief Get the value.
    double get_value() const noexcept { return value_; }
    /// \brief Set the value.
    void set_value(double v) noexcept { value_ = v; }

    /// \brief Get all variable IDs with non-zero partial derivatives.
    std::vector<int> get_ids() const {
        std::vector<int> ids;
        ids.reserve(duals_.size());
        for (const auto& p : duals_) ids.push_back(p.first);
        return ids;
    }

    /// \brief Get the partial derivative for variable \p id.
    /// \return the derivative, or 0 if not present.
    double get_derivative(int id) const {
        auto it = duals_.find(id);
        return (it != duals_.end()) ? it->second : 0.0;
    }

    /// \brief Set the partial derivative for variable \p id.
    void set_derivative(int id, double deriv) {
        duals_[id] = deriv;
    }

    /// \brief Set a partial derivative with auto-assigned ID.
    void set_derivative(double deriv) {
        duals_[new_id()] = deriv;
    }

    /// \brief Set partial derivative from a DualNumber's derivative.
    void set_derivative(int id, const DualNumber& dn) {
        duals_[id] = dn.get_derivative();
    }

    /// \brief Set partial derivative from a DualNumber, auto-assigning ID.
    void set_derivative(const DualNumber& dn) {
        duals_[new_id()] = dn.get_derivative();
    }

    /// \brief Extract a single-variable DualNumber for variable \p id.
    DualNumber get_part(int id) const {
        return DualNumber(value_, get_derivative(id));
    }

    /// \brief Compute the total uncertainty via quadrature sum.
    ///
    /// \f[ \sigma = \sqrt{ \sum_i d_i^2 + 2\sum_{i<j} d_i d_j \rho_{ij} } \f]
    /// where \f$ d_i \f$ are partial derivatives and \f$ \rho_{ij} \f$
    /// are correlation coefficients from the static correlation matrix.
    double get_uncertainty() const {
        double sum2 = 0.0;
        for (const auto& p : duals_) sum2 += p.second * p.second;
        if (correlations_.empty()) return std::sqrt(sum2);
        for (const auto& c : correlations_) {
            int i = c.first.first, j = c.first.second;
            double rho = c.second;
            auto it_i = duals_.find(i);
            auto it_j = duals_.find(j);
            if (it_i != duals_.end() && it_j != duals_.end()) {
                sum2 += 2.0 * it_i->second * it_j->second * rho;
            }
        }
        return std::sqrt(sum2);
    }

    // --- I/O ---

    /// \brief Convert to string.
    std::string to_string() const {
        std::ostringstream oss;
        oss << *this;
        return oss.str();
    }

    /// \brief LaTeX / ROOT formatted output via ValueErrors.
    std::string str_latex(const std::string& opt = "R") const {
        return ValueErrors(value_, get_uncertainty()).print_err1(opt);
    }

    /// \brief Stream output: "(value, [id:deriv]...)".
    friend std::ostream& operator<<(std::ostream& os, const DualMultiv& dm) {
        os << "(" << dm.value_;
        for (const auto& p : dm.duals_) {
            os << ", [" << p.first << ":" << p.second << "]";
        }
        os << ")";
        return os;
    }

    // --- unary operators ---

    /// \brief Unary plus (identity).
    friend DualMultiv operator+(const DualMultiv& dm) { return dm; }

    /// \brief Unary minus (negates value and all partials).
    friend DualMultiv operator-(const DualMultiv& dm) {
        DualMultiv result(-dm.value_);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = -p.second;
        }
        return result;
    }

    // --- dual-dual arithmetic ---

    /// \brief Addition: merges partial derivative maps.
    friend DualMultiv operator+(const DualMultiv& a, const DualMultiv& b) {
        DualMultiv result(a.value_ + b.value_);
        auto ia = a.duals_.begin(), ib = b.duals_.begin();
        while (ia != a.duals_.end() && ib != b.duals_.end()) {
            if (ia->first < ib->first) {
                result.duals_[ia->first] = ia->second;
                ++ia;
            } else if (ia->first > ib->first) {
                result.duals_[ib->first] = ib->second;
                ++ib;
            } else {
                result.duals_[ia->first] = ia->second + ib->second;
                ++ia; ++ib;
            }
        }
        for (; ia != a.duals_.end(); ++ia) result.duals_[ia->first] = ia->second;
        for (; ib != b.duals_.end(); ++ib) result.duals_[ib->first] = ib->second;
        return result;
    }

    /// \brief Subtraction.
    friend DualMultiv operator-(const DualMultiv& a, const DualMultiv& b) {
        return a + (-b);
    }

    /// \brief Multiplication (product rule for multi-variable case).
    friend DualMultiv operator*(const DualMultiv& a, const DualMultiv& b) {
        DualMultiv result(a.value_ * b.value_);
        auto ia = a.duals_.begin(), ib = b.duals_.begin();
        while (ia != a.duals_.end() && ib != b.duals_.end()) {
            if (ia->first < ib->first) {
                result.duals_[ia->first] = ia->second * b.value_;
                ++ia;
            } else if (ia->first > ib->first) {
                result.duals_[ib->first] = a.value_ * ib->second;
                ++ib;
            } else {
                result.duals_[ia->first] = ia->second * b.value_ + a.value_ * ib->second;
                ++ia; ++ib;
            }
        }
        for (; ia != a.duals_.end(); ++ia) result.duals_[ia->first] = ia->second * b.value_;
        for (; ib != b.duals_.end(); ++ib) result.duals_[ib->first] = a.value_ * ib->second;
        return result;
    }

    /// \brief Division (delegates to multiplication by reciprocal).
    friend DualMultiv operator/(const DualMultiv& a, const DualMultiv& b) {
        return a * (1.0 / b);
    }

    // --- dual-scalar arithmetic ---

    /// \brief DualMultiv + scalar.
    friend DualMultiv operator+(const DualMultiv& dm, double c) {
        DualMultiv result(dm.value_ + c);
        result.duals_ = dm.duals_;
        return result;
    }
    /// \brief Scalar + DualMultiv.
    friend DualMultiv operator+(double c, const DualMultiv& dm) { return dm + c; }

    /// \brief DualMultiv - scalar.
    friend DualMultiv operator-(const DualMultiv& dm, double c) {
        DualMultiv result(dm.value_ - c);
        result.duals_ = dm.duals_;
        return result;
    }
    /// \brief Scalar - DualMultiv.
    friend DualMultiv operator-(double c, const DualMultiv& dm) {
        DualMultiv result(c - dm.value_);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = -p.second;
        }
        return result;
    }

    /// \brief DualMultiv * scalar (scales all partials).
    friend DualMultiv operator*(const DualMultiv& dm, double c) {
        DualMultiv result(dm.value_ * c);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = p.second * c;
        }
        return result;
    }
    /// \brief Scalar * DualMultiv.
    friend DualMultiv operator*(double c, const DualMultiv& dm) { return dm * c; }

    /// \brief DualMultiv / scalar.
    friend DualMultiv operator/(const DualMultiv& dm, double c) {
        double inv = 1.0 / c;
        return dm * inv;
    }

    /// \brief Scalar / DualMultiv.
    friend DualMultiv operator/(double c, const DualMultiv& dm) {
        DualMultiv result(c / dm.value_);
        double coeff = -c / (dm.value_ * dm.value_);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = coeff * p.second;
        }
        return result;
    }

    // --- math functions ---

    /// \brief Power: \f$ (v^c,\ c v^{c-1} \cdot \partial_i v) \f$ for each variable.
    friend DualMultiv pow(const DualMultiv& dm, double c) {
        DualMultiv result(std::pow(dm.value_, c));
        double coeff = c * std::pow(dm.value_, c - 1.0);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = coeff * p.second;
        }
        return result;
    }

    /// \brief Square root (delegates to pow with c = 0.5).
    friend DualMultiv sqrt(const DualMultiv& dm) { return pow(dm, 0.5); }

    /// \brief Absolute value.
    friend DualMultiv abs(const DualMultiv& dm) {
        return (dm.value_ < 0.0) ? -dm : dm;
    }

    /// \brief Exponential.
    friend DualMultiv exp(const DualMultiv& dm) {
        double ev = std::exp(dm.value_);
        DualMultiv result(ev);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = ev * p.second;
        }
        return result;
    }

    /// \brief Natural logarithm.
    friend DualMultiv log(const DualMultiv& dm) {
        DualMultiv result(std::log(dm.value_));
        double inv = 1.0 / dm.value_;
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = p.second * inv;
        }
        return result;
    }

    /// \brief Sine.
    friend DualMultiv sin(const DualMultiv& dm) {
        DualMultiv result(std::sin(dm.value_));
        double cv = std::cos(dm.value_);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = cv * p.second;
        }
        return result;
    }

    /// \brief Cosine.
    friend DualMultiv cos(const DualMultiv& dm) {
        DualMultiv result(std::cos(dm.value_));
        double sv = -std::sin(dm.value_);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = sv * p.second;
        }
        return result;
    }

    /// \brief Tangent.
    friend DualMultiv tan(const DualMultiv& dm) {
        DualMultiv result(std::tan(dm.value_));
        double cv = std::cos(dm.value_);
        double coeff = 1.0 / (cv * cv);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = coeff * p.second;
        }
        return result;
    }

    /// \brief Arctangent.
    friend DualMultiv atan(const DualMultiv& dm) {
        DualMultiv result(std::atan(dm.value_));
        double coeff = 1.0 / (1.0 + dm.value_ * dm.value_);
        for (const auto& p : dm.duals_) {
            result.duals_[p.first] = coeff * p.second;
        }
        return result;
    }

    /// \brief Two-argument arctangent for multi-variable duals.
    friend DualMultiv atan2(const DualMultiv& y, const DualMultiv& x) {
        double denom = x.value_ * x.value_ + y.value_ * y.value_;
        DualMultiv result(std::atan2(y.value_, x.value_));
        auto iy = y.duals_.begin(), ix = x.duals_.begin();
        while (iy != y.duals_.end() && ix != x.duals_.end()) {
            if (iy->first < ix->first) {
                result.duals_[iy->first] = (x.value_ * iy->second) / denom;
                ++iy;
            } else if (iy->first > ix->first) {
                result.duals_[ix->first] = (-y.value_ * ix->second) / denom;
                ++ix;
            } else {
                result.duals_[iy->first] = (x.value_ * iy->second - y.value_ * ix->second) / denom;
                ++iy; ++ix;
            }
        }
        for (; iy != y.duals_.end(); ++iy)
            result.duals_[iy->first] = (x.value_ * iy->second) / denom;
        for (; ix != x.duals_.end(); ++ix)
            result.duals_[ix->first] = (-y.value_ * ix->second) / denom;
        return result;
    }
};


#endif
