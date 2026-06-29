/**
 * \file DualNumber.h
 * \brief Single-variable dual number for first-order automatic differentiation.
 *
 * \details A dual number \f$ (v, d) \f$ represents a value \f$ v \f$ and its
 * first derivative \f$ dv/dx \f$ with respect to a single variable.
 * All common arithmetic operators and math functions are overloaded
 * to propagate derivatives through the chain rule.
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   - Default constructor now represents constant 0 (was: variable x)
 *   - Removed trivial copy ctor/dtor (compiler-generated is optimal)
 *   - Removed const from return-by-value (enables move semantics)
 *   - Added operator<<, operator==, operator!=
 *   - Operators are friend functions with direct member access
 *   - Fixed pow name shadowing (std::pow qualified)
 *   - Fixed atan2 dual computation for numerical stability
 *   - Added noexcept qualifiers on trivial operations
 *   - Removed dead commented-out code
 *   - Qualified std::string explicitly
 *   - Renamed class to DualNumber with snake_case get/set accessors
 */

#ifndef DualNumber_H
#define DualNumber_H

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "ValueErrors.h"

/// \class DualNumber
/// \brief Single-variable dual number for first-order automatic differentiation.
///
/// Stores a value and its derivative \f$ dv/dx \f$. All arithmetic and
/// elementary math functions are overloaded via friend functions to
/// correctly propagate derivatives.
class DualNumber
{
private:
    double value_;  ///< Function value
    double deriv_;  ///< First derivative \f$ dv/dx \f$

public:
    // --- constructors ---

    /// \brief Default constructor — represents constant 0.
    DualNumber() noexcept : value_(0.0), deriv_(0.0) {}

    /// \brief Construct with value and optional derivative.
    /// \param value the function value
    /// \param deriv the derivative (default 1.0 for variable x)
    DualNumber(double value, double deriv = 1.0) noexcept : value_(value), deriv_(deriv) {}

    // --- accessors ---

    /// \brief Get the value.
    double get_value()      const noexcept { return value_; }
    /// \brief Get the derivative.
    double get_derivative() const noexcept { return deriv_; }
    /// \brief Set the value.
    void set_value(double v)      noexcept { value_ = v; }
    /// \brief Set the derivative.
    void set_derivative(double d) noexcept { deriv_ = d; }

    // --- I/O ---

    /// \brief Print to stdout.
    void print() const { std::cout << *this << std::endl; }

    /// \brief Stream output as "(value, deriv)".
    friend std::ostream& operator<<(std::ostream& os, const DualNumber& dn) {
        return os << "(" << dn.value_ << ", " << dn.deriv_ << ")";
    }

    /// \brief Convert to string.
    std::string to_string() const {
        std::ostringstream oss;
        oss << *this;
        return oss.str();
    }

    // --- comparison ---

    /// \brief Equality comparison (exact floating-point match).
    friend bool operator==(const DualNumber& a, const DualNumber& b) noexcept {
        return a.value_ == b.value_ && a.deriv_ == b.deriv_;
    }
    /// \brief Inequality comparison.
    friend bool operator!=(const DualNumber& a, const DualNumber& b) noexcept {
        return !(a == b);
    }

    // --- LaTeX / ROOT formatted output ---

    /// \brief Static helper: format value and error via ValueErrors.
    /// \param val the value
    /// \param err the error
    /// \param opt format option ("R" for ROOT, "L" for LaTeX)
    static std::string str_latex(double val, double err, const std::string& opt = "R") {
        return ValueErrors(val, err).print_err1(opt);
    }

    /// \brief Format this dual number via ValueErrors.
    /// \param opt format option ("R" for ROOT, "L" for LaTeX)
    std::string str_latex(const std::string& opt = "R") const {
        return str_latex(value_, deriv_, opt);
    }

    // --- unary operators ---

    /// \brief Unary plus (identity).
    friend DualNumber operator+(const DualNumber& dn) {
        return dn;
    }
    /// \brief Unary minus (negate both value and derivative).
    friend DualNumber operator-(const DualNumber& dn) {
        return DualNumber(-dn.value_, -dn.deriv_);
    }

    // --- dual-dual arithmetic ---

    /// \brief Addition: \f$ (v_a+v_b,\ d_a+d_b) \f$.
    friend DualNumber operator+(const DualNumber& a, const DualNumber& b) {
        return DualNumber(a.value_ + b.value_, a.deriv_ + b.deriv_);
    }
    /// \brief Subtraction: \f$ (v_a-v_b,\ d_a-d_b) \f$.
    friend DualNumber operator-(const DualNumber& a, const DualNumber& b) {
        return DualNumber(a.value_ - b.value_, a.deriv_ - b.deriv_);
    }
    /// \brief Multiplication (product rule): \f$ (v_a v_b,\ v_a d_b + v_b d_a) \f$.
    friend DualNumber operator*(const DualNumber& a, const DualNumber& b) {
        return DualNumber(a.value_ * b.value_,
                          a.value_ * b.deriv_ + b.value_ * a.deriv_);
    }
    /// \brief Division (quotient rule).
    friend DualNumber operator/(const DualNumber& a, const DualNumber& b) {
        double inv = 1.0 / b.value_;
        return DualNumber(a.value_ * inv,
                          (a.deriv_ * b.value_ - a.value_ * b.deriv_) / (b.value_ * b.value_));
    }

    // --- dual-scalar arithmetic ---

    /// \brief Dual + scalar: \f$ (v + c,\ d) \f$.
    friend DualNumber operator+(const DualNumber& dn, double c) {
        return DualNumber(dn.value_ + c, dn.deriv_);
    }
    /// \brief Scalar + dual.
    friend DualNumber operator+(double c, const DualNumber& dn) {
        return dn + c;
    }
    /// \brief Dual - scalar.
    friend DualNumber operator-(const DualNumber& dn, double c) {
        return DualNumber(dn.value_ - c, dn.deriv_);
    }
    /// \brief Scalar - dual.
    friend DualNumber operator-(double c, const DualNumber& dn) {
        return DualNumber(c - dn.value_, -dn.deriv_);
    }
    /// \brief Dual * scalar: \f$ (v c,\ d c) \f$.
    friend DualNumber operator*(const DualNumber& dn, double c) {
        return DualNumber(dn.value_ * c, dn.deriv_ * c);
    }
    /// \brief Scalar * dual.
    friend DualNumber operator*(double c, const DualNumber& dn) {
        return dn * c;
    }
    /// \brief Dual / scalar.
    friend DualNumber operator/(const DualNumber& dn, double c) {
        double inv = 1.0 / c;
        return DualNumber(dn.value_ * inv, dn.deriv_ * inv);
    }
    /// \brief Scalar / dual: \f$ (c/v,\ -c d / v^2) \f$.
    friend DualNumber operator/(double c, const DualNumber& dn) {
        double inv = 1.0 / dn.value_;
        return DualNumber(c * inv, -c * dn.deriv_ * inv * inv);
    }

    // --- math functions ---

    /// \brief Power: \f$ (v^c,\ c\, v^{c-1}\, d) \f$.
    friend DualNumber pow(const DualNumber& dn, double c) {
        double vp = std::pow(dn.value_, c);
        return DualNumber(vp, c * std::pow(dn.value_, c - 1.0) * dn.deriv_);
    }

    /// \brief Square root (delegates to pow with c = 0.5).
    friend DualNumber sqrt(const DualNumber& dn) {
        return pow(dn, 0.5);
    }

    /// \brief Absolute value: negates if value < 0.
    friend DualNumber abs(const DualNumber& dn) {
        return (dn.value_ < 0.0) ? -dn : dn;
    }

    /// \brief Exponential: \f$ (e^v,\ e^v\, d) \f$.
    friend DualNumber exp(const DualNumber& dn) {
        double ev = std::exp(dn.value_);
        return DualNumber(ev, ev * dn.deriv_);
    }

    /// \brief Natural logarithm: \f$ (\ln v,\ d/v) \f$.
    friend DualNumber log(const DualNumber& dn) {
        return DualNumber(std::log(dn.value_), dn.deriv_ / dn.value_);
    }

    /// \brief Sine: \f$ (\sin v,\ \cos v \cdot d) \f$.
    friend DualNumber sin(const DualNumber& dn) {
        return DualNumber(std::sin(dn.value_), std::cos(dn.value_) * dn.deriv_);
    }

    /// \brief Cosine: \f$ (\cos v,\ -\sin v \cdot d) \f$.
    friend DualNumber cos(const DualNumber& dn) {
        return DualNumber(std::cos(dn.value_), -std::sin(dn.value_) * dn.deriv_);
    }

    /// \brief Tangent: \f$ (\tan v,\ d / \cos^2 v) \f$.
    friend DualNumber tan(const DualNumber& dn) {
        double cv = std::cos(dn.value_);
        return DualNumber(std::tan(dn.value_), dn.deriv_ / (cv * cv));
    }

    /// \brief Arctangent: \f$ (\arctan v,\ d / (1 + v^2)) \f$.
    friend DualNumber atan(const DualNumber& dn) {
        return DualNumber(std::atan(dn.value_),
                          dn.deriv_ / (1.0 + dn.value_ * dn.value_));
    }

    /// \brief Two-argument arctangent: \f$ \atan2(y, x) \f$.
    /// Derivative: \f$ (x\,dy - y\,dx) / (x^2 + y^2) \f$.
    friend DualNumber atan2(const DualNumber& y, const DualNumber& x) {
        double denom = x.value_ * x.value_ + y.value_ * y.value_;
        return DualNumber(std::atan2(y.value_, x.value_),
                          (x.value_ * y.deriv_ - y.value_ * x.deriv_) / denom);
    }
};


#endif
