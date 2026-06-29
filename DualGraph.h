/**
 * \file DualGraph.h
 * \brief Core graph classes for dual-number uncertainty propagation (ROOT-independent).
 *
 * \details Defines DualPoint (a single (x, y) point with dual-number uncertainties)
 * and DualGraph (a collection of DualPoints with arithmetic operators, merge
 * operations, and math functions). All classes in this file are ROOT-independent.
 * For ROOT I/O and TGraphErrors mirroring, use TDualGraph instead.
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   - Renamed MyDualPoint -> DualPoint, MyDualGraph -> DualGraph
 *   - Data members made private; snake_case get_/set_ prefix throughout
 *   - Dropped "My" prefix from includes (DualNumber.h, DualMultiv.h)
 *   - Inlined inverse-variance weighting (previously AvePlus/average_plus)
 *   - Removed const from operator return-by-value (enables move semantics)
 *   - Replaced VLAs (double vx[n]) with std::vector<double>
 *   - std::string qualified everywhere; noexcept on trivial methods
 *   - Added operator<< for DualPoint
 *   - Split into DualGraph (base, ROOT-free) + TDualGraph (ROOT I/O)
 *   - points_ made protected; mutators virtual for derived-class override
 */

#ifndef DualGraph_H
#define DualGraph_H

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "DualNumber.h"
#include "DualMultiv.h"

// ===========================================================================
//  DualPoint — a single (x, y) point with dual-number uncertainties
// ===========================================================================

/// \class DualPoint
/// \brief A single (x, y) data point with dual-number uncertainties.
///
/// The x-coordinate is a single-variable DualNumber; the y-coordinate is a
/// multi-variable DualMultiv supporting multiple systematic sources.
class DualPoint
{
private:
    DualNumber   px_;  ///< x-coordinate with its uncertainty
    DualMultiv   py_;  ///< y-coordinate with multi-variable uncertainties

public:
    // --- constructors ---

    /// \brief Default constructor.
    DualPoint() noexcept = default;

    /// \brief Construct from x only (y defaults to 0).
    explicit DualPoint(const DualNumber& px) noexcept : px_(px) {}

    /// \brief Construct from x and y dual numbers.
    DualPoint(const DualNumber& px, const DualMultiv& py) noexcept
        : px_(px), py_(py) {}

    // --- accessors ---

    /// \brief Get the x-coordinate (const).
    const DualNumber&   get_px() const noexcept { return px_; }
    /// \brief Get the y-coordinate (const).
    const DualMultiv&   get_py() const noexcept { return py_; }
    /// \brief Get a mutable reference to the y-coordinate.
    DualMultiv&         get_py()       noexcept { return py_; }

    /// \brief Set the x-coordinate.
    void set_px(const DualNumber& px) noexcept { px_ = px; }
    /// \brief Set the y-coordinate.
    void set_py(const DualMultiv& py) noexcept { py_ = py; }

    // --- I/O ---

    /// \brief Stream output.
    friend std::ostream& operator<<(std::ostream& os, const DualPoint& dp) {
        return os << "{px=" << dp.px_ << ", py=" << dp.py_ << "}";
    }
};

// ===========================================================================
//  DualGraph — graph of DualPoint objects (ROOT-independent base class)
// ===========================================================================

/// \class DualGraph
/// \brief A graph of DualPoint objects supporting uncertainty propagation.
///
/// Provides arithmetic operators (+ - * /), math functions (pow, sqrt, exp,
/// log, sin, cos, tan, atan, atan2), weighted/unweighted bin averaging,
/// and bin selection. This is the ROOT-independent base class; for ROOT I/O
/// use TDualGraph.
class DualGraph
{
protected:
    std::vector<DualPoint> points_;  ///< The data points

    /// \brief Check that \p n matches the number of points.
    /// Prints a warning on mismatch.
    bool check_size(int n) const noexcept {
        if (n != static_cast<int>(points_.size())) {
            std::cerr << "DualGraph: warning: size mismatch! "
                      << n << " vs " << points_.size() << std::endl;
            return false;
        }
        return true;
    }

public:
    // --- constructors ---

    /// \brief Default constructor — empty graph.
    DualGraph() = default;

    /// \brief Construct with \p n default-initialized points.
    explicit DualGraph(int n) {
        reset(n);
    }

    /// \brief Construct from an existing vector of DualPoints.
    explicit DualGraph(const std::vector<DualPoint>& points)
        : points_(points) {}

    /// \brief Virtual destructor (enables derived-class polymorphism).
    virtual ~DualGraph() = default;

    // --- accessors ---

    /// \brief Get the number of points.
    int get_n() const noexcept { return static_cast<int>(points_.size()); }

    /// \brief Check that index \p i is in range.
    /// Prints a warning if out of range.
    bool check_index(int i) const noexcept {
        int n = static_cast<int>(points_.size());
        if (i < 0 || i >= n) {
            std::cerr << "DualGraph: warning: index out of range! "
                      << i << "/" << n << std::endl;
            return false;
        }
        return true;
    }

    /// \brief Get point at index \p i (bounds-checked).
    /// \return a copy of the point, or a default DualPoint if out of range.
    DualPoint get_point(int i) const {
        if (!check_index(i)) return DualPoint();
        return points_[i];
    }

    // --- mutators ---

    /// \brief Replace all points (virtual — overridden by TDualGraph).
    virtual void set_points(const std::vector<DualPoint>& points) {
        points_ = points;
    }

    /// \brief Set a single point at index \p i (virtual).
    virtual void set_point(int i, const DualPoint& dp) {
        if (!check_index(i)) return;
        points_[i] = dp;
    }

    // --- reset ---

    /// \brief Reset to \p n default-initialized points (virtual).
    virtual void reset(int n = 0) {
        points_.clear();
        points_.reserve(n);
        for (int i = 0; i < n; ++i) {
            points_.emplace_back();
        }
    }

    // --- LaTeX / ROOT formatted output ---

    /// \brief LaTeX-formatted Y-value at index \p i.
    /// \param i point index
    /// \param opt format option ("R" for ROOT, "L" for LaTeX)
    std::string str_latex(int i = 0,
                          const std::string& opt = "R") const {
        if (!check_index(i)) return std::string();
        return points_[i].get_py().str_latex(opt);
    }

    // --- merge operations ---

    /// \brief Weighted-average merge of specified bins (inverse-variance weights).
    ///
    /// Combines points at the given indices into a single bin using
    /// \f$ w_i = 1 / \sigma_i^2 \f$ weighting.
    /// \param range list of bin indices to merge
    /// \return a single-bin DualGraph with the weighted average
    DualGraph average_bin(const std::vector<int>& range) const {
        std::vector<int> unique_range = range;
        std::sort(unique_range.begin(), unique_range.end());
        unique_range.erase(
            std::unique(unique_range.begin(), unique_range.end()),
            unique_range.end());
        int n = static_cast<int>(unique_range.size());
        if (n <= 0) {
            std::cerr << "DualGraph::average_bin(): invalid range!" << std::endl;
            return DualGraph();
        }
        if (!check_index(unique_range[0]) || !check_index(unique_range[n - 1])) {
            std::cerr << "DualGraph::average_bin(): invalid range!" << std::endl;
            return DualGraph();
        }
        DualNumber   avex(0.0, 0.0);
        DualMultiv   avey;
        for (int k = 0; k < n; ++k) {
            int i = unique_range[k];
            avex = avex + points_[i].get_px();
            DualMultiv dm;
            dm.set_value(points_[i].get_py().get_value());
            std::vector<int> ids = points_[i].get_py().get_ids();
            for (std::size_t j = 0; j < ids.size(); ++j) {
                int idx = ids[j];
                dm.set_derivative(idx, points_[i].get_py().get_derivative(idx));
            }
            if (k == 0) {
                avey = dm;
            } else {
                double e1 = avey.get_uncertainty();
                double e2 = dm.get_uncertainty();
                double w1 = 1.0 / (e1 * e1);
                double w2 = 1.0 / (e2 * e2);
                avey = (w1 * avey + w2 * dm) / (w1 + w2);
            }
        }
        avex = avex / static_cast<double>(n);
        std::vector<DualPoint> out;
        out.emplace_back(avex, avey);
        return DualGraph(out);
    }

    /// \brief Simple (equal-weight) average merge of specified bins.
    /// \param range list of bin indices to merge
    /// \return a single-bin DualGraph with the arithmetic mean
    DualGraph average_bin_simple(const std::vector<int>& range) const {
        std::vector<int> unique_range = range;
        std::sort(unique_range.begin(), unique_range.end());
        unique_range.erase(
            std::unique(unique_range.begin(), unique_range.end()),
            unique_range.end());
        int n = static_cast<int>(unique_range.size());
        if (n <= 0) {
            std::cerr << "DualGraph::average_bin_simple(): invalid range!" << std::endl;
            return DualGraph();
        }
        if (!check_index(unique_range[0]) || !check_index(unique_range[n - 1])) {
            std::cerr << "DualGraph::average_bin_simple(): invalid range!" << std::endl;
            return DualGraph();
        }
        DualNumber   avex(0.0, 0.0);
        DualMultiv   avey;
        double inv_n = 1.0 / n;
        for (int k = 0; k < n; ++k) {
            int i = unique_range[k];
            avex = avex + points_[i].get_px();
            DualMultiv dm;
            dm.set_value(points_[i].get_py().get_value());
            std::vector<int> ids = points_[i].get_py().get_ids();
            for (std::size_t j = 0; j < ids.size(); ++j) {
                int idx = ids[j];
                dm.set_derivative(idx, points_[i].get_py().get_derivative(idx));
            }
            avey = avey + dm * inv_n;
        }
        avex = avex * inv_n;
        std::vector<DualPoint> out;
        out.emplace_back(avex, avey);
        return DualGraph(out);
    }

    // --- select subset of bins ---

    /// \brief Select a subset of bins by index.
    /// \param range list of bin indices to keep
    /// \return a new DualGraph with only the selected bins
    DualGraph select_bin(const std::vector<int>& range) const {
        std::vector<int> unique_range = range;
        std::sort(unique_range.begin(), unique_range.end());
        unique_range.erase(
            std::unique(unique_range.begin(), unique_range.end()),
            unique_range.end());
        int n = static_cast<int>(unique_range.size());
        if (n <= 0) {
            std::cerr << "DualGraph::select_bin(): invalid range!" << std::endl;
            return DualGraph();
        }
        if (!check_index(unique_range[0]) || !check_index(unique_range[n - 1])) {
            std::cerr << "DualGraph::select_bin(): invalid range!" << std::endl;
            return DualGraph();
        }
        std::vector<DualPoint> out;
        out.reserve(n);
        for (int k = 0; k < n; ++k) {
            out.push_back(points_[unique_range[k]]);
        }
        return DualGraph(out);
    }

};

// ===========================================================================
//  Unary operators
// ===========================================================================

/// \brief Unary plus (identity).
inline DualGraph operator+(const DualGraph& dg) {
    return dg;
}

/// \brief Unary minus (negates all y-values).
inline DualGraph operator-(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), -dp.get_py()));
    }
    return result;
}

// ===========================================================================
//  Free operators — DualGraph vs scalar
// ===========================================================================

/// \brief DualGraph + scalar (adds to all y-values).
inline DualGraph operator+(const DualGraph& dg, double c) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), dp.get_py() + c));
    }
    return result;
}
inline DualGraph operator+(double c, const DualGraph& dg) { return dg + c; }

/// \brief DualGraph - scalar.
inline DualGraph operator-(const DualGraph& dg, double c) { return dg + (-c); }
/// \brief Scalar - DualGraph.
inline DualGraph operator-(double c, const DualGraph& dg) { return -dg + c; }

/// \brief DualGraph * scalar.
inline DualGraph operator*(const DualGraph& dg, double c) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), dp.get_py() * c));
    }
    return result;
}
inline DualGraph operator*(double c, const DualGraph& dg) { return dg * c; }

/// \brief DualGraph / scalar.
inline DualGraph operator/(const DualGraph& dg, double c) { return dg * (1.0 / c); }
/// \brief Scalar / DualGraph.
inline DualGraph operator/(double c, const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), c / dp.get_py()));
    }
    return result;
}

// ===========================================================================
//  Free operators — DualGraph vs DualMultiv
// ===========================================================================

/// \brief DualGraph + DualMultiv (broadcast to all points).
inline DualGraph operator+(const DualGraph& dg, const DualMultiv& c) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), dp.get_py() + c));
    }
    return result;
}
inline DualGraph operator+(const DualMultiv& c, const DualGraph& dg) { return dg + c; }
inline DualGraph operator-(const DualGraph& dg, const DualMultiv& c) { return dg + (-c); }
inline DualGraph operator-(const DualMultiv& c, const DualGraph& dg) { return -dg + c; }
inline DualGraph operator*(const DualGraph& dg, const DualMultiv& c) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), dp.get_py() * c));
    }
    return result;
}
inline DualGraph operator*(const DualMultiv& c, const DualGraph& dg) { return dg * c; }
inline DualGraph operator/(const DualGraph& dg, const DualMultiv& c) { return dg * (1.0 / c); }
inline DualGraph operator/(const DualMultiv& c, const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), c / dp.get_py()));
    }
    return result;
}

// ===========================================================================
//  Free operators — DualGraph vs DualGraph (pointwise)
// ===========================================================================

/// \brief Pointwise DualGraph + DualGraph.
///
/// Broadcasting rules: if one side has a single point, it is broadcast to
/// all points on the other side. Otherwise, the two graphs must have the
/// same number of points.
inline DualGraph operator+(const DualGraph& dg1, const DualGraph& dg2) {
    int n1 = dg1.get_n();
    int n2 = dg2.get_n();
    if (n1 == 0 || n2 == 0) {
        std::cerr << "DualGraph operator: 0 points!" << std::endl;
        return DualGraph();
    }
    if (n1 == 1) {
        DualGraph result(n2);
        DualPoint dp1 = dg1.get_point(0);
        for (int i = 0; i < n2; ++i) {
            DualPoint dp2 = dg2.get_point(i);
            result.set_point(i, DualPoint(dp2.get_px(), dp1.get_py() + dp2.get_py()));
        }
        return result;
    }
    if (n2 == 1) {
        DualGraph result(n1);
        DualPoint dp2 = dg2.get_point(0);
        for (int i = 0; i < n1; ++i) {
            DualPoint dp1 = dg1.get_point(i);
            result.set_point(i, DualPoint(dp1.get_px(), dp1.get_py() + dp2.get_py()));
        }
        return result;
    }
    int n = n1;
    if (n != n2) {
        std::cerr << "DualGraph operator+: size mismatch!" << std::endl;
        return DualGraph();
    }
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp1 = dg1.get_point(i);
        DualPoint dp2 = dg2.get_point(i);
        result.set_point(i, DualPoint(dp1.get_px(), dp1.get_py() + dp2.get_py()));
    }
    return result;
}

/// \brief Pointwise DualGraph - DualGraph.
inline DualGraph operator-(const DualGraph& dg1, const DualGraph& dg2) {
    return dg1 + (-dg2);
}

/// \brief Pointwise DualGraph * DualGraph.
inline DualGraph operator*(const DualGraph& dg1, const DualGraph& dg2) {
    int n1 = dg1.get_n();
    int n2 = dg2.get_n();
    if (n1 == 0 || n2 == 0) {
        std::cerr << "DualGraph operator: 0 points!" << std::endl;
        return DualGraph();
    }
    if (n1 == 1) {
        DualGraph result(n2);
        DualPoint dp1 = dg1.get_point(0);
        for (int i = 0; i < n2; ++i) {
            DualPoint dp2 = dg2.get_point(i);
            result.set_point(i, DualPoint(dp2.get_px(), dp1.get_py() * dp2.get_py()));
        }
        return result;
    }
    if (n2 == 1) {
        DualGraph result(n1);
        DualPoint dp2 = dg2.get_point(0);
        for (int i = 0; i < n1; ++i) {
            DualPoint dp1 = dg1.get_point(i);
            result.set_point(i, DualPoint(dp1.get_px(), dp1.get_py() * dp2.get_py()));
        }
        return result;
    }
    int n = n1;
    if (n != n2) {
        std::cerr << "DualGraph operator*: size mismatch!" << std::endl;
        return DualGraph();
    }
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp1 = dg1.get_point(i);
        DualPoint dp2 = dg2.get_point(i);
        result.set_point(i, DualPoint(dp1.get_px(), dp1.get_py() * dp2.get_py()));
    }
    return result;
}

/// \brief Pointwise DualGraph / DualGraph.
inline DualGraph operator/(const DualGraph& dg1, const DualGraph& dg2) {
    return dg1 * (1.0 / dg2);
}

// ===========================================================================
//  Math functions on DualGraph (pointwise application)
// ===========================================================================

/// \brief Pointwise power: \f$ y_i^c \f$.
inline DualGraph pow(const DualGraph& dg, double c) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), pow(dp.get_py(), c)));
    }
    return result;
}
inline DualGraph sqrt(const DualGraph& dg) { return pow(dg, 0.5); }
inline DualGraph abs(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), abs(dp.get_py())));
    }
    return result;
}
inline DualGraph exp(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), exp(dp.get_py())));
    }
    return result;
}
inline DualGraph log(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), log(dp.get_py())));
    }
    return result;
}
inline DualGraph sin(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), sin(dp.get_py())));
    }
    return result;
}
inline DualGraph cos(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), cos(dp.get_py())));
    }
    return result;
}
inline DualGraph tan(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), tan(dp.get_py())));
    }
    return result;
}
inline DualGraph atan(const DualGraph& dg) {
    int n = dg.get_n();
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp = dg.get_point(i);
        result.set_point(i, DualPoint(dp.get_px(), atan(dp.get_py())));
    }
    return result;
}
/// \brief Pointwise two-argument arctangent.
inline DualGraph atan2(const DualGraph& dg1, const DualGraph& dg2) {
    int n = dg1.get_n();
    if (n != dg2.get_n()) {
        std::cerr << "DualGraph atan2: size mismatch!" << std::endl;
        return DualGraph();
    }
    DualGraph result(n);
    for (int i = 0; i < n; ++i) {
        DualPoint dp1 = dg1.get_point(i);
        DualPoint dp2 = dg2.get_point(i);
        result.set_point(i, DualPoint(dp1.get_px(),
                           atan2(dp1.get_py(), dp2.get_py())));
    }
    return result;
}

#endif  // DualGraph_H
