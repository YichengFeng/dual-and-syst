/**
 * \file SystGraph.h
 * \brief Systematic-uncertainty graph class (ROOT-independent base).
 *
 * \details Holds one DualGraph as the default and a map of DualGraph
 * variations keyed by integer variation ID (vid). Supports systematic
 * uncertainty evaluation using Simple and Barlow (subset / independent /
 * fully-correlated) prescriptions for both one-sided and two-sided
 * systematics. Arithmetic operators align variations by vid across two
 * SystGraphs.
 *
 * Barlow's prescription reference: https://arxiv.org/abs/hep-ex/0207026
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * \see TSystGraph (ROOT-aware counterpart), DualGraph, SystMergeMode
 */

#ifndef SystGraph_H
#define SystGraph_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "DualGraph.h"
#include "ValueErrors.h"

// ===========================================================================
//  SystMergeMode — how to combine a variation with the default
// ===========================================================================

/// \enum SystMergeMode
/// \brief Specifies how a systematic variation is combined with the default.
///
/// \par One-sided vs Two-sided
/// **One-sided**: the variation contributes to only one side of the error
/// (sign determines direction). **Two-sided**: the variation contributes
/// symmetrically to both +/- sides.
///
/// \par Simple vs Barlow
/// **Simple**: systematic = v - d (difference only).
/// **Barlow**: statistical uncertainty is subtracted from the difference.
/// Let \f$ v \f$ = variation value, \f$ d \f$ = default value,
/// \f$ e_v \f$ = stat. error on \f$ v \f$, \f$ e_d \f$ = stat. error on \f$ d \f$.
///
///   - **Subset**: \f$ \mathrm{sign}(v-d) \cdot \sqrt{(v-d)^2 - |e_v^2 - e_d^2|} \f$ (->0 if arg < 0)
///   - **Independent**: \f$ \mathrm{sign}(v-d) \cdot \sqrt{(v-d)^2 - |e_v^2 + e_d^2|} \f$ (->0 if arg < 0)
///   - **FullyCorrelated**: \f$ \mathrm{sign}(v-d) \cdot \sqrt{(v-d)^2 - |e_v - e_d|^2} \f$ (->0 if arg < 0)
enum class SystMergeMode {
    // ---- one-sided ----
    kOneSided_Simple,                    ///< One-sided, simple difference
    kOneSided_Barlow_Subset,             ///< One-sided, Barlow subset prescription
    kOneSided_Barlow_Independent,        ///< One-sided, Barlow independent prescription
    kOneSided_Barlow_FullyCorrelated,    ///< One-sided, Barlow fully-correlated prescription

    // ---- two-sided ----
    kTwoSided_Simple,                    ///< Two-sided, simple difference
    kTwoSided_Barlow_Subset,             ///< Two-sided, Barlow subset prescription
    kTwoSided_Barlow_Independent,        ///< Two-sided, Barlow independent prescription
    kTwoSided_Barlow_FullyCorrelated,    ///< Two-sided, Barlow fully-correlated prescription

    kNull  ///< Sentinel: no mode / end-of-valid-range marker
};

/// \brief Check whether a SystMergeMode is one-sided.
/// \return true for one-sided modes, false for two-sided or kNull.
inline bool is_one_sided(SystMergeMode sm) noexcept {
    return (sm == SystMergeMode::kOneSided_Simple ||
            sm == SystMergeMode::kOneSided_Barlow_Subset ||
            sm == SystMergeMode::kOneSided_Barlow_Independent ||
            sm == SystMergeMode::kOneSided_Barlow_FullyCorrelated);
}

// ===========================================================================
//  SystVariation — name, graph, weight, and merge mode for one variation
// ===========================================================================

/// \struct SystVariation
/// \brief Describes a single systematic variation.
///
/// Each variation has a human-readable name, a DualGraph, a weight factor,
/// and a SystMergeMode specifying how it is combined with the default.
struct SystVariation {
    std::string      name;    ///< Human-readable variation name
    DualGraph        graph;   ///< The variation graph
    double           weight   = 1.0;  ///< Weight factor (default 1.0)
    SystMergeMode    mode     = SystMergeMode::kOneSided_Simple;  ///< Merge mode
};

// ===========================================================================
//  SystResult — output of evaluate_systematics()
// ===========================================================================

/// \struct SystResult
/// \brief Result of a systematic uncertainty evaluation.
///
/// Contains per-variation raw differences, Barlow-corrected signed
/// systematics, per-bin statistical components, and quadrature-summed
/// total lower/upper uncertainties.
struct SystResult {
    std::vector<int>                  vids;          ///< Variation IDs evaluated
    std::vector<std::vector<double>>  var_diff;      ///< [v][bin] raw v - d (before Barlow correction)
    std::vector<std::vector<double>>  var_syst;      ///< [v][bin] signed systematic after Barlow
    /// one-sided: >0 -> upper only, <0 -> lower only (magnitude = |var_syst|)
    /// two-sided: always >=0 (symmetric)
    std::vector<std::vector<double>>  var_stat;      ///< [v][bin] statistical component of Barlow correction
    std::vector<double>               total_lower;   ///< [bin] total lower uncertainty (quadrature sum)
    std::vector<double>               total_upper;   ///< [bin] total upper uncertainty (quadrature sum)
};

/// \class SystGraph
/// \brief A graph with systematic uncertainty variations (ROOT-independent).
///
/// Holds a default DualGraph and a map of systematic variation DualGraphs
/// keyed by integer vids. Provides systematic evaluation via evaluate_systematics(),
/// LaTeX table generation via to_latex(), and full arithmetic operator support
/// with vid-aligned variation propagation.
///
/// \see TSystGraph (ROOT-aware counterpart), SystMergeMode, SystVariation, SystResult
class SystGraph
{
protected:
    DualGraph                       def_;          ///< Default (central) graph
    std::map<int, SystVariation>    var_;          ///< Map: vid -> variation
    SystMergeMode                   merge_mode_ = SystMergeMode::kNull;  ///< Global merge mode override

    // --- vid bookkeeping ---
    static int          next_vid_;   ///< Auto-incrementing variation ID counter
    static const int    vid_max_ = INT_MAX;  ///< Maximum allowed variation ID

public:
    // --- constructors ---

    /// \brief Default constructor — empty graph.
    SystGraph() = default;

    /// \brief Construct with a default graph only (no variations).
    explicit SystGraph(const DualGraph& def) : def_(def) {}

    /// \brief Construct with default graph and variation map.
    SystGraph(const DualGraph& def,
              const std::map<int, SystVariation>& var)
        : def_(def), var_(var) {}

    /// \brief Virtual destructor.
    virtual ~SystGraph() = default;

    // --- vid management ---

    /// \brief Allocate a new unique variation ID.
    /// \return the new vid.
    static int auto_new_vid() {
        if (next_vid_ >= vid_max_ - 1) {
            std::cerr << "SystGraph::auto_new_vid(): vid out of range! "
                      << next_vid_ << "/" << vid_max_ - 1 << std::endl;
        } else {
            ++next_vid_;
        }
        return next_vid_;
    }

    /// \brief Get the current variation ID counter.
    static int get_current_vid() noexcept { return next_vid_; }

    // --- default ---

    /// \brief Get the default graph (const).
    const DualGraph& get_def() const noexcept { return def_; }
    /// \brief Get a mutable reference to the default graph.
    DualGraph&       get_def()       noexcept { return def_; }

    /// \brief Set the default graph.
    void set_def(const DualGraph& def) { def_ = def; }

    // --- global merge mode override ---

    /// \brief Get the global merge mode override.
    ///
    /// When kNull (default), each variation uses its own mode.
    /// When set to anything else, all variations use this mode instead.
    SystMergeMode get_merge_mode() const noexcept { return merge_mode_; }

    /// \brief Set the global merge mode override.
    void set_merge_mode(SystMergeMode s) noexcept { merge_mode_ = s; }

    // --- variation map ---

    /// \brief Get the full variation map (const).
    const std::map<int, SystVariation>& get_var_all() const noexcept { return var_; }

    /// \brief Replace the entire variation map.
    void set_var_all(const std::map<int, SystVariation>& var) { var_ = var; }

    /// \brief Get the list of all variation IDs.
    std::vector<int> get_list() const {
        std::vector<int> list;
        list.reserve(var_.size());
        for (const auto& kv : var_) list.push_back(kv.first);
        return list;
    }

    /// \brief Check whether a variation ID exists.
    bool check_var(int vid) const {
        return var_.count(vid) != 0;
    }

    // --- variation accessors / mutators ---

    /// \brief Get a variation by ID (const).
    /// \return the SystVariation, or a dummy if not found.
    const SystVariation& get_var(int vid) const {
        auto it = var_.find(vid);
        if (it == var_.end()) {
            std::cerr << "SystGraph::get_var(): vid " << vid
                      << " not found!" << std::endl;
            static const SystVariation dummy;
            return dummy;
        }
        return it->second;
    }

    /// \brief Add a variation with a specific vid.
    void add_var(int vid, const SystVariation& sv) {
        if (var_.count(vid) != 0) {
            std::cerr << "SystGraph::add_var(): vid " << vid
                      << " already exists!" << std::endl;
            return;
        }
        var_.emplace(vid, sv);
    }

    /// \brief Add a variation with an auto-assigned vid.
    /// \return the assigned vid.
    int add_var(const SystVariation& sv) {
        int vid = auto_new_vid();
        add_var(vid, sv);
        return vid;
    }

    /// \brief Set (overwrite or insert) a variation.
    void set_var(int vid, const SystVariation& sv) {
        var_[vid] = sv;
    }

    /// \brief Remove a single variation by vid.
    void erase_var(int vid) {
        var_.erase(vid);
    }

    /// \brief Remove multiple variations by vid.
    void erase_var(const std::vector<int>& vids) {
        for (int v : vids) erase_var(v);
    }

    /// \brief Remove all variations.
    void clear_var() { var_.clear(); }

    // --- extract a subset of variations ---

    /// \brief Extract a subset of variations by vid.
    /// \param vidlist list of vids to keep
    /// \return a new SystGraph with only the specified variations
    SystGraph get_sub(const std::vector<int>& vidlist) const {
        SystGraph sg(def_);
        for (int vid : vidlist) {
            auto it = var_.find(vid);
            if (it != var_.end()) sg.add_var(vid, it->second);
        }
        return sg;
    }

    // --- systematic evaluation ---

    /// \brief Evaluate systematics for all variations.
    /// \return a SystResult with per-variation and total uncertainties.
    SystResult evaluate_systematics() const;

    /// \brief Evaluate systematics for a specified subset of vids.
    /// \param vids list of variation IDs to evaluate
    /// \return a SystResult (skips non-existent vids).
    SystResult evaluate_systematics(const std::vector<int>& vids) const;

    // --- LaTeX table ---

    /// \brief Generate a LaTeX tabular of systematic uncertainties.
    ///
    /// Default layout: rows = variations (+ Total), columns = bins.
    /// Pass transpose = true to swap axes (rows = bins, columns = variations).
    /// \param transpose if true, bins are rows and variations are columns
    /// \return vector of strings, one per line of the LaTeX tabular
    std::vector<std::string> to_latex(bool transpose = false) const;

    /// \brief Print the LaTeX table to stdout.
    void print(bool transpose = false) const {
        for (const auto& line : to_latex(transpose))
            std::cout << line << std::endl;
    }

    // --- internal helper: pick variation graph if present, else default ---
    /// \brief Pick the variation graph for vid if it exists, otherwise the default.
    const DualGraph& pick(int vid) const {
        auto it = var_.find(vid);
        return (it == var_.end()) ? def_ : it->second.graph;
    }

    // --- average over bins (weighted, inverse-variance) ---
    /// \brief Weighted average over specified bins for default and all variations.
    SystGraph average_bin(const std::vector<int>& range) const {
        SystGraph sg(*this);
        sg.set_def(def_.average_bin(range));
        for (auto& kv : sg.var_)
            kv.second.graph = kv.second.graph.average_bin(range);
        return sg;
    }

    // --- average over bins (simple, equal-weight) ---
    /// \brief Equal-weight average over specified bins.
    SystGraph average_bin_simple(const std::vector<int>& range) const {
        SystGraph sg(*this);
        sg.set_def(def_.average_bin_simple(range));
        for (auto& kv : sg.var_)
            kv.second.graph = kv.second.graph.average_bin_simple(range);
        return sg;
    }

    // --- select subset of bins ---
    /// \brief Select a subset of bins for default and all variations.
    SystGraph select_bin(const std::vector<int>& range) const {
        SystGraph sg(*this);
        sg.set_def(def_.select_bin(range));
        for (auto& kv : sg.var_)
            kv.second.graph = kv.second.graph.select_bin(range);
        return sg;
    }

};

// initialize the static counter (header-only: inline avoids ODR clash)
inline int SystGraph::next_vid_ = 0;

// ===========================================================================
//  SystGraph::evaluate_systematics  —  core systematic evaluation
// ===========================================================================

/// \brief Evaluate systematics for all variations.
inline SystResult SystGraph::evaluate_systematics() const
{
    std::vector<int> all_vids;
    for (const auto& kv : var_) all_vids.push_back(kv.first);
    return evaluate_systematics(all_vids);
}

/// \brief Evaluate systematics for specified vids.
///
/// For each variation, computes the raw difference \f$ v - d \f$,
/// applies the Barlow correction if configured, and accumulates
/// quadrature-summed totals. One-sided systematics contribute to
/// either the upper or lower total based on the sign of the systematic.
inline SystResult SystGraph::evaluate_systematics(const std::vector<int>& vids) const
{
    SystResult result;

    std::vector<int> eval_vids;
    for (int vid : vids) {
        if (var_.count(vid)) {
            eval_vids.push_back(vid);
        } else {
            std::cerr << "SystGraph::evaluate_systematics(): vid "
                      << vid << " not found, skipped." << std::endl;
        }
    }

    const int nvar = static_cast<int>(eval_vids.size());
    const int nbin = def_.get_n();
    if (nvar == 0 || nbin == 0) return result;

    result.vids = eval_vids;
    result.var_syst.assign(nvar, std::vector<double>(nbin, 0.0));
    result.var_stat.assign(nvar, std::vector<double>(nbin, 0.0));
    result.var_diff.assign(nvar, std::vector<double>(nbin, 0.0));
    result.total_lower.assign(nbin, 0.0);
    result.total_upper.assign(nbin, 0.0);

    for (int iv = 0; iv < nvar; ++iv) {
        int vid = eval_vids[iv];
        const SystVariation& sv = var_.find(vid)->second;
        const DualGraph& vg = sv.graph;
        SystMergeMode sm = (merge_mode_ != SystMergeMode::kNull) ? merge_mode_ : sv.mode;
        double w = sv.weight;

        if (sm == SystMergeMode::kNull) continue;

        const bool one_sided = ::is_one_sided(sm);

        for (int ib = 0; ib < nbin; ++ib) {
            double d  = def_.get_point(ib).get_py().get_value();
            double ed = def_.get_point(ib).get_py().get_uncertainty();
            double v  = vg.get_point(ib).get_py().get_value();
            double ev = vg.get_point(ib).get_py().get_uncertainty();

            double diff = v - d;
            result.var_diff[iv][ib] = diff;
            double raw  = diff;
            double stat = 0.0;

            switch (sm) {
            case SystMergeMode::kOneSided_Barlow_Subset:
            case SystMergeMode::kTwoSided_Barlow_Subset: {
                double dvar = std::fabs(ev * ev - ed * ed);
                stat = std::sqrt(dvar);
                double arg = diff * diff - dvar;
                raw = (arg > 0.0) ? std::copysign(std::sqrt(arg), diff) : 0.0;
                break;
            }
            case SystMergeMode::kOneSided_Barlow_Independent:
            case SystMergeMode::kTwoSided_Barlow_Independent: {
                double dvar = ev * ev + ed * ed;
                stat = std::sqrt(dvar);
                double arg = diff * diff - dvar;
                raw = (arg > 0.0) ? std::copysign(std::sqrt(arg), diff) : 0.0;
                break;
            }
            case SystMergeMode::kOneSided_Barlow_FullyCorrelated:
            case SystMergeMode::kTwoSided_Barlow_FullyCorrelated: {
                stat = std::fabs(ev - ed);
                double arg = diff * diff - stat * stat;
                raw = (arg > 0.0) ? std::copysign(std::sqrt(arg), diff) : 0.0;
                break;
            }
            default:
                break;
            }

            double syst = (one_sided) ? raw : std::fabs(raw);
            result.var_syst[iv][ib] = syst;
            result.var_stat[iv][ib] = stat;

            if (one_sided) {
                if (syst > 0.0) result.total_upper[ib] += w * syst * syst;
                else            result.total_lower[ib] += w * syst * syst;
            } else {
                result.total_lower[ib] += w * syst * syst;
                result.total_upper[ib] += w * syst * syst;
            }
        }
    }

    for (int ib = 0; ib < nbin; ++ib) {
        result.total_lower[ib] = std::sqrt(result.total_lower[ib]);
        result.total_upper[ib] = std::sqrt(result.total_upper[ib]);
    }

    return result;
}

// ===========================================================================
//  SystGraph::to_latex  —  format SystResult as LaTeX tabular
// ===========================================================================

inline std::vector<std::string> SystGraph::to_latex(bool transpose) const
{
    SystResult result = evaluate_systematics();
    const int nvar = static_cast<int>(result.vids.size());
    const int nbin = static_cast<int>(result.total_lower.size());

    std::vector<std::string> lines;
    if (nvar == 0 || nbin == 0) return lines;

    auto fmt_num = [](double v) -> std::string {
        if (v == 0.0) return "0";
        ValueErrors ve(v);
        return ve.print_err0("L");
    };

    auto fmt_val = [](double v) -> std::string {
        if (v == 0.0) return "0";
        ValueErrors ve(v);
        if (ve.get_precision() < 3) ve.set_precision(3);
        return ve.print_err0("L");
    };

    auto fmt_diff = [&](double diff) -> std::string {
        if (diff == 0.0) return "0";
        double ad = std::fabs(diff);
        if (diff > 0.0)
            return "$+" + fmt_num(ad) + "$";
        else
            return "$-" + fmt_num(ad) + "$";
    };

    auto fmt_cell = [&](double lo, double hi) -> std::string {
        if (lo == 0.0 && hi == 0.0) return "0";
        double max_err = std::max(std::fabs(lo), std::fabs(hi));
        ValueErrors ve_ref(0.0, max_err);
        int prec = ve_ref.get_precision();

        auto fmt_side = [prec](double v) -> std::string {
            ValueErrors ve(v);
            ve.set_precision(prec);
            return ve.print_err0("L");
        };

        if (lo == hi) {
            return "$\\pm " + fmt_side(lo) + "$";
        } else if (lo == 0.0) {
            return "$+" + fmt_side(hi) + "$";
        } else if (hi == 0.0) {
            return "$-" + fmt_side(lo) + "$";
        } else {
            return "$_{-" + fmt_side(lo)
                 + "}^{+" + fmt_side(hi) + "}$";
        }
    };

    std::vector<std::string> names(nvar);
    std::vector<std::string> weights(nvar);
    std::vector<std::string> sides(nvar);
    for (int iv = 0; iv < nvar; ++iv) {
        auto it = var_.find(result.vids[iv]);
        if (it == var_.end()) {
            names[iv]   = "vid_" + std::to_string(result.vids[iv]);
            weights[iv] = "1";
            sides[iv]   = "?";
        } else {
            SystMergeMode sm = (merge_mode_ != SystMergeMode::kNull)
                               ? merge_mode_ : it->second.mode;
            names[iv] = it->second.name;
            {
                std::ostringstream ws;
                ws << it->second.weight;
                weights[iv] = ws.str();
            }
            bool one_sided = ::is_one_sided(sm);
            sides[iv] = one_sided ? "1-sided" : "2-sided";
        }
    }

    auto bin_label = [&](int ib) -> std::string {
        double x = def_.get_point(ib).get_px().get_value();
        ValueErrors ve(x);
        return ve.print_err0("L", true);
    };

    if (!transpose) {
        std::string coldef = "l|cc|";
        for (int ib = 0; ib < nbin; ++ib) coldef += "c";
        lines.push_back("\\begin{tabular}{" + coldef + "}");
        lines.push_back("\\hline");

        std::string hdr = "Variation & Weight & Side";
        for (int ib = 0; ib < nbin; ++ib)
            hdr += " & " + bin_label(ib);
        hdr += " \\\\";
        lines.push_back(hdr);
        lines.push_back("\\hline");

        std::string def_row = "Default & & & ";
        for (int ib = 0; ib < nbin; ++ib) {
            if (ib > 0) def_row += " & ";
            def_row += fmt_val(def_.get_point(ib).get_py().get_value());
        }
        def_row += " \\\\";
        lines.push_back(def_row);

        std::string stat_row = "Stat. unc. & & & ";
        for (int ib = 0; ib < nbin; ++ib) {
            if (ib > 0) stat_row += " & ";
            stat_row += fmt_num(def_.get_point(ib).get_py().get_uncertainty());
        }
        stat_row += " \\\\";
        lines.push_back(stat_row);

        lines.push_back("\\hline");
        for (int iv = 0; iv < nvar; ++iv) {
            std::string row = names[iv] + " & " + weights[iv]
                            + " & " + sides[iv];
            for (int ib = 0; ib < nbin; ++ib)
                row += " & " + fmt_diff(result.var_diff[iv][ib]);
            row += " \\\\";
            lines.push_back(row);
        }

        lines.push_back("\\hline");
        std::string tot = "Total & & & ";
        for (int ib = 0; ib < nbin; ++ib) {
            if (ib > 0) tot += " & ";
            tot += fmt_cell(result.total_lower[ib],
                            result.total_upper[ib]);
        }
        tot += " \\\\";
        lines.push_back(tot);
        lines.push_back("\\hline");

        lines.push_back("\\end{tabular}");

    } else {
        std::string coldef = "l|cc|";
        for (int iv = 0; iv < nvar; ++iv) coldef += "c";
        coldef += "c";
        lines.push_back("\\begin{tabular}{" + coldef + "}");
        lines.push_back("\\hline");

        std::string hdr1 = "Bin & Default & Stat. unc.";
        for (int iv = 0; iv < nvar; ++iv)
            hdr1 += " & " + names[iv];
        hdr1 += " & Total \\\\";
        lines.push_back(hdr1);
        lines.push_back("\\hline");

        std::string hdr2 = " & & ";
        for (int iv = 0; iv < nvar; ++iv)
            hdr2 += " & " + weights[iv];
        hdr2 += " & \\\\";
        lines.push_back(hdr2);

        std::string hdr3 = " & & ";
        for (int iv = 0; iv < nvar; ++iv)
            hdr3 += " & " + sides[iv];
        hdr3 += " & \\\\";
        lines.push_back(hdr3);

        lines.push_back("\\hline");

        for (int ib = 0; ib < nbin; ++ib) {
            std::string row = bin_label(ib);
            row += " & " + fmt_val(def_.get_point(ib).get_py().get_value());
            row += " & " + fmt_num(def_.get_point(ib).get_py().get_uncertainty());
            for (int iv = 0; iv < nvar; ++iv)
                row += " & " + fmt_diff(result.var_diff[iv][ib]);
            row += " & " + fmt_cell(result.total_lower[ib],
                                    result.total_upper[ib]);
            row += " \\\\";
            lines.push_back(row);
        }

        lines.push_back("\\hline");
        lines.push_back("\\end{tabular}");
    }

    return lines;
}

// ===========================================================================
//  Helpers — apply a binary op across two SystGraphs, aligning by vid
// ===========================================================================

namespace syst_detail {

/// \brief Apply a binary operator to two SystGraphs, aligning variations by vid.
///
/// Vids absent on one side fall back to that side's default graph.
template <typename Op>
inline SystGraph apply_binary(const SystGraph& a, const SystGraph& b, Op op) {
    SystGraph out(op(a.get_def(), b.get_def()));
    const auto& va = a.get_var_all();
    const auto& vb = b.get_var_all();
    auto ia = va.begin();
    auto ib = vb.begin();
    while (ia != va.end() && ib != vb.end()) {
        if (ia->first < ib->first) {
            SystVariation sv = ia->second;
            sv.graph = op(ia->second.graph, b.get_def());
            out.add_var(ia->first, sv);
            ++ia;
        } else if (ia->first > ib->first) {
            SystVariation sv = ib->second;
            sv.graph = op(a.get_def(), ib->second.graph);
            out.add_var(ib->first, sv);
            ++ib;
        } else {
            SystVariation sv = ia->second;
            sv.graph = op(ia->second.graph, ib->second.graph);
            out.add_var(ia->first, sv);
            ++ia; ++ib;
        }
    }
    for (; ia != va.end(); ++ia) {
        SystVariation sv = ia->second;
        sv.graph = op(ia->second.graph, b.get_def());
        out.add_var(ia->first, sv);
    }
    for (; ib != vb.end(); ++ib) {
        SystVariation sv = ib->second;
        sv.graph = op(a.get_def(), ib->second.graph);
        out.add_var(ib->first, sv);
    }
    return out;
}

/// \brief Apply a unary operator to a SystGraph (default + all variations).
template <typename Op>
inline SystGraph apply_unary(const SystGraph& a, Op op) {
    SystGraph out(op(a.get_def()));
    for (const auto& kv : a.get_var_all()) {
        SystVariation sv = kv.second;
        sv.graph = op(kv.second.graph);
        out.add_var(kv.first, sv);
    }
    return out;
}

}  // namespace syst_detail

// ===========================================================================
//  Free operators — SystGraph vs scalar
// ===========================================================================

inline SystGraph operator+(const SystGraph& sg, double c) {
    return syst_detail::apply_unary(sg, [c](const DualGraph& g){ return g + c; });
}
inline SystGraph operator+(double c, const SystGraph& sg) { return sg + c; }
inline SystGraph operator-(const SystGraph& sg, double c) { return sg + (-c); }
inline SystGraph operator-(double c, const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [c](const DualGraph& g){ return c - g; });
}
inline SystGraph operator*(const SystGraph& sg, double c) {
    return syst_detail::apply_unary(sg, [c](const DualGraph& g){ return g * c; });
}
inline SystGraph operator*(double c, const SystGraph& sg) { return sg * c; }
inline SystGraph operator/(const SystGraph& sg, double c) { return sg * (1.0 / c); }
inline SystGraph operator/(double c, const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [c](const DualGraph& g){ return c / g; });
}

// ===========================================================================
//  Free operators — SystGraph vs DualMultiv
// ===========================================================================

inline SystGraph operator+(const SystGraph& sg, const DualMultiv& c) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return g + c; });
}
inline SystGraph operator+(const DualMultiv& c, const SystGraph& sg) { return sg + c; }
inline SystGraph operator-(const SystGraph& sg, const DualMultiv& c) { return sg + (-c); }
inline SystGraph operator-(const DualMultiv& c, const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return c - g; });
}
inline SystGraph operator*(const SystGraph& sg, const DualMultiv& c) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return g * c; });
}
inline SystGraph operator*(const DualMultiv& c, const SystGraph& sg) { return sg * c; }
inline SystGraph operator/(const SystGraph& sg, const DualMultiv& c) { return sg * (1.0 / c); }
inline SystGraph operator/(const DualMultiv& c, const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return c / g; });
}

// ===========================================================================
//  Free operators — SystGraph vs DualGraph
// ===========================================================================

inline SystGraph operator+(const SystGraph& sg, const DualGraph& c) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return g + c; });
}
inline SystGraph operator+(const DualGraph& c, const SystGraph& sg) { return sg + c; }
inline SystGraph operator-(const SystGraph& sg, const DualGraph& c) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return g - c; });
}
inline SystGraph operator-(const DualGraph& c, const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return c - g; });
}
inline SystGraph operator*(const SystGraph& sg, const DualGraph& c) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return g * c; });
}
inline SystGraph operator*(const DualGraph& c, const SystGraph& sg) { return sg * c; }
inline SystGraph operator/(const SystGraph& sg, const DualGraph& c) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return g / c; });
}
inline SystGraph operator/(const DualGraph& c, const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [&c](const DualGraph& g){ return c / g; });
}

// ===========================================================================
//  Unary operators
// ===========================================================================

/// \brief Unary plus (identity).
inline SystGraph operator+(const SystGraph& sg) { return sg; }

/// \brief Unary minus (negates default and all variations).
inline SystGraph operator-(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return -g; });
}

// ===========================================================================
//  Free operators — SystGraph vs SystGraph (vid-aligned)
// ===========================================================================

/// \brief vid-aligned SystGraph + SystGraph.
inline SystGraph operator+(const SystGraph& a, const SystGraph& b) {
    return syst_detail::apply_binary(a, b,
        [](const DualGraph& x, const DualGraph& y){ return x + y; });
}
/// \brief vid-aligned SystGraph - SystGraph.
inline SystGraph operator-(const SystGraph& a, const SystGraph& b) {
    return syst_detail::apply_binary(a, b,
        [](const DualGraph& x, const DualGraph& y){ return x - y; });
}
/// \brief vid-aligned SystGraph * SystGraph.
inline SystGraph operator*(const SystGraph& a, const SystGraph& b) {
    return syst_detail::apply_binary(a, b,
        [](const DualGraph& x, const DualGraph& y){ return x * y; });
}
/// \brief vid-aligned SystGraph / SystGraph.
inline SystGraph operator/(const SystGraph& a, const SystGraph& b) {
    return syst_detail::apply_binary(a, b,
        [](const DualGraph& x, const DualGraph& y){ return x / y; });
}

// ===========================================================================
//  Math functions on SystGraph
// ===========================================================================

/// \brief Pointwise power on default and all variations.
inline SystGraph pow(const SystGraph& sg, double c) {
    return syst_detail::apply_unary(sg, [c](const DualGraph& g){ return pow(g, c); });
}
inline SystGraph sqrt(const SystGraph& sg) { return pow(sg, 0.5); }
inline SystGraph abs(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return abs(g); });
}
inline SystGraph exp(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return exp(g); });
}
inline SystGraph log(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return log(g); });
}
inline SystGraph sin(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return sin(g); });
}
inline SystGraph cos(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return cos(g); });
}
inline SystGraph tan(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return tan(g); });
}
inline SystGraph atan(const SystGraph& sg) {
    return syst_detail::apply_unary(sg, [](const DualGraph& g){ return atan(g); });
}
/// \brief Pointwise vid-aligned atan2 on two SystGraphs.
inline SystGraph atan2(const SystGraph& a, const SystGraph& b) {
    return syst_detail::apply_binary(a, b,
        [](const DualGraph& x, const DualGraph& y){ return atan2(x, y); });
}

#endif  // SystGraph_H
