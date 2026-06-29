/**
 * \file TDualGraph.h
 * \brief ROOT-aware DualGraph with TGraphErrors mirror, I/O, and plotting.
 *
 * \details Extends DualGraph with ROOT integration: constructors from TH1
 * and TGraphErrors, automatic synchronization of a TGraphErrors mirror,
 * and TFile persistence via TObject inheritance. The "T" prefix follows
 * ROOT naming convention (TGraphErrors, TH1, …).
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   - Inherits from DualGraph (ROOT-free base) + TObject, enabling Write()
 *   - is_updated_ tracks whether graph_ is in sync with points_
 *   - Overrides set_points / set_point / reset to maintain is_updated_ flag
 *   - ROOT constructors forward to base reset(), then populate via set_derivative()
 *
 * \see DualGraph, TSystGraph
 */

#ifndef TDualGraph_H
#define TDualGraph_H

#include "DualGraph.h"

#include "TBuffer.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"


// ===========================================================================
//  TDualGraph — DualGraph with ROOT TGraphErrors mirror and I/O
// ===========================================================================

/// \class TDualGraph
/// \brief DualGraph with a ROOT TGraphErrors mirror for plotting and I/O.
///
/// The internal TGraphErrors (`graph_`) is lazily synced from the DualGraph
/// points via calc(). The is_updated_ flag tracks sync state. Supports
/// construction from TH1 and TGraphErrors, ROOT TFile I/O, and x-axis shifts.
///
/// \see DualGraph (ROOT-independent base), TSystGraph (systematic-uncertainty graphs)
class TDualGraph : public DualGraph, public TObject
{
private:
    bool          is_updated_ = false;  ///< Whether graph_ is in sync with points_
    TGraphErrors  graph_;              ///< ROOT TGraphErrors mirror

public:
    // --- constructors ---

    /// \brief Default constructor.
    TDualGraph() = default;

    /// \brief Construct with \p n default points.
    explicit TDualGraph(int n) : DualGraph(n) {}

    /// \brief Construct from a TH1 histogram.
    /// \param h the input histogram
    /// \param xerror whether to set x-errors from bin widths
    TDualGraph(const TH1& h, bool xerror = false) {
        reset(h.GetXaxis()->GetNbins());
        set_derivative(h, true, xerror, true);
    }

    /// \brief Construct from a TGraphErrors.
    /// \param g the input graph
    /// \param xerror whether to propagate x-errors
    TDualGraph(const TGraphErrors& g, bool xerror = false) {
        reset(g.GetN());
        set_derivative(g, true, xerror, true);
    }

    /// \brief Construct from a TH1 with explicit systematic ID.
    /// \param index the systematic variation ID
    /// \param h the input histogram
    /// \param xerror whether to set x-errors from bin widths
    TDualGraph(int index, const TH1& h, bool xerror = false) {
        reset(h.GetXaxis()->GetNbins());
        set_derivative(index, h, true, xerror, true);
    }

    /// \brief Construct from a TGraphErrors with explicit systematic ID.
    /// \param index the systematic variation ID
    /// \param g the input graph
    /// \param xerror whether to propagate x-errors
    TDualGraph(int index, const TGraphErrors& g, bool xerror = false) {
        reset(g.GetN());
        set_derivative(index, g, true, xerror, true);
    }

    /// \brief Construct from a vector of DualPoints.
    explicit TDualGraph(const std::vector<DualPoint>& points)
        : DualGraph(points) {}

    // --- copy constructor (fast, no ID sync) ---
    TDualGraph(const TDualGraph&) = default;

    // --- converting constructor from DualGraph ---
    /// \brief Convert a DualGraph to a TDualGraph.
    TDualGraph(const DualGraph& dg) : DualGraph(dg) {}

    // --- construct from pointer (e.g. from TFile::Get) and sync IDs ---
    /// \brief Construct from a pointer (e.g. from TFile::Get) and sync
    ///        DualMultiv ID counters to prevent collisions.
    explicit TDualGraph(const TDualGraph* p) : TDualGraph(*p) {
        sync_id_counter();
    }

    // --- overrides (maintain is_updated_ flag) ---

    /// \brief Replace all points — invalidates the ROOT mirror.
    void set_points(const std::vector<DualPoint>& points) override {
        DualGraph::set_points(points);
        is_updated_ = false;
    }

    /// \brief Set a single point — invalidates the ROOT mirror.
    void set_point(int i, const DualPoint& dp) override {
        DualGraph::set_point(i, dp);
        is_updated_ = false;
    }

    /// \brief Reset to \p n default points — invalidates the ROOT mirror.
    void reset(int n = 0) override {
        DualGraph::reset(n);
        is_updated_ = false;
    }

    // --- accessors ---

    /// \brief Check whether the ROOT graph mirror is up to date.
    bool get_is_updated() const noexcept { return is_updated_; }

    /// \brief Get the ROOT TGraphErrors mirror (const).
    const TGraphErrors& get_graph() const noexcept { return graph_; }
    /// \brief Get a mutable reference to the ROOT TGraphErrors mirror.
    TGraphErrors&       get_graph()       noexcept { return graph_; }

    // --- set_derivative from TH1 ---

    /// \brief Set derivatives from a TH1 histogram (specific systematic ID).
    /// \param index the systematic variation ID
    /// \param h the histogram
    /// \param xreset reset x-values from bin centers
    /// \param xerror set x-errors from bin widths
    /// \param yreset reset y-values from bin contents
    void set_derivative(int index, const TH1& h,
                        bool xreset, bool xerror, bool yreset) {
        int n = h.GetXaxis()->GetNbins();
        if (!check_size(n)) return;
        for (int i = 0; i < n; ++i) {
            if (xreset) points_[i].set_px(DualNumber(
                h.GetXaxis()->GetBinCenter(i + 1),
                xerror ? h.GetXaxis()->GetBinWidth(i + 1) * 0.5 : 0.0));
            if (yreset) points_[i].get_py().set_value(h.GetBinContent(i + 1));
            points_[i].get_py().set_derivative(index, h.GetBinError(i + 1));
        }
        is_updated_ = false;
    }

    /// \brief Set derivatives from a TH1 (specific ID, preserve x/y values).
    void set_derivative(int index, const TH1& h) {
        set_derivative(index, h, false, false, false);
    }

    /// \brief Set derivatives from a TH1 (auto-assigned systematic ID).
    void set_derivative(const TH1& h,
                        bool xreset, bool xerror, bool yreset) {
        int n = h.GetXaxis()->GetNbins();
        if (!check_size(n)) return;
        for (int i = 0; i < n; ++i) {
            if (xreset) points_[i].set_px(DualNumber(
                h.GetXaxis()->GetBinCenter(i + 1),
                xerror ? h.GetXaxis()->GetBinWidth(i + 1) * 0.5 : 0.0));
            if (yreset) points_[i].get_py().set_value(h.GetBinContent(i + 1));
            points_[i].get_py().set_derivative(h.GetBinError(i + 1));
        }
        is_updated_ = false;
    }

    /// \brief Set derivatives from a TH1 (auto ID, preserve x/y).
    void set_derivative(const TH1& h) {
        set_derivative(h, false, false, false);
    }

    // --- set_derivative from TGraphErrors ---

    /// \brief Set derivatives from a TGraphErrors (specific systematic ID).
    void set_derivative(int index, const TGraphErrors& g,
                        bool xreset, bool xerror, bool yreset) {
        int n = g.GetN();
        if (!check_size(n)) return;
        for (int i = 0; i < n; ++i) {
            if (xreset) points_[i].set_px(DualNumber(
                g.GetX()[i],
                xerror ? g.GetEX()[i] : 0.0));
            if (yreset) points_[i].get_py().set_value(g.GetY()[i]);
            points_[i].get_py().set_derivative(index, g.GetEY()[i]);
        }
        is_updated_ = false;
    }

    void set_derivative(int index, const TGraphErrors& g) {
        set_derivative(index, g, false, false, false);
    }

    /// \brief Set derivatives from a TGraphErrors (auto-assigned systematic ID).
    void set_derivative(const TGraphErrors& g,
                        bool xreset, bool xerror, bool yreset) {
        int n = g.GetN();
        if (!check_size(n)) return;
        for (int i = 0; i < n; ++i) {
            if (xreset) points_[i].set_px(DualNumber(
                g.GetX()[i],
                xerror ? g.GetEX()[i] : 0.0));
            if (yreset) points_[i].get_py().set_value(g.GetY()[i]);
            points_[i].get_py().set_derivative(g.GetEY()[i]);
        }
        is_updated_ = false;
    }

    void set_derivative(const TGraphErrors& g) {
        set_derivative(g, false, false, false);
    }

    // --- sync TGraphErrors from internal points ---

    /// \brief Synchronize the ROOT TGraphErrors mirror from the internal
    ///        DualGraph points. Sets is_updated_ to true.
    void calc() {
        int n = static_cast<int>(points_.size());
        std::vector<double> vx(n), ex(n), vy(n), ey(n);
        for (int i = 0; i < n; ++i) {
            vx[i] = points_[i].get_px().get_value();
            ex[i] = points_[i].get_px().get_derivative();
            vy[i] = points_[i].get_py().get_value();
            ey[i] = points_[i].get_py().get_uncertainty();
        }
        graph_ = TGraphErrors(n, vx.data(), vy.data(), ex.data(), ey.data());
        is_updated_ = true;
    }

    // --- shift x-values ---

    /// \brief Shift graph x-values to match the internal DualGraph points.
    void shift_x() {
        if (!is_updated_) calc();
        for (int i = 0; i < graph_.GetN(); ++i) {
            graph_.GetX()[i] = points_[i].get_px().get_value();
        }
    }

    /// \brief Shift graph x-values by \p dx (additive offset).
    void shift_x(double dx) {
        if (!is_updated_) calc();
        for (int i = 0; i < graph_.GetN(); ++i) {
            graph_.GetX()[i] += dx;
        }
    }

    // --- merge uncertainties (quadrature-sum x and y errors from graph) ---
    /// \brief Merge uncertainties: reconstruct from the TGraphErrors
    ///        including both x and y errors.
    TDualGraph merge_uncertainty() {
        if (!is_updated_) calc();
        return TDualGraph(graph_, true);
    }

    // --- ROOT I/O ---

    /// \brief Write this object to a ROOT TFile.
    /// \param f the open TFile
    /// \param name the key name in the file
    void write_to(TFile* f, const char* name) const {
        f->WriteTObject(this, name);
    }

    /// \brief Sync DualMultiv ID counters after reading from file.
    ///
    /// Prevents DualMultiv::new_id() from producing IDs that collide
    /// with those already stored in the read-back DualMultiv objects.
    /// Call this after loading a TDualGraph from a TFile.
    void sync_id_counter() const {
        int max_id = 0;
        for (const auto& dp : points_) {
            for (int id : dp.get_py().get_ids()) {
                if (id > max_id) max_id = id;
            }
        }
        DualMultiv::set_id_counter(max_id);
    }

    ClassDefOverride(TDualGraph, 1)  ///< ROOT streamer info (version 1)
};
#endif  // TDualGraph_H
