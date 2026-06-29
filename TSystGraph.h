/**
 * \file TSystGraph.h
 * \brief ROOT-aware systematic-uncertainty graph with plotting and I/O.
 *
 * \details Extends SystGraph with ROOT TGraphErrors / TGraphAsymmErrors
 * output graphs, drawing, overlay plotting, stacked histogram generation,
 * LaTeX table formatting, and TFile I/O. It is the ROOT-aware counterpart
 * of SystGraph, analogous to how TDualGraph extends DualGraph.
 *
 * The "T" prefix follows ROOT naming convention (TGraphErrors, TH1, …).
 *
 * Systematic uncertainties are evaluated with Barlow's prescription:
 * https://arxiv.org/abs/hep-ex/0207026
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code/MySystGraph.h:
 *   - Inherits from SystGraph (ROOT-free base) + TObject, enabling Write()
 *   - Uses SystMergeMode enum instead of magic-number modes (1–6)
 *   - is_updated_ flag tracks whether ROOT output graphs are synced
 *   - Replaced VLAs with std::vector<double>
 *   - Uses ValueErrors instead of MyValuErrs
 *   - TDualGraph for ROOT-enabled default / variation graphs
 *   - std::string for names; const char* / TString at ROOT boundary
 *   - snake_case naming, noexcept on trivial methods, const-correct
 *   - Removed const from operator return-by-value (move semantics)
 *
 * \see SystGraph (ROOT-independent base), TDualGraph, SystMergeMode
 */

#ifndef TSystGraph_H
#define TSystGraph_H

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "SystGraph.h"
#include "TDualGraph.h"
#include "ValueErrors.h"

#include "TBox.h"
#include "TBuffer.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TObject.h"
#include "TString.h"


// ===========================================================================
//  TSystVariation — variation with ROOT-enabled TDualGraph
// ===========================================================================

/// \struct TSystVariation
/// \brief A systematic variation with a ROOT-enabled TDualGraph.
///
/// Analogous to SystVariation but uses TDualGraph for the variation graph,
/// enabling ROOT I/O. Sliced to SystVariation when stored in the base class.
///
/// \see SystVariation (ROOT-independent counterpart)
struct TSystVariation {
    std::string    name;    ///< Human-readable variation name
    TDualGraph     graph;   ///< The variation graph (ROOT-enabled)
    double         weight = 1.0;  ///< Weight factor
    SystMergeMode  mode   = SystMergeMode::kOneSided_Simple;  ///< Merge mode
};

// ===========================================================================
//  TSystGraph — ROOT-enabled systematic-graph class
// ===========================================================================

/// \class TSystGraph
/// \brief ROOT-enabled systematic-uncertainty graph.
///
/// Provides ROOT output graphs (stat_, syst_, comb_) synced via calc(),
/// canvas drawing (draw_stat, draw_syst, draw_comb), overlay plots
/// (make_plot), stacked systematic^2 histograms (make_hist), LaTeX
/// formatting, and TFile I/O via TObject inheritance.
///
/// \see SystGraph (ROOT-independent base), TDualGraph, TSystVariation
class TSystGraph : public SystGraph, public TObject
{
private:
    bool  is_updated_ = false;  ///< Whether ROOT output graphs are synced
    double width_      = 0.8;   ///< Systematic box half-width in x (for syst/comb graphs)

    /// x-axis limit filter (points outside [xlim_low_, xlim_high_) are excluded from calc)
    double xlim_low_  = -std::numeric_limits<double>::infinity();
    double xlim_high_ =  std::numeric_limits<double>::infinity();

    SystResult calc_result_;  ///< Cached result of last systematic evaluation

    // ROOT output graphs (synced by calc)
    TGraphErrors        stat_;  ///< Statistical-uncertainty graph
    TGraphAsymmErrors   syst_;  ///< Systematic-uncertainty graph
    TGraphAsymmErrors   comb_;  ///< Combined (stat + syst in quadrature) uncertainty graph

public:
    // --- ROOT output graph accessors (trigger calc if needed) ---

    /// \brief Get the statistical-uncertainty TGraphErrors (syncs if needed).
    TGraphErrors&       stat() { calc_if_not_updated(); return stat_; }
    /// \brief Get the systematic-uncertainty TGraphAsymmErrors (syncs if needed).
    TGraphAsymmErrors&  syst() { calc_if_not_updated(); return syst_; }
    /// \brief Get the combined-uncertainty TGraphAsymmErrors (syncs if needed).
    TGraphAsymmErrors&  comb() { calc_if_not_updated(); return comb_; }

    // --- constructors ---

    TSystGraph() = default;

    /// \brief Construct with a default TDualGraph.
    explicit TSystGraph(const TDualGraph& def) {
        set_def(def);
    }

    /// \brief Construct with a default TGraphErrors.
    explicit TSystGraph(const TGraphErrors& def) {
        set_def(TDualGraph(def));
    }

    /// \brief Construct with default and variation map.
    TSystGraph(const TDualGraph& def,
               const std::map<int, TSystVariation>& var) {
        set_def(def);
        for (const auto& kv : var)
            add_var(kv.first, kv.second);
    }

    // --- copy constructor (fast, no ID sync) ---
    TSystGraph(const TSystGraph&) = default;

    // --- construct from pointer (e.g. from TFile::Get) and sync IDs ---
    explicit TSystGraph(const TSystGraph* p) : TSystGraph(*p) {
        sync_id_counter();
    }

    ~TSystGraph() override = default;

    // --- update flag ---

    /// \brief Check whether ROOT output graphs are up to date.
    bool get_is_updated() const noexcept { return is_updated_; }
    void set_is_updated(bool v) noexcept { is_updated_ = v; }

    // --- reset ---

    /// \brief Reset to empty state.
    void reset() {
        SystGraph::clear_var();
        def_ = DualGraph();
        width_      = 0.8;
        is_updated_ = false;
        xlim_low_   = -std::numeric_limits<double>::infinity();
        xlim_high_  =  std::numeric_limits<double>::infinity();
        calc_result_ = SystResult();
    }

    // --- width ---

    /// \brief Set the systematic box half-width.
    ///
    /// Also updates the ROOT output graphs (syst_ and comb_) if they are
    /// already synced.
    void set_width(double w) {
        width_ = std::fabs(w);
        if (!is_updated_) return;
        for (int i = 0; i < syst_.GetN(); ++i) {
            syst_.GetEXlow()[i]   = width_;
            syst_.GetEXhigh()[i]  = width_;
            comb_.GetEXlow()[i]   = width_;
            comb_.GetEXhigh()[i]  = width_;
        }
    }

    /// \brief Get the systematic box half-width.
    double get_width() const noexcept { return width_; }

    // --- default graph ---

    /// \brief Set the default graph from a TDualGraph.
    /// Auto-computes width_ from the first two x-values.
    void set_def(const TDualGraph& def) {
        def_ = def;
        if (def_.get_n() >= 2) {
            width_ = 0.1 * std::fabs(
                def_.get_point(0).get_px().get_value() -
                def_.get_point(1).get_px().get_value());
        }
        is_updated_ = false;
    }

    /// \brief Set the default graph from a DualGraph.
    void set_def(const DualGraph& def) {
        def_ = def;
        if (def_.get_n() >= 2) {
            width_ = 0.1 * std::fabs(
                def_.get_point(0).get_px().get_value() -
                def_.get_point(1).get_px().get_value());
        }
        is_updated_ = false;
    }

    using SystGraph::get_def;   // const DualGraph& get_def()

    // --- mode override (adds is_updated_ tracking) ---

    /// \brief Set the global merge mode, invalidating the ROOT sync.
    void set_merge_mode(SystMergeMode sm) noexcept {
        merge_mode_ = sm;
        is_updated_ = false;
    }

    using SystGraph::get_merge_mode;

    // --- variation management (TDualGraph-aware) ---

    /// \brief Add a TSystVariation with auto-assigned vid.
    /// \return the assigned vid.
    int add_var(const TSystVariation& sv) {
        int vid = auto_new_vid();
        add_var(vid, sv);
        return vid;
    }

    /// \brief Add a TSystVariation with a specific vid.
    void add_var(int vid, const TSystVariation& sv) {
        SystVariation base_sv;
        base_sv.name   = sv.name;
        base_sv.graph  = sv.graph;   // slice TDualGraph -> DualGraph
        base_sv.weight = sv.weight;
        base_sv.mode   = sv.mode;
        SystGraph::add_var(vid, base_sv);
        is_updated_ = false;
    }

    /// \brief Set (overwrite) a TSystVariation.
    void set_var(int vid, const TSystVariation& sv) {
        SystVariation base_sv;
        base_sv.name   = sv.name;
        base_sv.graph  = sv.graph;
        base_sv.weight = sv.weight;
        base_sv.mode   = sv.mode;
        SystGraph::set_var(vid, base_sv);
        is_updated_ = false;
    }

    // --- extract a subset of variations ---

    /// \brief Extract a subset of variations, preserving TDualGraph types.
    TSystGraph get_sub(const std::vector<int>& vidlist) const {
        TSystGraph sg;
        sg.set_def(TDualGraph(def_));
        for (int vid : vidlist) {
            auto it = var_.find(vid);
            if (it == var_.end()) continue;
            TSystVariation sv;
            sv.name   = it->second.name;
            sv.graph  = TDualGraph(it->second.graph);
            sv.weight = it->second.weight;
            sv.mode   = it->second.mode;
            sg.add_var(vid, sv);
        }
        return sg;
    }

    /// \brief Get the number of variations used in the last calc().
    int get_var_calc() const noexcept {
        return static_cast<int>(calc_result_.vids.size());
    }

    // =====================================================================
    //  calc — evaluate systematics and sync ROOT output graphs
    // =====================================================================

    /// \brief Evaluate systematics for specified vids and sync ROOT graphs.
    ///
    /// Applies x-axis limit filter, computes stat/syst/comb TGraphs,
    /// and sets is_updated_ to true. Points outside [xlim_low_, xlim_high_)
    /// are excluded from the ROOT output graphs.
    void calc(const std::vector<int>& vidlist) {
        int n = def_.get_n();
        if (n <= 0) return;

        calc_result_ = evaluate_systematics(vidlist);

        std::vector<double> vx, vex, vy, vey, syl, syh;
        vx.reserve(n); vex.reserve(n); vy.reserve(n);
        vey.reserve(n); syl.reserve(n); syh.reserve(n);

        for (int i = 0; i < n; ++i) {
            double xv = def_.get_point(i).get_px().get_value();
            if (xv < xlim_low_ || xv >= xlim_high_) continue;
            vx.push_back(xv);
            vex.push_back(width_);
            vy.push_back(def_.get_point(i).get_py().get_value());
            vey.push_back(def_.get_point(i).get_py().get_uncertainty());
            syl.push_back(calc_result_.total_lower[i]);
            syh.push_back(calc_result_.total_upper[i]);
        }

        int np = static_cast<int>(vx.size());

        std::vector<double> sylc(np), syhc(np);
        for (int i = 0; i < np; ++i) {
            sylc[i] = std::sqrt(vey[i] * vey[i] + syl[i] * syl[i]);
            syhc[i] = std::sqrt(vey[i] * vey[i] + syh[i] * syh[i]);
        }

        if (stat_.GetN() != np) stat_.Set(np);
        for (int i = 0; i < np; ++i) {
            stat_.GetX()[i]  = vx[i];
            stat_.GetY()[i]  = vy[i];
            stat_.GetEX()[i] = vex[i];
            stat_.GetEY()[i] = vey[i];
        }
        stat_.SetHistogram(nullptr);

        if (syst_.GetN() != np) syst_.Set(np);
        for (int i = 0; i < np; ++i) {
            syst_.GetX()[i]       = vx[i];
            syst_.GetY()[i]       = vy[i];
            syst_.GetEXlow()[i]   = vex[i];
            syst_.GetEXhigh()[i]  = vex[i];
            syst_.GetEYlow()[i]   = syl[i];
            syst_.GetEYhigh()[i]  = syh[i];
        }
        syst_.SetHistogram(nullptr);

        if (comb_.GetN() != np) comb_.Set(np);
        for (int i = 0; i < np; ++i) {
            comb_.GetX()[i]       = vx[i];
            comb_.GetY()[i]       = vy[i];
            comb_.GetEXlow()[i]   = vex[i];
            comb_.GetEXhigh()[i]  = vex[i];
            comb_.GetEYlow()[i]   = sylc[i];
            comb_.GetEYhigh()[i]  = syhc[i];
        }
        comb_.SetHistogram(nullptr);

        is_updated_ = true;
    }

    /// \brief Evaluate systematics for all variations.
    void calc() {
        calc(get_list());
    }

    /// \brief Trigger calc() only if ROOT graphs are not up to date.
    void calc_if_not_updated() {
        if (!is_updated_) calc();
    }

    // =====================================================================
    //  Value / LaTeX formatting
    // =====================================================================

    /// \brief LaTeX-formatted value with stat+syst uncertainties at index i.
    std::string str_latex(int i = 0,
                          const std::string& opt = "R") {
        if (!def_.check_index(i)) return std::string();
        calc_if_not_updated();
        ValueErrors ve(def_.get_point(i).get_py().get_value(),
                       def_.get_point(i).get_py().get_uncertainty(),
                       calc_result_.total_lower[i],
                       calc_result_.total_upper[i]);
        return ve.print_err3(opt);
    }

    /// \brief Get a ValueErrors object for bin i with stat and syst.
    ValueErrors valu_stat_syst(int i = 0) {
        if (!def_.check_index(i)) return ValueErrors();
        calc_if_not_updated();
        return ValueErrors(def_.get_point(i).get_py().get_value(),
                           def_.get_point(i).get_py().get_uncertainty(),
                           calc_result_.total_lower[i],
                           calc_result_.total_upper[i]);
    }

    // =====================================================================
    //  Shift x-values
    // =====================================================================

    /// \brief Shift ROOT graph x-values to match internal DualGraph x-values.
    void shift_x() {
        calc_if_not_updated();
        int ig = 0;
        int nd = def_.get_n();
        for (int i = 0; i < nd && ig < stat_.GetN(); ++i) {
            double xv = def_.get_point(i).get_px().get_value();
            if (xv < xlim_low_ || xv >= xlim_high_) continue;
            stat_.GetX()[ig] = xv;
            syst_.GetX()[ig] = xv;
            comb_.GetX()[ig] = xv;
            ++ig;
        }
    }

    /// \brief Shift ROOT graph x-values by an additive offset \p dx.
    void shift_x(double dx) {
        shift_x();
        for (int i = 0; i < stat_.GetN(); ++i) {
            stat_.GetX()[i] += dx;
            syst_.GetX()[i] += dx;
            comb_.GetX()[i] += dx;
        }
    }

    // =====================================================================
    //  Style / axis control
    // =====================================================================

    /// \brief Set marker style and color for all three ROOT output graphs.
    void set_style_color(Int_t style, Int_t color) {
        calc_if_not_updated();
        syst_.SetFillStyle(0);
        comb_.SetFillStyle(0);
        syst_.SetFillColor(color);
        comb_.SetFillColor(color);
        stat_.SetMarkerStyle(style);
        syst_.SetMarkerStyle(style);
        comb_.SetMarkerStyle(style);
        stat_.SetMarkerColor(color);
        syst_.SetMarkerColor(color);
        comb_.SetMarkerColor(color);
        stat_.SetLineColor(color);
        syst_.SetLineColor(color);
        comb_.SetLineColor(color);
    }

    /// \brief Set the axis range for all ROOT output graphs.
    void set_xaxis_range(double xl, double xh) {
        calc_if_not_updated();
        comb_.GetXaxis()->SetLimits(xl, xh);
        syst_.GetXaxis()->SetLimits(xl, xh);
        stat_.GetXaxis()->SetLimits(xl, xh);
    }

    /// \brief Set the x-axis limit filter applied during calc().
    ///
    /// Points outside [xl, xh) are excluded from the ROOT output graphs.
    void set_xaxis_limit(double xl, double xh) {
        xlim_low_  = xl;
        xlim_high_ = xh;
        is_updated_ = false;
    }

    /// \brief Clear the x-axis limit filter.
    void clear_xaxis_limit() {
        xlim_low_  = -std::numeric_limits<double>::infinity();
        xlim_high_ =  std::numeric_limits<double>::infinity();
        is_updated_ = false;
    }

    // =====================================================================
    //  Drawing
    // =====================================================================

    /// \brief Draw the statistical-uncertainty graph (stat_).
    void draw_stat(const char* opt = "PL same") {
        calc_if_not_updated();
        stat_.DrawClone(opt);
    }

    /// \brief Draw the combined-uncertainty graph (comb_).
    void draw_comb(const char* opt = "PL same") {
        calc_if_not_updated();
        comb_.DrawClone(opt);
    }

    /// \brief Draw systematic-uncertainty boxes with stat overlay.
    void draw_syst(const char* opt = "") {
        calc_if_not_updated();
        syst_.DrawClone(TString("P5") + opt + " same");
        stat_.DrawClone(TString("PL") + opt + " same");
    }

    /// \brief Draw systematic boxes and stat graph with separate options.
    void draw_syst(const char* optstat, const char* optsyst) {
        calc_if_not_updated();
        syst_.DrawClone(optsyst);
        stat_.DrawClone(optstat);
    }

    // =====================================================================
    //  MakePlot — systematic variation overlay plot
    // =====================================================================

    /// \brief Create an overlay plot showing all variation graphs.
    ///
    /// Each variation is drawn with a distinct marker and color,
    /// overlaid on the default graph. A legend is automatically generated.
    /// \param name base name for the output file
    /// \param hFrame a TH2D defining the plot frame (axis ranges and titles)
    /// \param path output directory path
    /// \param file_types output file formats (e.g. {"pdf", "png"})
    void make_plot(const char* name, TH2D* hFrame,
                   const char* path = "systplot",
                   std::vector<TString> file_types = {"pdf"}) {
        shift_x();

        int vcol[10] = {2, 3, 4, kOrange, 6, 7, 8, 9, 16, 28};
        int vmrk[10] = {20, 24, 21, 25, 22, 26, 23, 32, 29, 30};

        int nvar = static_cast<int>(var_.size());
        int ncln = nvar / 20 + 1;
        if (nvar == 0) ncln = 0;
        int wpxl = 700 + ncln * 175;
        int hpxl = 500;

        TCanvas* c = new TCanvas(
            TString("SystPlot") + name,
            TString("SystPlot") + name, wpxl, hpxl);
        gPad->SetLeftMargin(0.18 * 700 / wpxl);
        gPad->SetRightMargin((0.04 * 700 + ncln * 175) / wpxl);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.16);

        double xl = hFrame->GetXaxis()->GetXmin();
        double xh = hFrame->GetXaxis()->GetXmax();
        double yl = hFrame->GetYaxis()->GetXmin();
        double yh = hFrame->GetYaxis()->GetXmax();

        TH2D* hFrameTmp = new TH2D(
            TString("SystPlotFrame") + name,
            TString("SystPlotFrame") + name,
            100, xl, xh, 100, yl, yh);
        hFrameTmp->GetXaxis()->SetTitle(
            hFrame->GetXaxis()->GetTitle());
        hFrameTmp->GetYaxis()->SetTitle(
            hFrame->GetYaxis()->GetTitle());
        hFrameTmp->GetYaxis()->SetTitleOffset(1.0 / wpxl * 700);
        hFrameTmp->Draw();

        TLine* line = nullptr;
        if (yl < 0 && yh > 0) {
            line = new TLine(xl, 0, xh, 0);
            line->SetLineColor(kGray);
            line->Draw("same");
        }
        if (yl < 1 && yh > 1) {
            line = new TLine(xl, 1, xh, 1);
            line->SetLineColor(kGray);
            line->Draw("same");
        }

        TGraphErrors gTmp = stat_;
        gTmp.SetLineColor(kBlack);
        gTmp.SetMarkerColor(kBlack);
        gTmp.SetMarkerStyle(20);
        gTmp.Draw("PL same");

        double bw = std::fabs(xh - xl) / (gTmp.GetN() + 1.0);
        for (int i = 1; i < gTmp.GetN(); ++i) {
            double tmpbw = std::fabs(
                gTmp.GetX()[i] - gTmp.GetX()[i - 1]);
            if (i == 1) bw = tmpbw;
            else        bw = std::min(bw, tmpbw);
        }

        TLegend* leg = new TLegend(
            700.0 / wpxl, 0.01, 1 - 0.01 * hpxl / wpxl, 0.99);
        leg->SetFillStyle(0);
        leg->SetBorderSize(1);
        leg->SetNColumns(ncln);
        leg->SetTextSize(0.04);
        leg->AddEntry(&gTmp, "default", "lp");

        int nsx = nvar + 2;
        int isx = 1;
        int iv  = 0;
        double sw = 0.06 * nsx;
        if (sw < 0.6) sw = 0.6;
        if (sw > 1.0) sw = 1.0;

        std::vector<TGraphErrors> gVarStore;
        gVarStore.reserve(nvar);

        for (auto& kv : var_) {
            const SystVariation& sv = kv.second;
            const DualGraph& vg = sv.graph;
            int vn = vg.get_n();
            std::vector<double> vx(vn), vy(vn);
            for (int j = 0; j < vn; ++j) {
                vx[j] = vg.get_point(j).get_px().get_value()
                      + sw * bw * isx / nsx;
                vy[j] = vg.get_point(j).get_py().get_value();
            }
            gVarStore.emplace_back(vn, vx.data(), vy.data());
            TGraphErrors& gVar = gVarStore.back();
            isx++;

            int col = vcol[iv % 10];
            int imrk = iv / 10;
            while (imrk >= 10) imrk /= 10;
            int mrk = vmrk[imrk];
            gVar.SetMarkerStyle(mrk);
            gVar.SetMarkerColor(col);
            gVar.SetLineColor(col);
            gVar.Draw("P same");
            leg->AddEntry(&gVar, sv.name.c_str(), "lp");
            iv++;
        }

        leg->Draw("same");
        gPad->SetTicks();

        for (const auto& ft : file_types)
            c->SaveAs(TString(path) + "/Syst" + name + "." + ft);

        delete hFrameTmp;
        delete leg;
        if (line) delete line;
        delete c;
    }

    /// \brief Create an overlay plot with auto-computed frame.
    ///
    /// The frame is automatically generated from the stat_ graph's range.
    void make_plot(const char* name = "",
                   const char* path = "systplot",
                   std::vector<TString> file_types = {"pdf"}) {
        TString nm(name);
        if (nm == "") nm = stat_.GetName();
        int n = stat_.GetN();
        double xl = TMath::MinElement(n, stat_.GetX());
        double xh = TMath::MaxElement(n, stat_.GetX());
        double yl = TMath::MinElement(n, stat_.GetY());
        double yh = TMath::MaxElement(n, stat_.GetY());
        double xw = xh - xl;
        double yw = yh - yl;
        if (n == 1) {
            xw = 10;
            yw = 40 * stat_.GetEY()[0];
        }
        xl = xl - 0.05 * xw;
        xh = xh + 0.05 * xw;
        yl = yl - 0.05 * yw;
        yh = yh + 0.15 * yw;

        TH2D* hFrame = new TH2D(
            TString("FramePlot") + nm,
            TString("FramePlot") + nm,
            100, xl, xh, 100, yl, yh);
        hFrame->GetXaxis()->SetTitle(stat_.GetXaxis()->GetTitle());
        hFrame->GetYaxis()->SetTitle(stat_.GetYaxis()->GetTitle());
        make_plot(nm, hFrame, path, file_types);
        delete hFrame;
    }

    // =====================================================================
    //  MakeHist — stacked histogram of systematic^2 breakdown
    // =====================================================================

    /// \brief Create a stacked histogram showing the systematic^2 breakdown.
    ///
    /// Each variation's contribution \f$ w \cdot \sigma_{\mathrm{syst}}^2 \f$
    /// is drawn as a stacked box. Upper uncertainties stack upward (positive y);
    /// lower uncertainties stack downward (negative y).
    ///
    /// \param name base name for the output file
    /// \param hFrame a TH2D defining the plot frame
    /// \param is_percent if true, normalize to percent of total systematic^2
    /// \param path output directory path
    /// \param file_types output file formats
    void make_hist(const char* name, TH2D* hFrame,
                   bool is_percent = false,
                   const char* path = "systplot",
                   std::vector<TString> file_types = {"pdf"}) {
        shift_x();

        TString nm(name);
        if (is_percent) nm = nm + "Percent";

        int vcol[10] = {2, 3, 4, kOrange, 6, 7, 8, 9, 16, 28};
        int vfls[10] = {1001, 3001, 3002, 3004, 3244,
                        3003, 3344, 3345, 3444, 3354};

        int nvar = static_cast<int>(var_.size());
        int ncln = nvar / 20 + 1;
        int wpxl = 700 + ncln * 175;
        int hpxl = 500;

        TCanvas* c = new TCanvas(
            TString("SystHist") + nm,
            TString("SystHist") + nm, wpxl, hpxl);
        gPad->SetLeftMargin(0.18 * 700 / wpxl);
        gPad->SetRightMargin((0.04 * 700 + ncln * 175) / wpxl);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.16);

        double xl = hFrame->GetXaxis()->GetXmin();
        double xh = hFrame->GetXaxis()->GetXmax();
        TString xtitle = hFrame->GetXaxis()->GetTitle();
        TString ytitle = TString(hFrame->GetYaxis()->GetTitle())
                       + ": syst^{2}";
        if (is_percent) ytitle = ytitle + " (%)";

        TGraphErrors gTmp = stat_;
        gTmp.SetLineColor(kBlack);
        gTmp.SetMarkerColor(kBlack);
        gTmp.SetMarkerStyle(20);

        int nBin = gTmp.GetN();
        if (nBin <= 0) { delete c; return; }

        double bw = std::fabs(xh - xl) / (nBin + 1.0);
        for (int i = 1; i < nBin; ++i) {
            double tmpbw = std::fabs(
                gTmp.GetX()[i] - gTmp.GetX()[i - 1]);
            if (i == 1) bw = tmpbw;
            else        bw = std::min(bw, tmpbw);
        }

        bool is_good_xedge = true;
        std::vector<double> xcenter(nBin);
        for (int i = 0; i < nBin; ++i)
            xcenter[i] = gTmp.GetX()[i];
        std::sort(xcenter.begin(), xcenter.end());

        std::vector<double> xedge(nBin + 1);
        xedge[0] = xl;
        for (int i = 0; i < nBin; ++i) {
            xedge[i + 1] = xedge[i]
                         + (xcenter[i] - xedge[i]) * 2.0;
            if (xedge[i + 1] <= xedge[i]) is_good_xedge = false;
        }
        if (!is_good_xedge) {
            xedge[0] = xcenter[0] - 0.5 * bw;
            for (int i = 0; i < nBin; ++i)
                xedge[i + 1] = xcenter[i] + 0.5 * bw;
        }

        std::vector<int> xbinedge(nBin);
        for (int i = 0; i < nBin; ++i) {
            double xc = gTmp.GetX()[i];
            for (int j = 0; j < nBin; ++j) {
                if (xedge[j] <= xc && xedge[j + 1] >= xc) {
                    xbinedge[i] = j;
                    break;
                }
            }
        }

        TLegend* leg = new TLegend(
            700.0 / wpxl, 0.01,
            1 - 0.01 * hpxl / wpxl, 0.99);
        leg->SetFillStyle(0);
        leg->SetBorderSize(1);
        leg->SetNColumns(ncln);
        leg->SetTextSize(0.04);

        int nsy = static_cast<int>(calc_result_.vids.size());
        if (nsy == 0) { delete c; delete leg; return; }

        std::vector<std::vector<TBox>> syBoxUpper(nsy);
        std::vector<std::vector<TBox>> syBoxLower(nsy);
        double syMaxUpper = 0.0;
        double syMaxLower = 0.0;

        for (int j = 0; j < nBin; ++j) {
            double s2_total_upper = calc_result_.total_upper[j]
                                  * calc_result_.total_upper[j];
            double s2_total_lower = calc_result_.total_lower[j]
                                  * calc_result_.total_lower[j];

            double syBinUpper = 0.0;
            double syBinLower = 0.0;

            for (int i = nsy - 1; i >= 0; --i) {
                int vid = calc_result_.vids[i];
                const SystVariation& sv = var_.find(vid)->second;
                double w = sv.weight;
                double syst = calc_result_.var_syst[i][j];

                SystMergeMode sm = (merge_mode_ != SystMergeMode::kNull)
                                   ? merge_mode_ : sv.mode;
                bool one_sided = ::is_one_sided(sm);

                double contrib = w * syst * syst;

                bool to_upper = !one_sided || syst > 0.0;
                bool to_lower = !one_sided || syst < 0.0;

                int col = vcol[i % 10];
                int ifls = i / 10;
                while (ifls >= 10) ifls /= 10;
                int fls = vfls[ifls];

                if (to_upper && contrib > 0.0) {
                    double upval = contrib;
                    if (is_percent && s2_total_upper != 0.0)
                        upval = 100.0 * upval / s2_total_upper;

                    TBox box(xedge[xbinedge[j]], syBinUpper,
                             xedge[xbinedge[j] + 1],
                             syBinUpper + upval);
                    box.SetLineColor(col);
                    box.SetFillStyle(fls);
                    box.SetFillColor(col);
                    syBoxUpper[i].push_back(box);
                    syBinUpper += upval;
                }

                if (to_lower && contrib > 0.0) {
                    double lowval = contrib;
                    if (is_percent && s2_total_lower != 0.0)
                        lowval = 100.0 * lowval / s2_total_lower;

                    TBox box(xedge[xbinedge[j]],
                             -(syBinLower + lowval),
                             xedge[xbinedge[j] + 1],
                             -syBinLower);
                    box.SetLineColor(col);
                    box.SetFillStyle(fls);
                    box.SetFillColor(col);
                    syBoxLower[i].push_back(box);
                    syBinLower += lowval;
                }
            }
            if (syBinUpper > syMaxUpper) syMaxUpper = syBinUpper;
            if (syBinLower > syMaxLower) syMaxLower = syBinLower;
        }

        double yMin = (syMaxLower > 0.0) ? -syMaxLower : 0.0;
        double yMax = (syMaxUpper > 0.0) ?  syMaxUpper : 0.0;
        if (yMin == 0.0 && yMax == 0.0) yMax = 1.0;

        TH2D* hFrameTmp = new TH2D(
            TString("SystHistFrame") + nm,
            TString("SystHistFrame") + nm,
            100, xl, xh, 100, yMin, yMax);
        hFrameTmp->GetXaxis()->SetTitle(xtitle);
        hFrameTmp->GetYaxis()->SetTitle(ytitle);
        hFrameTmp->GetYaxis()->SetTitleOffset(1.0 / wpxl * 700);
        hFrameTmp->Draw();

        for (int i = 0; i < nsy; ++i) {
            int vid = calc_result_.vids[i];
            const SystVariation& sv = var_.find(vid)->second;
            if (!syBoxUpper[i].empty())
                leg->AddEntry(&syBoxUpper[i][0], sv.name.c_str(), "f");
            else if (!syBoxLower[i].empty())
                leg->AddEntry(&syBoxLower[i][0], sv.name.c_str(), "f");
        }

        for (int i = 0; i < nsy; ++i) {
            for (auto& box : syBoxUpper[i]) box.Draw();
            for (auto& box : syBoxLower[i]) box.Draw();
        }

        TLine *tmp0 = new TLine(xl,0, xh,0);
        tmp0->Draw("same");

        leg->Draw("same");
        gPad->SetTicks();

        for (const auto& ft : file_types)
            c->SaveAs(TString(path) + "/Syst" + nm
                      + "Hist." + ft);

        delete hFrameTmp;
        delete leg;
        delete c;
    }

    /// \brief Create a stacked histogram with auto-computed frame.
    void make_hist(const char* name = "",
                   bool is_percent = false,
                   const char* path = "systplot",
                   std::vector<TString> file_types = {"pdf"}) {
        TString nm(name);
        if (nm == "") nm = stat_.GetName();
        int n = stat_.GetN();
        double xl = TMath::MinElement(n, stat_.GetX());
        double xh = TMath::MaxElement(n, stat_.GetX());
        double yl = TMath::MinElement(n, stat_.GetY());
        double yh = TMath::MaxElement(n, stat_.GetY());
        double xw = xh - xl;
        double yw = yh - yl;
        if (n == 1) {
            xw = 10;
            yw = 40 * stat_.GetEY()[0];
        }
        xl = xl - 0.05 * xw;
        xh = xh + 0.05 * xw;
        yl = yl - 0.05 * yw;
        yh = yh + 0.15 * yw;

        TH2D* hFrame = new TH2D(
            TString("FramePlot") + nm,
            TString("FramePlot") + nm,
            100, xl, xh, 100, yl, yh);
        hFrame->GetXaxis()->SetTitle(
            stat_.GetXaxis()->GetTitle());
        hFrame->GetYaxis()->SetTitle(
            stat_.GetYaxis()->GetTitle());
        make_hist(nm, hFrame, is_percent, path, file_types);
        delete hFrame;
    }

    // =====================================================================
    //  ROOT TFile I/O
    // =====================================================================

    /// \brief Sync DualMultiv and SystGraph ID counters after reading from file.
    ///
    /// Prevents ID collisions when new DualMultiv or SystGraph objects are
    /// created after reading this object from a TFile.
    /// Call this after: `auto* sg = f->Get<TSystGraph>("name");`
    void sync_id_counter() const {
        int max_dm_id = 0;
        int max_vid   = 0;
        auto scan = [&](const DualGraph& dg) {
            for (int i = 0; i < dg.get_n(); ++i) {
                for (int id : dg.get_point(i).get_py().get_ids()) {
                    if (id > max_dm_id) max_dm_id = id;
                }
            }
        };
        scan(def_);
        for (const auto& kv : var_) {
            scan(kv.second.graph);
            if (kv.first > max_vid) max_vid = kv.first;
        }
        DualMultiv::set_id_counter(max_dm_id);
        if (max_vid >= next_vid_) next_vid_ = max_vid;
    }

    /// \brief Write individual ROOT output graphs for quick plotting.
    /// \param name base name for the output objects
    /// \param is_stat if true, write the stat_ graph
    /// \param is_asym if true, write the syst_ (asymmetric) graph
    void write_graph(const char* name = "",
                     bool is_stat = true,
                     bool is_asym = true) {
        shift_x();
        TString nm(name);
        if (is_stat) stat_.Write(nm + "Stat");
        if (is_asym) syst_.Write(nm + "Asym");
    }

    /// \brief Write default and variation graphs as plain TGraphErrors.
    void write_def_var(const char* name = "") {
        TString nm(name);
        {
            TDualGraph tdg(def_);
            tdg.calc();
            tdg.get_graph().Write(nm + "Def");
        }
        for (auto& kv : var_) {
            TString snm = TString::Format(
                "Syst%d", kv.first);
            TDualGraph tdg(kv.second.graph);
            tdg.calc();
            tdg.get_graph().Write(nm + snm);
        }
    }

    // =====================================================================
    //  Bin-level operations (override SystGraph for TDualGraph output)
    // =====================================================================

    /// \brief Weighted average over bins, preserving TDualGraph type.
    TSystGraph average_bin(const std::vector<int>& range) const {
        TSystGraph sg(*this);
        sg.set_def(TDualGraph(def_.average_bin(range)));
        for (auto& kv : sg.var_)
            kv.second.graph = kv.second.graph.average_bin(range);
        return sg;
    }

    /// \brief Equal-weight average over bins, preserving TDualGraph type.
    TSystGraph average_bin_simple(const std::vector<int>& range) const {
        TSystGraph sg(*this);
        sg.set_def(TDualGraph(def_.average_bin_simple(range)));
        for (auto& kv : sg.var_)
            kv.second.graph = kv.second.graph.average_bin_simple(range);
        return sg;
    }

    /// \brief Select a subset of bins, preserving TDualGraph type.
    TSystGraph select_bin(const std::vector<int>& range) const {
        TSystGraph sg(*this);
        sg.set_def(TDualGraph(def_.select_bin(range)));
        for (auto& kv : sg.var_)
            kv.second.graph = kv.second.graph.select_bin(range);
        return sg;
    }

    /// \brief Select a contiguous range of bins [bl, bh].
    TSystGraph select_bin(int bl, int bh) const {
        std::vector<int> range;
        for (int i = bl; i <= bh; ++i) range.push_back(i);
        return select_bin(range);
    }

    /// \brief Select a single bin.
    TSystGraph select_bin(int bl) const {
        return select_bin(bl, bl);
    }

    ClassDefOverride(TSystGraph, 1)  ///< ROOT streamer info (version 1)
};
#endif  // TSystGraph_H
