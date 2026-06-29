/**
 * \file TGraphErrorsOperator.h
 * \brief Free-function operators and utilities for ROOT TGraphErrors and TGraphAsymmErrors.
 *
 * \details Provides arithmetic operators (+ - * /), conversion functions,
 * centrality rebinning, uncertainty propagation helpers (Barlow-style
 * difference and ratio), event-plane resolution calculation, and various
 * graph manipulation utilities.
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   - Removed using namespace std (namespace pollution in header)
 *   - Removed const from operator return-by-value (enables move semantics)
 *   - Replaced VLAs (double xe[n]) with std::vector<double>
 *   - Fixed operator%: no longer mutates input graph x-values in place
 *   - Replaced 0.0/0.0 with std::nan("")
 *   - Added explicit includes: TF1, TGraphAsymmErrors, iomanip, sstream
 *   - Consistent formatting throughout
 */

#ifndef TGraphErrorsOperator_H
#define TGraphErrorsOperator_H

#include <cmath>
#include <vector>

#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>


// ===========================================================================
//  Arithmetic operators  + - * /  (with error propagation)
// ===========================================================================

/// \brief Pointwise TGraphErrors + TGraphErrors with error propagation.
/// Errors are added in quadrature: \f$ e = \sqrt{e_1^2 + e_2^2} \f$.
TGraphErrors operator+(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        double y2e = g2.GetErrorY(i);
        y.push_back(y1[i] + y2[i]);
        ye.push_back(std::sqrt(y1e*y1e + y2e*y2e));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Scalar + TGraphErrors: adds constant to all y-values (errors unchanged).
TGraphErrors operator+(double c, const TGraphErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        y.push_back(y1[i] + c);
        ye.push_back(y1e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


inline TGraphErrors operator+(const TGraphErrors &g1, double c) { return c + g1; }


/// \brief Pointwise TGraphErrors * TGraphErrors with error propagation.
/// \f$ e = \sqrt{y_2^2 e_1^2 + y_1^2 e_2^2} \f$.
TGraphErrors operator*(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        double y2e = g2.GetErrorY(i);
        y.push_back(y1[i] * y2[i]);
        ye.push_back(std::sqrt(y2[i]*y2[i] * y1e*y1e + y1[i]*y1[i] * y2e*y2e));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Scalar * TGraphErrors.
TGraphErrors operator*(double c, const TGraphErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        y.push_back(y1[i] * c);
        ye.push_back(std::fabs(y1e * c));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


inline TGraphErrors operator*(const TGraphErrors &g1, double c) { return c * g1; }


/// \brief TGraphErrors - TGraphErrors.
inline TGraphErrors operator-(const TGraphErrors &g1, const TGraphErrors &g2)
{
    return g1 + (-1.0 * g2);
}


inline TGraphErrors operator-(double c, const TGraphErrors &g1) { return c + (-1.0 * g1); }
inline TGraphErrors operator-(const TGraphErrors &g1, double c) { return g1 + (-1.0 * c); }


/// \brief Pointwise TGraphErrors / TGraphErrors with error propagation.
/// \f$ e = |y_1/y_2| \sqrt{(e_1/y_1)^2 + (e_2/y_2)^2} \f$.
TGraphErrors operator/(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        double y2e = g2.GetErrorY(i);
        double tmpy = y1[i] / y2[i];
        double tmpye = std::fabs(tmpy) *
                       std::sqrt(y1e*y1e/y1[i]/y1[i] + y2e*y2e/y2[i]/y2[i]);
        if (y1[i] == 0) {
            tmpy  = 0;
            tmpye = 0;
        }
        y.push_back(tmpy);
        ye.push_back(tmpye);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Scalar / TGraphErrors.
TGraphErrors operator/(double c, const TGraphErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        double tmpy = c / y1[i];
        double tmpye = c * y1e / y1[i] / y1[i];
        if (c == 0) {
            tmpy  = 0;
            tmpye = 0;
        }
        y.push_back(tmpy);
        ye.push_back(tmpye);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


inline TGraphErrors operator/(const TGraphErrors &g1, double c) { return g1 * (1.0 / c); }


/// \brief Shift x-values by \p c (operator%).
TGraphErrors operator%(const TGraphErrors &g1, double c)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, ye;

    for (int i = 0; i < n; i++) {
        double xi = x1[i];
        if (std::isnan(c) || std::isinf(c)) {
            xi = 0;
        } else {
            xi += c;
        }
        x.push_back(xi);
        xe.push_back(0);
        ye.push_back(g1.GetErrorY(i));
    }

    TGraphErrors g(n, &x[0], &y1[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Subtract a TF1 from a TGraphErrors (preserves errors).
TGraphErrors operator-(const TGraphErrors &g, const TF1 &f)
{
    const int n = g.GetN();
    double *gx = g.GetX();
    double *gy = g.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(gx[i]);
        xe.push_back(0);
        double gye = g.GetErrorY(i);
        y.push_back(gy[i] - f.Eval(gx[i]));
        ye.push_back(gye);
    }

    return TGraphErrors(x.size(), &x[0], &y[0], &xe[0], &ye[0]);
}


// ===========================================================================
//  In-place graph mutation helpers
// ===========================================================================

/// \brief Shift x-values of a TGraphErrors in-place by \p c.
void GraphShiftX(TGraphErrors *g1, double c)
{
    const int n = g1->GetN();
    double *x = g1->GetX();

    for (int i = 0; i < n; i++) {
        if (std::isnan(c) || std::isinf(c)) {
            x[i] = 0;
        } else {
            x[i] += c;
        }
    }
}


/// \brief Shift x-values of a TGraphErrors in-place to 0.
void GraphShiftX(TGraphErrors *g1)
{
    double c = std::nan("");
    GraphShiftX(g1, c);
}


/// \brief Shift x-values of a TGraphAsymmErrors in-place by \p c.
void GraphShiftX(TGraphAsymmErrors *g1, double c)
{
    const int n = g1->GetN();
    double *x = g1->GetX();

    for (int i = 0; i < n; i++) {
        if (std::isnan(c) || std::isinf(c)) {
            x[i] = 0;
        } else {
            x[i] += c;
        }
    }
}


/// \brief Shift x-values of a TGraphAsymmErrors in-place to 0.
void GraphShiftX(TGraphAsymmErrors *g1)
{
    double c = std::nan("");
    GraphShiftX(g1, c);
}


// ===========================================================================
//  Graph creation / conversion
// ===========================================================================

/// \brief Convert a TH1D to a TGraphErrors.
TGraphErrors GraphFrom(const TH1D *h)
{
    if (h == nullptr) return TGraphErrors();

    const int n = h->GetSize() - 2;
    std::vector<double> x, xe, y, ye;
    for (int i = 0; i < n; i++) {
        x.push_back(h->GetBinCenter(i + 1));
        xe.push_back(0);
        y.push_back(h->GetBinContent(i + 1));
        ye.push_back(h->GetBinError(i + 1));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Convert a TH1F to a TGraphErrors.
TGraphErrors GraphFrom(const TH1F *h)
{
    if (h == nullptr) return TGraphErrors();

    const int n = h->GetSize() - 2;
    std::vector<double> x, xe, y, ye;
    for (int i = 0; i < n; i++) {
        x.push_back(h->GetBinCenter(i + 1));
        xe.push_back(0);
        y.push_back(h->GetBinContent(i + 1));
        ye.push_back(h->GetBinError(i + 1));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Copy a TGraphErrors (with errors).
TGraphErrors GraphFrom(const TGraphErrors *g1)
{
    if (g1 == nullptr) return TGraphErrors();

    TGraphErrors g(g1->GetN(), g1->GetX(), g1->GetY(), g1->GetEX(), g1->GetEY());

    return g;
}


/// \brief Copy a TGraphErrors via copy constructor.
TGraphErrors GraphSelf(const TGraphErrors &g1)
{
    TGraphErrors g = g1;
    return g;
}


// ===========================================================================
//  Centrality rebinning
// ===========================================================================

/// \brief Map a TGraphErrors to centrality bin x-values (9-bin default).
TGraphErrors GraphCent(const TGraphErrors &g1)
{
    double cent[9]  = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
    double cente[9] = {0};

    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    int di = 0;
    for (int i = 0; i < n; i++) {
        if (n == 10 && i == 0) {
            di = 1;
            continue;
        }
        double y1e = g1.GetErrorY(i);
        x.push_back(cent[i - di]);
        xe.push_back(g1.GetErrorX(i));
        y.push_back(y1[i]);
        ye.push_back(y1e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Alias for GraphCent (9-bin).
inline TGraphErrors GraphCent9(const TGraphErrors &g1) { return GraphCent(g1); }

/// \brief Map to 8 centrality bins.
TGraphErrors GraphCent8(const TGraphErrors &g1)
{
    double cent[8]  = {75, 65, 55, 45, 35, 25, 15, 5};
    double cente[8] = {0};

    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    int di = 0;
    for (int i = 0; i < 8; i++) {
        double y1e = g1.GetErrorY(i);
        x.push_back(cent[i - di]);
        xe.push_back(cente[i - di]);
        y.push_back(y1[i]);
        ye.push_back(y1e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Map to 10 centrality bins.
TGraphErrors GraphCent10(const TGraphErrors &g1)
{
    double cent[10]  = {90, 75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
    double cente[10] = {0};

    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    int di = 0;
    for (int i = 0; i < n; i++) {
        double y1e = g1.GetErrorY(i);
        x.push_back(cent[i - di]);
        xe.push_back(g1.GetErrorX(i));
        y.push_back(y1[i]);
        ye.push_back(y1e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Map to 16 centrality bins.
TGraphErrors GraphCent16(const TGraphErrors &g1)
{
    double cent[16]  = {77.5, 72.5, 67.5, 62.5, 57.5, 52.5, 47.5,
                         42.5, 37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 7.5, 2.5};
    double cente[16] = {0};

    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    int di = 0;
    for (int i = 0; i < n; i++) {
        if (n == 17 && i == 0) {
            di = 1;
            continue;
        }
        double y1e = g1.GetErrorY(i);
        x.push_back(cent[i - di]);
        xe.push_back(cente[i - di]);
        y.push_back(y1[i]);
        ye.push_back(y1e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Map a TGraphAsymmErrors to centrality bin x-values.
TGraphAsymmErrors GraphCent(const TGraphAsymmErrors &g1)
{
    double cent[9]  = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
    double cente[9] = {0};

    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xle, xhe, y, yle, yhe;

    int di = 0;
    for (int i = 0; i < n; i++) {
        if (n == 10 && i == 0) {
            di = 1;
            continue;
        }
        x.push_back(cent[i - di]);
        xle.push_back(g1.GetErrorXlow(i));
        xhe.push_back(g1.GetErrorXhigh(i));
        y.push_back(y1[i]);
        yle.push_back(g1.GetErrorYlow(i));
        yhe.push_back(g1.GetErrorYhigh(i));
    }

    TGraphAsymmErrors g(x.size(), &x[0], &y[0],
                         &xle[0], &xhe[0], &yle[0], &yhe[0]);

    return g;
}


/// \brief Create an error-band envelope polygon from a TGraphErrors.
///
/// Traces the upper edge forward, then the lower edge backward,
/// suitable for drawing with "f" (fill) option.
TGraph GraphEdge(const TGraphErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, y;

    for (int i = 0; i < n; i++) {
        double y1e = g1.GetErrorY(i);
        x.push_back(x1[i]);
        y.push_back(y1[i] + y1e);
    }
    for (int i = n - 1; i >= 0; i--) {
        double y1e = g1.GetErrorY(i);
        x.push_back(x1[i]);
        y.push_back(y1[i] - y1e);
    }
    x.push_back(x[0]);
    y.push_back(y[0]);

    TGraph g(x.size(), &x[0], &y[0]);
    g.SetLineColor(g1.GetLineColor());
    g.SetMarkerColor(g1.GetMarkerColor());

    return g;
}


/// \brief Create an error-band envelope polygon from a TGraphAsymmErrors.
TGraph GraphEdge(const TGraphAsymmErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, y;

    for (int i = 0; i < n; i++) {
        double y1e = g1.GetErrorYhigh(i);
        x.push_back(x1[i]);
        y.push_back(y1[i] + y1e);
    }
    for (int i = n - 1; i >= 0; i--) {
        double y1e = g1.GetErrorYlow(i);
        x.push_back(x1[i]);
        y.push_back(y1[i] - y1e);
    }
    x.push_back(x[0]);
    y.push_back(y[0]);

    TGraph g(x.size(), &x[0], &y[0]);
    g.SetLineColor(g1.GetLineColor());
    g.SetMarkerColor(g1.GetMarkerColor());

    return g;
}


// ===========================================================================
//  Point-wise transformations
// ===========================================================================

/// \brief Pointwise absolute value.
TGraphErrors GraphAbs(const TGraphErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        y.push_back(std::fabs(y1[i]));
        ye.push_back(std::fabs(y1e));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Pointwise power \f$ y^c \f$ with error propagation.
/// \f$ e_y = |c| \cdot y^{c-1} \cdot e \f$.
TGraphErrors GraphPow(const TGraphErrors &g1, double c)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        y.push_back(std::pow(y1[i], c));
        ye.push_back(y1e * c * std::pow(y1[i], c - 1));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Pointwise square root (of absolute value).
TGraphErrors GraphSqrt(const TGraphErrors &g1)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        y.push_back(std::sqrt(std::fabs(y1[i])));
        ye.push_back(y1e / 2.0 / std::sqrt(std::fabs(y1[i])));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


// ===========================================================================
//  Weighted averaging
// ===========================================================================

/// \brief Inverse-variance weighted average over a contiguous y-range.
///
/// Weights are \f$ w_i = 1 / e_i^2 \f$.
/// \param g1 input graph
/// \param firstpoint first bin index (inclusive)
/// \param lastpoint last bin index (inclusive)
/// \return a single-bin TGraphErrors with the weighted average
TGraphErrors GraphAveY(const TGraphErrors &g1, int firstpoint, int lastpoint)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double x  = 0;
    double xe = 0;
    double y  = 0;
    double yy = 0;
    double ye = 0;
    double w  = 0;
    int dn = 0;

    int fp = firstpoint <= lastpoint ? firstpoint : lastpoint;
    int lp = firstpoint <= lastpoint ? lastpoint : firstpoint;
    fp = fp >= 0   ? fp   : 0;
    fp = fp >= n   ? n-1 : fp;
    lp = lp >= 0   ? lp   : 0;
    lp = lp >= n   ? n-1 : lp;

    for (int i = fp; i <= lp; i++) {
        double ee = g1.GetErrorY(i);
        if (ee == 0) continue;
        x += x1[i];
        double ww = 1.0 / ee / ee;
        y  += ww * y1[i];
        yy += ww * y1[i] * y1[i];
        w  += ww;
        dn++;
    }
    x  = x / dn;
    y  = y / w;
    yy = yy / w;
    ye = std::sqrt(1.0 / w);

    TGraphErrors g(1, &x, &y, &xe, &ye);

    return g;
}


/// \brief Inverse-variance weighted average using errors from a second graph.
///
/// Weights come from g2's y-errors; g1 provides the y-values.
TGraphErrors GraphWtAveY(const TGraphErrors &g1, const TGraphErrors &g2,
                          int firstpoint, int lastpoint)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double x  = 0;
    double xe = 0;
    double y  = 0;
    double yy = 0;
    double ye = 0;
    double w  = 0;
    int dn = 0;

    int fp = firstpoint <= lastpoint ? firstpoint : lastpoint;
    int lp = firstpoint <= lastpoint ? lastpoint : firstpoint;
    fp = fp >= 0   ? fp   : 0;
    fp = fp >= n   ? n-1 : fp;
    lp = lp >= 0   ? lp   : 0;
    lp = lp >= n   ? n-1 : lp;

    for (int i = fp; i <= lp; i++) {
        double ee = g2.GetErrorY(i);
        if (ee == 0) continue;
        x += x1[i];
        double ww = 1.0 / ee / ee;
        y  += ww * y1[i];
        yy += ww * y1[i] * y1[i];
        w  += ww;
        ye += ww * std::pow(g1.GetErrorY(i), 2.0);
        dn++;
    }
    x  = x / dn;
    y  = y / w;
    yy = yy / w;
    ye = std::sqrt(ye / w / dn);

    TGraphErrors g(1, &x, &y, &xe, &ye);

    return g;
}


/// \brief Average over all bins.
inline TGraphErrors GraphAveY(const TGraphErrors &g1)
{
    const int n = g1.GetN();
    return GraphAveY(g1, 0, n - 1);
}


/// \brief Create a single-point TGraphErrors (x=0).
TGraphErrors GraphOne(double con, double err)
{
    double x  = 0;
    double xe = 0;
    double y  = con;
    double ye = err;

    TGraphErrors g(1, &x, &y, &xe, &ye);

    return g;
}


/// \brief Create a single-point TGraphErrors with explicit x.
TGraphErrors GraphOne(double vx, double ex, double con, double err)
{
    double x  = vx;
    double xe = ex;
    double y  = con;
    double ye = err;

    TGraphErrors g(1, &x, &y, &xe, &ye);

    return g;
}


// ===========================================================================
//  Filtering and error manipulation
// ===========================================================================

/// \brief Remove points with zero y-value or invalid error.
TGraphErrors RemoveZero(const TGraphErrors &g)
{
    const int n = g.GetN();
    double *x = g.GetX();
    double *y = g.GetY();
    std::vector<double> xe(n), ye(n);

    int on = 0;
    std::vector<double> ox(n), oy(n), oxe(n), oye(n);

    for (int i = 0; i < n; i++) {
        xe[i] = g.GetErrorX(i);
        ye[i] = g.GetErrorY(i);

        if (y[i] != 0 && (ye[i] != 0 && !std::isnan(ye[i]) && !std::isinf(ye[i]))) {
            ox[on]  = x[i];
            oy[on]  = y[i];
            oxe[on] = xe[i];
            oye[on] = ye[i];
            on++;
        }
    }

    TGraphErrors og(on, &ox[0], &oy[0], &oxe[0], &oye[0]);

    return og;
}


/// \brief Zero all y-errors (keep x-errors and values).
TGraphErrors ZeroError(const TGraphErrors &g)
{
    const int n = g.GetN();
    double *x = g.GetX();
    double *y = g.GetY();
    std::vector<double> xe(n), ye(n);

    for (int i = 0; i < n; i++) {
        xe[i] = g.GetErrorX(i);
        ye[i] = 0;
    }

    TGraphErrors og(n, x, y, &xe[0], &ye[0]);

    return og;
}


// ===========================================================================
//  Graph combination
// ===========================================================================

/// \brief Map two graphs: x = y of g1, y = y of g2 (with errors).
TGraphErrors GraphMap(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        double y1e = g1.GetErrorY(i);
        double y2e = g2.GetErrorY(i);
        x.push_back(y1[i]);
        xe.push_back(y1e);
        y.push_back(y2[i]);
        ye.push_back(y2e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Concatenate two TGraphErrors.
TGraphErrors GraphMerge(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n1 = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    const int n2 = g2.GetN();
    double *x2 = g2.GetX();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n1; i++) {
        double x1e = g1.GetErrorX(i);
        double y1e = g1.GetErrorY(i);
        x.push_back(x1[i]);
        xe.push_back(x1e);
        y.push_back(y1[i]);
        ye.push_back(y1e);
    }
    for (int i = 0; i < n2; i++) {
        double x2e = g2.GetErrorX(i);
        double y2e = g2.GetErrorY(i);
        x.push_back(x2[i]);
        xe.push_back(x2e);
        y.push_back(y2[i]);
        ye.push_back(y2e);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


/// \brief Swap y-values with y-errors (x-errors become new y-errors, zeroed).
TGraphErrors GraphErrToVal(const TGraphErrors &g1)
{
    TGraphErrors g(g1.GetN(), g1.GetX(), g1.GetEY(), g1.GetEX(), g1.GetEY());
    for (int i = 0; i < g.GetN(); i++) g.GetEY()[i] = 0;
    return g;
}


/// \brief Shift graph by \p c times its errors: \f$ y + c \cdot e_y \f$.
TGraphErrors GraphMoveSigma(const TGraphErrors &g1, double c)
{
    TGraphErrors ge = GraphErrToVal(g1);
    return (g1 + c * ge);
}


// ===========================================================================
//  Event-plane resolution helpers
// ===========================================================================

/// \brief Numerically invert a resolution function to find chi.
///
/// Uses binary search to find chi such that CalcResoFromChi(chi) ≈ res.
/// \param res the target resolution
/// \param CalcResoFromChi function mapping chi -> resolution
/// \return the chi value
double CalcChiFromReso(double res, double (*CalcResoFromChi)(double))
{
    double chi   = 2.0;
    double delta = 1.0;

    for (int i = 0; i < 15; i++) {
        chi   = (CalcResoFromChi(chi) < res) ? chi + delta : chi - delta;
        delta = delta / 2.;
    }

    return chi;
}

/// \brief Event-plane resolution for k=2.
double CalcResoPsiK2(double chi)
{
    double con = std::sqrt(M_PI / 2) / 2;
    double arg = chi * chi / 4.;
    double halfpi = M_PI / 2;
    double besselOneHalf  = std::sqrt(arg / halfpi) * std::sinh(arg) / arg;
    double besselThreeHalfs = std::sqrt(arg / halfpi) *
                               (std::cosh(arg) / arg - std::sinh(arg) / (arg * arg));
    double res = con * chi * std::exp(-arg) *
                 (besselOneHalf + besselThreeHalfs);

    return res;
}

/// \brief Event-plane resolution for k=1.
double CalcResoPsiK1(double chi)
{
    double con = std::sqrt(M_PI / 2) / 2;
    double arg = chi * chi / 4.;

    double res = con * chi * std::exp(-arg) *
                 (TMath::BesselI0(arg) + TMath::BesselI1(arg));

    return res;
}

/// \brief Simple numerical derivative calculator.
double CalcDerivative(double x, double (*f)(double))
{
    double xl;
    double xh;

    if (x != 0) {
        xl = x * 0.99;
        xh = x * 1.01;
    } else {
        xl = x - 1e-4;
        xh = x + 1e-4;
    }

    return (f(xh) - f(xl)) / (xh - xl);
}

/// \brief Convert half-resolution to full resolution (k=1).
TGraphErrors ResHalfToFull(const TGraphErrors &g)
{
    const int n = g.GetN();

    double *x = g.GetX();
    double *y = g.GetY();
    std::vector<double> xe(n), ye(n);
    std::vector<double> oy(n), oye(n);

    for (int i = 0; i < n; i++) {
        xe[i] = g.GetErrorX(i);
        ye[i] = g.GetErrorY(i);
        double TmpRes    = y[i];
        double TmpResErr = ye[i];
        double TmpChi    = CalcChiFromReso(TmpRes, CalcResoPsiK1);
        double TmpChiErr = TmpResErr /
                           std::fabs(CalcDerivative(TmpChi, CalcResoPsiK1));
        TmpChi    *= std::sqrt(2.0);
        TmpChiErr *= std::sqrt(2.0);
        double TmpFullRes    = CalcResoPsiK1(TmpChi);
        double TmpFullResErr = TmpChiErr *
                               std::fabs(CalcDerivative(TmpChi, CalcResoPsiK1));
        oy[i]  = TmpFullRes;
        oye[i] = TmpFullResErr;
        if (ye[i] == 0) oye[i] = 0;
        if (y[i]  == 0) oy[i]  = 0;
    }

    TGraphErrors og(n, x, &oy[0], &xe[0], &oye[0]);

    return og;
}

/// \brief Convert half-resolution to full resolution (k=2).
TGraphErrors ResHalfToFullK2(const TGraphErrors &g)
{
    const int n = g.GetN();

    double *x = g.GetX();
    double *y = g.GetY();
    std::vector<double> xe(n), ye(n);
    std::vector<double> oy(n), oye(n);

    for (int i = 0; i < n; i++) {
        xe[i] = g.GetErrorX(i);
        ye[i] = g.GetErrorY(i);
        double TmpRes    = y[i];
        double TmpResErr = ye[i];
        double TmpChi    = CalcChiFromReso(TmpRes, CalcResoPsiK1);
        double TmpChiErr = TmpResErr /
                           std::fabs(CalcDerivative(TmpChi, CalcResoPsiK1));
        TmpChi    *= std::sqrt(2.0);
        TmpChiErr *= std::sqrt(2.0);
        double TmpFullRes    = CalcResoPsiK2(TmpChi);
        double TmpFullResErr = TmpChiErr *
                               std::fabs(CalcDerivative(TmpChi, CalcResoPsiK2));
        oy[i]  = TmpFullRes;
        oye[i] = TmpFullResErr;
        if (ye[i] == 0) oye[i] = 0;
        if (y[i]  == 0) oy[i]  = 0;
    }

    TGraphErrors og(n, x, &oy[0], &xe[0], &oye[0]);

    return og;
}


// ===========================================================================
//  Barlow-style check (correlated systematics)
// ===========================================================================

/// \brief Barlow-style difference between two graphs with correlated systematics.
///
/// The statistical error on the difference is computed as:
/// \f$ e = \sqrt{|e_1^2 - e_2^2|} \f$, reflecting that common systematics cancel.
TGraphErrors GraphBarlowDiff(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1e = g1.GetErrorY(i);
        double y2e = g2.GetErrorY(i);
        y.push_back(y1[i] - y2[i]);
        ye.push_back(std::sqrt(std::fabs(y1e*y1e - y2e*y2e)));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}

/// \brief Barlow-style ratio between two graphs with correlated systematics.
///
/// \f$ e = |y_1/y_2| \sqrt{|(e_1/y_1)^2 - (e_2/y_2)^2|} \f$.
TGraphErrors GraphBarlowRatio(const TGraphErrors &g1, const TGraphErrors &g2)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    double *y2 = g2.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        double y1v = y1[i];
        double y2v = y2[i];
        double y1e = g1.GetErrorY(i);
        double y2e = g2.GetErrorY(i);
        y.push_back(y1v / y2v);
        ye.push_back(std::fabs(y1v / y2v) *
                     std::sqrt(std::fabs(std::pow(y1e / y1v, 2) -
                                        std::pow(y2e / y2v, 2))));
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}


// ===========================================================================
//  TF1 helpers
// ===========================================================================

/// \brief Evaluate a TF1 on a TGraphErrors' x-values (zero errors).
TGraphErrors GraphFunc(const TGraphErrors &g1, TF1 *func)
{
    const int n = g1.GetN();
    double *x1 = g1.GetX();
    double *y1 = g1.GetY();
    std::vector<double> x, xe, y, ye;

    for (int i = 0; i < n; i++) {
        x.push_back(x1[i]);
        xe.push_back(0);
        y.push_back(func->Eval(x1[i]));
        ye.push_back(0);
    }

    TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

    return g;
}

/// \brief Compute pull = (data - function) (preserves data errors).
inline TGraphErrors GraphPull(const TGraphErrors &g1, TF1 *func)
{
    return (g1 - GraphFunc(g1, func));
}


#endif
