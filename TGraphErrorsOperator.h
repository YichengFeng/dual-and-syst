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
TGraphErrors operator+(const TGraphErrors &g1, const TGraphErrors &g2);


/// \brief Scalar + TGraphErrors: adds constant to all y-values (errors unchanged).
TGraphErrors operator+(double c, const TGraphErrors &g1);


TGraphErrors operator+(const TGraphErrors &g1, double c);


/// \brief Pointwise TGraphErrors * TGraphErrors with error propagation.
/// \f$ e = \sqrt{y_2^2 e_1^2 + y_1^2 e_2^2} \f$.
TGraphErrors operator*(const TGraphErrors &g1, const TGraphErrors &g2);


/// \brief Scalar * TGraphErrors.
TGraphErrors operator*(double c, const TGraphErrors &g1);


TGraphErrors operator*(const TGraphErrors &g1, double c);


/// \brief TGraphErrors - TGraphErrors.
TGraphErrors operator-(const TGraphErrors &g1, const TGraphErrors &g2);


TGraphErrors operator-(double c, const TGraphErrors &g1);
TGraphErrors operator-(const TGraphErrors &g1, double c);


/// \brief Pointwise TGraphErrors / TGraphErrors with error propagation.
/// \f$ e = |y_1/y_2| \sqrt{(e_1/y_1)^2 + (e_2/y_2)^2} \f$.
TGraphErrors operator/(const TGraphErrors &g1, const TGraphErrors &g2);


/// \brief Scalar / TGraphErrors.
TGraphErrors operator/(double c, const TGraphErrors &g1);


TGraphErrors operator/(const TGraphErrors &g1, double c);


/// \brief Shift x-values by \p c (operator%).
TGraphErrors operator%(const TGraphErrors &g1, double c);


/// \brief Subtract a TF1 from a TGraphErrors (preserves errors).
TGraphErrors operator-(const TGraphErrors &g, const TF1 &f);


// ===========================================================================
//  In-place graph mutation helpers
// ===========================================================================

/// \brief Shift x-values of a TGraphErrors in-place by \p c.
void GraphShiftX(TGraphErrors *g1, double c);


/// \brief Shift x-values of a TGraphErrors in-place to 0.
void GraphShiftX(TGraphErrors *g1);


/// \brief Shift x-values of a TGraphErrors in-place by \p c (reference overload).
void GraphShiftX(TGraphErrors &g1, double c);


/// \brief Shift x-values of a TGraphErrors in-place to 0 (reference overload).
void GraphShiftX(TGraphErrors &g1);


/// \brief Shift x-values of a TGraphAsymmErrors in-place by \p c.
void GraphShiftX(TGraphAsymmErrors *g1, double c);


/// \brief Shift x-values of a TGraphAsymmErrors in-place to 0.
void GraphShiftX(TGraphAsymmErrors *g1);


// ===========================================================================
//  Graph creation / conversion
// ===========================================================================

/// \brief Convert a TH1D to a TGraphErrors.
TGraphErrors GraphFrom(const TH1D *h);


/// \brief Convert a TH1F to a TGraphErrors.
TGraphErrors GraphFrom(const TH1F *h);


/// \brief Copy a TGraphErrors (with errors).
TGraphErrors GraphFrom(const TGraphErrors *g1);


/// \brief Copy a TGraphErrors via copy constructor.
TGraphErrors GraphSelf(const TGraphErrors &g1);


// ===========================================================================
//  Centrality rebinning
// ===========================================================================

/// \brief Map a TGraphErrors to centrality bin x-values (9-bin default).
TGraphErrors GraphCent(const TGraphErrors &g1);


/// \brief Alias for GraphCent (9-bin).
TGraphErrors GraphCent9(const TGraphErrors &g1);

/// \brief Map to 8 centrality bins.
TGraphErrors GraphCent8(const TGraphErrors &g1);


/// \brief Map to 10 centrality bins.
TGraphErrors GraphCent10(const TGraphErrors &g1);


/// \brief Map to 16 centrality bins.
TGraphErrors GraphCent16(const TGraphErrors &g1);


/// \brief Map a TGraphAsymmErrors to centrality bin x-values.
TGraphAsymmErrors GraphCent(const TGraphAsymmErrors &g1);


/// \brief Create an error-band envelope polygon from a TGraphErrors.
///
/// Traces the upper edge forward, then the lower edge backward,
/// suitable for drawing with "f" (fill) option.
TGraph GraphEdge(const TGraphErrors &g1);


/// \brief Create an error-band envelope polygon from a TGraphAsymmErrors.
TGraph GraphEdge(const TGraphAsymmErrors &g1);


// ===========================================================================
//  Point-wise transformations
// ===========================================================================

/// \brief Pointwise absolute value.
TGraphErrors GraphAbs(const TGraphErrors &g1);


/// \brief Pointwise power \f$ y^c \f$ with error propagation.
/// \f$ e_y = |c| \cdot y^{c-1} \cdot e \f$.
TGraphErrors GraphPow(const TGraphErrors &g1, double c);


/// \brief Pointwise square root (of absolute value).
TGraphErrors GraphSqrt(const TGraphErrors &g1);


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
TGraphErrors GraphAveY(const TGraphErrors &g1, int firstpoint, int lastpoint);


/// \brief Inverse-variance weighted average using errors from a second graph.
///
/// Weights come from g2's y-errors; g1 provides the y-values.
TGraphErrors GraphWtAveY(const TGraphErrors &g1, const TGraphErrors &g2,
                          int firstpoint, int lastpoint);


/// \brief Average over all bins.
TGraphErrors GraphAveY(const TGraphErrors &g1);


/// \brief Create a single-point TGraphErrors (x=0).
TGraphErrors GraphOne(double con, double err);


/// \brief Create a single-point TGraphErrors with explicit x.
TGraphErrors GraphOne(double vx, double ex, double con, double err);


// ===========================================================================
//  Filtering and error manipulation
// ===========================================================================

/// \brief Remove points with zero y-value or invalid error.
TGraphErrors RemoveZero(const TGraphErrors &g);


/// \brief Zero all y-errors (keep x-errors and values).
TGraphErrors ZeroError(const TGraphErrors &g);


// ===========================================================================
//  Graph combination
// ===========================================================================

/// \brief Map two graphs: x = y of g1, y = y of g2 (with errors).
TGraphErrors GraphMap(const TGraphErrors &g1, const TGraphErrors &g2);


/// \brief Concatenate two TGraphErrors.
TGraphErrors GraphMerge(const TGraphErrors &g1, const TGraphErrors &g2);


/// \brief Swap y-values with y-errors (x-errors become new y-errors, zeroed).
TGraphErrors GraphErrToVal(const TGraphErrors &g1);


/// \brief Shift graph by \p c times its errors: \f$ y + c \cdot e_y \f$.
TGraphErrors GraphMoveSigma(const TGraphErrors &g1, double c);


// ===========================================================================
//  Event-plane resolution helpers
// ===========================================================================

/// \brief Numerically invert a resolution function to find chi.
///
/// Uses binary search to find chi such that CalcResoFromChi(chi) ≈ res.
/// \param res the target resolution
/// \param CalcResoFromChi function mapping chi -> resolution
/// \return the chi value
double CalcChiFromReso(double res, double (*CalcResoFromChi)(double));

/// \brief Event-plane resolution for k=2.
double CalcResoPsiK2(double chi);

/// \brief Event-plane resolution for k=1.
double CalcResoPsiK1(double chi);

/// \brief Simple numerical derivative calculator.
double CalcDerivative(double x, double (*f)(double));

/// \brief Convert half-resolution to full resolution (k=1).
TGraphErrors ResHalfToFull(const TGraphErrors &g);

/// \brief Convert half-resolution to full resolution (k=2).
TGraphErrors ResHalfToFullK2(const TGraphErrors &g);


// ===========================================================================
//  Barlow-style check (correlated systematics)
// ===========================================================================

/// \brief Barlow-style difference between two graphs with correlated systematics.
///
/// The statistical error on the difference is computed as:
/// \f$ e = \sqrt{|e_1^2 - e_2^2|} \f$, reflecting that common systematics cancel.
TGraphErrors GraphBarlowDiff(const TGraphErrors &g1, const TGraphErrors &g2);

/// \brief Barlow-style ratio between two graphs with correlated systematics.
///
/// \f$ e = |y_1/y_2| \sqrt{|(e_1/y_1)^2 - (e_2/y_2)^2|} \f$.
TGraphErrors GraphBarlowRatio(const TGraphErrors &g1, const TGraphErrors &g2);


// ===========================================================================
//  TF1 helpers
// ===========================================================================

/// \brief Evaluate a TF1 on a TGraphErrors' x-values (zero errors).
TGraphErrors GraphFunc(const TGraphErrors &g1, TF1 *func);

/// \brief Compute pull = (data - function) (preserves data errors).
TGraphErrors GraphPull(const TGraphErrors &g1, TF1 *func);


#endif
