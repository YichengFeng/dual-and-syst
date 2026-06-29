/**
 * \file ValueErrors.h
 * \brief Formatted value-with-errors output for LaTeX and ROOT.
 *
 * \details Provides precision-aware formatting of values with statistical
 * and systematic uncertainties. Supports multiple output modes:
 * value-only, value±stat, value±stat±syst, and value±stat_{-systL}^{+systH}.
 * Output can be formatted for LaTeX ("L") or ROOT text ("R").
 *
 * \author Yicheng Feng
 * Email:  fengyich@outlook.com
 *
 * Improvements over old_code/MyValuErrs.h:
 *   - Renamed class to ValueErrors
 *   - Data members made private with snake_case getters
 *   - Replaced std::vector<double> returns with a private result struct
 *   - Added enum Mode with sentinel MODE_COUNT for extensibility
 *   - Const-correct print methods (no longer mutate internal state)
 *   - snake_case naming throughout
 *   - Removed redundant AddPrecision / SetPrecision (use recalc instead)
 *   - Efficient mode check: recalculates only when mode changes
 *   - Added noexcept where appropriate
 */

#ifndef ValueErrors_H
#define ValueErrors_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

/// \class ValueErrors
/// \brief Formats a value with its statistical and systematic uncertainties.
///
/// Normalizes the value and errors to scientific notation and determines
/// the appropriate number of significant digits for the uncertainty.
/// Supports LaTeX and ROOT text output formats.
class ValueErrors
{
public:
    /// \enum Mode
    /// \brief Output format mode.
    enum Mode {
        VALUE_ONLY      = 0,  ///< Value only, no errors
        VALUE_STAT      = 1,  ///< Value ± statistical error
        VALUE_STAT_SYST = 2,  ///< Value ± stat ± symmetric syst
        VALUE_STAT_ASYM = 3,  ///< Value ± stat _{-sysL}^{+sysH}
        MODE_COUNT            ///< Sentinel — always last
    };

    // --- constructors ---

    /// \brief Default constructor: value 0, stat mode.
    ValueErrors() { recalc(VALUE_STAT); }

    /// \brief Value only (no errors).
    explicit ValueErrors(double val) : value_(val)                       { recalc(VALUE_ONLY); }

    /// \brief Value with symmetric statistical error.
    /// \param val the central value
    /// \param stat the statistical uncertainty
    ValueErrors(double val, double stat) : value_(val), stat_(stat) { recalc(VALUE_STAT); }

    /// \brief Value with statistical and symmetric systematic errors.
    /// \param val the central value
    /// \param stat the statistical uncertainty
    /// \param syst the symmetric systematic uncertainty
    ValueErrors(double val, double stat, double syst)
        : value_(val), stat_(stat), syst_(syst), sysl_(syst)    { recalc(VALUE_STAT_SYST); }

    /// \brief Value with statistical and asymmetric systematic errors.
    /// \param val the central value
    /// \param stat the statistical uncertainty
    /// \param sysl the lower systematic uncertainty
    /// \param sysh the upper systematic uncertainty
    ValueErrors(double val, double stat, double sysl, double sysh)
        : value_(val), stat_(stat), sysl_(sysl), sysh_(sysh)    { recalc(VALUE_STAT_ASYM); }

    // --- accessors ---

    /// \brief Get the central value.
    double get_value() const noexcept { return value_; }
    /// \brief Get the statistical error.
    double get_stat()  const noexcept { return stat_; }
    /// \brief Get the symmetric systematic error.
    double get_syst()  const noexcept { return syst_; }
    /// \brief Get the lower asymmetric systematic error.
    double get_sysl()  const noexcept { return sysl_; }
    /// \brief Get the upper asymmetric systematic error.
    double get_sysh()  const noexcept { return sysh_; }

    /// \brief Get the computed output precision (decimal places).
    int    get_precision() const noexcept { return precision_; }
    /// \brief Get the scientific notation exponent.
    int    get_exponent()  const noexcept { return exponent_; }
    /// \brief Get the normalized value coefficient.
    double get_value_coeff() const noexcept { return value_coeff_; }
    /// \brief Get the normalized stat error coefficient.
    double get_stat_coeff()  const noexcept { return stat_coeff_; }
    /// \brief Get the normalized syst error coefficient.
    double get_syst_coeff()  const noexcept { return syst_coeff_; }
    /// \brief Get the normalized lower syst coefficient.
    double get_sysl_coeff()  const noexcept { return sysl_coeff_; }
    /// \brief Get the normalized upper syst coefficient.
    double get_sysh_coeff()  const noexcept { return sysh_coeff_; }

    // --- precision control ---

    /// \brief Increase the output precision by \p p decimal places.
    void add_precision(int p) noexcept { precision_ += p; }

    /// \brief Set the output precision to \p p decimal places.
    void set_precision(int p) noexcept { precision_ = p; }

    // --- formatting ---

    /// \brief Print the value with errors in the requested format.
    /// \param opt format option: "R" for ROOT text, "L" for LaTeX
    /// \param mode output mode (default -1: use current mode)
    /// \param force_normal force non-scientific notation
    /// \return formatted string
    std::string print(const std::string& opt = "R", int mode = -1,
                       bool force_normal = false) {
        if (mode >= VALUE_ONLY && mode < MODE_COUNT) recalc(static_cast<Mode>(mode));
        switch (mode_) {
            case VALUE_ONLY:      return print_err0(opt, force_normal);
            case VALUE_STAT:      return print_err1(opt, force_normal);
            case VALUE_STAT_SYST: return print_err2(opt, force_normal);
            case VALUE_STAT_ASYM: return print_err3(opt, force_normal);
            default:              return std::string();
        }
    }

    /// \brief Print value only.
    std::string print_err0(const std::string& opt = "R",
                           bool force_normal = false) {
        if (mode_ != VALUE_ONLY) recalc(VALUE_ONLY);
        int prec = precision_;
        int expo = exponent_;
        double vc = value_coeff_;
        int de = 0;
        int np = prec;
        bool isnormal = force_normal || (expo >= -1 && expo <= 0);
        if (isnormal) { de = expo; np = prec - expo; if (np < 0) np = 0; }
        std::ostringstream sv;
        sv << std::fixed << std::setprecision(np) << vc * std::pow(10.0, de);
        if (isnormal) return sv.str();
        return format_sci(sv.str(), "", expo, opt);
    }

    /// \brief Print value ± statistical error.
    std::string print_err1(const std::string& opt = "R",
                           bool force_normal = false) {
        if (mode_ != VALUE_STAT) recalc(VALUE_STAT);
        return print_err_impl(value_coeff_, stat_coeff_, precision_, exponent_, opt, force_normal);
    }

    /// \brief Print value ± stat ± symmetric syst.
    std::string print_err2(const std::string& opt = "R",
                           bool force_normal = false) {
        if (mode_ != VALUE_STAT_SYST) recalc(VALUE_STAT_SYST);
        return print_err2_impl(value_coeff_, stat_coeff_, syst_coeff_, precision_, exponent_, opt, force_normal);
    }

    /// \brief Print value ± stat _{-sysL}^{+sysH}.
    std::string print_err3(const std::string& opt = "R",
                           bool force_normal = false) {
        if (mode_ != VALUE_STAT_ASYM) recalc(VALUE_STAT_ASYM);
        return print_err3_impl(value_coeff_, stat_coeff_, sysl_coeff_, sysh_coeff_, precision_, exponent_, opt, force_normal);
    }

    /// \brief Print raw values (space-separated internal fields).
    std::string print_raw() const {
        std::ostringstream ss;
        ss << value_ << "   " << stat_ << "   " << syst_ << "   " << sysl_ << "   " << sysh_;
        return ss.str();
    }

private:
    // --- raw data ---
    double value_ = 0.0;   ///< Central value
    double stat_  = 0.0;   ///< Statistical uncertainty
    double syst_  = 0.0;   ///< Symmetric systematic uncertainty
    double sysl_  = 0.0;   ///< Lower asymmetric systematic uncertainty
    double sysh_  = 0.0;   ///< Upper asymmetric systematic uncertainty

    // --- formatted output parameters ---
    int    mode_      = VALUE_STAT;  ///< Current output mode
    int    precision_ = 0;           ///< Decimal places of the error
    int    exponent_  = 0;           ///< Scientific notation exponent
    double value_coeff_ = 0.0;       ///< Normalized value coefficient
    double stat_coeff_  = 0.0;       ///< Normalized stat error coefficient
    double syst_coeff_  = 0.0;       ///< Normalized syst error coefficient
    double sysl_coeff_  = 0.0;       ///< Normalized lower syst coefficient
    double sysh_coeff_  = 0.0;       ///< Normalized upper syst coefficient

    /// \struct FormatResult
    /// \brief Internal struct holding normalized formatting parameters.
    struct FormatResult {
        int precision;
        int exponent;
        double value_coeff;
        double stat_coeff;
        double syst_coeff;
        double sysl_coeff;
        double sysh_coeff;
    };

    // --- static helpers ---

    /// \brief Normalize a number to the range [1, 10).
    static double calc_coefficient(double x) {
        if (x == 0.0) return 0.0;
        double ax = std::fabs(x);
        while (ax < 1.0)  ax *= 10.0;
        while (ax >= 10.0) ax *= 0.1;
        return (x > 0.0) ? ax : -ax;
    }

    /// \brief Compute the scientific notation exponent.
    static int calc_exponent(double x) {
        if (x == 0.0) return 0;
        int n = 0;
        double ax = std::fabs(x);
        while (ax < 1.0)  { ax *= 10.0; --n; }
        while (ax >= 10.0) { ax *= 0.1;  ++n; }
        return n;
    }

    // --- core formatting logic ---

    /// \brief Compute formatting parameters for value ± error.
    static FormatResult valu_err1(double val, double err) {
        FormatResult r = {};
        if (val == 0.0 && err == 0.0) return r;

        double val_c = calc_coefficient(val);
        int    val_e = calc_exponent(val);
        double err_c = calc_coefficient(err);
        int    err_e = calc_exponent(err);

        if (val == 0.0) val_e = err_e;
        if (err == 0.0) err_e = val_e;

        double value, error;
        int exponent, precision;

        if (std::fabs(val) >= std::fabs(err) || val_e >= err_e) {
            value    = val_c;
            exponent = val_e;
            int err_de = err_e - val_e;
            error    = err_c * std::pow(10.0, err_de);
            precision = (std::fabs(err_c) < 3.5) ? (1 - err_de) : (0 - err_de);
        } else if (val_e + 1 == err_e && std::fabs(err_c) < 3.5) {
            value    = val_c / 10.0;
            exponent = err_e;
            error    = err_c;
            precision = 1;
        } else if (std::fabs(err_c) < 3.5) {
            value    = val_c * std::pow(10.0, val_e - err_e);
            exponent = err_e;
            error    = err_c;
            precision = 1;
        } else {
            value    = val_c * std::pow(10.0, val_e - err_e);
            exponent = err_e;
            error    = err_c;
            precision = 0;
        }

        r.precision   = precision;
        r.exponent    = exponent;
        r.value_coeff = value;
        r.stat_coeff  = error;
        return r;
    }

    /// \brief Compute formatting parameters for value only.
    static FormatResult valu_err0(double val) {
        FormatResult r = valu_err1(val, val);
        r.stat_coeff = 0.0;
        return r;
    }

    /// \brief Compute formatting parameters for value ± stat ± syst.
    static FormatResult valu_err2(double val, double stat, double syst) {
        double err = std::fabs(stat) > std::fabs(syst) ? std::fabs(stat) : std::fabs(syst);

        double val_c  = calc_coefficient(val);
        int    val_e  = calc_exponent(val);
        double stat_c = calc_coefficient(stat);
        int    stat_e = calc_exponent(stat);
        double syst_c = calc_coefficient(syst);
        int    syst_e = calc_exponent(syst);

        FormatResult r1 = valu_err1(val, err);
        int exponent = r1.exponent;

        FormatResult r = {};
        r.precision   = r1.precision;
        r.exponent    = exponent;
        r.value_coeff = val_c  * std::pow(10.0, val_e  - exponent);
        r.stat_coeff  = stat_c * std::pow(10.0, stat_e - exponent);
        r.syst_coeff  = syst_c * std::pow(10.0, syst_e - exponent);
        r.sysl_coeff  = r.syst_coeff;
        return r;
    }

    /// \brief Compute formatting parameters for asymmetric errors.
    static FormatResult valu_err3(double val, double stat, double sysl, double sysh) {
        double err = std::fabs(stat);
        if (std::fabs(sysl) > err) err = std::fabs(sysl);
        if (std::fabs(sysh) > err) err = std::fabs(sysh);

        double val_c  = calc_coefficient(val);
        int    val_e  = calc_exponent(val);
        double stat_c = calc_coefficient(stat);
        int    stat_e = calc_exponent(stat);
        double sysl_c = calc_coefficient(sysl);
        int    sysl_e = calc_exponent(sysl);
        double sysh_c = calc_coefficient(sysh);
        int    sysh_e = calc_exponent(sysh);

        FormatResult r1 = valu_err1(val, err);
        int exponent = r1.exponent;

        double sysl_val = sysl_c * std::pow(10.0, sysl_e - exponent);
        double sysh_val = sysh_c * std::pow(10.0, sysh_e - exponent);

        FormatResult r = {};
        r.precision   = r1.precision;
        r.exponent    = exponent;
        r.value_coeff = val_c  * std::pow(10.0, val_e  - exponent);
        r.stat_coeff  = stat_c * std::pow(10.0, stat_e - exponent);
        r.syst_coeff  = std::sqrt(sysl_val * sysl_val + sysh_val * sysh_val);
        r.sysl_coeff  = sysl_val;
        r.sysh_coeff  = sysh_val;
        return r;
    }

    /// \brief Recalculate all formatted coefficients for a given mode.
    void recalc(Mode mode) {
        mode_ = mode;
        FormatResult r;
        switch (mode_) {
            case VALUE_ONLY:      r = valu_err0(value_);                        break;
            case VALUE_STAT:      r = valu_err1(value_, stat_);                 break;
            case VALUE_STAT_SYST: r = valu_err2(value_, stat_, syst_);          break;
            case VALUE_STAT_ASYM: r = valu_err3(value_, stat_, sysl_, sysh_);   break;
            default: return;
        }
        precision_   = r.precision;
        exponent_    = r.exponent;
        value_coeff_ = r.value_coeff;
        stat_coeff_  = r.stat_coeff;
        syst_coeff_  = r.syst_coeff;
        sysl_coeff_  = r.sysl_coeff;
        sysh_coeff_  = r.sysh_coeff;
    }

    // --- print helpers ---

    /// \brief Format a number in scientific notation.
    /// \param val_str the mantissa string
    /// \param err_str the error mantissa string (empty for value-only)
    /// \param exponent the power-of-10 exponent
    /// \param opt format option ("L" for LaTeX, "R" for ROOT)
    static std::string format_sci(const std::string& val_str, const std::string& err_str,
                                   int exponent, const std::string& opt) {
        if (err_str.empty()) {
            if (opt == "L") return val_str + "\\times10^{" + std::to_string(exponent) + "}";
            return val_str + "#times10^{" + std::to_string(exponent) + "}";
        } else {
            if (opt == "L") return "(" + val_str + "\\pm" + err_str + ")\\times10^{" + std::to_string(exponent) + "}";
            return "(" + val_str + "#pm" + err_str + ")#times10^{" + std::to_string(exponent) + "}";
        }
    }

    /// \brief Implementation of value ± error printing.
    std::string print_err_impl(double vc, double ec, int prec, int expo,
                                const std::string& opt, bool force_normal = false) {
        int de = 0, np = prec;
        bool isnormal = force_normal || (expo >= -1 && expo <= 0);
        if (isnormal) { de = expo; np = prec - expo; if (np < 0) np = 0; }
        double factor = std::pow(10.0, de);
        std::ostringstream sv, se;
        sv << std::fixed << std::setprecision(np) << vc * factor;
        se << std::fixed << std::setprecision(np) << ec * factor;
        if (isnormal) {
            if (opt == "L") return sv.str() + "\\pm" + se.str();
            return sv.str() + "#pm" + se.str();
        }
        return format_sci(sv.str(), se.str(), expo, opt);
    }

    /// \brief Implementation of value ± stat ± syst printing.
    std::string print_err2_impl(double vc, double sc, double syc, int prec, int expo,
                                 const std::string& opt, bool force_normal = false) {
        int de = 0, np = prec;
        bool isnormal = force_normal || (expo >= -1 && expo <= 0);
        if (isnormal) { de = expo; np = prec - expo; if (np < 0) np = 0; }
        double factor = std::pow(10.0, de);
        std::ostringstream sv, ss, sy;
        sv << std::fixed << std::setprecision(np) << vc * factor;
        ss << std::fixed << std::setprecision(np) << sc * factor;
        sy << std::fixed << std::setprecision(np) << syc * factor;
        if (isnormal) {
            if (opt == "L") return sv.str() + "\\pm" + ss.str() + "\\pm" + sy.str();
            return sv.str() + "#pm" + ss.str() + "#pm" + sy.str();
        }
        if (opt == "L") return "(" + sv.str() + "\\pm" + ss.str() + "\\pm" + sy.str() + ")\\times10^{" + std::to_string(expo) + "}";
        return "(" + sv.str() + "#pm" + ss.str() + "#pm" + sy.str() + ")#times10^{" + std::to_string(expo) + "}";
    }

    /// \brief Implementation of asymmetric error printing.
    std::string print_err3_impl(double vc, double sc, double slc, double shc, int prec, int expo,
                                 const std::string& opt, bool force_normal = false) {
        int de = 0, np = prec;
        bool isnormal = force_normal || (expo >= -1 && expo <= 0);
        if (isnormal) { de = expo; np = prec - expo; if (np < 0) np = 0; }
        double factor = std::pow(10.0, de);
        std::ostringstream sv, ss, sl, sh;
        sv << std::fixed << std::setprecision(np) << vc * factor;
        ss << std::fixed << std::setprecision(np) << sc * factor;
        sl << std::fixed << std::setprecision(np) << slc * factor;
        sh << std::fixed << std::setprecision(np) << shc * factor;
        if (isnormal) {
            if (opt == "L") return sv.str() + "\\pm" + ss.str() + "_{-" + sl.str() + "}^{+" + sh.str() + "}";
            return sv.str() + "#pm" + ss.str() + "_{-" + sl.str() + "}^{+" + sh.str() + "}";
        }
        if (opt == "L") return "(" + sv.str() + "\\pm" + ss.str() + "_{-" + sl.str() + "}^{+" + sh.str() + "})\\times10^{" + std::to_string(expo) + "}";
        return "(" + sv.str() + "#pm" + ss.str() + "_{-" + sl.str() + "}^{+" + sh.str() + "})#times10^{" + std::to_string(expo) + "}";
    }
};


#endif
