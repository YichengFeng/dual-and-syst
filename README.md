# DualAndSyst

**Dual-number graph toolkit with systematic uncertainty propagation for ROOT.**

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)](https://en.cppreference.com/w/cpp/17)
[![ROOT](https://img.shields.io/badge/ROOT-6.34-green)](https://root.cern/)

---

## Overview

DualAndSyst is a C++17 header-only library for high-energy physics data analysis. It provides:

- **Automatic differentiation** via dual numbers (`DualNumber`, `DualMultiv`) — compute first-order derivatives of arbitrary expressions without manual calculus.
- **Uncertainty propagation** through graph arithmetic (`DualGraph`) — add, subtract, multiply, and divide graphs while correctly propagating statistical and systematic uncertainties.
- **Systematic uncertainty evaluation** (`SystGraph`) — manage multiple systematic variation graphs per observable, evaluate them under Simple or Barlow (subset / independent / fully-correlated) prescriptions, and produce publication-ready LaTeX tables.
- **ROOT integration** (`TDualGraph`, `TSystGraph`) — ROOT-aware counterparts with `TGraphErrors` mirrors, `TFile` I/O, canvas drawing, overlay plots, and stacked systematic<sup>2</sup>-breakdown histograms.
- **TGraphErrors utilities** (`TGraphErrorsOperator`) — free functions for arithmetic, rebinning, Barlow difference/ratio, event-plane resolution, and more.

Systematic uncertainties follow [Barlow's prescription](https://arxiv.org/abs/hep-ex/0207026) (hep-ex/0207026).

---

## Class hierarchy

```
                    DualNumber ── ValueErrors
                         │
                    DualMultiv
                    ╱        ╲
             DualPoint    (correlation matrix)
                │
           DualGraph ──────────── TDualGraph
                │                      │
          SystVariation          TSystVariation
                │                      │
           SystGraph ──────────── TSystGraph
                                     │
                            TGraphErrorsOperator
```

| Layer | ROOT-free | ROOT-aware |
|---|---|---|
| **Math core** | `DualNumber`, `DualMultiv`, `ValueErrors` | — |
| **Graph core** | `DualPoint`, `DualGraph` | `TDualGraph` |
| **Systematics** | `SystGraph`, `SystVariation`, `SystResult` | `TSystGraph`, `TSystVariation` |
| **Utilities** | — | `TGraphErrorsOperator` |

---

## Quick start

### Prerequisites

- **ROOT** ≥ 6.24 (tested with 6.34)
- **g++** with C++17 support
- **Doxygen** + **Graphviz** (for documentation)

### Build the ROOT dictionaries

```bash
# Edit Makefile to point ROOT_BINDIR to your ROOT installation, then:
make
```

This produces `TDualGraph_cxx.so` and `TSystGraph_cxx.so` — loadable ROOT shared libraries.

### Use in a ROOT macro

```cpp
// Load the library
gSystem->Load("TSystGraph_cxx.so");

// Create a default graph from a histogram
TFile f("data.root");
TH1D* h = (TH1D*)f.Get("myHist");
TSystGraph sg(TDualGraph(*h));

// Add systematic variations
TSystVariation sv;
sv.name   = "Jet Energy Scale";
sv.graph  = TDualGraph(*h_jes_up);
sv.mode   = SystMergeMode::kOneSided_Barlow_Subset;
sv.weight = 1.0;
sg.add_var(sv);

// Evaluate and draw
sg.calc();
sg.set_style_color(20, kBlue);
sg.draw_syst();
sg.make_plot("myPlot", "plots/");
sg.make_hist("myBreakdown", true, "plots/");

// Print LaTeX table
sg.print();
```

### Use header-only (no ROOT dictionaries needed)

Most classes are header-only. If you only need the ROOT-independent classes:

```cpp
#include "SystGraph.h"

// Build a graph point-by-point
DualGraph dg(5);
for (int i = 0; i < 5; i++) {
    dg.set_point(i, DualPoint(
        DualNumber(i * 10.0, 5.0),        // x ± xerr
        DualMultiv(100.0 - i * 5.0, 3.0)  // y ± yerr (auto-assigned ID)
    ));
}

// Arithmetic with full uncertainty propagation
DualGraph ratio = dg / dg;  // ratio to self = 1.0 everywhere
```

### Systematic uncertainty evaluation

The library supports 8 merge modes via `SystMergeMode`:

| Mode | Type | Description |
|---|---|---|
| `kOneSided_Simple` | 1-sided | `σ = (v − d)` |
| `kOneSided_Barlow_Subset` | 1-sided | σ² = (v−d)² − &#124;e<sub>v</sub>² − e<sub>d</sub>²&#124; |
| `kOneSided_Barlow_Independent` | 1-sided | σ² = (v−d)² − &#124;e<sub>v</sub>² + e<sub>d</sub>²&#124; |
| `kOneSided_Barlow_FullyCorrelated` | 1-sided | σ² = (v−d)² − &#124;e<sub>v</sub> − e<sub>d</sub>&#124;² |
| `kTwoSided_*` | 2-sided | Same formulas, symmetric ± |

One-sided systematics contribute to either upper or lower error (sign-dependent); two-sided contribute symmetrically.

---

## Documentation

**Online:** [yichengfeng.github.io/dual-and-syst](https://yichengfeng.github.io/dual-and-syst/)

**Generate locally:**

```bash
doxygen Doxyfile
xdg-open html/index.html
```

The documentation includes class references, method descriptions, call/caller graphs, and LaTeX-rendered formulas.

---

## Project structure

```
.
├── DualNumber.h              # Single-variable dual number
├── DualMultiv.h              # Multi-variable dual number + correlation matrix
├── ValueErrors.h             # Formatted value±error output (LaTeX / ROOT)
├── DualGraph.h               # DualPoint + DualGraph (ROOT-free)
├── TDualGraph.h / .cxx       # ROOT-aware DualGraph
├── SystGraph.h               # Systematic uncertainty graph (ROOT-free)
├── TSystGraph.h / .cxx       # ROOT-aware SystGraph
├── TGraphErrorsOperator.h    # TGraphErrors utility functions
├── LinkDef_TDualGraph.h      # ROOT dictionary linkdef
├── LinkDef_TSystGraph.h      # ROOT dictionary linkdef
├── Makefile                  # Build ROOT shared libraries
├── Doxyfile                  # Doxygen configuration
└── test/                     # Test macros and scripts
```

---

## Author

**Yicheng Feng** — [fengyich@outlook.com](mailto:fengyich@outlook.com)

---

## License

[MIT](LICENSE) © 2026 Yicheng Feng
