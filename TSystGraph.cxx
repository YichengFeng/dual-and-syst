/**************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   — Added custom Streamer for ROOT I/O that serialises the SystGraph
 *     base-class members (def_, var_, merge_mode_) that the auto-generated
 *     Streamer must skip (SystGraph has no ROOT dictionary).
 **************************************************************************/

#include "TSystGraph.h"

#include <utility>   // std::move

ClassImp(TSystGraph)

// ===========================================================================
//  Internal helpers — serialise a DualGraph (points_ array) via TBuffer
//  primitives.  Same wire format as TDualGraph::Streamer so that the two
//  are interchangeable on disk.
// ===========================================================================

namespace {

// ── helper: write one DualGraph's points_ to buffer ──────────────────────
void write_dualgraph(TBuffer &b, const DualGraph &dg)
{
    Int_t npts = static_cast<Int_t>(dg.get_n());
    b << npts;

    for (Int_t i = 0; i < npts; ++i) {
        auto pt = dg.get_point(i);
        b << pt.get_px().get_value();
        b << pt.get_px().get_derivative();
        b << pt.get_py().get_value();

        auto ids   = pt.get_py().get_ids();
        Int_t nd   = static_cast<Int_t>(ids.size());
        b << nd;
        for (Int_t id : ids) {
            b << id;
            b << pt.get_py().get_derivative(id);
        }
    }
}

// ── helper: read one DualGraph's points_ from buffer ────────────────────
void read_dualgraph(TBuffer &b, DualGraph &dg)
{
    Int_t npts;
    b >> npts;

    // Reset dg.  We need write-access to the underlying points_ vector.
    // Since DualGraph exposes reset(n) (virtual) and set_point() we can
    // use the public API, but it is simpler to access points_ directly
    // from the derived class.  Here we operate on a standalone DualGraph
    // via reset() + set_point().
    dg.reset(npts);
    for (Int_t i = 0; i < npts; ++i) {
        Double_t px_val, px_der, py_val;
        Int_t    nduals;
        b >> px_val >> px_der >> py_val >> nduals;

        DualPoint dp(DualNumber(px_val, px_der), DualMultiv(py_val));
        for (Int_t j = 0; j < nduals; ++j) {
            Int_t    id;
            Double_t deriv;
            b >> id >> deriv;
            dp.get_py().set_derivative(id, deriv);
        }
        // Use set_point() which is virtual and checks bounds.
        // After reset(npts) the indices 0..npts-1 are valid.
        dg.set_point(i, dp);
    }
}

}  // anonymous namespace

// ===========================================================================
//  Custom ROOT I/O Streamer
//
//  The auto-generated Streamer (suppressed via LinkDef.h's '-' flag) only
//  knows about TSystGraph's own ROOT-type data members (stat_, syst_,
//  comb_) and the TObject base.  It skips the SystGraph base entirely
//  (def_, var_, merge_mode_) because SystGraph has no ClassDef/dictionary.
//
//  This hand-written Streamer replaces it and manually serialises the
//  SystGraph base members using TBuffer primitive operations, preserving
//  the full dual-number and variation information.
// ===========================================================================

void TSystGraph::Streamer(TBuffer &R__b)
{
    UInt_t R__s, R__c;
    if (R__b.IsReading()) {
        Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
        if (R__v) { }   // (silence unused-variable warning)

        TObject::Streamer(R__b);

        R__b >> is_updated_;
        R__b >> width_;
        R__b >> xlim_low_ >> xlim_high_;

        // calc_result_ is a transient cache — skip it (rebuilt by calc())
        stat_.Streamer(R__b);
        syst_.Streamer(R__b);
        comb_.Streamer(R__b);

        // ---- SystGraph base members ----
        {
            Int_t mode;
            R__b >> mode;
            merge_mode_ = static_cast<SystMergeMode>(mode);
        }

        read_dualgraph(R__b, def_);

        {
            Int_t nvars;
            R__b >> nvars;
            var_.clear();
            for (Int_t i = 0; i < nvars; ++i) {
                Int_t vid;
                R__b >> vid;

                SystVariation sv;
                R__b.ReadStdString(&sv.name);
                read_dualgraph(R__b, sv.graph);
                R__b >> sv.weight;
                {
                    Int_t m;
                    R__b >> m;
                    sv.mode = static_cast<SystMergeMode>(m);
                }
                var_.emplace(vid, std::move(sv));
            }
        }

        R__b.CheckByteCount(R__s, R__c, TSystGraph::IsA());

    } else {
        R__c = R__b.WriteVersion(TSystGraph::IsA(), kTRUE);

        TObject::Streamer(R__b);

        R__b << is_updated_;
        R__b << width_;
        R__b << xlim_low_ << xlim_high_;

        stat_.Streamer(R__b);
        syst_.Streamer(R__b);
        comb_.Streamer(R__b);

        // ---- SystGraph base members ----
        R__b << static_cast<Int_t>(merge_mode_);

        write_dualgraph(R__b, def_);

        {
            Int_t nvars = static_cast<Int_t>(var_.size());
            R__b << nvars;

            for (const auto &kv : var_) {
                R__b << kv.first;     // vid
                R__b.WriteStdString(&kv.second.name);
                write_dualgraph(R__b, kv.second.graph);
                R__b << kv.second.weight;
                R__b << static_cast<Int_t>(kv.second.mode);
            }
        }

        R__b.SetByteCount(R__c, kTRUE);
    }
}
