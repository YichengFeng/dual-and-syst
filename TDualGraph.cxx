/**************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 *
 * Improvements over old_code version:
 *   — Added custom Streamer for ROOT I/O that serialises the DualGraph
 *     base-class points_ vector that the auto-generated Streamer must
 *     skip (DualGraph has no ROOT dictionary).
 **************************************************************************/

#include "TDualGraph.h"

#include <utility>   // std::move

ClassImp(TDualGraph)

// ===========================================================================
//  Custom ROOT I/O Streamer
//
//  The auto-generated Streamer (suppressed via LinkDef.h's '-' flag) only
//  knows about TDualGraph's own data members (is_updated_, graph_) and
//  the TObject base.  It skips the DualGraph base entirely because
//  DualGraph has no ClassDef/dictionary.
//
//  This hand-written Streamer replaces it and manually serialises the
//  points_ vector element-by-element using TBuffer primitive operations,
//  preserving the full dual-number information (IDs and partial
//  derivatives).
// ===========================================================================

void TDualGraph::Streamer(TBuffer &R__b)
{
    UInt_t R__s, R__c;
    if (R__b.IsReading()) {
        Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
        if (R__v) { }   // (silence unused-variable warning in older ROOT)

        TObject::Streamer(R__b);
        R__b >> is_updated_;
        graph_.Streamer(R__b);

        // ---- read points_ manually ----
        Int_t npts;
        R__b >> npts;
        points_.clear();
        points_.reserve(static_cast<std::size_t>(npts));

        for (Int_t i = 0; i < npts; ++i) {
            Double_t px_val, px_der, py_val;
            Int_t    nduals;

            R__b >> px_val >> px_der >> py_val >> nduals;

            DualPoint dp(DualNumber(px_val, px_der), DualMultiv(py_val));
            for (Int_t j = 0; j < nduals; ++j) {
                Int_t    id;
                Double_t deriv;
                R__b >> id >> deriv;
                dp.get_py().set_derivative(id, deriv);
            }
            points_.push_back(std::move(dp));
        }

        R__b.CheckByteCount(R__s, R__c, TDualGraph::IsA());

    } else {
        R__c = R__b.WriteVersion(TDualGraph::IsA(), kTRUE);

        TObject::Streamer(R__b);
        R__b << is_updated_;
        graph_.Streamer(R__b);

        // ---- write points_ manually ----
        Int_t npts = static_cast<Int_t>(points_.size());
        R__b << npts;

        for (const auto &pt : points_) {
            R__b << pt.get_px().get_value();
            R__b << pt.get_px().get_derivative();
            R__b << pt.get_py().get_value();

            auto ids  = pt.get_py().get_ids();
            Int_t nduals = static_cast<Int_t>(ids.size());
            R__b << nduals;

            for (Int_t id : ids) {
                R__b << id;
                R__b << pt.get_py().get_derivative(id);
            }
        }

        R__b.SetByteCount(R__c, kTRUE);
    }
}
