/**************************************************************************
 * LinkDef_TSystGraph.h — rootcling selection file (TSystGraph only)
 *
 * The trailing '-' tells rootcling to NOT generate an auto-Streamer so
 * that we can provide a hand-written one that manually serialises the
 * SystGraph base-class members (def_, var_, merge_mode_) which are
 * otherwise invisible to ROOT I/O.
 *
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 **************************************************************************/

#ifdef __CLING__
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ class TSystGraph-;
#endif
