/**************************************************************************
 * LinkDef.h — rootcling selection file (TDualGraph only)
 *
 * The trailing '-' tells rootcling to NOT generate an auto-Streamer so
 * that we can provide a hand-written one that manually serialises the
 * DualGraph base-class points_ (otherwise invisible to ROOT I/O).
 *
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 **************************************************************************/

#ifdef __CLING__
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ class TDualGraph-;
#endif
