#ifndef AnExCacheSim_h
#define AnExCacheSim_h

#include "AnalExtractor.h"

#include "SXrdClasses.h"

class TH1;

//==============================================================================

struct CacheState
{
  long64            f_fs;  // file size
  long64            f_bs;  // block size
  int               f_nb;  // num blocks
  double            f_pfr; // prefetch rate in B/s

  // Performance vars / accumulators
  long64            f_bytes_needed = 0; // asked for by client
  long64            f_trips_needed = 0;
  long64            f_bytes_done   = 0; // downloaded via read requests
  long64            f_trips_done   = 0;
  long64            f_bytes_extra  = 0; // extra bytes in issued read requests;
  long64            f_trips_extra  = 0; // these can be reused later, see saved.
  long64            f_bytes_saved  = 0; // bytes that were served from cache
  long64            f_trips_saved  = 0;
  long64            f_bytes_pref   = 0; // prefetched
  long64            f_trips_pref   = 0;

  double done_o_req()     { return (double)  f_bytes_done / f_bytes_needed; }
  double gotten_o_req()   { return (double) (f_bytes_done + f_bytes_pref) / f_bytes_needed; }
  double gotten_o_fs()    { return (double) (f_bytes_done + f_bytes_pref) / f_fs; }
  double saved_o_req()    { return (double)  f_bytes_saved / f_bytes_needed; }
  double trpsaved_o_req() { return (double)  f_trips_saved / f_trips_needed; }
  double extra_o_req()    { return (double)  f_bytes_extra / f_bytes_needed; }
  double extra_o_fs()     { return (double)  f_bytes_extra / f_fs; }
  double unused_o_req()   { return (double) (f_bytes_extra - f_bytes_saved) / f_bytes_needed; }
  double unused_o_fs()    { return (double) (f_bytes_extra - f_bytes_saved) / f_fs; }
  double saved_o_done()   { return (double)  f_bytes_saved / f_bytes_done; }
  double extra_o_done()   { return (double)  f_bytes_extra / f_bytes_done; }


  // State
  std::vector<bool> f_blocks;
  std::set<int>     f_blocks_to_get;

  long64            f_prev_bytes_needed = 0;
  double            f_pref_carry   = 0;
  int               f_curr_time    = 0;
  int               f_prev_time    = 0;

  int               f_pf_block     = 0;


  CacheState(long64 fs, long64 bs, double pfr);

  void BeginRequest(int t);
  void Read(long64 req_off, int req_len);
  void EndRequest();

  void Finish();

  void Print();
};

//==============================================================================

struct ZHistos
{
  static const int NN = 7;

  TH1 *done_req[NN], *donewpf_req[NN], *donewpf_fs[NN], *saved_req[NN], *trpsaved_req[NN],
      *extr_req[NN], *extr_fs[NN], *unused_req[NN], *unused_fs[NN],
      *saved_done[NN], *extra_done[NN];

  void BookSet(TH1 *h[NN], const TString &nnn, const TString &ttt,
               int Xn, double Xl, double Xh);
  void Book(const TString &nnn, const TString &ttt);


  void FillInter(CacheState &cs, double t);

  void Fill(CacheState &cs, double frac_read, double frac_vread, double data_rate,
            double duration, double file_size);

  static double shl_vread_min;
};

//==============================================================================

class AnExCacheSim : public AnalExtractor
{
  static const int NN = 8, PP = 6;

  // cache sim block sizes in kB, prefetch rates in MB/s
  int     csbs[NN] = { 64, 128, 256, 512, 1024, 2048, 4096, 8192 };
  double  cspf[PP] = { 0, 0.25, 1, 4, 16, 64 }; 

  const char *nmbs[NN] = { "64kB", "128kB", "256kB", "512kB", "1MB", "2MB", "4MB", "8MB" };
  const char *nmpf[PP] = { "0", "256kBps", "1MBps", "4MBps", "16MBps", "64MBps" }; 

  // cache sim stats histos
  ZHistos cshi[PP][NN];

  double  f_min_vread;
  bool    f_print_and_wait = false;


public:

  AnExCacheSim(const TString& name, AnalManager &mgr, const TString &out_file="",
               double min_vread=0);

  // ----------------------------------------------------------------

  void BookHistos();

  void Process();

  void WriteHistos();

  // ----------------------------------------------------------------

  void SimulateCache(ZHistos &zhis, int blk_size, double pref_rate);
};

#endif
