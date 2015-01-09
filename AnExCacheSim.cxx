#include "AnExCacheSim.h"
#include "AnalManager.h"

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

namespace
{
  const double OneMB = 1024 * 1024;

  bool overlap(int     blk,      // block to query
               long64  blk_size, //
               long64  req_off,  // offset of user request
               int     req_size, // size of user request
               // output:
               long64 &off,     // offset in user buffer
               long64 &blk_off, // offset in block
               long64 &size)    // size to copy
  {
    const long64 beg     = blk * blk_size;
    const long64 end     = beg + blk_size;
    const long64 req_end = req_off + req_size;

    if (req_off < end && req_end > beg)
    {
      const long64 ovlp_beg = std::max(beg, req_off);
      const long64 ovlp_end = std::min(end, req_end);

      off     = ovlp_beg - req_off;
      blk_off = ovlp_beg - beg;
      size    = ovlp_end - ovlp_beg;

      return true;
    }
    else
    {
      return false;
    }
  }
}


//==============================================================================
// CacheState
//==============================================================================

CacheState::CacheState(long64 fs, long64 bs, double pfr) :
  f_fs(fs), f_bs(bs), f_nb((fs - 1) / bs + 1), f_pfr(pfr),
  f_blocks(f_nb)
{}

void CacheState::BeginRequest(int t)
{
  // Mark blocks that got prefetched since the last request.

  f_prev_bytes_needed = f_bytes_needed;

  f_prev_time = f_curr_time;
  f_curr_time = t;

  double pf = f_pfr * (t - f_prev_time) + f_pref_carry;
  if (pf > 0 && f_pf_block < f_nb)
  {
    int nb_to_mark = TMath::CeilNint(pf / f_bs);
    f_pref_carry = pf - f_bs * nb_to_mark;

    while (nb_to_mark > 0)
    {
      if ( ! f_blocks[f_pf_block])
      {
        f_blocks[f_pf_block] = true;
        --nb_to_mark;

        f_bytes_pref += f_bs;
        f_trips_pref += 1;
      }

      if (++f_pf_block >= f_nb) break;
    }
  }
  else
  {
    f_pref_carry = pf;
  }
}

void CacheState::Read(long64 req_off, int req_len)
{
  // Single read request or a sub-request from a vec read.
  // Chop into blocks and process them individually.

  if (req_off + req_len > f_fs)
  {
    printf("CacheState::Read Achtung, req over the end of the file: f_fs=%lld req_off=%lld req_len=%d\n",
           f_fs, req_off, req_len);
    return;
  }

  int b_min = req_off / f_bs;
  int b_max = (req_off + req_len - 1) / f_bs + 1;

  long64 off, blk_off, size;
  for (int bi = b_min; bi < b_max; ++bi)
  {
    if (overlap(bi, f_bs, req_off, req_len, off, blk_off, size))
    {
      f_bytes_needed += size;

      if (f_blocks[bi])
      {
        // The block is already here
        f_bytes_saved += size;
      }
      else
      {
        if (f_blocks_to_get.insert(bi).second)
        {
          // Block newly added
          f_bytes_done  += f_bs;
          f_bytes_extra += f_bs - size;
        }
        else
        {
          // Block was already there ... adjust counters
          f_bytes_extra -= size;
        }
      }
    }
    else
    {
      fprintf(stderr, "Verdamt ... block range calc im CacheState::Read() ist krappen.\n");
    }
  }
}

void CacheState::EndRequest()
{
  // This is the end of current read request.
  // For single read this is rather trivial ... but we had
  // to accumulate things over subreqs of vector reads.
  // Do we need to extract some stats here?
  // It would be nice to be able to see how stuff that we loaded
  // as extra can be reused in later requests. If we indeed can ...
  // or it's all in vain. It's really non trivial ... as it's hard
  // to introduce some measure how far into file access we are.
  // Could have "checkpoints" at 25%, 50%, 75%, 100% of data read.
  // Not bad at all! But don't want to be the guy implementing this ... :)

  ++f_trips_needed;
      
  int     nbl = f_blocks_to_get.size();
  long64  nby = nbl * f_bs;
  int     ntr = nby > 0 ? (nby - 1) / (128 * 1024 * 1024) + 1 : 0;

  if (ntr > 0)
  {
    f_trips_done  += ntr;
    f_trips_extra += ntr - 1;
  }
  else
  {
    ++f_trips_saved;
  }

  for (auto i : f_blocks_to_get) f_blocks[i] = true;

  f_blocks_to_get.clear();
}

void CacheState::Finish()
{}

void CacheState::Print()
{
  printf("  Needed  %5lld   %'lld\n", f_trips_needed, f_bytes_needed);
  printf("  Done    %5lld   %'lld\n", f_trips_done,   f_bytes_done);
  printf("  Extra   %5lld   %'lld\n", f_trips_extra,  f_bytes_extra);
  printf("  Saved   %5lld   %'lld\n", f_trips_saved,  f_bytes_saved);
  printf("  Pfetchd %5lld   %'lld\n", f_trips_pref,   f_bytes_pref);
}


//==============================================================================
// ZHistos
//==============================================================================

double ZHistos::shl_vread_min = 0;

void ZHistos::BookSet(TH1 *h[4],
                      const TString &nnn, const TString &ttt,
                      int Xn, double Xl, double Xh)
{
  h[0] = new TH1D(nnn, ttt, Xn, Xl, Xh);

  h[1] = new TH2D(nnn + "_t", ttt + " vs. job progress",
                  Xn, Xl, Xh, 20, 0.025, 1.025);

  h[2] = new TH2D(nnn + "__vs__frac_read", ttt + " vs. fraction of file read",
                  Xn, Xl, Xh, 200, 0, 2);

  h[3] = new TH2D(nnn + "__vs__frac_vread", ttt + " vs. fraction of data vector-read",
                  Xn, Xl, Xh, 200, shl_vread_min, 1);

  h[4] = new TH2D(nnn + "__vs__read_rate", ttt + " vs. log10 of data read rate",
                  Xn, Xl, Xh, 200, 3, 8);

  h[5] = new TH2D(nnn + "__vs__open_duration", ttt + " vs. log10 of open duration",
                  Xn, Xl, Xh, 200, 1, 6);

  h[6] = new TH2D(nnn + "__vs__file_size", ttt + " vs. log10 of file size",
                  Xn, Xl, Xh, 200, 6, 11);
}

void ZHistos::Book(const TString &nnn, const TString &ttt)
{
  BookSet(done_req,     nnn + "done_o_req",    ttt + "Bytes transferred / requested",
          210, -0.5, 20.5);
  BookSet(donewpf_req,  nnn + "donewpf_o_req", ttt + "Bytes transferred (with prefetch) / requested",
          210, -0.5, 20.5);
  BookSet(donewpf_fs,   nnn + "donewpf_o_fs",  ttt + "Bytes transferred (with prefetch) / file-size",
          101, -0.005, 1.005);
  BookSet(saved_req,    nnn + "saved_o_req",   ttt + "Bytes saved / requested",
          101, -0.005, 1.005);
  BookSet(trpsaved_req, nnn + "trpsave_o_req", ttt + "Trips saved / requested",
          101, -0.005, 1.005);
  BookSet(extr_req,     nnn + "extra_o_req",   ttt + "Extra bytes transferred / requested",
          210, -0.5, 20.5);
  BookSet(extr_fs,      nnn + "extra_o_fs",    ttt + "Extra bytes transferred / file-size",
          101, -0.005, 1.005);
  BookSet(unused_req,   nnn + "unused_o_req",  ttt + "Unused bytes transferred / requested",
          210, -0.5, 20.5);
  BookSet(unused_fs,    nnn + "unused_o_fs",   ttt + "Unused bytes transferred / file-size",
          101, -0.005, 1.005);
  BookSet(saved_done,   nnn + "saved_o_done",  ttt + "Bytes saved / done",
          101, -0.005, 1.005);
  BookSet(extra_done,   nnn + "extra_o_done",  ttt + "Extra bytes transferred / done",
          101, -0.005, 1.005);
}

void ZHistos::FillInter(CacheState &cs, double t)
{
  done_req    [1]->Fill(cs.done_o_req(),     t);
  donewpf_req [1]->Fill(cs.gotten_o_req(),   t);
  donewpf_fs  [1]->Fill(cs.gotten_o_fs(),    t);
  saved_req   [1]->Fill(cs.saved_o_req(),    t);
  trpsaved_req[1]->Fill(cs.trpsaved_o_req(), t);
  extr_req    [1]->Fill(cs.extra_o_req(),    t);
  extr_fs     [1]->Fill(cs.extra_o_fs(),     t);
  unused_req  [1]->Fill(cs.unused_o_req(),   t);
  unused_fs   [1]->Fill(cs.unused_o_fs(),    t);
  saved_done  [1]->Fill(cs.saved_o_done(),   t);
  extra_done  [1]->Fill(cs.extra_o_done(),   t);
}

void ZHistos::Fill(CacheState &cs, double frac_read, double frac_vread, double data_rate,
                   double duration, double file_size)
{
  double oo[NN] = { 0, 0, frac_read, frac_vread, data_rate, duration, file_size };

  done_req    [0]->Fill(cs.done_o_req());
  donewpf_req [0]->Fill(cs.gotten_o_req());
  donewpf_fs  [0]->Fill(cs.gotten_o_fs());
  saved_req   [0]->Fill(cs.saved_o_req());
  trpsaved_req[0]->Fill(cs.trpsaved_o_req());
  extr_req    [0]->Fill(cs.extra_o_req());
  extr_fs     [0]->Fill(cs.extra_o_fs());
  unused_req  [0]->Fill(cs.unused_o_req());
  unused_fs   [0]->Fill(cs.unused_o_fs());
  saved_done  [0]->Fill(cs.saved_o_done());
  extra_done  [0]->Fill(cs.extra_o_done());

  for (int i = 2; i < NN; ++i)
  {
    done_req    [i]->Fill(cs.done_o_req(),     oo[i]);
    donewpf_req [i]->Fill(cs.gotten_o_req(),   oo[i]);
    donewpf_fs  [i]->Fill(cs.gotten_o_fs(),    oo[i]);
    saved_req   [i]->Fill(cs.saved_o_req(),    oo[i]);
    trpsaved_req[i]->Fill(cs.trpsaved_o_req(), oo[i]);
    extr_req    [i]->Fill(cs.extra_o_req(),    oo[i]);
    extr_fs     [i]->Fill(cs.extra_o_fs(),     oo[i]);
    unused_req  [i]->Fill(cs.unused_o_req(),   oo[i]);
    unused_fs   [i]->Fill(cs.unused_o_fs(),    oo[i]);
    saved_done  [i]->Fill(cs.saved_o_done(),   oo[i]);
    extra_done  [i]->Fill(cs.extra_o_done(),   oo[i]);
  }
}


//==============================================================================
// AnExCacheSim
//==============================================================================

AnExCacheSim::AnExCacheSim(const TString &name, AnalManager &mgr,
                           const TString &out_file,
                           double min_vread) :
  AnalExtractor(name, mgr, out_file),
  f_min_vread(min_vread)
{}

//==============================================================================

void AnExCacheSim::BookHistos()
{
  DirHolder xx0;

  OpenFile();

  ZHistos::shl_vread_min = f_min_vread;

  for (int p = 0; p < PP; ++p)
  {
    TString pdn; pdn.Form("PF_%s", nmpf[p]);

    DirHolder xx1(gDirectory->mkdir(pdn));

    for (int i = 0; i < NN; ++i)
    {
      TString dn; dn.Form("BS_%s", nmbs[i]);
      TString tt; tt.Form("PF=%s, BS=%s, ", nmpf[p], nmbs[i]);

      DirHolder xx2(gDirectory->mkdir(dn));
    
      cshi[p][i].Book(TString(), tt);
    }
  }
}

//------------------------------------------------------------------------------

void AnExCacheSim::WriteHistos()
{
  CloseFile();
}

//==============================================================================

void AnExCacheSim::Process()
{
  // Call simulate cache with various parameters.

  for (int p = 0; p < PP; ++p)
  {
    for (int i = 0; i < NN; ++i)
    {
      SimulateCache(cshi[p][i], 1024 * csbs[i], OneMB * cspf[p]);
    }
  }
}

//==============================================================================

void AnExCacheSim::SimulateCache(ZHistos &zhis, int blk_size, double pref_rate)
{
  CacheState cs(M.F.mSizeMB * OneMB, blk_size, pref_rate);

  for (vSXrdReq_i i = M.I.mReqs.begin(); i != M.I.mReqs.end(); ++i)
  {
    if (i->Type() == SXrdReq::R_Write) continue;

    cs.BeginRequest(i->Time());

    if (i->Type() == SXrdReq::R_Read)
    {
      cs.Read(i->Offset(), i->Length());
    }
    else
    {
      int sr_idx = i->SubReqIndex();

      if (sr_idx >= 0)
      {
        const int max = sr_idx + i->SubReqsStored();
        for (int si = sr_idx; si < max; ++si)
        {
          cs.Read(M.I.mOffsetVec[si], M.I.mLengthVec[si]);
        }
      }
      else
      {
        // fprintf(stderr, "Argh, sub-req details not available, this will be all krap!\n");
      }
    }

    cs.EndRequest();
    
    // Fill if frac over some threshold
    int pperc = TMath::FloorNint(20 * cs.f_prev_bytes_needed / M.F.mReadStats.mSumX / OneMB);
    int cperc = TMath::FloorNint(20 * cs.f_bytes_needed / M.F.mReadStats.mSumX / OneMB);
    if (cperc > pperc)
    {
      zhis.FillInter(cs, 0.05 * cperc);
    }
  }

  cs.Finish();


  // Fill histos

  zhis.Fill(cs,
            M.F.mReadStats.mSumX / M.F.mSizeMB,
            M.F.mVecReadStats.mSumX / M.F.mReadStats.mSumX,
            TMath::Log10(OneMB * M.F.mReadStats.mSumX / M.mDt),
            TMath::Log10(M.mDt),
            TMath::Log10(OneMB * M.F.mSizeMB));

  //------------------------------------------------------------------------------

  if (f_print_and_wait)
  {
    printf("%s\n  %s.%s   %s@%s.%s %s\n"
           "  S=%.3f MB   Nr=%lld Ns=%lld Nv=%lld Vavg=%.1f Nerr=%d\n"
           "  Frac_file %.3f - (%.3f / %.3f)   Frac_size (%.3f / %.3f)   Frac_count (%.3f / %.3f)\n",
           M.F.mName.Data(),
           M.S.mHost.Data(), M.S.mDomain.Data(),
           M.U.mName.Data(), M.U.mFromHost.Data(), M.U.mFromDomain.Data(), M.U.mRealName.Data(),
           M.F.mSizeMB,
           M.F.mReadStats.mN, M.F.mSingleReadStats.mN, M.F.mVecReadStats.mN,
           M.F.mVecReadCntStats.GetAverage(), M.I.mNErrors,
           M.F.mReadStats.mSumX / M.F.mSizeMB, 
           M.F.mSingleReadStats.mSumX / M.F.mSizeMB, M.F.mVecReadStats.mSumX / M.F.mSizeMB,
           M.F.mSingleReadStats.mSumX / M.F.mReadStats.mSumX, M.F.mVecReadStats.mSumX / M.F.mReadStats.mSumX,
           (double) M.F.mSingleReadStats.mN / M.F.mReadStats.mN, (double) M.F.mVecReadStats.mN / M.F.mReadStats.mN
    );

    cs.Print();

    printf("\n");
    char uline[256];
    fgets(uline, 256, stdin);
  }
}
