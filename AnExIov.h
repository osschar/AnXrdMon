#ifndef AnExIov_h
#define AnExIov_h

#include "AnalExtractor.h"

#include "SXrdClasses.h"

class TH1;

//==============================================================================

struct YHistos
{
  SRange   srng;
  TH1     *hcum[2], *hN[2], *hsum[2], *hmin[2], *hmax[2], *havg[2], *hsgm[2];
  Double_t fsize;
  Bool_t   is_count;

  void Book(const TString &name, const TString &title, int maxN, Double_t logB=10);
  void BookCount(const TString &name, const TString &title, int maxN, int MaxTotN);

  // Count or Fraction
  Double_t CoF(Double_t x) { return is_count ? x : x / fsize; }
  // Standard fill
  void     FillStat(TH1 *h[2], Double_t x);

  void BeginFile(Double_t file_size);
  void AddSample(Double_t x);
  void EndFile();
};

struct YYHistos
{
  YHistos *G, *H;
  TH1     *hcum[2], *hN[2], *hsum[2], *havg[2];
  Double_t fsize;

  TH1* Make2D(TH1 *g, TH1 *h);

  void Book(YHistos &x, YHistos &y);

  void FillStat(TH1 *h[2], Double_t x, Double_t y);

  void BeginFile(Double_t file_size);
  void AddSample(Double_t x, Double_t y);
  void EndFile();
};

//==============================================================================

class AnExIov : public AnalExtractor
{
  YHistos   rall, rsin, rvec;            // sin
  YHistos   prr, psrs, psrv, pvrv, pvrs; // pos offsets
  YHistos   nrr, nsrs, nsrv, nvrv, nvrs; // neg offsets
  YHistos   vn, vl, vir, vtl, vfr;       // vec

  YYHistos  rvec_vn, vir_vl;

  std::vector<YHistos*>  allYs;
  std::vector<YYHistos*> allYYs;

  // Histos for:
  // 1. cumulants for individual reqs / subreqs
  // 2. per file stats:
  //    - N, min, max, avg, sigma
  // ? !! Total backwards seek in bytes, in file size fraction

  // Can I join pos/neg stats in histos?
  // Then, there will be 2 entries per file ... hmmh.
  // Separate them first, can hack in a joint one later.

public:

  AnExIov(const TString& name, AnalManager &mgr, const TString& out_file="");

  // ----------------------------------------------------------------

  void BookHistos();

  void Process();

  void WriteHistos();

  // ----------------------------------------------------------------

  void BeginFileForYHistos();

  void EndFileForYHistos();

};

#endif
