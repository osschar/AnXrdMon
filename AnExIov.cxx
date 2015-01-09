#include "AnExIov.h"
#include "AnalManager.h"

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

namespace
{
  const double OneMB = 1024 * 1024;

  const TString H_pre_titles[]  = { "S", "F" };

  const TString H_name_a[]  = { "full", "N", "min", "max", "avg", "sgm" };
  const TString H_name_b[]  = { "S", "F" };

  const TString H_title_a[] = { "for all requests", "N", "min", "max", "average", "sigma" };
  const TString H_title_b[] = { "in bytes", "in file size fraction" };

  const TString HL10 = "log 10 of ";

  Double_t PLog10(Double_t x)
  {
    // Positive log10, 0 is the lowest return value.

    return x < 1 ? 0 : TMath::Log10(x);
  }

  Long64_t NearestLong(double x)
  {
    // Round to nearest integer. Rounds half integers to the nearest
    // even integer.
    Long64_t i;
    if (x >= 0) {
      i = Long64_t(x + 0.5);
      if ( i & 1 && x + 0.5 == double(i) ) i--;
    } else {
      i = Long64_t(x - 0.5);
      if ( i & 1 && x - 0.5 == double(i) ) i++;
    }
    return i;
  }
}


void YHistos::Book(const TString &name, const TString &title,
                   int maxN, Double_t logB)
{
  // !!! N also goes from: 0 - 0.5 to Max + 0.5 with Nbins N + 1
  // !!! Checked with req_min which is log 10 of minimum request:
  //       it does have underflow for 0. So, either:
  //       - show underflow bin printed
  //       - add it to first bin in WriteHistos Yes, do that
  //       - fake values ?

  // Hmmh, N could also have two histos: one in absolute, another in log10 x-axis.
  // Still, it will be special.

  // Aha .. and plot also SumX ... this is important

  is_count = false;

  hcum[0] = new TH1D(name + "_cum_B", HL10 + title + " - individual reqs in bytes",
                     200, 0, logB);
  hcum[1] = new TH1D(name + "_cum_F", title + " - individual reqs in file size",
                     200, 0, 1);

  hN[0] = new TH1D(name + "_num_log", HL10 + title + " - number of entries",
                   200, 0, TMath::Ceil(TMath::Log10(maxN)));
  hN[1] = new TH1D(name + "_num",  title + " - number of entries",
                   200, 0, maxN);

  hsum[0] = new TH1D(name + "_sum_B", HL10 + title + " - sum in bytes",
                     200, 0, logB + 3);
  hsum[1] = new TH1D(name + "_sum_F", HL10 + title + " - sum in file size",
                     200, -5, 5);

  hmin[0] = new TH1D(name + "_min_B", HL10 + title + " - min in bytes",
                     200, 0, logB);
  hmin[1] = new TH1D(name + "_min_F", title + " - min in file size",
                     200, 0, 1);

  hmax[0] = new TH1D(name + "_max_B", HL10 + title + " - max in bytes",
                     200, 0, logB);
  hmax[1] = new TH1D(name + "_max_F", title + " - max in file size",
                     200, 0, 1);

  havg[0] = new TH1D(name + "_avg_B", HL10 + title + " - average in bytes",
                     200, 0, logB);
  havg[1] = new TH1D(name + "_avg_F", title + " - average in file size",
                     200, 0, 1);

  hsgm[0] = new TH1D(name + "_sgm_B", HL10 + title + " - sigma in bytes",
                     200, 0, logB);
  hsgm[1] = new TH1D(name + "_sgm_F", title + " - sigma in file size",
                     200, 0, 1);
}

void YHistos::BookCount(const TString &name, const TString &title,
                        int maxN, int maxTotN)
{
  // See comments for Book()

  is_count = true;

  Double_t maxNLog    = TMath::Ceil(TMath::Log10(maxN));
  Double_t maxTotNLog = TMath::Ceil(TMath::Log10(maxTotN));

  hcum[0] = new TH1D(name + "_cum_log", HL10 + title + " - individual counts",
                     200, 0, 3);
  hcum[1] = new TH1D(name + "_cum_N", title + " - individual counts",
                     200, 0, maxN);

  hN[0] = new TH1D(name + "_num_log", HL10 + title + " - number of entries",
                   200, 0, maxNLog);
  hN[1] = new TH1D(name + "_num",  title + " - number of entries",
                   200, 0, maxN);

  hsum[0] = new TH1D(name + "_sum_log", HL10 + title + " - sum per file",
                     200, 0, maxTotNLog);
  hsum[1] = new TH1D(name + "_sum_N", title + " - sum per file",
                     200, 0, maxTotN);

  hmin[0] = new TH1D(name + "_min_log", HL10 + title + " - min number",
                     200, 0, maxNLog);
  hmin[1] = new TH1D(name + "_min_N", title + " - min number",
                     200, 0, maxN);

  hmax[0] = new TH1D(name + "_max_log", HL10 + title + " - max number",
                     200, 0, maxNLog);
  hmax[1] = new TH1D(name + "_max_N", title + " - max number",
                     200, 0, maxN);

  havg[0] = new TH1D(name + "_avg_log", HL10 + title + " - average number",
                     200, 0, maxNLog);
  havg[1] = new TH1D(name + "_avg_N", title + " - average number",
                     200, 0, maxN);

  hsgm[0] = new TH1D(name + "_sgm_log", HL10 + title + " - sigma of counts",
                     200, 0, maxNLog);
  hsgm[1] = new TH1D(name + "_sgm_N", title + " - sigma of counts",
                     200, 0, maxN);
}

//------------------------------------------------------------------------------

void YHistos::FillStat(TH1 *h[2], Double_t x)
{
  h[0]->Fill(PLog10(x));
  h[1]->Fill(CoF(x));
}

//------------------------------------------------------------------------------

void YHistos::BeginFile(Double_t file_size)
{
  fsize = file_size;
  srng.Reset();
}

void YHistos::AddSample(Double_t x)
{
  srng.AddSample(x);

  FillStat(hcum, x);
}

void YHistos::EndFile()
{
  hN[0]->Fill(PLog10(srng.mN));
  hN[1]->Fill(srng.mN);

  hsum[0]->Fill(PLog10(srng.mSumX));
  hsum[1]->Fill(is_count ? srng.mSumX : TMath::Log10(srng.mSumX / fsize));

  FillStat(hmin, srng.mMin);
  FillStat(hmax, srng.mMax);
  FillStat(havg, srng.GetAverage());
  FillStat(hsgm, srng.GetSigma());

  fsize = 0;
}


//==============================================================================
// YYHistos
//==============================================================================

TH1* YYHistos::Make2D(TH1 *g, TH1 *h)
{
  return new TH2D
    (
      TString::Format("%s__vs__%s", g->GetName(), h->GetName()),
      TString::Format("%s  versus  %s", g->GetTitle(), h->GetTitle()),
      g->GetNbinsX(), g->GetXaxis()->GetXmin(), g->GetXaxis()->GetXmax(),
      h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()
    );
}

void YYHistos::Book(YHistos &g, YHistos &h)
{
  G = &g; H = &h;
  YHistos &X = *G, &Y = *H;

  for (int i : { 0, 1 })
  {
    hcum[i] = Make2D(X.hcum[i], Y.hcum[i]);
    hN  [i] = Make2D(X.hN  [i], Y.hN  [i]);
    hsum[i] = Make2D(X.hsum[i], Y.hsum[i]);
    havg[i] = Make2D(X.havg[i], Y.havg[i]);
  }
}

//------------------------------------------------------------------------------

void YYHistos::FillStat(TH1 *h[2], Double_t x, Double_t y)
{
  h[0]->Fill(PLog10(x), PLog10(y));
  h[1]->Fill(G->CoF(x), H->CoF(y));
}


//------------------------------------------------------------------------------

void YYHistos::BeginFile(Double_t file_size)
{
  fsize = file_size;
}

void YYHistos::AddSample(Double_t x, Double_t y)
{
  FillStat(hcum, x, y);
}

void YYHistos::EndFile()
{
  YHistos &X = *G, &Y = *H;

  hN[0]->Fill(PLog10(X.srng.mN), PLog10(Y.srng.mN));
  hN[1]->Fill(X.srng.mN, Y.srng.mN);

  hsum[0]->Fill(PLog10(X.srng.mSumX), PLog10(Y.srng.mSumX));
  hsum[1]->Fill(X.is_count ? X.srng.mSumX : TMath::Log10(X.srng.mSumX / fsize),
                Y.is_count ? Y.srng.mSumX : TMath::Log10(Y.srng.mSumX / fsize));

  FillStat(havg, X.srng.GetAverage(), Y.srng.GetAverage());

  fsize = 0;
}

//==============================================================================
// AnExIov
//==============================================================================

AnExIov::AnExIov(const TString &name, AnalManager &mgr,
                 const TString &out_file) :
  AnalExtractor(name, mgr, out_file)
{
  allYs.push_back(&rall);
  allYs.push_back(&rsin);
  allYs.push_back(&rvec);

  allYs.push_back(&prr);
  allYs.push_back(&psrs);
  allYs.push_back(&psrv);
  allYs.push_back(&pvrv);
  allYs.push_back(&pvrs);

  allYs.push_back(&nrr);
  allYs.push_back(&nsrs);
  allYs.push_back(&nsrv);
  allYs.push_back(&nvrv);
  allYs.push_back(&nvrs);

  allYs.push_back(&vn);
  allYs.push_back(&vl);
  allYs.push_back(&vir);
  allYs.push_back(&vtl);
  allYs.push_back(&vfr);


  allYYs.push_back(&rvec_vn);
  allYYs.push_back(&vir_vl);
}

//==============================================================================

// ranges:
// for cumulant: 0, 10  ... 0, 1 
// for N       : wait, N actually is special -- no Bytes / Fraction and range is different.
//               It's also diffrent for read reqs (~100) and vec-sub-reqs (~1000).
//               Probably not much sense in making it log10;

void AnExIov::BookHistos()
{
  static const int Ns   = 1000;
  static const int Nl   = 5000;
  static const int Nvec = 50000;

  DirHolder xxx;

  OpenFile();

  rall.Book("read_all", "all read request sizes",    20000);
  rsin.Book("read_sin", "single read request sizes", 20000);
  rvec.Book("read_vec", "vector read request sizes", 2000);

  prr .Book("pos_all",     "positive offsets of all read reqs", Ns);
  psrs.Book("pos_sin_sin", "positive offsets of single reads preceeded by single reads", Ns);
  psrv.Book("pos_sin_vec", "positive offsets of single reads preceeded by vector reads", Ns);
  pvrv.Book("pos_vec_vec", "positive offsets of vector reads preceeded by vector reads", Ns);
  pvrs.Book("pos_vec_sin", "positive offsets of vector reads preceeded by single reads", Ns);

  nrr .Book("neg_all",     "negative offsets of all read reqs", Ns);
  nsrs.Book("neg_sin_sin", "negative offsets of single reads preceeded by single reads", Ns);
  nsrv.Book("neg_sin_vec", "negative offsets of single reads preceeded by vector reads", Ns);
  nvrv.Book("neg_vec_vec", "negative offsets of vector reads preceeded by vector reads", Ns);
  nvrs.Book("neg_vec_sin", "negative offsets of vector reads preceeded by single reads", Ns);

  vl  .Book("vrd_extent",     "extent of all sub-reqs in vec reads", Nvec, 8);
  vir .Book("vrd_inner_offs", "offsets within vector reads", Nvec);
  vtl .Book("vrd_tot_extent", "total extent of a vec req, from start to last byte of last sub req", Ns);

  vfr .BookCount("vrd_tot_frac",   "1000 * fraction of total extent actually requested", 1000, 1000);

  vn  .BookCount("vrd_num_subreq", "number of sub-reqs in vec reads", 600, 100000);

  rvec_vn.Book(rvec, vn);
  vir_vl .Book(vir, vl);
}

//------------------------------------------------------------------------------

void AnExIov::WriteHistos()
{
  CloseFile();
}


//==============================================================================

void AnExIov::Process()
{
  // This used to be 1.0 / oneMB for deep_dump -- wanted reports in megs.
  // Here we do transformation into log or file size range.

  static const double OM = 1.0;

  const Long64_t fsize = NearestLong(M.F.mSizeMB * OneMB);

  Long64_t ro = 0, so = 0, vo = 0, wo = 0, ddo;

  bool prev_read_was_vector = false;

  BeginFileForYHistos();

  for (vSXrdReq_i i = M.I.mReqs.begin(); i != M.I.mReqs.end(); ++i)
  {
    switch (i->Type())
    {
      case SXrdReq::R_Write:
      {
        ddo = i->Offset() - wo; wo = i->Offset() + i->Length();
        // wr.AddSample(OM * ddo);
	break;
      }
      case SXrdReq::R_Read:
      {
        ddo = i->Offset() - ro;
        if ( ! prev_read_was_vector)
        {
          if (ddo >= 0) psrs.AddSample(OM * ddo); else nsrs.AddSample(-OM * ddo);
        }
        else
        {
          if (ddo >= 0) psrv.AddSample(OM * ddo); else nsrv.AddSample(-OM * ddo);
        }
        if (ddo >= 0) prr.AddSample(OM * ddo); else nrr.AddSample(-OM * ddo);

        ro = so = i->Offset() + i->Length();

        rall.AddSample(i->Length());
        rsin.AddSample(i->Length());

        prev_read_was_vector = false;
	break;
      }
      case SXrdReq::R_VecRead:
      {
	Int_t sr_idx = i->SubReqIndex();

        Int_t sum_len = 0;

	if (sr_idx >= 0)
	{
          ddo = M.I.mOffsetVec[sr_idx] - ro;
          if (prev_read_was_vector)
          {
            if (ddo >= 0) pvrv.AddSample(OM * ddo); else nvrv.AddSample(-OM * ddo);
          }
          else
          {
            if (ddo >= 0) pvrs.AddSample(OM * ddo); else nvrs.AddSample(-OM * ddo);
          }
          if (ddo >= 0) prr.AddSample(OM * ddo); else nrr.AddSample(-OM * ddo);

	  const Int_t max = sr_idx + i->SubReqsStored();
	  for (Int_t  si  = sr_idx; si < max; ++si)
	  {
            const Long64_t off = M.I.mOffsetVec[si];
            const Int_t    len = M.I.mLengthVec[si];

            ddo = off - vo; vo = off + len;

            if (si == sr_idx) {  }

            if (si >  sr_idx)
            {
              vir   .AddSample(OM * ddo);
              vir_vl.AddSample(OM * ddo, OM * len);
            }

            vl.AddSample(OM * len);

            sum_len += len;
	  }
          // Total length and fraction of total lenght actually read.
          vtl.AddSample(OM * (vo - M.I.mOffsetVec[sr_idx]));
          vfr.AddSample(1000.0 * sum_len / (vo - M.I.mOffsetVec[sr_idx]));

          rall   .AddSample(sum_len);
          rvec   .AddSample(sum_len);
          vn     .AddSample(i->SubReqsStored());
          rvec_vn.AddSample(sum_len, i->SubReqsStored());
	}

        ro = vo;

        prev_read_was_vector = true;
	break;
      }
    }

    if (ro > fsize || ro < 0)
    {
      printf("Shitty offset: Event %lld, Off = %lld (fsize=%lld)\n",
             M.GetChainI(), ro, fsize);
      return;
    }

    if (prr.srng.mSumX > nrr.srng.mSumX + fsize)
    {
      printf("Pos offsets larger than neg + fsize: Event %lld, pos = %lld, neg = %lld  (fsize=%lld)\n",
             M.GetChainI(), NearestLong(prr.srng.mSumX), NearestLong(nrr.srng.mSumX), fsize);
      return;
    }
  }

  EndFileForYHistos();
}


//==============================================================================

void AnExIov::BeginFileForYHistos()
{
  for (auto i : allYs )  { i->BeginFile(M.F.mSizeMB * OneMB); }
  for (auto i : allYYs)  { i->BeginFile(M.F.mSizeMB * OneMB); }
}

void AnExIov::EndFileForYHistos()
{
  for (auto i : allYs )  { i->EndFile(); }
  for (auto i : allYYs)  { i->EndFile(); }
}
