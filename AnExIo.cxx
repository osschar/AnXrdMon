#include "AnExIo.h"
#include "AnalManager.h"

#include "SXrdClasses.h"

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>

namespace
{
  const double OneMB = 1024 * 1024;
  const double Sec_to_Week  = 3600 * 24 * 7;
  const double Sec_to_Hour  = 3600;
}


//==============================================================================

AnExIo::AnExIo(const TString& name, AnalManager& mgr,
               const TString& out_file) :
  AnalExtractor(name, mgr, out_file)
{}

//==============================================================================

void AnExIo::AddDir(const char *name, int accum_idx)
{
  A_dirs.push_back(XDir(name, accum_idx));
}

void AnExIo::SetupAaaDirs()
{
  AddDir("all", -1);
  AddDir("all_no_grpusr", -1);
  AddDir("all_grp_usr", -1);

  AddDir("data", 1);
  AddDir("generator", 1);
  AddDir("mc", 1 );
  AddDir("unmerged", 1);
  AddDir("hidata", 1);
  AddDir("himc", 1);

  AddDir("group", 2);
  AddDir("user", 2);

  // "group", // nothing read from 2012-06 to 2013-08
  // "temp",  // nothing read in 2013-03, 29 reads from 2012-06 to 2013-08
  // "test",  // test is only Andrea Sciaba getting stuff from "/store/test/xrootd/T2_US_Wisconsin/"
}

//==============================================================================

void AnExIo::AddHisto(const char *name, const char *title, int nbx,
                             double xl, double xh, val_foo_t foo)
{
  A_histos.push_back(XHisto(name, title, nbx, xl, xh, foo));
}

void AnExIo::AddHistoCum(const char *name, const char *title, int nbx,
                              double xl, double xh, val_foo_t foo)
{
  C_histos.push_back(XHistoCum(name, title, nbx, xl, xh, foo));
  C_histos.back().f_cum.resize(M.mTotalDtHour);
}

void AnExIo::SetupAaaHistos()
{
#define XX(a) &AnExIo::a

  AddHisto("bytes_read",    "log 10 of bytes read",             140,  0, 14, XX(foo_bytes_read));
  AddHisto("data_rate",     "log 10 of avg data rate",          140,  0, 10, XX(foo_data_rate));
  AddHisto("frac_read",     "fraction of file read",            200,  0,  2, XX(foo_frac_read));
  AddHisto("frac_vread",    "fraction of bytes read via vread", 100,  0,  1, XX(foo_frac_vread));
  AddHisto("num_reqs",      "log 10 of number of requests",     160,  0,  6, XX(foo_num_reqs));
  AddHisto("num_sin_reqs",  "log 10 of number of single reads", 160,  0,  6, XX(foo_num_sin_reqs));
  AddHisto("num_vec_reqs",  "log 10 of number of vector reads", 160,  0,  4, XX(foo_num_vec_reqs));
  AddHisto("num_avg_vec_subreqs", "log 10 of average number of vector subreqs", 160,  0,  3, XX(foo_avg_vec_subreqs));
  AddHisto("open_duration", "log 10 of seconds file was open",  120,  0,  6, XX(foo_open_duration));
  AddHisto("req_avg",       "log 10 of average request size",   100,  0,  9, XX(foo_req_avg));
  AddHisto("req_max",       "log 10 of max request size",       100,  0,  9, XX(foo_req_max));
  AddHisto("req_min",       "log 10 of min request size",       100,  0,  9, XX(foo_req_min));
  AddHisto("req_rate",      "log 10 of request rate",           100, -4,  4, XX(foo_req_rate));

  // Bottom limit is respected in histogram, upper llimit is calculated.
  AddHistoCum("open_events_ph",     "log 10 of open events per hour",            200,  0,  5, XX(foo_zero));
  AddHistoCum("close_events_ph",    "log 10 of close events per hour",           200,  0,  5, XX(foo_zero));
  AddHistoCum("open_file_hours_ph", "log 10 of number of file-hours per hour",   200,  0,  5, XX(foo_zero));
  AddHistoCum("read_rate_ph",       "log 10 of read rate averaged over an hour", 200,  4, 11, XX(foo_zero));
  AddHistoCum("bytes_read_ph",      "log 10 of all data read in an hour",        200,  7, 14, XX(foo_zero));

#undef XX
}


//==============================================================================
// Init, book & write functions
//==============================================================================

void AnExIo::BookHistos()
{
  DirHolder xxx;

  OpenFile();

  N_dirs     = A_dirs.size();
  N_A_histos = A_histos.size();
  N_C_histos = C_histos.size();

  const Int_t N_weeks = TMath::CeilNint(M.mTotalDtWeek);

  for (int k = 0; k < N_dirs; ++k)
  {
    XDir &d = A_dirs[k];

    if (d.f_accum_idx > 0) f_dir_to_idx_map[d.f_name] = k;

    TDirectory *dir = mFile->mkdir(d.f_name);
    dir->cd();

    for (int i = 0; i < N_A_histos; ++i)
    {
      XHisto &h = A_histos[i];

      // printf("%2d %-14s %-34s %3d %4.0f %4.0f\n", i, h.f_name, h.f_title, h.f_nbx, h.f_xl, h.f_xh);

      d.f_hvec.push_back(new TH1D(h.f_name, h.f_title, h.f_nbx, h.f_xl, h.f_xh));

      TString nameT (h.f_name);  nameT  += "_t";
      TString titleT(h.f_title); titleT += " vs. time";
      d.f_hvec_t.push_back (new TH2D(nameT, titleT, h.f_nbx, h.f_xl, h.f_xh, N_weeks, 0, N_weeks));

      for (int j = i + 1; j < N_A_histos; ++j)
      {
        XHisto &k = A_histos[j];

        TString name2 (h.f_name);  name2  += "_vs_";  name2  += k.f_name;
        TString title2(h.f_title); title2 += " vs. "; title2 += k.f_title;
        d.f_hvec_2d.push_back (new TH2D(name2, title2, h.f_nbx, h.f_xl, h.f_xh, k.f_nbx, k.f_xl, k.f_xh));
      }
    }
  }
}

//------------------------------------------------------------------------------

void AnExIo::WriteHistos()
{
  // Convert counting vectors to top-level histos
  {
    DirHolder xxx(mFile);

    std::vector<TH1D*>   hvec(N_C_histos);
    std::vector<TGraph*> gvec(N_C_histos);
    std::vector<TH1D*>   qvec(N_C_histos);

    std::vector<double> vTime(M.mTotalDtHour);
    for (int i = 0; i < M.mTotalDtHour; ++i)
    {
      vTime[i] = i;
    }

    // Calc min, max, create histos, create & fill graphs

    TString time_axis_title;
    time_axis_title.Form("t in hours, %s - %s", M.mMinDate.Data(), M.mMaxDate.Data());

    for (int j = 0; j < N_C_histos; ++j)
    {
      XHistoCum           &xh = C_histos[j];
      std::vector<double> &v  = xh.f_cum;

      double min, max;
      min = max = v[0];
      for (int i = 1; i < M.mTotalDtHour; ++i)
      {
        if (v[i] < min) min = v[i];
        if (v[i] > max) max = v[i];
      }

      // Conclusion here is: low is zero, log options are good. Which to take depends on whether
      // we take log10 of the entry itself.

      double l_logfix = 0, l_frcfix = 0, l_cf = 0, l_medfix = 0;
      if (min > 0)
      {
        l_logfix = TMath::Power(10, TMath::Floor(TMath::Log10(min)));
        l_frcfix = 1000 * TMath::Floor(0.1 * min);
        l_cf     = TMath::Power(10, TMath::Floor(TMath::Log10(min)));
        // Round on first most significant digit
        l_medfix = l_cf * TMath::Floor(min / l_cf);
      }

      printf("Minimum for '%s' = %f, lower bound estimates: %f, %f, %f\n"
             "    this will not be applied, limit from setup %f will be kept.\n",
             xh.f_name, min, l_logfix, l_medfix, l_frcfix,xh. f_xl);

      double h_logfix = TMath::Power(10, TMath::Ceil(TMath::Log10(max)));
      double h_frcfix = 10 * TMath::Ceil(0.1 * max);
      double h_cf     = TMath::Power(10, TMath::Floor(TMath::Log10(max)));
      // Round on second most significant digit
      double h_medfix = h_cf + h_cf/10 * TMath::Ceil((max - h_cf) / h_cf * 10);

      printf("Maximum for '%s' = %f, upper bound estimates: %f, %f, %f\n",
             xh.f_name, max, h_logfix, h_medfix, h_frcfix);

      // xh.f_xl = min; // Don't touch, respect what is in config.
      xh.f_xh = h_medfix;

      hvec[j] = new TH1D(xh.f_name, xh.f_title, xh.f_nbx, xh.f_xl, TMath::Log10(xh.f_xh));

      gvec[j] = new TGraph(M.mTotalDtHour, &vTime[0], &v[0]);
      TString g_name (xh.f_name);  g_name += "_grf";
      TString g_title(xh.f_title); g_title.ReplaceAll("log 10 of ", "");
      gvec[j]->SetNameTitle(g_name, g_title);
      gvec[j]->GetXaxis()->SetTitle(time_axis_title);
      gDirectory->Add(gvec[j]);

      TString q_name (xh.f_name);  q_name += "_qhist";
      TString q_title(xh.f_title); q_title.ReplaceAll("log 10 of ", "");
      qvec[j] = new TH1D(q_name, q_title, M.mTotalDtHour, 0, M.mTotalDtHour);
      qvec[j]->GetXaxis()->SetTitle(time_axis_title);
      for (int q = 0; q < M.mTotalDtHour; ++q)
      {
        qvec[j]->SetBinContent(q+1, v[q]);
      }
    }

    // Fill histos

    for (int j = 0; j < N_C_histos; ++j)
    {
      XHistoCum           &xh = C_histos[j];
      std::vector<double> &v  = xh.f_cum;
      TH1D                *h  = hvec[j];

      for (int i = 0; i < M.mTotalDtHour; ++i)
      {
        h->Fill(v[i] > 1 ? TMath::Log10(v[i]) : xh.f_xl);
      }
    }
  }

  CloseFile();
}


//==============================================================================
// Process
//==============================================================================

void AnExIo::Process()
{
  std::vector<int> d_idcs;

  // Select what sub-dirs to fill

  d_idcs.push_back(0); // always fill "all"

  int dir_idx = -1;
  {
    auto xx = f_dir_to_idx_map.find(M.mSlashRe[2]);
    if (xx != f_dir_to_idx_map.end())
    {
      dir_idx = xx->second;
      d_idcs.push_back(A_dirs[dir_idx].f_accum_idx);
      d_idcs.push_back(dir_idx);
    }
  }

  // Now fill the histos

  double week = (M.F.mOpenTime - M.mMinT) / Sec_to_Week;

  int i2d = 0; // index for 2d histos; we fill in same order as during creation

  for (int i = 0; i < N_A_histos; ++i)
  {
    XHisto &xh = A_histos[i];

    double val = (this->*xh.f_foo)();

    for (unsigned int k = 0; k < d_idcs.size(); ++k)
    {
      XDir &d = A_dirs[d_idcs[k]];

      d.f_hvec[i]  ->Fill(val);
      d.f_hvec_t[i]->Fill(val, week);
    }

    for (int j = i + 1; j < N_A_histos; ++j)
    {
      XHisto &xk = A_histos[j];

      double kval = (this->*xk.f_foo)();

      for (unsigned int k = 0; k < d_idcs.size(); ++k)
      {
        XDir &d = A_dirs[d_idcs[k]];

        d.f_hvec_2d[i2d]->Fill(val, kval);
      }

      ++i2d;
    }
  }

  // Per hour statistics, cumulative

  // Nasty hack ?
  std::vector<double> &cvOpenEvents    = C_histos[0].f_cum;
  std::vector<double> &cvCloseEvents   = C_histos[1].f_cum;
  std::vector<double> &cvOpenFileHours = C_histos[2].f_cum;
  std::vector<double> &cvReadRate      = C_histos[3].f_cum;
  std::vector<double> &cvBytesRead     = C_histos[4].f_cum;

  double open_hour  = (M.F.mOpenTime - M.mMinT) / Sec_to_Hour;
  double duration   = M.mDt / Sec_to_Hour;
  double close_hour = open_hour + duration;

  //printf("M.F.mOpenTime = %lld, open_hour = %f     oe = %p\n", M.F.mOpenTime, open_hour, & cvOpenEvents[0]);
  //fflush(stdout);

  ++cvOpenEvents [(int) open_hour ];
  ++cvCloseEvents[(int) close_hour];

  double yebo;
  double time_left     = duration;
  double hour_frac_max = TMath::Min(1 - modf(open_hour, &yebo), duration);
  int    hour          = (int) open_hour;
  while (time_left > 0)
  {
    double hour_frac = TMath::Min(hour_frac_max, time_left);

    cvOpenFileHours[hour] += hour_frac;
    cvReadRate     [hour] += OneMB * M.F.mReadStats.mSumX * hour_frac / duration / 3600;
    cvBytesRead    [hour] += OneMB * M.F.mReadStats.mSumX * hour_frac / duration;

    time_left    -= hour_frac;
    hour_frac_max = 1;
    ++hour;
  }
}



//==============================================================================
// Value extracting functions
//==============================================================================

double AnExIo::foo_frac_read()     { return M.F.mReadStats.mSumX / M.F.mSizeMB; }
double AnExIo::foo_frac_vread()    { return M.F.mVecReadStats.mSumX / M.F.mReadStats.mSumX; }
double AnExIo::foo_bytes_read()    { return TMath::Log10(OneMB * M.F.mReadStats.mSumX); }
double AnExIo::foo_open_duration() { return TMath::Log10(M.mDt); }
double AnExIo::foo_num_reqs()      { return TMath::Log10(M.F.mReadStats.mN); }
double AnExIo::foo_num_sin_reqs()  { return TMath::Log10(M.F.mSingleReadStats.mN); }
double AnExIo::foo_num_vec_reqs()  { return TMath::Log10(M.F.mVecReadStats.mN); }
double AnExIo::foo_avg_vec_subreqs() { return TMath::Log10(M.F.mVecReadCntStats.GetAverage()); }
double AnExIo::foo_req_avg()       { return TMath::Log10(OneMB * M.F.mReadStats.GetAverage()); }
double AnExIo::foo_req_min()       { return TMath::Log10(OneMB * M.F.mReadStats.mMin); }
double AnExIo::foo_req_max()       { return TMath::Log10(OneMB * M.F.mReadStats.mMax); }
double AnExIo::foo_req_rate()      { return TMath::Log10(M.F.mReadStats.mN / M.mDt); }
double AnExIo::foo_data_rate()     { return TMath::Log10(OneMB * M.F.mReadStats.mSumX / M.mDt); }

double AnExIo::foo_zero()          { return 0; }
