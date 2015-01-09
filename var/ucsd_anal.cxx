#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TPRegexp.h"

#include "TEntryList.h"
#include "TTree.h"
#include "TChain.h"

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

#include "TCanvas.h"

#include "SXrdClasses.h"


void deep_dump();


//==============================================================================
// Notes
//==============================================================================
//
// F.mReadStats.mSumX == 0 or even F.mSizeMB == 0
// XrdFar->Draw("F.mReadStats.mSumX", "F.mReadStats.mSumX/F.mSizeMB < 0.01 && ! F.mName.BeginsWith(\"/store/test/\") && F.mName != \"/store/user/dan/test_file\"")
// XrdFar->Scan("F.mReadStats.mSumX:F.mSizeMB:F.mName", "F.mSizeMB == 0", "col=.7f:.7f:180s")

// What's with overreading? Filter it out?
//
// Brian / JobRobot is not here because of wisc-wisc restriction

// min open-time = 1340521900
// max open-time = 1378018800 (also close time)
// delta = 37496900 s = 624949 min = 10416 h = 433 days = 61.86 weeks = 14.2 months

// All records = 107,586,853
//
//            run-1        run-2
// W->W  56,528,302   17,108,892
// W->X   1,009,351      694,696
// X->W     855,761      702,403
//

//==============================================================================
// Config
//==============================================================================

enum EAccessType
{
  AT_any,       // all records pass
  AT_local,     // server and client at _local_domain_
  AT_nonlocal,  // server _local_domain_, client elsewhere
  AT_remote,    // server elsewhere, client at _local_domain_
  AT_US,        // remote access in US, .edu + .gov
};

enum ETimeCut
{
  TC_none,
  TC_keep_long,
  TC_keep_short
};

#ifndef __CINT__

struct XConfig
{
  const char           *f_name;
  const char           *f_out_file_name;
  const char           *f_local_domain;
  const EAccessType     f_access_type;
  const ETimeCut        f_time_cut;
  const bool            f_cut_lazy_download;
};

const XConfig A_Configs[] =
{
  { "r1_full_loc",    "uw1/full_local.root",    ".wisc.edu", AT_local,    TC_none, false }, // W->W  56,528,302
  { "r1_full_nonloc", "uw1/full_nonlocal.root", ".wisc.edu", AT_nonlocal, TC_none, false }, // W->X   1,009,351
  { "r1_full_remote", "uw1/full_remote.root",   ".wisc.edu", AT_remote,   TC_none, false }, // X->W     855,761

  { "r2_full_loc",    "uw2/full_local.root",    ".wisc.edu", AT_local,    TC_keep_long,  false }, // W->W  17,108,892
  { "r2_full_nonloc", "uw2/full_nonlocal.root", ".wisc.edu", AT_nonlocal, TC_keep_long,  false }, // W->X     694,696
  { "r2_full_remote", "uw2/full_remote.root",   ".wisc.edu", AT_remote,   TC_keep_long,  false }, // X->W     702,403

  { "us1_all_n_lazy", "us1/all_n_lazy.root", ".wisc.edu", AT_US, TC_none,       false }, // all w/ preload   6,310,605
  { "us1_short",      "us1/short.root",      ".wisc.edu", AT_US, TC_keep_short, true  }, // short, t < 100s    703,793
  { "us1_long",       "us1/long.root",       ".wisc.edu", AT_US, TC_keep_long,  true  }, // long,  t > 100s  4,534,612

  { "uw3_all_n_lazy", "uw3/all_n_lazy.root", ".wisc.edu", AT_local, TC_none,       false }, // all w/ preload   56,528,278
  { "uw3_short",      "uw3/short.root",      ".wisc.edu", AT_local, TC_keep_short, true  }, // short, t < 100s  35,672,599
  { "uw3_long",       "uw3/long.root",       ".wisc.edu", AT_local, TC_keep_long,  true  }, // long,  t > 100s   7,839,134


  { "ucsd_test", "ucsd_test.root", ".ucsd.edu", AT_any, TC_none, false },
};

const int N_A_Configs = sizeof(A_Configs) / sizeof(XConfig);

const XConfig* theConfig = 0;

#endif


//==============================================================================
// Globals
//==============================================================================

TTree  *theTree    = 0;

TFile  *theOutFile = 0;

SXrdFileInfo   F, *_fp = &F;
SXrdUserInfo   U, *_up = &U;
SXrdServerInfo S, *_sp = &S;
SXrdIoInfo     I, *_ip = &I;

// Counting vectors per hour
const Int_t N_hours = 10416 + 1;
std::vector<double> cvOpenEvents(N_hours);
std::vector<double> cvCloseEvents(N_hours);
std::vector<double> cvOpenFileHours(N_hours);
std::vector<double> cvBytesRead(N_hours);

#ifndef __CINT__
std::vector<double> *cvVecs[] =
{
  &cvOpenEvents, &cvCloseEvents, &cvOpenFileHours, &cvBytesRead
};
#endif

// ------------------------------------------------------------------------

struct XDir
{
  const char *f_name;
  int         f_accum_idx;

  std::vector<TH1*> f_hvec, f_hvec_t;
};

XDir A_dirs[] =
{
  { "all",  -1 }, { "all_no_grpusr", -1 }, { "all_grp_usr", -1 },
  { "data",  1 }, { "generator", 1 }, { "mc", 1 }, { "results", 1 }, { "unmerged", 1 },
  { "group", 2 }, { "user", 2 }
  // "group", // nothing read from 2012-06 to 2013-08
  // "temp",  // nothing read in 2013-03, 29 reads from 2012-06 to 2013-08
  // "test",  // test is only Andrea Sciaba getting stuff from "/store/test/xrootd/T2_US_Wisconsin/"
};

const int N_dirs = sizeof(A_dirs) / sizeof(XDir);

std::map<TString, int> M_dir_to_idx;

// ------------------------------------------------------------------------

const double oneMB = 1024 * 1024;

double delta_t;

double foo_frac_read()     { return F.mReadStats.mSumX / F.mSizeMB; }
double foo_frac_vread()    { return F.mVecReadStats.mSumX / F.mReadStats.mSumX; }
double foo_bytes_read()    { return TMath::Log10(oneMB * F.mReadStats.mSumX); }
double foo_open_duration() { return TMath::Log10(delta_t); }
double foo_num_reqs()      { return TMath::Log10(F.mReadStats.mN); }
double foo_req_avg()       { return TMath::Log10(oneMB * F.mReadStats.GetAverage()); }
double foo_req_min()       { return TMath::Log10(oneMB * F.mReadStats.mMin); }
double foo_req_max()       { return TMath::Log10(oneMB * F.mReadStats.mMax); }
double foo_req_rate()      { return TMath::Log10(F.mReadStats.mN / delta_t); }
double foo_data_rate()     { return TMath::Log10(oneMB * F.mReadStats.mSumX / delta_t); }

double foo_zero()          { return 0; }

// ------------------------------------------------------------------------

typedef double (*valfoo)();

struct XHisto
{
  const char *f_name, *f_title;
  int         f_nbx;
  double      f_xl, f_xh;
  valfoo      f_foo; 
};

XHisto A_histos[] =
{
  { "frac_read",     "fraction of file read",            200,  0,  2, foo_frac_read },
  { "frac_vread",    "fraction of bytes read via vread", 100,  0,  1, foo_frac_vread },
  { "bytes_read",    "log 10 of bytes read",             140,  0, 14, foo_bytes_read },
  { "open_duration", "log 10 of seconds file was open",  120,  0,  6, foo_open_duration },
  { "num_reqs",      "log 10 of number of requests",     160,  0,  8, foo_num_reqs},
  { "req_avg",       "log 10 of average request size",   100,  0, 10, foo_req_avg },
  { "req_min",       "log 10 of min request size",       100,  0, 10, foo_req_min },
  { "req_max",       "log 10 of max request size",       100,  0, 10, foo_req_max },
  { "req_rate",      "log 10 of request rate",           100, -5,  5, foo_req_rate },
  { "data_rate",     "log 10 of read rate",              140, -4, 10, foo_data_rate },
};

const int N_A_histos = sizeof(A_histos) / sizeof(XHisto);

XHisto C_histos[] =
{
  { "open_events_ph",     "log 10 of open events per hour",            200,  0, 0, foo_zero },
  { "close_events_ph",    "log 10 of close events per hour",           200,  0, 0, foo_zero },
  { "open_file_hours_ph", "log 10 of number of file-hours per hour",   200,  0, 0, foo_zero },
  { "bytes_read_ph",      "log 10 of read rate averaged over an hour", 200,  3, 0, foo_zero },
};

const int N_C_histos = sizeof(C_histos) / sizeof(XHisto);


//==============================================================================
// Init, book & write functions
//==============================================================================

void init_input(int argc, char **argv)
{
  TChain *t = new TChain(argv[0]);
  for (int i = 1; i < argc; ++i)
  {
    t->Add(argv[i]);
  }
  theTree = t;

  theTree->SetBranchAddress("F.", &_fp);
  theTree->SetBranchAddress("U.", &_up);
  theTree->SetBranchAddress("S.", &_sp);
  theTree->SetBranchAddress("I.", &_ip);
}

//------------------------------------------------------------------------------

void book_histos()
{
  TDirectory *xxx = gDirectory;

  theOutFile = TFile::Open(theConfig->f_out_file_name, "create");
  if ( ! theOutFile)
  {
    fprintf(stderr, "Opening of output file '%s' failed. Probably it exists already.\n",
            theConfig->f_out_file_name);
    exit(2);
  }

  for (int k = 0; k < N_dirs; ++k)
  {
    XDir &d = A_dirs[k];

    if (k > 2) M_dir_to_idx[d.f_name] = k;

    TDirectory *dir = theOutFile->mkdir(d.f_name);
    dir->cd();

    for (int i = 0; i < N_A_histos; ++i)
    {
      XHisto &h = A_histos[i];

      // printf("%2d %-14s %-34s %3d %4.0f %4.0f\n", i, h.f_name, h.f_title, h.f_nbx, h.f_xl, h.f_xh);

      d.f_hvec.push_back(new TH1D(h.f_name, h.f_title, h.f_nbx, h.f_xl, h.f_xh));

      TString name2 (h.f_name);  name2  += "_t";
      TString title2(h.f_title); title2 += " vs. time";
      d.f_hvec_t.push_back (new TH2D(name2, title2, h.f_nbx, h.f_xl, h.f_xh, 62, 0, 62));
    }
  }

  xxx->cd();
}

//------------------------------------------------------------------------------

void write_histos()
{
  // Convert counting vectors to top-level histos
  {
    TDirectory *xxx = gDirectory;

    theOutFile->cd();

    std::vector<TH1D*>   hvec(N_C_histos);
    std::vector<TGraph*> gvec(N_C_histos);
    std::vector<TH1D*>   qvec(N_C_histos);

    std::vector<double> vTime(N_hours);
    for (int i = 0; i < N_hours; ++i)
    {
      vTime[i] = i;
    }

    // Calc min, max, create histos, create & fill graphs

    for (int j = 0; j < N_C_histos; ++j)
    {
      std::vector<double> &v  = * cvVecs[j];
      XHisto              &xh = C_histos[j];

      double min, max;
      min = max = v[0];
      for (int i = 1; i < N_hours; ++i)
      {
        if (v[i] < min) min = v[i];
        if (v[i] > max) max = v[i];
      }

      // Conclusion here is: low is zero, log options are good. Which to take depends on whether
      // we take log10 of the entry itself.

      double l_logfix = TMath::Power(10, TMath::Floor(TMath::Log10(min)));
      double l_frcfix = 1000 * TMath::Floor(0.1 * min);
      double l_cf = TMath::Power(10, TMath::Floor(TMath::Log10(min)));
      double l_medfix = l_cf * TMath::Floor(min / l_cf);

      printf("Minimum for '%s' = %f, lower bound estimates: %f, %f, %f\n"
             "    this will not be applied, limit from setup %f will be kept.\n",
             xh.f_name, min, l_logfix, l_medfix, l_frcfix,xh. f_xl);

      double h_logfix = TMath::Power(10, TMath::Ceil(TMath::Log10(max)));
      double h_frcfix = 10 * TMath::Ceil(0.1 * max);

      double h_cf = TMath::Power(10, TMath::Floor(TMath::Log10(max)));
      double h_medfix = h_cf * TMath::Ceil(max / h_cf);

      printf("Maximum for '%s' = %f, upper bound estimates: %f, %f, %f\n",
             xh.f_name, max, h_logfix, h_medfix, h_frcfix);

      // xh.f_xl = min; // Don't touch, respect what is in config.
      xh.f_xh = h_medfix;

      hvec[j] = new TH1D(xh.f_name, xh.f_title, xh.f_nbx, xh.f_xl, TMath::Log10(xh.f_xh));

      gvec[j] = new TGraph(N_hours, &vTime[0], &v[0]);
      TString g_name (xh.f_name);  g_name += "_grf";
      TString g_title(xh.f_title); g_title.ReplaceAll("log 10 of ", "");
      gvec[j]->SetNameTitle(g_name, g_title);
      gvec[j]->GetXaxis()->SetTitle("t in hours, June 2012 - Aug 2013");
      gDirectory->Add(gvec[j]);

      TString q_name (xh.f_name);  q_name += "_qhist";
      TString q_title(xh.f_title); q_title.ReplaceAll("log 10 of ", "");
      qvec[j] = new TH1D(q_name, q_title, N_hours, 0, N_hours);
      for (int q = 0; q < N_hours; ++q)
      {
        qvec[j]->SetBinContent(q+1, v[q]);
      }
    }

    // Fill histos

    for (int j = 0; j < N_C_histos; ++j)
    {
      std::vector<double> &v  = * cvVecs[j];
      XHisto              &xh = C_histos[j];
      TH1D                *h  = hvec[j];

      for (int i = 0; i < N_hours; ++i)
      {
        h->Fill(v[i] > 1 ? TMath::Log10(v[i]) : xh.f_xl);
      }
    }

    xxx->cd();
  }

  theOutFile->Write();
  theOutFile->Close();
  delete theOutFile;
  theOutFile = 0;
}


//==============================================================================
// Fill & loop functions
//==============================================================================

std::vector<int> d_idcs;

const Long64_t min_opentime = 1340521900;
const double   sec_to_week  = 3600 * 24 * 7;
const double   sec_to_hour  = 3600;

void fill_histos()
{
   double week = (F.mOpenTime - min_opentime) / sec_to_week;

   for (int i = 0; i < N_A_histos; ++i)
   {
     XHisto &xh = A_histos[i];

     double val = xh.f_foo();

     for (unsigned int k = 0; k < d_idcs.size(); ++k)
     {
       XDir &d = A_dirs[d_idcs[k]];

       d.f_hvec[i]  ->Fill(val);
       d.f_hvec_t[i]->Fill(val, week);
     }
   }

   // Per hour statistics, cumulative

   double open_hour  = (F.mOpenTime - min_opentime) / sec_to_hour;
   double duration   = delta_t / sec_to_hour;
   double close_hour = open_hour + duration;

   //printf("F.mOpenTime = %lld, open_hour = %f     oe = %p\n", F.mOpenTime, open_hour, & cvOpenEvents[0]);
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
     cvBytesRead    [hour] += oneMB * F.mReadStats.mSumX * hour_frac / duration / 3600;

     time_left    -= hour_frac;
     hour_frac_max = 1;
     ++hour;
   }
}

//------------------------------------------------------------------------------

void loop_etc()
{
  Long64_t N = theTree->GetEntries(), N_pass = 0;
  printf("loop_etc(), going over %lld entries ...\n", N);

  TPMERegexp slash_re("/", "o");

  TPMERegexp domain_re("[^.]+\\.[^.]+$", "o");

  const Int_t NDiv = TMath::Power(10, TMath::Floor(TMath::Log10(N) - 4));

  // Crash at 76605732, valgrind error at 69.07%, ~74300000
  for (Long64_t i = 0; i < N; ++i)
  {
    // Progress report
    if (i % NDiv == 0)
    {
      // printf("%lld ", i);
      // printf("\x1b[2K\x1b[0E\x1b[31mProgress: %5.2f%%\x1b[0m", 100*(double)i/N);
      fflush(stdout);
    }

    theTree->GetEntry(i);

    switch (theConfig->f_access_type)
    {
      case AT_any:
        break;
      case AT_local:
        if ( ! S.mDomain.EndsWith(theConfig->f_local_domain) || ! U.mFromDomain.EndsWith(theConfig->f_local_domain))
          continue;
        break;
      case AT_nonlocal:
        if ( ! S.mDomain.EndsWith(theConfig->f_local_domain) || U.mFromDomain.EndsWith(theConfig->f_local_domain))
          continue;
        break;
      case AT_remote:
        if ( S.mDomain.EndsWith(theConfig->f_local_domain) || ! U.mFromDomain.EndsWith(theConfig->f_local_domain))
          continue;
        break;
      case AT_US:
        if (( ! S.mDomain    .EndsWith(".edu") && ! S.mDomain    .EndsWith(".gov")) ||
            ( ! U.mFromDomain.EndsWith(".edu") && ! U.mFromDomain.EndsWith(".gov")))
        {
          continue;
        }
        else
        {
          TString sdom, udom;
          if (domain_re.Match(S.mDomain))     sdom = domain_re[0]; else continue;
          if (domain_re.Match(U.mFromDomain)) udom = domain_re[0]; else continue;
          if (sdom == udom) continue;
        }
        break;
      default:
        exit(2);
    }
    if (F.mName.BeginsWith("/store/test/") || F.mName.BeginsWith("/store/user/dan/") ||
        F.mName.BeginsWith("/store/temp/") ||
        F.mName.BeginsWith("/store/mc/JobRobot/RelValProdTTbar/") ||
        F.mName.BeginsWith("/store/mc/SAM/GenericTTbar/"))
    {
      continue;
    }
    // Filter out Alja and me.
    if (U.mRealName.Contains("Tadel"))
    {
      continue;
    }
    // Filter out events with a one single read request only. Brian and some other tests
    if (F.mReadStats.mN <= 1 && F.mReadStats.mN == F.mSingleReadStats.mN)
    {
      continue;
    }

    delta_t = TMath::Max(1ll, F.mCloseTime - F.mOpenTime);

    // Pathological time open / close / duration.
    if (F.mCloseTime < min_opentime || F.mOpenTime < min_opentime || delta_t > 1e6)
    {
      printf("\x1b[2K\x1b[0EYEBO Event=%lld CloseTime=%lld OpenTime=%lld delta_t=%f ... skipping.\n",
             i, F.mCloseTime, F.mOpenTime, delta_t);
      continue;
    }

    // Filter out small and short accesses
    // Q. Should be small and short  or  small or short ? How correlated is this?
    // A. Just short, weird things happen there, analyze separately.
    switch (theConfig->f_time_cut)
    {
      case TC_none:       break;
      case TC_keep_long:  if (delta_t <= 100) continue; break;
      case TC_keep_short: if (delta_t >  100) continue; break;
      default: exit(3);
    }

    // Filter out lazy-preload ... is there additional thing to check?
    // XXXX Collect who does that and from, to domain pairs.
    // XXXX Also, check distribution of fraction of file read.
    if (theConfig->f_cut_lazy_download && F.mReadStats.mMax == 128)
    {
      continue;
    }


    // -------------------
    // End of top cuts ...
    // -------------------

    ++N_pass;

    slash_re.Split(F.mName);
    TString dname = slash_re[2], tname = slash_re[5];

    // Now fill the histos

    d_idcs.push_back(0); // always fill "all"

    int dir_idx = -1;
    {
      std::map<TString, int>::iterator xx = M_dir_to_idx.find(dname);
      if (xx != M_dir_to_idx.end())
      {
        dir_idx = xx->second;
        d_idcs.push_back(A_dirs[dir_idx].f_accum_idx);
        d_idcs.push_back(dir_idx);
      }
    }

    fill_histos();

    d_idcs.clear();

    deep_dump();
  }

  printf("\x1b[2K\x1b[0E... done! Pass count = %lld\n", N_pass);
}


//==============================================================================
// main
//==============================================================================

int main(int argc, char **argv)
{
  if (argc < 4)
  {
    fprintf(stderr, "Usage: %s config-name tree-name file-name-or-pattern ...\n", argv[0]);
    exit(1);
  }
  for (int i = 0; i < N_A_Configs; ++i)
  {
    if (strcmp(A_Configs[i].f_name, argv[1]) == 0)
    {
      theConfig = & A_Configs[i];
      break;
    }
  }
  if (theConfig == 0)
  {
    fprintf(stderr, "Configuration with name '%s' not found, check the code. Exiting.\n", argv[1]);
    exit(1);
  }

  gSystem->Load("libSXrdClasses.so");


  init_input(argc - 2, & argv[2]);

  book_histos();

  loop_etc();

  write_histos();

  return 0;
}
