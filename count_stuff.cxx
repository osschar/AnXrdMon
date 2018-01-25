//Quick script to get the names of all the domains used during a particular
//timeframe, in order to set up smarter domain name parsing

/*
  g++ `root-config --cflags --libs` -Wl,-rpath=. count_stuff.cxx -o count_stuff libSXrdClasses.so
*/

#include "SXrdClasses.h"

#include "TTree.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TChain.h"
#include "TSystem.h"
#include "TPRegexp.h"
#include "TMath.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRint.h"

#include <map>
#include <iostream>
#include <cassert>

//==============================================================================
// Countor
//==============================================================================

const int  N_DUMP       = 50;
const int  N_DUMP_FILES = 500;

const bool DO_FILES = false;


//==============================================================================
// Countor
//==============================================================================

struct Countor
{
  struct CountDuration
  {
    Int_t    f_count;
    Double_t f_duration;
    Double_t f_data;

    CountDuration() : f_count(0), f_duration(0), f_data(0) {}

    void inc(double dur, double dat)
    {
      ++f_count;
      f_duration += dur;
      f_data     += dat;
    }
  };

  typedef std::map<TString, CountDuration>  CDMap_t;
  typedef CDMap_t::iterator                 CDMap_i;

  CDMap_t srv_domains, cli_domains, files;

  CDMap_t d_dir, d_tier, d_dir_tier;

  CDMap_t users;

  std::vector<Long64_t> bad_hours;

  // ----------------------------------------------------------------

  TString merge(TPMERegexp& slash, int n)
  {
    TString s;
    for (int i = 1; i <= n; ++i)
    {
      if (slash[i].IsNull()) break;
      s += "/";
      s += slash[i];
    }
    return s;
  }

  void inc_dir_tier(TPMERegexp& slash, int n, double dur, double dat)
  {
    TString d = merge(slash, n);

    d_dir[d].inc(dur, dat);
    d_tier[slash[5]].inc(dur, dat);
    d_dir_tier[slash[2] + "___" + slash[5]].inc(dur, dat);
  }

  void inc_dir(TPMERegexp& slash, int n, double dur, double dat)
  {
    TString d = merge(slash, n);
    d_dir[d].inc(dur, dat);
  }

  void print_sep(char c='=')
  {
    printf("\n\n");
    for (int i = 0; i < 80; ++i)
      fputc(c, stdout);
    printf("\n\n\n");
  }

  void print_cdmap(CDMap_t& m, const TString& title, Int_t NMax=-1, Int_t NWidth=-1)
  {
    Int_t nn = m.size();

    if (nn <= 0)
    {
      printf("%s -- Nothing passed your selections\n\n", title.Data());
      return;
    }

    std::vector<TString>  names(nn);
    std::vector<Int_t>    cnts(nn);
    std::vector<Double_t> durs(nn);
    std::vector<Double_t> dats(nn);
    std::vector<Int_t>    idcs(nn);

    Int_t max_len = 0;

    {
      Int_t i = 0;
      for (CDMap_i mi = m.begin(); mi != m.end(); ++mi, ++i)
      {
        names[i] = mi->first;
        cnts [i] = mi->second.f_count;
        durs [i] = mi->second.f_duration;
        dats [i] = mi->second.f_data;

        Int_t len = mi->first.Length();
        if (len > max_len) max_len = len;
      }
    }

    if (NMax < 0 || NMax > nn) NMax   = nn;
    if (NWidth < 0)            NWidth = max_len;

    TMath::Sort(nn, &cnts[0], &idcs[0], true);
    print_sep('=');
    printf("%s -- Sort by open count\n\n", title.Data());
    for (Int_t i = 0; i < NMax; ++i)
    {
      Int_t ii = idcs[i];
      printf("%*s %8d %8.2f %8.2f\n", -NWidth, names[ii].Data(),
             cnts[ii], durs[ii], dats[ii]);
    }

    TMath::Sort(nn, &durs[0], &idcs[0], true);
    print_sep('-');
    printf("%s -- Sort by open duration [hours]\n\n", title.Data());
    for (Int_t i = 0; i < NMax; ++i)
    {
      Int_t ii = idcs[i];
      printf("%*s %8d %8.2f %8.2f\n", -NWidth, names[ii].Data(),
             cnts[ii], durs[ii], dats[ii]);
    }

    TMath::Sort(nn, &dats[0], &idcs[0], true);
    print_sep('-');
    printf("%s -- Sort by open data transffered [MB]\n\n", title.Data());
    for (Int_t i = 0; i < NMax; ++i)
    {
      Int_t ii = idcs[i];
      printf("%*s %8d %8.2f %8.2f\n", -NWidth, names[ii].Data(),
             cnts[ii], durs[ii], dats[ii]);
    }
  }
};

bool is_usa(const TString& d)
{
  // d is truncated domain, last two names only.
  // also, this is quite pitiful. sigh.

  return (d.EndsWith(".edu") || d.EndsWith(".gov") ||
          d == "ultralight.org" || d == "batlab.org" || d == "aglt2.org" ||
          d == "rr.com" || d == "amazonaws.com" || d == "akamaitechnologies.com"
         );
}


//==============================================================================
// main
//==============================================================================

int main(int argc, char *argv[])
{
  setlocale(LC_NUMERIC, "en_US");

  // Chain up several reports
  TChain mychain("XrdFar");

  // mychain.Add("/net/xrootd.t2/data/xrdmon/far/xmfar-2012-06-*.root");
  // mychain.Add("/net/xrootd.t2/data/xrdmon/far/xmfar-2014-*-*.root");
  // mychain.Add("/bar/xrdmon-far-merged/xmfar-2017-01.root");

  // mychain.Add("/bar/xrdmon-far-merged/xmfar-*.root");
  // mychain.Add("/bar/xrdmon-xxx-merged/xmxxx-*.root");
  mychain.Add("/bar/xrdmon-xxx-merged/xmxxx-2017-06.root");

  // Get set up to read domain data from the chain
  SXrdFileInfo   F, *fp = &F;
  SXrdUserInfo   U, *up = &U;
  SXrdServerInfo S, *sp = &S;

  //mychain.SetBranchStatus("*",1);
  //mychain.SetBranchAddress("S.", &sp);
  //mychain.SetBranchAddress("U.", &up);

  //mychain.SetBranchStatus("*", 0);
  //mychain.SetBranchStatus("S.mDomain",  1);
  //mychain.SetBranchStatus("U.mFromDomain", 1);

  mychain.SetBranchAddress("F.", &fp);
  mychain.SetBranchAddress("U.", &up);
  mychain.SetBranchAddress("S.", &sp);


  Countor C;

  Long64_t N = mychain.GetEntries();
  std::cout << N << " entries found.\n";

  TH1I *hp  = new TH1I("Proc",    "", 20, 0, 0);
  TH2I *hps = new TH2I("ProcSiz", "", 20, 0, 0, 20, 0, 0);
  TH1I *hd  = new TH1I("Duration",    "", 20, 0, 0);
  TH2I *hds = new TH2I("DurationSiz", "", 20, 0, 0, 20, 0, 0);

  Long64_t tot_count = 0;
  Double_t tot_size = 0, tot_size_read = 0, tot_hours = 0;

  TPMERegexp slash("/", "o");
  TPMERegexp domain_re("[^.]+\\.[^.]+$", "o");

  Long64_t acc_count = 0, rej_time = 0, rej_domain = 0, rej_pref = 0, rej_user = 0,
    rej_not_cmsrun = 0, rej_not_miniaod = 0;


  for (Long64_t i = 0; i < N; ++i)
  {
    if (i % 10000 == 0 && i > 0)
    {
      printf("  At entry %lld, so far accepted %lld\n", i, acc_count);
    }

    mychain.GetEntry(i);

    Double_t hours  = TMath::Max(1ll, F.mCloseTime - F.mOpenTime) / 3600.0;
    Double_t mbytes = F.mReadStats.mSumX;


    // ------------------------------------------------------------------------
    // Sanity checks
    // ------------------------------------------------------------------------

    const Long64_t min_open_time = 1340521949;

    if (F.mCloseTime < min_open_time || F.mOpenTime < min_open_time || hours > 240)
    {
      ++rej_time;
      C.bad_hours.push_back(i);
      continue;
    }

    // ------------------------------------------------------------------------
    // Conditions
    // ------------------------------------------------------------------------

    TString s_domain = (domain_re.Match(S.mDomain))     ? domain_re[0] : "";
    TString u_domain = (domain_re.Match(U.mFromDomain)) ? domain_re[0] : "";


    // {
    //   bool s_usa = is_usa(s_domain);
    //   bool u_usa = is_usa(u_domain);

    //   if ( ! s_usa || ! u_usa) continue;

    //   // if (s_domain == u_domain) continue;
    // }

    // debugging transfers from fnal to ral
    // if ( ! (U.mFromDomain.EndsWith("rl.ac.uk") && S.mDomain.EndsWith("fnal.gov")))
    // {
    //   ++rej_domain;
    //   continue;
    // }

    if (F.mName.BeginsWith("/store/test/") ||
        F.mName.BeginsWith("/store/user/dan/") ||
        F.mName.BeginsWith("/store/temp/") ||
        F.mName.BeginsWith("/store/mc/JobRobot/RelValProdTTbar/") ||
        F.mName.BeginsWith("/store/mc/SAM/GenericTTbar/")
        )
    {
      ++rej_pref;
      continue;
    }

    // Take MINIAOD and MINIAODSIM only
    if ( ! F.mName.Contains("/MINIAOD/") && ! F.mName.Contains("/MINIAODSIM/"))
    {
      ++rej_not_miniaod;
      continue;
    }

    // People who only do monitoring / development / testing.
    if (U.mRealName.Contains("xrootd-proxy.t2.ucsd.edu") ||
        U.mRealName.Contains("cms nanoAOD integration bot") ||
        U.mRealName.Contains("Bockelman") ||
        U.mRealName.Contains("Andrea Sciaba") ||
        U.mRealName.Contains("Vuosalo") ||
        U.mRealName.Contains("Daniel Charles Bradley") ||
        U.mRealName.Contains("Tadel") ||
        U.mRealName.IsNull()
        )
    {
      ++rej_user;
      continue;
    }

    // Require cmsRun access
    if ( ! U.mRealName.Contains("&x=cmsRun"))
    {
      ++rej_not_cmsrun;
      continue;
    }

    ++acc_count;

    // ------------------------------------------------------------------------
    // 
    // ------------------------------------------------------------------------

    slash.Split(F.mName);
    TString d = slash[2], t = slash[5];

    if (slash[1] == "store")
    {
      if (d == "data" || d == "mc" || d == "himc" || d == "hidata" || d == "generator" ||
          d == "unmerged" || d == "relval")
      {
        C.inc_dir_tier(slash, 5, hours, mbytes);
      }
      else if (d == "user" || d == "group")
      {
        C.inc_dir(slash, 3, hours, mbytes);
      }
      else
      {
        C.inc_dir(slash, 2, hours, mbytes);
      }
    }
    else
    {
      C.inc_dir(slash, 1, hours, mbytes);
    }

    int proc = F.mSizeMB > 0 ? TMath::Nint(100 * F.mReadStats.mSumX / F.mSizeMB) : 0;

    hp->Fill(proc);
    hps->Fill(proc, F.mSizeMB);
    hd->Fill(hours);
    hds->Fill(hours, F.mSizeMB);

    // printf("%8.3f %8.3f %2d%% | %5llu %6.3f %6.3f -- %5llds -- %-.30s :: %-.50s\n",
    //        F.mReadStats.mSumX, F.mSizeMB, proc,
    //        F.mReadStats.mN, F.mReadStats.GetAverage(), F.mReadStats.GetSigma(),
    //        F.mCloseTime - F.mOpenTime,
    //        U.mRealName.Data(), F.mName.Data());

    C.users[U.mRealName].inc(hours, mbytes);

    C.srv_domains[S.mDomain].inc(hours, mbytes);
    C.cli_domains[U.mFromDomain].inc(hours, mbytes);

    if (DO_FILES)
    {
      C.files[F.mName].inc(hours, mbytes);
    }

    tot_size      += F.mSizeMB;
    tot_size_read += F.mReadStats.mSumX;
    tot_hours     += hours;

    ++tot_count;
  }

  // ------------------------------------------------------------------------

  printf("acc_count = %lld, rej_time = %lld, rej_domain = %lld, rej_not_miniaod = %lld, rej_pref = %lld, rej_user = %lld, rej_not_cmsrun = %lld\n",
         acc_count, rej_time, rej_domain, rej_not_miniaod, rej_pref, rej_user, rej_not_cmsrun);
  // ------------------------------------------------------------------------

  C.print_cdmap(C.users,       "Users",            N_DUMP);

  C.print_cdmap(C.srv_domains, "Server domains",   N_DUMP);
  C.print_cdmap(C.cli_domains, "Client domains",   N_DUMP);

  C.print_cdmap(C.d_dir,       "Directories",      N_DUMP);
  C.print_cdmap(C.d_tier,      "Data-tiers",       N_DUMP);
  C.print_cdmap(C.d_dir_tier,  "Data-type___tier", N_DUMP);

  // ----------------------------------------------------------------

  if (DO_FILES)
  {
    C.print_cdmap(C.files, "Files", N_DUMP_FILES, 30);

    std::map<int, int> cnt_map;
    int sumus = 0;
    for (Countor::CDMap_i j = C.files.begin(); j != C.files.end(); ++j)
    {
      cnt_map[j->second.f_count]++;
      sumus += j->second.f_count;
    }

    C.print_sep('=');
    printf("Frequencies for multi-access:\n\n");

    int sus = 0;
    for (std::map<int,int>::iterator j = cnt_map.begin(); j != cnt_map.end(); ++j)
    {
      printf(" N_acc = %3d | Count = %6d --- %8d %f\n",
             j->first, j->second, j->first * j->second, j->first*((double)j->second)/tot_count);
      sus += j->first * j->second;
    }

    printf("\n  DO_FILES counters: sumus = %d,  sus = %d\n", sumus, sus);
  }

  printf("\nN = %lld tot_size=%'fTB, tot_size_read = %'fPB\n",
         tot_count, tot_size/1024/1024, tot_size_read/1024/1024/1024);

  // ----------------------------------------------------------------

  new TRint("count_stuff", &argc, argv);
  gApplication->InitializeGraphics();

  TCanvas *cvs = new TCanvas("Lojze", "", 1280, 1024);
  cvs->Divide(2,2);
  cvs->cd(1);
  hp->Draw();
  cvs->cd(2);
  hps->Draw("colz");
  cvs->cd(3);
  hd->Draw();
  cvs->cd(4);
  hds->Draw("colz");

  gApplication->Run();

  return 0;
}
