#include "AnalManager.h"

// Needed for init functions ... should eventually go elsewhere
#include "AnExIo.h"
#include "AnExIov.h"
#include "AnExCacheSim.h"

#include <TChain.h>
#include <TMath.h>
#include <TSystem.h>
#include <TH1.h>
#include <TDatime.h>

//==============================================================================

AnalManager::AnalManager(const TString& name,      const TString& out_dir,
                         const TString& tree_name, const TString& pfx,
                         bool setup_I_branch) :
  AnalFilter(name, *this),
  mChn(0), mChnN(-1), mChnI(-1),
  mInFilePrefix(pfx),
  mOutDirName(out_dir),
  _fp(&F), _up(&U), _sp(&S), _ip(&I),
  mBranchIActive(setup_I_branch),
  mSDomainRe("[^.]+\\.[^.]+$", "o"),
  mUDomainRe("[^.]+\\.[^.]+$", "o"),
  mSlashRe("/", "o")
{
  // Make sure out-dir does not exist and then create it.
  if (gSystem->AccessPathName(mOutDirName) == false)
  {
    fprintf(stderr, "Output directory '%s' already exists. Cowardly refusing to use it :(\n",
            mOutDirName.Data());
    exit(1);
  }
  if (gSystem->mkdir(mOutDirName, true) == -1)
  {
    fprintf(stderr, "Creation of output directory '%s' failed. Dying ...\n");
    exit(1);
  }

  mOnTty = isatty(fileno(stdout));

  mChn = new TChain(tree_name);

  // ?? Can this be done before adding of the files ?
  mChn->SetBranchAddress("F.", &_fp);
  mChn->SetBranchAddress("U.", &_up);
  mChn->SetBranchAddress("S.", &_sp);
  if (mBranchIActive)
    mChn->SetBranchAddress("I.", &_ip);
}

//------------------------------------------------------------------------------

AnalManager::~AnalManager()
{
  // In principle should delete all prefilters, filters and extractors.
}

//==============================================================================

void AnalManager::AddFile(const TString& files)
{
  if ( ! mAnalExs.empty())
  {
    fprintf(stderr, "Adding files after analyses have been registered is wrong! Dying ...\n");
    exit(1);
  }

  mChn->Add(mInFilePrefix + files);
}

void AnalManager::AddPreFilter(AnalFilter* flt)
{
  mPreFilters.push_back(flt);
}

void AnalManager::AddExtractor(AnalExtractor* ext)
{
  mAnalExs.push_back(ext);
  for (auto f : ext->mFilters)     mAnalFis.insert(f);
  for (auto f : ext->mAntiFilters) mAnalFis.insert(f);
}

//==============================================================================

void AnalManager::ScanEdgeTimes(Long64_t scan_entries)
{
  printf("AnalManager::ScanEdgeTimes entered ...\n");

  Long64_t N = mChn->GetEntries();

  printf("Chain has %lld entries.\n", N);

  if (scan_entries > N) scan_entries = N;

  Long64_t mo, mc, Mo, Mc;
  Long64_t imo, imc, iMo, iMc;

  mChn->GetEntry(0);
  mo = Mo = F.mOpenTime;
  mc = Mc = F.mCloseTime;
  imo = imc = iMo = iMc = 0;

  for (Long64_t i = 1; i < scan_entries; ++i)
  {
    mChn->GetEntry(i);
    if (F.mOpenTime  < mo) { mo = F.mOpenTime;  imo = i; }
    if (F.mCloseTime < mc) { mc = F.mCloseTime; imc = i; }
  }

  for (Long64_t i = N - scan_entries; i < N; ++i)
  {
    mChn->GetEntry(i);
    if (F.mOpenTime  > Mo) { Mo = F.mOpenTime;  iMo = i; }
    if (F.mCloseTime > Mc) { Mc = F.mCloseTime; iMc = i; }
  }

  printf("Min time open %lld (%lld), close %lld (%lld)\n", mo, imo, mc, imc);
  printf("Max time open %lld (%lld), close %lld (%lld)\n", Mo, iMo, Mc, iMc);

  SetEdgeTimes(TMath::Min(mo, mc), TMath::Max(Mo, Mc), true);

  printf("  // mgr.SetEdgeTimes(%lld, %lld);\n", mMinT, mMaxT);

  printf("AnalManager::ScanEdgeTimes finished.\n");
}

void AnalManager::SetEdgeTimes(Long64_t min, Long64_t max, bool verbose)
{
  printf("Time Min %lld -- %lld Max\n", min, max);

  // round to nearest hour
  mMinT = min - min % 3600;
  mMaxT = max + (3600 - max % 3600);

  printf("  Rounded Min %lld -- %lld Max\n", mMinT, mMaxT);

  mTotalDtSec  = mMaxT - mMinT;
  mTotalDtMin  = mTotalDtSec  / 60;
  mTotalDtHour = mTotalDtMin  / 60;
  mTotalDtDay  = mTotalDtHour / 24.0;
  mTotalDtWeek = mTotalDtDay  / 7;
  mTotalDtMonth= mTotalDtDay  / 30.4375; // = 365.25 / 12

  TDatime dt;
  dt.Set(mMinT);
  mMinDate.Form("%d-%02d-%02d %02d:%02d PST",
                dt.GetYear(), dt.GetMonth(), dt.GetDay(), dt.GetHour(), dt.GetMinute());
  dt.Set(mMaxT);
  mMaxDate.Form("%d-%02d-%02d %02d:%02d PST",
                dt.GetYear(), dt.GetMonth(), dt.GetDay(), dt.GetHour(), dt.GetMinute());

  printf("  Duration = %lld s, %lld min, %lld hours,\n"
         "             %.2f days, %.2f weeks, %.2f months\n"
         "  From %s to %s\n",
         mTotalDtSec, mTotalDtMin,  mTotalDtHour,
         mTotalDtDay, mTotalDtWeek, mTotalDtMonth,
         mMinDate.Data(), mMaxDate.Data());
}

//==============================================================================

bool AnalManager::Filter()
{
  // Extract commonly used data & filter out crap

  mDt = TMath::Max(1ll, F.mCloseTime - F.mOpenTime);

  // Pathological time open / close / duration.
  if (F.mCloseTime < mMinT || F.mOpenTime < mMinT || mDt > 1e6)
  {
    printf("YEBO Event=%lld CloseTime=%lld OpenTime=%lld delta_t=%f ... skipping.\n",
           mChnI, F.mCloseTime, F.mOpenTime, mDt);
    return false;
  }

  mSDomain = (mSDomainRe.Match(S.mDomain))     ? mSDomainRe[0] : "";
  mUDomain = (mUDomainRe.Match(U.mFromDomain)) ? mUDomainRe[0] : "";

  mSlashRe.Split(F.mName);

  for (auto flt : mPreFilters)
  {
    if ( ! flt->Filter())  return false;
  }

  return true;
}

//==============================================================================

void AnalManager::Process()
{
  for (auto ext : mAnalExs) ext->BookHistos();

  mChnN = mChn->GetEntries();

  const Int_t NDiv = TMath::Power(10, TMath::Floor(TMath::Log10(mChnN) - 4));

  printf("AnalManager::Process(), going over %lld entries ...\n", mChnN);

  for (mChnI = 0; mChnI < mChnN; ++mChnI)
  {
    // Progress report
    if (mChnI % NDiv == 0)
    {
      // printf("%lld ", mChnI);
      if (mOnTty)
      {
        printf("\x1b[2K\x1b[31mProgress: %5.2f%%\x1b[0m\x1b[0E", 100*(double)mChnI/mChnN);
        fflush(stdout);
      }
    }

    mChn->GetEntry(mChnI);


    // Extract commonly used data & filter out crap
    if ( ! FilterAndStore())
    {
      continue;
    }

    // Call filters
    for (auto flt : mAnalFis)
    {
      flt->FilterAndStore();
    }

    // Call extractors
    for (auto ext : mAnalExs)
    {
      if (ext->FilterAndStore())
      {
        ext->Process();
      }
    }
  }

  printf("%sDone!\n\n", mOnTty ? "\n" : "");

  for (auto ext : mAnalExs) ext->WriteHistos();

  // XXXX Output entry lists

  // Print pass counts from all extractors and filters.
  printf("\n");
  printf("Manager pass count                 = %'12lld\n", GetPassCount());
  for (auto ext : mAnalExs)
    printf("Extractor %-24s = %'12lld\n", ext->RefName().Data(), ext->GetPassCount());

  for (auto fil : mAnalFis)
    printf("Filter    %-24s = %'12lld\n", fil->RefName().Data(), fil->GetPassCount());

}


//==============================================================================
// Stuff that should really go elsewhere ... setup functions specific to AAA
// and main().
//==============================================================================

void AnalManager::SetupAaaStuffonAllExtractors()
{
  for (auto ext : mAnalExs)
  {
    auto io = dynamic_cast<AnExIo*>(ext);
    if (io)
    {
      io->SetupAaaDirs();
      io->SetupAaaHistos();
    }
  }
}

void SetupAaaTest(AnalManager& M)
{
  auto pf_AaaMon = new AnFiAaaMoniTest("AaaMonitoringAndTests", M);

  M.AddPreFilter(pf_AaaMon);


  auto fi_InUsa  = new AnFiUsa("InUsa", M, AT_local);

  auto fi_Remote = new AnFiDomain("RemoteAccess", M, "", AT_remote);

  auto fi_ND_All = new AnFiAnyFoo("ReadFromND",  M, [&M]() {
      return M.mUDomain.EndsWith("nd.edu"); });

  auto fi_UNL = new AnFiAnyFoo("ServeFromUNL",  M, [&M]() {
      return M.mSDomain.EndsWith("unl.edu"); });

  auto fi_UCSD = new AnFiAnyFoo("ServeFromUCSD",  M, [&M]() {
      return M.mSDomain.EndsWith("ucsd.edu"); });

  auto ex_All = new AnExIo("AllND", M);
  ex_All->AddFilter(fi_ND_All);
  ex_All->AddFilter(fi_Remote);
  ex_All->AddFilter(fi_InUsa);

  auto ex_UNL = new AnExIo("UNL", M);
  ex_UNL->AddFilter(fi_ND_All);
  ex_UNL->AddFilter(fi_UNL);

  auto ex_UCSD = new AnExIo("UCSD", M);
  ex_UCSD->AddFilter(fi_ND_All);
  ex_UCSD->AddFilter(fi_UCSD);

  M.AddExtractor(ex_All);
  M.AddExtractor(ex_UNL);
  M.AddExtractor(ex_UCSD);

  M.SetupAaaStuffonAllExtractors();
}

void SetupAaaUsa1(AnalManager& M)
{
  auto pf_AaaMon = new AnFiAaaMoniTest("AaaMonitoringAndTests", M);

  // Let all krappe in !!!
  // M.AddPreFilter(pf_AaaMon);


  auto fi_InUsa  = new AnFiUsa("InUsa", M, AT_local);
  auto fi_Remote = new AnFiDomain("RemoteAccess", M, "", AT_remote);
  auto fi_AodSim = new AnFiAodAodsim("AodAodSim", M);

  auto fi_Dur100 = new AnFiDuration("MoreThan100s", M, VC_greater_than, 100);
  // lazy download / xrdcp: F.mReadStats.mMax == 128

  auto fi_Frac10 = new AnFiValueCut<double>
    ("FractionMoreThan10p", M, VC_greater_than, 0.1,
     [&M]() { return M.F.mReadStats.mSumX / M.F.mSizeMB; }
    );

  auto ex_All = new AnExIo("All", M);

  auto ex_UsaAll = new AnExIo("UsaAll", M);
  ex_UsaAll->AddFilter(fi_InUsa);

  auto ex_UsaLoc = new AnExIo("UsaLoc", M);
  ex_UsaLoc->AddFilter(fi_InUsa);
  ex_UsaLoc->AddAntiFilter(fi_Remote);

  auto ex_UsaRem = new AnExIo("UsaRem", M);
  ex_UsaRem->AddFilter(fi_InUsa);
  ex_UsaRem->AddFilter(fi_Remote);

  auto ex_UsaRemAod = new AnExIo("UsaRemAod", M);
  ex_UsaRemAod->AddFilter(fi_InUsa);
  ex_UsaRemAod->AddFilter(fi_Remote);
  ex_UsaRemAod->AddFilter(fi_AodSim);

  auto ex_UsaRemNoAod = new AnExIo("UsaRemNoAod", M);
  ex_UsaRemNoAod->AddFilter(fi_InUsa);
  ex_UsaRemNoAod->AddFilter(fi_Remote);
  ex_UsaRemNoAod->AddAntiFilter(fi_AodSim);

  auto ex_UsaRemAodFrac10Dur100 = new AnExIo("UsaRemAodFrac10Dur100", M);
  ex_UsaRemAodFrac10Dur100->AddFilter(fi_InUsa);
  ex_UsaRemAodFrac10Dur100->AddFilter(fi_Remote);
  ex_UsaRemAodFrac10Dur100->AddFilter(fi_AodSim);
  ex_UsaRemAodFrac10Dur100->AddFilter(fi_Frac10);
  ex_UsaRemAodFrac10Dur100->AddFilter(fi_Dur100);


  M.AddExtractor(ex_All);
  M.AddExtractor(ex_UsaAll);
  M.AddExtractor(ex_UsaLoc);

  M.AddExtractor(ex_UsaRem);
  M.AddExtractor(ex_UsaRemAod);
  M.AddExtractor(ex_UsaRemNoAod);
  M.AddExtractor(ex_UsaRemAodFrac10Dur100);

  M.SetupAaaStuffonAllExtractors();
}

//------------------------------------------------------------------------------

void SetupAaaIov(AnalManager& M)
{
  auto pf_AaaMon = new AnFiAaaMoniTest("AaaMonitoringAndTests", M);

  M.AddPreFilter(pf_AaaMon);

  auto fi_IovLoc    = new AnFiDomain("IovLoc",    M, "", AT_local);
  auto fi_IovNonLoc = new AnFiDomain("IovNonLoc", M, "", AT_nonlocal);

  auto fi_InUsa  = new AnFiUsa("InUsa", M, AT_local);
  auto fi_AodSim = new AnFiAodAodsim("AodAodSim", M);

  auto fi_YesVread  = new AnFiAnyFoo("YesVread",  M, [&M]() {
      return M.F.mVecReadStats.mN > 0; });

  auto fi_Vread60p  = new AnFiAnyFoo("Vread60p",  M, [&M]() {
      return M.F.mVecReadStats.mSumX / M.F.mReadStats.mSumX >= 0.6; });

  auto fi_CrappyIov = new AnFiCrappyIov("CrappyIov", M);


  auto ex_IovAny = new AnExIo("IovAnyAod", M);
  ex_IovAny->AddFilter(fi_InUsa);
  ex_IovAny->AddFilter(fi_AodSim);

  auto ex_IovLoc = new AnExIo("IovLocAod", M);
  ex_IovLoc->AddFilter(fi_InUsa);
  ex_IovLoc->AddFilter(fi_AodSim);
  ex_IovLoc->AddFilter(fi_IovLoc);

  auto ex_IovNonLoc = new AnExIo("IovNonLocAod", M);
  ex_IovNonLoc->AddFilter(fi_InUsa);
  ex_IovNonLoc->AddFilter(fi_AodSim);
  ex_IovNonLoc->AddFilter(fi_IovNonLoc);

  auto ex_IovAnyNonAod = new AnExIo("IovAnyNoAod", M);
  ex_IovAnyNonAod->AddFilter(fi_InUsa);
  ex_IovAnyNonAod->AddAntiFilter(fi_AodSim);

  auto ex_SeekNonLoc = new AnExIov("SeekNonLocAod", M);
  ex_SeekNonLoc->AddFilter(fi_InUsa);
  ex_SeekNonLoc->AddFilter(fi_AodSim);
  ex_SeekNonLoc->AddFilter(fi_IovNonLoc);
  ex_SeekNonLoc->AddFilter(fi_YesVread);
  ex_SeekNonLoc->AddFilter(fi_CrappyIov);


  //M.AddExtractor(ex_IovAny);
  //M.AddExtractor(ex_IovLoc);
  //M.AddExtractor(ex_IovNonLoc);
  //M.AddExtractor(ex_IovAnyNonAod);

  M.AddExtractor(ex_SeekNonLoc);


  M.SetupAaaStuffonAllExtractors();
}

void SetupAaaCacheSim(AnalManager& M)
{
  auto pf_AaaMon = new AnFiAaaMoniTest("AaaMonitoringAndTests", M);

  M.AddPreFilter(pf_AaaMon);

  auto fi_IovLoc    = new AnFiDomain("IovLoc",    M, "", AT_local);
  auto fi_IovNonLoc = new AnFiDomain("IovNonLoc", M, "", AT_nonlocal);

  auto fi_InUsa  = new AnFiUsa("InUsa", M, AT_local);
  auto fi_AodSim = new AnFiAodAodsim("AodAodSim", M);

  auto fi_YesVread  = new AnFiAnyFoo("YesVread",  M, [&M]() {
      return M.F.mVecReadStats.mN > 0; });

  auto fi_Vread60p  = new AnFiAnyFoo("Vread60p",  M, [&M]() {
      return M.F.mVecReadStats.mSumX / M.F.mReadStats.mSumX >= 0.6; });

  auto fi_CrappyIov = new AnFiCrappyIov("CrappyIov", M);


  auto ex_CacheSim = new AnExCacheSim("CacheSim", M);
  ex_CacheSim->AddFilter(fi_InUsa);
  ex_CacheSim->AddFilter(fi_AodSim);
  ex_CacheSim->AddFilter(fi_IovNonLoc);
  ex_CacheSim->AddFilter(fi_YesVread);
  ex_CacheSim->AddFilter(fi_CrappyIov);

  auto ex_CacheSimVec60 = new AnExCacheSim("CacheSimVec60", M, "", 0.6);
  ex_CacheSimVec60->AddFilter(fi_InUsa);
  ex_CacheSimVec60->AddFilter(fi_AodSim);
  ex_CacheSimVec60->AddFilter(fi_IovNonLoc);
  ex_CacheSimVec60->AddFilter(fi_YesVread);
  ex_CacheSimVec60->AddFilter(fi_Vread60p);
  ex_CacheSimVec60->AddFilter(fi_CrappyIov);


  M.AddExtractor(ex_CacheSim);
  M.AddExtractor(ex_CacheSimVec60);


  M.SetupAaaStuffonAllExtractors();
}

void SetupAaaFnalRal(AnalManager& M)
{
  auto pf_AaaMon = new AnFiAaaMoniTest("AaaMonitoringAndTests", M);

  M.AddPreFilter(pf_AaaMon);

  auto fi_FnalToRal  = new AnFiAnyFoo("FnalToRal",  M, [&M]() {
      return M.U.mFromDomain.EndsWith("rl.ac.uk") && M.S.mDomain.EndsWith("fnal.gov");
    });

  auto ex_FnalToRal = new AnExIo("FnalToRal", M);
  ex_FnalToRal->AddFilter(fi_FnalToRal);

  M.AddExtractor(ex_FnalToRal);

  M.SetupAaaStuffonAllExtractors();
}

//==============================================================================

AnalManager* setup_iov()
{
  AnalManager *mgp = new AnalManager("Mgr", "Iov-3", "XrdFar", "/bar/xrdmon-xxx-merged/");
  AnalManager &mgr = *mgp;

  mgr.AddFile("*.root");
  mgr.SetEdgeTimes(1359748800, 1393661701);
  //   mgr.ScanEdgeTimes();

  //mgr.AddFile("*-2014-*.root");
  //mgr.SetEdgeTimes(1388563336, 1393661701);
  //   mgr.ScanEdgeTimes();

  SetupAaaIov(mgr);

  return mgp;
}

AnalManager* setup_cache_sim()
{
  AnalManager *mgp = new AnalManager("Mgr", "CacheSim-3", "XrdFar", "/bar/xrdmon-xxx-merged/");
  AnalManager &mgr = *mgp;

  mgr.AddFile("*.root");
  mgr.SetEdgeTimes(1359748800, 1393661701);
  //   mgr.ScanEdgeTimes();

  //mgr.AddFile("*-2014-*.root");
  //mgr.SetEdgeTimes(1388563336, 1393661701);
  //   mgr.ScanEdgeTimes();

  SetupAaaCacheSim(mgr);

  return mgp;
}

AnalManager* setup_all_krappe()
{
  AnalManager *mgp = new AnalManager("Mgr", "AllKrappe-1", "XrdFar", "/bar/xrdmon-xxx-merged/");
  AnalManager &mgr = *mgp;

  // BEWARE: There are entries with stoopid open/close times and durations.
  // ScanEdgeTimes() only looks at 100k first / last entries.
  // Make sure the values returned make sense !!!

  /*
  // TEST RUN:
  //
  // mgr.AddFile("xmfar-2012-06.root");
  // mgr.AddFile("xmfar-2013-01.root");
  mgr.AddFile("xmfar-2013-02.root");
  mgr.ScanEdgeTimes();
  */

  // TRUE RUN:
  //
  mgr.AddFile("*.root");
  // mgr.ScanEdgeTimes();
  // These get rounded up to a full hour (1340521200 -- 1393664400)
  mgr.SetEdgeTimes(1340521948, 1393661701);


  SetupAaaUsa1(mgr);

  return mgp;
}

AnalManager* setup_aaa_test()
{
  AnalManager *mgp = new AnalManager("Mgr", "AaaTest-1", "XrdFar", "/net/xrootd.t2/data/xrdmon/far/", false);
  AnalManager &mgr = *mgp;

  // BEWARE: There are entries with stoopid open/close times and durations.
  // ScanEdgeTimes() only looks at 100k first / last entries.
  // Make sure the values returned make sense !!!

  /*
  // TEST RUN:
  //
  // mgr.AddFile("xmfar-2012-06.root");
  // mgr.AddFile("xmfar-2013-01.root");
  mgr.AddFile("xmfar-2013-02.root");
  mgr.ScanEdgeTimes();
  */

  // TRUE RUN:
  //
  //mgr.AddFile("xmfar-2014-05-*.root");
  //mgr.AddFile("xmfar-2014-06-*.root");
  //mgr.AddFile("xmfar-2014-07-*.root");
  // mgr.SetEdgeTimes(1399010400, 1406642400);

  mgr.AddFile("xmfar-2014-07-23-*.root");
  mgr.AddFile("xmfar-2014-07-24-*.root");
  mgr.ScanEdgeTimes();

  SetupAaaTest(mgr);

  return mgp;
}

/*
  g++ `root-config --cflags --libs` -Wl,-rpath=. AnalManager.cxx libSXrdClasses.so
*/

int main()
{
  setlocale(LC_NUMERIC, "en_US");

  // AnalManager &mgr = * setup_iov();
  // AnalManager &mgr = * setup_cache_sim();
  // AnalManager &mgr = * setup_all_krappe();

  AnalManager &mgr = * setup_aaa_test();

  mgr.Process();

  delete &mgr;

  return 0;
}
