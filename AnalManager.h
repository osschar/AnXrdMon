#ifndef AnalManager_h
#define AnalManager_h

#include "AnalFilter.h"
#include "AnalExtractor.h"

#include "SXrdClasses.h"

#include "TPRegexp.h"

#include <vector>
#include <set>


class TChain;
class TFile;

class AnalManager : private AnalFilter
{
  TChain           *mChn;
  Long64_t          mChnN;
  Long64_t          mChnI;
  // XXX
  // For entry-lists need current tree and current index.
  // Then, need to detect tree change and notify all filters.
  // TTree    *mTree
  // Long64_t  mTreeI;

  TString           mInFilePrefix;
  TString           mOutDirName;

  vpAnalFilter_t    mPreFilters;   // Results not stored.

  vpAnalExtractor_t mAnalExs;
  spAnalFilter_t    mAnalFis;

public:
  TChain*        GetChain()           { return mChn;  }
  Long64_t       GetChainN()          { return mChnN; }
  Long64_t       GetChainI()          { return mChnI; }
  // TTree*         GetTree()        { return ; }

  const TString& RefOutDirName() const { return mOutDirName; }

  SXrdFileInfo      F, *_fp;
  SXrdUserInfo      U, *_up;
  SXrdServerInfo    S, *_sp;
  SXrdIoInfo        I, *_ip;

  Bool_t      mBranchIActive;
  Bool_t      mOnTty;

  // Whole data-set constants
  Long64_t    mMinT, mMaxT;
  Long64_t    mTotalDtSec, mTotalDtMin,  mTotalDtHour;
  Double_t    mTotalDtDay, mTotalDtWeek, mTotalDtMonth;
  TString     mMinDate, mMaxDate;

  // Per event variables
  Double_t    mDt;

  TPMERegexp  mSDomainRe;
  TPMERegexp  mUDomainRe;
  TString     mSDomain, mUDomain;

  TPMERegexp  mSlashRe;

public:

  AnalManager(const TString& name,      const TString& out_dir,
              const TString& tree_name, const TString& prefix="",
              bool setup_I_branch = true);
  virtual ~AnalManager();

  void AddFile(const TString& files);

  void AddPreFilter(AnalFilter*    flt);
  void AddExtractor(AnalExtractor* ext);

  void ScanEdgeTimes(Long64_t scan_entries=100000);
  void SetEdgeTimes(Long64_t min, Long64_t max, bool verbose=true);

  virtual bool Filter();

  void Process();

  // To get rid of ...
  void SetupAaaStuffonAllExtractors();
};

#endif
