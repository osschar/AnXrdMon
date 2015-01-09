//Quick script to get the names of all the domains used during a particular
//timeframe, in order to set up smarter domain name parsing

/*
  g++ `root-config --cflags --libs` -Wl,-rpath=. count_dn_nodn.cxx -o count_dn_nodn libSXrdClasses.so
*/

#include "SXrdClasses.h"

#include "TTree.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TChain.h"
#include "TSystem.h"
#include "TPRegexp.h"
#include "TMath.h"

#include <map>
#include <iostream>
#include <cassert>

//==============================================================================
// main
//==============================================================================

int main()
{
  setlocale(LC_NUMERIC, "en_US");

  // Chain up several reports
  TChain mychain("XrdFar");

  // mychain.Add("/net/xrootd.t2/data/xrdmon/far/xmfar-2012-06-*.root");
  mychain.Add("/net/xrootd.t2/data/xrdmon/far/xmfar-2014-09-1*.root");
  // mychain.Add("/bar/xrdmon-far-merged/xmfar-2014-02.root");

  // mychain.Add("/bar/xrdmon-far-merged/xmfar-*.root");
  // mychain.Add("/bar/xrdmon-xxx-merged/xmxxx-*.root");

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

  Long64_t N = mychain.GetEntries();
  std::cout << N << " entries found.\n";

  std::map<TString, int> dn_count;
  std::map<TString, int> no_dn_count;

  for (Long64_t i = 0; i < N; ++i)
  {
    mychain.GetEntry(i);

    // if (S.mDomain == "datarid.cea.fr" && U.mDN.IsNull())
    if (U.mDN.IsNull())
    {
      printf("'%-20s' %-10s %-40s %s\n", S.mDomain.Data(), U.mVO.Data(), U.mDN.Data(), F.mName.Data());
    }

    // TString host = S.mHost + S.mDomain;

    // if (U.mDN.IsNull() )
    // {
    //   ++no_dn_count[host];
    // }
    // else
    // {
    //   ++dn_count[host];
    // }
  }

  printf("==================================================================\n");
  printf("Yes DN\n");
  printf("==================================================================\n");

  for (auto i = dn_count.begin(); i != dn_count.end(); ++i)
  {
    printf("%-40s %d\n", i->first.Data(), i->second);
  }

  printf("==================================================================\n");
  printf("No DN\n");
  printf("==================================================================\n");

  for (auto i = no_dn_count.begin(); i != no_dn_count.end(); ++i)
  {
    printf("%-40s %d\n", i->first.Data(), i->second);
  }


  return 0;
}
