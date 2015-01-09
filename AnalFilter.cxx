#include "AnalFilter.h"
#include "AnalManager.h"

#include "SXrdClasses.h"

namespace
{
  const double OneMB = 1024 * 1024;

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

//==============================================================================

AnalFilter::AnalFilter(const TString& name, AnalManager& mgr) :
  mState(false), mName(name), M(mgr),
  mPassCount(0), mTotalCount(0),
  mEntryList(0)
{}

bool AnalFilter::FilterAndStore()
{
  mState = Filter();

  if (mState)
  {
    ++mPassCount;
    if (mEntryList)
    {
      // XXX Half cooked, need lots of crap to make this efficient.
      // mEntryList->Enter(M.GetTreeII(), M.GetTree());
    }
  }
  ++mTotalCount;

  return mState;
}

//==============================================================================
// User, domain, etc filters
//==============================================================================

bool AnFiUserRealName::Filter()
{
  return M.U.mRealName.Contains(mUName);
}

//------------------------------------------------------------------------------

bool AnFiDomain::Filter()
{
  if (mDomain.IsNull()) {
    switch (mType)
    {
      case AT_any:      return true;
      case AT_local:    return M.mSDomain == M.mUDomain;
      case AT_nonlocal: return M.mSDomain != M.mUDomain;
      case AT_remote:   return M.mSDomain != M.mUDomain;
    }
  } else {
    switch (mType)
    {
      case AT_any:      return M.mSDomain == mDomain || M.mUDomain == mDomain;
      case AT_local:    return M.mSDomain == mDomain && M.mUDomain == mDomain;
      case AT_nonlocal: return M.mSDomain == mDomain && M.mUDomain != mDomain;
      case AT_remote:   return M.mSDomain != mDomain && M.mUDomain == mDomain;
    }
  }
}

//------------------------------------------------------------------------------

bool is_usa(const TString& d)
{
  // d is truncated domain, last two names only.
  // this is pathetish, and it will only get worse :(

  return (d.EndsWith(".edu") || d.EndsWith(".gov") ||
          d == "ultralight.org" || d == "batlab.org" || d == "aglt2.org" ||
          d == "rr.com" || d == "amazonaws.com" || d == "akamaitechnologies.com"
         );
}

bool AnFiUsa::Filter()
{
  bool s_usa = is_usa(M.mSDomain);
  bool u_usa = is_usa(M.mUDomain);

  switch (mType)
  {
    case AT_any:      return   s_usa ||   u_usa;
    case AT_local:    return   s_usa &&   u_usa;
    case AT_nonlocal: return   s_usa && ! u_usa;
    case AT_remote:   return ! s_usa &&   u_usa;
  }
}

//------------------------------------------------------------------------------

bool AnFiAodAodsim::Filter()
{
  return M.mSlashRe[5] == "AOD" || M.mSlashRe[5] == "AODSIM";
}

//------------------------------------------------------------------------------

bool AnFiDuration::Filter()
{
  switch (mType)
  {
    case VC_greater_than: return M.mDt > mDuration;
    case VC_less_than:    return M.mDt < mDuration;
   }
}


//==============================================================================
// Crappy IOV data filter
//==============================================================================

bool AnFiCrappyIov::Filter()
{
  const Long64_t fsize = NearestLong(M.F.mSizeMB * OneMB);

  Long64_t ro = 0, so = 0, vo = 0, ddo;

  Long64_t sum_pos = 0, sum_neg = 0;

  bool prev_read_was_vector = false;

  for (vSXrdReq_i i = M.I.mReqs.begin(); i != M.I.mReqs.end(); ++i)
  {
    switch (i->Type())
    {
      case SXrdReq::R_Write:
      {
	break;
      }
      case SXrdReq::R_Read:
      {
        ddo = i->Offset() - ro;

        if (ddo >= 0) sum_pos += ddo; else sum_neg -= ddo;

        ro = so = i->Offset() + i->Length();
	break;
      }
      case SXrdReq::R_VecRead:
      {
	Int_t sr_idx = i->SubReqIndex();

	if (sr_idx >= 0)
	{
          ddo = M.I.mOffsetVec[sr_idx] - ro;

          if (ddo >= 0) sum_pos += ddo; else sum_neg -= ddo;

	  const Int_t max = sr_idx + i->SubReqsStored();
	  for (Int_t  si  = sr_idx; si < max; ++si)
	  {
            const Long64_t off = M.I.mOffsetVec[si];
            const Int_t    len = M.I.mLengthVec[si];

            ddo = off - vo; vo = off + len;
	  }
	}

        ro = vo;

	break;
      }
    }

    if (ro > fsize || ro < 0)
    {
      // printf("Shitty offset: Event %lld, Off = %lld (fsize=%lld)\n",
      //        M.GetChainI(), ro, fsize);
      return false;
    }

    if (sum_pos > sum_neg + fsize)
    {
      // printf("Pos offsets larger than neg + fsize: Event %lld, pos = %lld, neg = %lld  (fsize=%lld)\n",
      //        M.GetChainI(), sum_pos, sum_neg, fsize);
      return false;
    }
  }

  return true;
}

//==============================================================================
// Sanity / crap pre-filters
//==============================================================================

bool AnFiAaaMoniTest::Filter()
{
  // Assert prefix
  if ( ! M.F.mName.BeginsWith("/store/"))
  {
    return false;
  }
  // Paths used for monitoring / tests.
  if (M.F.mName.BeginsWith("/store/test/") ||
      M.F.mName.BeginsWith("/store/user/dan/") ||
      M.F.mName.BeginsWith("/store/temp/") ||
      M.F.mName.BeginsWith("/store/mc/JobRobot/RelValProdTTbar/") ||
      M.F.mName.BeginsWith("/store/mc/SAM/GenericTTbar/")
  )
  {
    return false;
  }

  // People who only do monitoring / development / testing.
  if (M.U.mRealName.Contains("Bockelman") ||
      M.U.mRealName.Contains("Andrea Sciaba") ||
      M.U.mRealName.Contains("Vuosalo") ||
      M.U.mRealName.Contains("Daniel Charles Bradley") ||
      M.U.mRealName.Contains("Tadel") ||
      M.U.mRealName.IsNull()
  )
  {
    return false;
  }

  return true;
}

//==============================================================================

/*
bool XXXX::Filter()
{
  return ;
}
*/
