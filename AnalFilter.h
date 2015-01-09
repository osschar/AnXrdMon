#ifndef AnalFilter_h
#define AnalFilter_h

#include <TString.h>

#include <vector>
#include <set>
#include <functional>

class TEntryList;

class SXrdFileInfo;
class SXrdUserInfo;
class SXrdServerInfo;

typedef unsigned int       uint;
typedef long long          long64;
typedef unsigned long long ulong64;

//==============================================================================

enum AccessType_e
{
  AT_any,
  AT_local,     // server and client at X
  AT_nonlocal,  // server at X, client elsewhere
  AT_remote,    // server elsewhere, client at X,
};

enum ValueCut_e
{
  VC_greater_than,
  VC_less_than
};

//==============================================================================

class AnalManager;
class AnalFilter;

typedef std::vector<AnalFilter*>  vpAnalFilter_t;
typedef std::set<AnalFilter*>     spAnalFilter_t;


class AnalFilter
{
private:
  bool            mState;

protected:
  TString         mName;

  AnalManager    &M;

  Long64_t        mPassCount;
  Long64_t        mTotalCount;

  TEntryList     *mEntryList;

public:
  AnalFilter(const TString& name, AnalManager& mgr);
  virtual ~AnalFilter() {}

  const TString& RefName() const { return mName; }

  Long64_t GetPassCount()  const { return mPassCount;  }
  Long64_t GetTotalCount() const { return mTotalCount; }
  

  void SetEntryList(TEntryList* el) { mEntryList = el; }

  bool Passed() const { return mState; }

  bool FilterAndStore();

  virtual bool Filter() = 0;


  bool AllFiltersPass(const vpAnalFilter_t& v)
  {
    for (auto const &filt : v)
    {
      if ( ! filt->mState)  return false;
    }
    return true;
  }

  bool AllFiltersFail(const vpAnalFilter_t& v)
  {
    for (auto const &filt : v)
    {
      if (filt->mState)  return false;
    }
    return true;
  }
};


class AnFiAnyFoo : public AnalFilter
{
public:
  typedef std::function<bool ()> Foo_t;

protected:
  Foo_t      mFunc;

public:
  AnFiAnyFoo(const TString& n, AnalManager& m, Foo_t f) :
    AnalFilter(n, m),
    mFunc(f)
  {}
  virtual ~AnFiAnyFoo() {}

  virtual bool Filter()
  {
    return mFunc();
  }
};


//==============================================================================
// User, domain, etc filters
//==============================================================================

class AnFiUserRealName : public AnalFilter
{
protected:
  TString mUName;

public:
  AnFiUserRealName(const TString& n, AnalManager& m, const TString& uname) :
    AnalFilter(n, m), mUName(uname) {}
  virtual ~AnFiUserRealName() {}

  virtual bool Filter();
};

//------------------------------------------------------------------------------

class AnFiDomain : public AnalFilter
{
protected:
  TString      mDomain;
  AccessType_e mType;

public:
  AnFiDomain(const TString& n, AnalManager& m, const TString& domain, AccessType_e at) :
    AnalFilter(n, m), mDomain(domain), mType(at) {}
  virtual ~AnFiDomain() {}

  virtual bool Filter();
};

//------------------------------------------------------------------------------

class AnFiUsa : public AnalFilter
{
protected:
  AccessType_e mType;

public:
  AnFiUsa(const TString& n, AnalManager& m, AccessType_e at) :
    AnalFilter(n, m), mType(at) {}
  virtual ~AnFiUsa() {}

  virtual bool Filter();
};

//------------------------------------------------------------------------------

class AnFiDuration : public AnalFilter
{
protected:
  Double_t   mDuration;
  ValueCut_e mType;

public:
  AnFiDuration(const TString& n, AnalManager& m, ValueCut_e t, Double_t d) :
    AnalFilter(n, m),
    mDuration(d), mType(t)
  {}
  virtual ~AnFiDuration() {}

  virtual bool Filter();
};

template<typename TT>
class AnFiValueCut : public AnalFilter
{
public:
  typedef std::function<TT ()> Foo_t;

protected:
  Foo_t      mFunc;
  TT         mCutValue;
  ValueCut_e mCutType;

public:
  AnFiValueCut(const TString& n, AnalManager& m, ValueCut_e t, TT v, Foo_t f) :
    AnalFilter(n, m),
    mFunc(f), mCutValue(v), mCutType(t)
  {}
  virtual ~AnFiValueCut() {}

  virtual bool Filter()
  {
    switch (mCutType)
    {
      case VC_greater_than: return mFunc() > mCutValue;
      case VC_less_than:    return mFunc() < mCutValue;
    }
  }
};


//==============================================================================
// Crappy IOV data filter
//==============================================================================

class AnFiCrappyIov : public AnalFilter
{
protected:

public:
  AnFiCrappyIov(const TString& n, AnalManager& m) : AnalFilter(n, m) {}
  virtual ~AnFiCrappyIov() {}

  virtual bool Filter();
};


//==============================================================================
// AAA/CMS specific filters
//==============================================================================

class AnFiAodAodsim : public AnalFilter
{
protected:

public:
  AnFiAodAodsim(const TString& n, AnalManager& m) : AnalFilter(n, m) {}
  virtual ~AnFiAodAodsim() {}

  virtual bool Filter();
};

//------------------------------------------------------------------------------


//==============================================================================
// Sanity / crap pre-filters. Also AAA specific
//==============================================================================

class AnFiAaaMoniTest : public AnalFilter
{
protected:

public:
  AnFiAaaMoniTest(const TString& n, AnalManager& m) : AnalFilter(n, m) {}
  virtual ~AnFiAaaMoniTest() {}

  virtual bool Filter();
};

//------------------------------------------------------------------------------

//==============================================================================


/*
class XXXX : public AnalFilter
{
protected:

public:
  XXXX(const TString& n, AnalManager& m) : AnalFilter(n, m) {}
  virtual ~XXXX() {}

  virtual bool Filter();
};
*/

#endif
