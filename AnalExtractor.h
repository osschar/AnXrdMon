#ifndef AnalExtractor_h
#define AnalExtractor_h

#include "AnalFilter.h"

class AnalManager;

class TFile;
class TDirectory;

//------------------------------------------------------------------------------

class AnalExtractor : public AnalFilter
{
  friend class AnalManager;

protected:

  class DirHolder
  {
    TDirectory *m_ex_dir, *m_cur_dir;
  public:
    DirHolder(TDirectory *dir = 0);
    ~DirHolder();

    void cd_cur();
    void cd_ex();
  };

  TString           mOutFileName;

  TFile            *mFile;

  vpAnalFilter_t    mFilters;     // Must all pass
  vpAnalFilter_t    mAntiFilters; // Must all fail

public:

  AnalExtractor(const TString& name, AnalManager &mgr, const TString& out_file="");

  void AddFilter(AnalFilter* f)     { mFilters    .push_back(f); }
  void AddAntiFilter(AnalFilter* f) { mAntiFilters.push_back(f); }

  void OpenFile();
  void CloseFile();

  // ----------------------------------------------------------------

  virtual void BookHistos()  {}

  virtual void WriteHistos() {}

  // ----------------------------------------------------------------

  virtual bool Filter();

  virtual void Process() = 0;
};


typedef std::vector<AnalExtractor*>  vpAnalExtractor_t;

#endif
