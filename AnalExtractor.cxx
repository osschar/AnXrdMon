#include "AnalExtractor.h"
#include "AnalManager.h"

#include <TFile.h>

//==============================================================================

AnalExtractor::DirHolder::DirHolder(TDirectory *dir)
{
  m_ex_dir  = gDirectory;
  m_cur_dir = dir;

  cd_cur();
}

AnalExtractor::DirHolder::~DirHolder()
{
  cd_ex();
}

void AnalExtractor::DirHolder::cd_cur()
{
  if (m_cur_dir)
    m_cur_dir->cd();
  else
    gDirectory = 0;
}

void AnalExtractor::DirHolder::cd_ex()
{
  if (m_ex_dir)
    m_ex_dir->cd();
  else
    gDirectory = 0;
}

//==============================================================================

AnalExtractor::AnalExtractor(const TString& name, AnalManager& mgr,
                             const TString& out_file) :
  AnalFilter(name, mgr),
  mFile(0)
{
  mOutFileName  = M.RefOutDirName() + "/";
  mOutFileName += out_file.IsNull() ? name : out_file;

  static const char *pstfx = ".root";
  if ( ! mOutFileName.EndsWith(pstfx))
  {
    mOutFileName += pstfx;
  }
}

void AnalExtractor::OpenFile()
{
  mFile = TFile::Open(mOutFileName, "create");
  if ( ! mFile)
  {
    fprintf(stderr, "Opening of output file '%s' failed. Probably it exists already.\n",
            mOutFileName.Data());
    exit(2);
  }
}

void AnalExtractor::CloseFile()
{
  mFile->Write();
  mFile->Close();
  delete mFile;
  mFile = 0;
}

//==============================================================================

bool AnalExtractor::Filter()
{
  if ( ! AllFiltersPass(mFilters)) return false;

  if ( ! AllFiltersFail(mAntiFilters)) return false;

  return true;
}
