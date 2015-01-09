#ifndef AnExIo_h
#define AnExIo_h

#include "AnalExtractor.h"

#include <TString.h>

#include <vector>
#include <map>


class TH1;

//------------------------------------------------------------------------------

class AnExIo;

typedef double (AnExIo::*val_foo_t) ();


struct XDir
{
  const char *f_name;
  int         f_accum_idx;

  std::vector<TH1*> f_hvec, f_hvec_t, f_hvec_2d;

  XDir(const char *name, int accum_idx) :
    f_name(name), f_accum_idx(accum_idx)
  {}
};


struct XHisto
{
  const char *f_name, *f_title;
  int         f_nbx;
  double      f_xl, f_xh;
  val_foo_t   f_foo; 

  XHisto(const char *name, const char *title, int nbx,
         double xl, double xh, val_foo_t foo) :
    f_name(name), f_title(title), f_nbx(nbx),
    f_xl(xl), f_xh(xh), f_foo(foo)
  {}
};


struct XHistoCum : public XHisto
{
  std::vector<double> f_cum;

  XHistoCum(const char *name, const char *title, int nbx,
            double xl, double xh, val_foo_t foo) :
    XHisto(name, title, nbx, xl, xh, foo)
  {}
};


//==============================================================================


class AnExIo : public AnalExtractor
{
protected:
  std::vector<XDir>      A_dirs;
  Int_t                  N_dirs;

  std::map<TString, int> f_dir_to_idx_map;

  // ----------------------------------------------------------------

  std::vector<XHisto>      A_histos;
  Int_t                    N_A_histos;

  std::vector<XHistoCum>   C_histos;
  Int_t                    N_C_histos;

public:

  AnExIo(const TString& name, AnalManager &mgr, const TString& out_file="");

  void AddDir(const char *name, int accum_idx);
  void SetupAaaDirs();

  void AddHisto(const char *name, const char *title, int nbx,
                double xl, double xh, val_foo_t foo);
  void AddHistoCum(const char *name, const char *title, int nbx,
                   double xl, double xh, val_foo_t foo);
  void SetupAaaHistos();


  // ----------------------------------------------------------------

  virtual void BookHistos();

  virtual void WriteHistos();

  // ----------------------------------------------------------------

  virtual void Process();

// ----------------------------------------------------------------

  double foo_frac_read();
  double foo_frac_vread();
  double foo_bytes_read();
  double foo_open_duration();
  double foo_num_reqs();

  double foo_num_sin_reqs();
  double foo_num_vec_reqs();
  double foo_avg_vec_subreqs();

  double foo_req_avg();
  double foo_req_min();
  double foo_req_max();
  double foo_req_rate();
  double foo_data_rate();

  double foo_zero();
};

#endif
