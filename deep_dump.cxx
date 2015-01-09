#include "SXrdClasses.h"

extern SXrdFileInfo   F;
extern SXrdUserInfo   U;
extern SXrdServerInfo S;
extern SXrdIoInfo     I;

static const double OM = 1.0/1024/1024;

// E.g. run:
// rm -f ucsd_test.root; ./ucsd_anal ucsd_test XrdFar /bar/xrdmon-xxx-merged/xmxxx-2013-07.root

void deep_dump()
{
  // Testing / inspection printouts
  printf("%s\n  %s.%s   %s@%s.%s %s\n"
         "  S=%.3f MB   Nr=%lld Ns=%lld Nv=%lld Vavg=%.1f Nerr=%d\n"
         "  Frac_file %.3f - (%.3f / %.3f)   Frac_size (%.3f / %.3f)   Frac_count (%.3f / %.3f)\n",
         F.mName.Data(),
         S.mHost.Data(), S.mDomain.Data(),
         U.mName.Data(), U.mFromHost.Data(), U.mFromDomain.Data(), U.mRealName.Data(),
         F.mSizeMB,
         F.mReadStats.mN, F.mSingleReadStats.mN, F.mVecReadStats.mN,
         F.mVecReadCntStats.GetAverage(), I.mNErrors,
         F.mReadStats.mSumX / F.mSizeMB, 
         F.mSingleReadStats.mSumX / F.mSizeMB, F.mVecReadStats.mSumX / F.mSizeMB,
         F.mSingleReadStats.mSumX / F.mReadStats.mSumX, F.mVecReadStats.mSumX / F.mReadStats.mSumX,
         (double) F.mSingleReadStats.mN / F.mReadStats.mN, (double) F.mVecReadStats.mN / F.mReadStats.mN
  );

  Long64_t ro = 0, so = 0, vo = 0, wo = 0, ddo;
  SRange   wr;
  SRange   prr, psrs, psrv, pvrv, pvrs; // pos
  SRange   nrr, nsrs, nsrv, nvrv, nvrs; // neg
  SRange   vl, vir, vtl, vfr;           // vec

  bool prev_read_was_vector = false;
  int  count = 0;

  for (vSXrdReq_i i = I.mReqs.begin(); i != I.mReqs.end(); ++i, ++count)
  {
    switch (i->Type())
    {
      case SXrdReq::R_Write:
      {
        ddo = i->Offset() - wo; wo = i->Offset() + i->Length();
        wr.AddSample(ddo);

	break;
      }
      case SXrdReq::R_Read:
      {
        ddo = i->Offset() - ro;
        if ( ! prev_read_was_vector)
        {
          if (ddo >= 0) psrs.AddSample(OM * ddo); else nsrs.AddSample(OM * ddo);
        }
        else
        {
          if (ddo >= 0) psrv.AddSample(OM * ddo); else nsrv.AddSample(OM * ddo);
        }
        if (ddo >= 0) prr.AddSample(OM * ddo); else nrr.AddSample(OM * ddo);

        ro = so = i->Offset() + i->Length();

        prev_read_was_vector = false;
	break;
      }
      case SXrdReq::R_VecRead:
      {
	Int_t sr_idx = i->SubReqIndex();

        Double_t sum_len = 0;

	if (sr_idx >= 0)
	{
          ddo = I.mOffsetVec[sr_idx] - ro;
          if (prev_read_was_vector)
          {
            if (ddo >= 0) pvrv.AddSample(OM * ddo); else nvrv.AddSample(OM * ddo);
          }
          else
          {
            if (ddo >= 0) pvrs.AddSample(OM * ddo); else nvrs.AddSample(OM * ddo);
          }
          if (ddo >= 0) prr.AddSample(OM * ddo); else nrr.AddSample(OM * ddo);

	  const Int_t max = sr_idx + i->SubReqsStored();
	  for (Int_t  si  = sr_idx; si < max; ++si)
	  {
            const Long64_t off = I.mOffsetVec[si];
            const Int_t    len = I.mLengthVec[si];

            ddo = off - vo; vo = off + len;

            if (si == sr_idx) {  }

            if (si >  sr_idx) { vir.AddSample(OM * ddo); }

            vl.AddSample(OM * len);

            sum_len += len;
	  }
          // Total length and fraction of total lenght actually read.
          vtl.AddSample(OM * (vo - I.mOffsetVec[sr_idx]));
          vfr.AddSample(sum_len / (vo - I.mOffsetVec[sr_idx]));
	}
        else
        {
          printf("Grr, sr_idx=%d for req num %d (out of %d)\n", sr_idx, count, I.mReqs.size());
        }

        ro = vo;

        prev_read_was_vector = true;
	break;
      }
    }
  }

  F.mReadStats      .Dump("RS:    ", "\t\tAll reads sizes\n");
  F.mSingleReadStats.Dump("SRS:   ", "\t\tSingle read sizes\n");
  F.mVecReadStats   .Dump("VRS:   ", "\t\tVector read sizes (summed for each vec req)\n");
  F.mVecReadCntStats.Dump("VRCS:  ", "\t\tVector request number of sub-requests\n");

  prr .Dump("O_R+:  ", "\t\tPositive offsets of read requests\n");
  psrs.Dump("O_Ss+: ", "\t\tPositive offsets of single reads (preceeded by single reads)\n");
  psrv.Dump("O_Sv+: ", "\t\tPositive offsets of single reads (preceeded by vector reads)\n");
  pvrv.Dump("O_Vv+: ", "\t\tPositive offsets of vector reads (preceeded by vector reads)\n");
  pvrs.Dump("O_Vs+: ", "\t\tPositive offsets of vector reads (preceeded by single reads)\n");
  vir .Dump("Vinof: ", "\t\tOffsets within vector reads\n");
  vl  .Dump("Vtrul: ", "\t\tExtent of all sub-reqs in vec reads\n");
  vtl .Dump("Vtotl: ", "\t\tTotal extent of a vec req, from start to last byte of last sub req\n");
  vfr .Dump("Vfr:   ", "\t\tFraction of above actually requested\n");

  nrr .Dump("O_R-:  ", "\t\tNegative offsets of read requests\n");
  nsrs.Dump("O_Ss-: ", "\t\tNegative offsets of single reads (preceeded by single reads)\n");
  nsrv.Dump("O_Sv-: ", "\t\tNegative offsets of single reads (preceeded by vector reads)\n");
  nvrv.Dump("O_Vv-: ", "\t\tNegative offsets of vector reads (preceeded by vector reads)\n");
  nvrs.Dump("O_Vs-: ", "\t\tNegative offsets of vector reads (preceeded by single reads)\n");

  char uline[256];
  while (true)
  {
    fgets(uline, 256, stdin);

    if      (uline[0] == '1') I.Dump(1);
    else if (uline[0] == '2') I.Dump(2);
    else    break;
  }
}
