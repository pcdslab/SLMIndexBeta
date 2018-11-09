/*
   Copyright 2017 Western Michigan University
*/

/* Include Required Headers */
#include "Common.h"
#include "CometSearch.h"
#include "ThreadPool.h"
#include "CometStatus.h"
#include "CometSearchManager.h"
#include "CometStatus.h"
#include "Threading.h"
#include "CometDataInternal.h"
#include "bwt.h"
#include <omp.h>

#define PP_MAX_FRAG_CHRG          3
#define PP_MAX_ION_SERIES         2
#define PP_ARRAYSIZE              PP_MAX_FRAG_CHRG * PP_MAX_ION_SERIES * MAX_PEPTIDE_LEN

/* Function Prototypes */
bool AllocateArraysForPeptides(UINT number);

int CheckMassTolerance(double dCalcPepMass);
