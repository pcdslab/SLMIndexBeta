/*
   Copyright 2012 University of Washington

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef _COMETSEARCH_H_
#define _COMETSEARCH_H_

#include "Common.h"
#include "CometDataInternal.h"

struct SearchThreadData
{
   sDBEntry dbEntry;
   UINT Query_num;
   UINT sdbEntryNum;
   bool *pbSearchMemoryPool;
   /* HM: Add the Search Query here */
   Query *iQuery;

   SearchThreadData()
   {

   }

   SearchThreadData(Query *_iQuery, UINT number)
   {
      Query_num = number;
      iQuery = _iQuery;
   }

   SearchThreadData(sDBEntry &dbEntry_in)
   {
      dbEntry.strName = dbEntry_in.strName;
      dbEntry.strSeq = dbEntry_in.strSeq;
      dbEntry.iSeqFilePosition = dbEntry_in.iSeqFilePosition;
      dbEntry.iLenPep = dbEntry_in.iLenPep;
   }

   ~SearchThreadData()
   {
      // Mark that the memory is no longer in use.
      // DO NOT FREE MEMORY HERE. Just release pointer.
      Threading::LockMutex(g_searchMemoryPoolMutex);

      if (pbSearchMemoryPool!=NULL)
      {
         *pbSearchMemoryPool=false;
         pbSearchMemoryPool=NULL;
      }

      Threading::UnlockMutex(g_searchMemoryPoolMutex);
   }
};

class CometSearch
{
public:
   CometSearch();
   ~CometSearch();

   // Manages memory in the search memory pool
   static bool AllocateMemory(int maxNumThreads);
   static bool DeallocateMemory(int maxNumThreads);
   static bool RunSearch(int minNumThreads,
                         int maxNumThreads,
                         int iPercentStart,
                         int iPercentEnd);
   static void SearchThreadProc(SearchThreadData *pSearchThreadData);
   static void PreprocessThreadProc(SearchThreadData *pSearchThreadData);
   bool DoSearch(Query *iQuery, bool *pbDuplFragment);
   bool PreprocessDatabase(sDBEntry dbe, bool *pbDuplFragment);
   bool Burrows_Wheeler_Transform();

public:
   /* HM: sdbEntry number */
   UINT sdbEntryNum;
   int ith;

private:

   // Core search functions
   int BinarySearchMass(int start,
                        int end,
                        double dCalcPepMass);
   void SubtractVarMods(int *piVarModCounts,
                        int cResidue,
                        int iResiduePosition);
   void CountVarMods(int *piVarModCounts,
                     int cResidue,
                     int iResiduePosition);
   bool HasVariableMod(int varModCounts[],
                       int iCVarModCount,
                       int iNVarModCount);
   int WithinMassTolerance(double dCalcPepMass,
                           char *szProteinSeq,
                           int iStartPos,
                           int iEndPos);
   void XcorrScore(char *szProteinSeq,
                   char *szProteinName,
                   int iStartPos,
                   int iEndPos,
                   bool bFoundVariableMod,
                   double dCalcPepMass,
                   bool bDecoyPep,
                   int iWhichQuery,
                   int iLenPeptide,
                   char *pcVarModSites);
   bool CheckEnzymeTermini(char *szProteinSeq,
                           int iStartPos,
                           int iEndPos);
   bool CheckMassMatch(int iWhichQuery,
                       double dCalcPepMass);
   double GetFragmentIonMass(int iWhichIonSeries,
                             int i,
                             int ctCharge,
                             double *pdAAforward,
                             double *pdAAreverse);
   /*************** HM: Adding my new and modified APIs ****************/
   bool GenerateAndPopulatePeptides(char *szProteinSeq,
                                    char *szProteinName,
                                    bool bNtermPeptideOnly,
                                    bool *pbDuplFragment);
   bool Pattern_Search(iMSQPeptide *iQuery,
                       short *list_ptr,
                       int star_peps,
                       int end_pep);

   bool Pattern_Mods_Search(iMSQPeptide *iQuery,
                            short *list_ptr,
                            USHORT index,
                            int star_peps,
                            int end_peps);

   bool SearchForModPeptides(char *dbe,
                             varModEntry *entry,
                             int iWhichQuery,
                             bool *pbDuplFragment);

   bool SearchForPeptides(char *dbe,
                          char *entry,
                          bool bNtermPeptideOnly,
                          UINT iWhichQuery,
                          bool *pbDuplFragment);

   bool GenerateAndPopulateModPeptide(char *szProteinSeq,
                                      int iWhichQuery,
                                      bool *pbDuplFragment);

   /***************************************************/
   int CheckDuplicate(int iWhichQuery,
                      int iStartPos,
                      int iEndPos,
                      bool bFoundVariableMod,
                      double dCalcPepMass,
                      char *szProteinSeq,
                      char *szProteinName,
                      bool bDecoyResults,
                      char *pcVarModSites);
   void StorePeptide(int iWhichQuery,
                     int iStartPos,
                     int iEndPos,
                     bool bFoundVariableMod,
                     char *szProteinSeq,
                     double dCalcPepMass,
                     double dScoreSp,
                     char *szProteinName,
                     bool bStoreSeparateDecoy,
                     char *pcVarModSites);
   bool VarModSearch(char *szProteinSeq,
                     char *szProteinName,
                     int varModCounts[],
                     int iStartPos,
                     int iEndPos,
                     bool *pbDuplFragment);
   double TotalVarModMass(int *pVarModCounts);
   bool PermuteMods(char *szProteinSeq,
                    int iWhichQuery,
                    int iWhichMod,
                    bool *pbDuplFragments);
   int  twiddle( int *x, int *y, int *z, int *p);
   void inittwiddle(int m, int n, int *p);
   bool TranslateNA2AA(int *frame,
                       int iDirection,
                       char *sDNASequence);
   char GetAA(int i,
              int iDirection,
              char *sDNASequence);

   // Processing results
   void CalculateEvalue(struct Results *_pResults);
   void GenerateXcorrDecoys(struct Results *_pResults,
                            bool bDecoy);

    // Displaying results
   void PrintResults(struct Results *_pResults,
                     bool bDecoySearch);
   void PrintOutputLine(struct Results *_pResults,
                        int iRankXcorr,
                        int iLenMaxDuplicates,
                        int iMaxWidthReference,
                        int iWhichResult,
                        bool bDecoySearch,
                        FILE *fpOut);
   void PrintIons(int iTmpCharge,
                  FILE *fpOut);

   // Cleaning up
   void CleanUp();

private:

   struct VarModStat
   {
       int iTotVarModCt;                     // # mod positions in peptide
       int iTotBinaryModCt;                  // # mod positions of all binary mods in group in the peptide
       int iMatchVarModCt;
       int iVarModSites[MAX_PEPTIDE_LEN+2];  // last 2 positions are for n- and c-term
   };

   struct VarModInfo
   {
       VarModStat varModStatList[VMODS];
       int        iStartPos;     // The start position of the peptide sequence
       int        iEndPos;       // The end position of the peptide sequence
       double     dCalcPepMass;  // Mass of peptide with mods
   };

   struct PepMassTolerance
   {
       double dPeptideMassTolerance;
       double dPeptideMassToleranceMinus;
       double dPeptideMassTolerancePlus;
   };

   struct PepMassInfo
   {
       double           dCalcPepMass;
       PepMassTolerance pepMassTol;
   };

   struct MatchedPeaksStruct
   {
      double dMass;
      double dIntensity;
   };

   struct ProteinInfo
   {
       char szProteinName[WIDTH_REFERENCE];
       char *pszProteinSeq;
       int  iProteinSeqLength;
       int  iAllocatedProtSeqLength;
       int  iSeqFilePosition;
   };

   struct SpectrumInfoInternal
   {
      Spectrum *pSpectrum;
      int      iChargeState;
      int      iArraySize;     // m/z versus intensity array
      int      iHighestIon;
      double   dExperimentalPeptideMass;
      double   dHighestIntensity;
      double   dTotalIntensity;
   };

   double             _pdAAforward[MAX_PEPTIDE_LEN];      // Stores fragment ion fragment ladder calc.; sum AA masses including mods
   double             _pdAAreverse[MAX_PEPTIDE_LEN];      // Stores n-term fragment ion fragment ladder calc.; sum AA masses including mods
   double             _pdAAforwardDecoy[MAX_PEPTIDE_LEN]; // Stores fragment ion fragment ladder calc.; sum AA masses including mods
   double             _pdAAreverseDecoy[MAX_PEPTIDE_LEN]; // Stores n-term fragment ion fragment ladder calc.; sum AA masses including mods
   int                _iSizepcVarModSites;
   VarModInfo         _varModInfo;
   ProteinInfo        _proteinInfo;

   unsigned int       _uiBinnedIonMasses[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN];
   unsigned int       _uiBinnedIonMassesDecoy[MAX_FRAGMENT_CHARGE+1][9][MAX_PEPTIDE_LEN];

   static bool *_pbSearchMemoryPool;    // Pool of memory to be shared by search threads
   static bool **_ppbDuplFragmentArr;      // Number of arrays equals number of threads
};

#endif // _COMETSEARCH_H_
