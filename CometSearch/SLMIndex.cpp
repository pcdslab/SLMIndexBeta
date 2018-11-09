#include "SLMIndex.h"
#include "libdivsufsort/divsufsort.h"


/* Global variables */
/* vector that holds the raw data arrays */
DataRaw raw_data;
DataRaw *raw_mods;

bool **_ppbDuplFragmentArray;
bool *_pbSearchMemoryPools;

bool *_pbSignalPool;
UINT **_ppbFirstRowOccurrences;
UINT **_ppbSecondRowOccurrences;
UINT **_ppbFirstRowModOccurrences;
UINT **_ppbSecondRowModOccurrences;

UINT iNumPeps = 0;
UINT iNumModPeps = 0;
UINT iNumSequences = 0;
int split = 0;
int raw_mods_size = 0;


#ifdef PTHREADING
Mutex g_raw_modsMutex;

#endif /* PTHREADING */

static bool AllocateMemory(int maxNumThreads);
static bool DeallocateMemory(int maxNumThreads);
static bool InitiatePopulateStructures(int minNumThreads, int maxNumThreads, int iPercentStart, int iPercentEnd);
void Qsort(int first, int last);
void Qsort2(int first, int last);
void Qsort3(int n, int first, int last);
void quicksort(int *arr, int low, int high);
void UpdateTable(int split, int ith);

#ifdef VARMOD_SEARCH
static bool HandleMods();
static void PopulateModifiedPeptides(unsigned int* spectrum, UINT currEntries, UINT splt, UINT ith);
#endif

static void PopulatePeptides(unsigned int* spectrum, UINT currEntries, UINT ith);
void merge(int a[], int left_low, int left_high, int right_low, int right_high);
void merge_sort(int a[], int low, int high);
void merge_sort(int a[], int length);

static bool AllocateMemory(int maxNumThreads)
{
   bool status = true;
   int i;

   // Must be equal to largest possible array
   int iArraySize = (int)((g_staticParams.options.dPeptideMassHigh + 100.0) * g_staticParams.dInverseBinWidth);

   // Initally mark all arrays as available (i.e. false == not in use)
   _pbSearchMemoryPools = new bool[maxNumThreads];

   for (i=0; i < maxNumThreads; i++)
   {
      _pbSearchMemoryPools[i] = false;
   }

   if (g_staticParams.variableModParameters.bVarModSearch)
   {
      // Initally mark all signals as zero (split not updated)
      _pbSignalPool = new bool[maxNumThreads];

      for (i = 0; i < maxNumThreads; i++)
      {
         _pbSignalPool[i] = false;
      }
   }

   // Allocate array
   _ppbDuplFragmentArray = new bool*[maxNumThreads];

   for (i=0; i < maxNumThreads; i++)
   {
      try
      {
         _ppbDuplFragmentArray[i] = new bool[iArraySize];
      }
      catch (std::bad_alloc& ba)
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg,  " Error - new(_ppbDuplFragmentArr[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         sprintf(szErrorMsg+strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
         logerr(szErrorMsg);
         status = false;
      }
   }

   // Allocate array
   _ppbFirstRowOccurrences = new UINT*[maxNumThreads];

   for (i=0; i < maxNumThreads; i++)
   {
      try
      {
         _ppbFirstRowOccurrences[i] = new UINT[FIRSTROW_CHARS];
         memset(_ppbFirstRowOccurrences[i], 0x0, sizeof(UINT) * FIRSTROW_CHARS);
      }
      catch (std::bad_alloc& ba)
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg,  " Error - new(_ppbFirstRowOccurrences[%d]). bad_alloc: %s.\n", FIRSTROW_CHARS, ba.what());
         sprintf(szErrorMsg+strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
         logerr(szErrorMsg);
         status = false;
      }
   }

   // Allocate array
   _ppbSecondRowOccurrences = new UINT*[maxNumThreads];

   for (i=0; i < maxNumThreads; i++)
   {
      try
      {
         _ppbSecondRowOccurrences[i] = new UINT[MAX_RANGE];
         memset(_ppbSecondRowOccurrences[i], 0x0, sizeof(UINT) * MAX_RANGE);
      }
      catch (std::bad_alloc& ba)
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg,  " Error - new(_ppbSecondRowOccurrences[%d]). bad_alloc: %s.\n", MAX_RANGE, ba.what());
         sprintf(szErrorMsg+strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
         logerr(szErrorMsg);
         status = false;
      }
   }

#if (PREPROCESS == 1)

   if (status == true)
   {
      status = AllocateArraysForPeptides(iNumSequences);
   }

   if (status == true)
   {
      try
      {
         pepEntries = new pepEntry[iNumSequences];
      }
      catch (std::bad_alloc& ba)
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg, " Error - new(pepEntries[%d]). bad_alloc: %s.\n", iNumSequences, ba.what());
         sprintf(szErrorMsg + strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
         logerr(szErrorMsg);
         status = false;
      }
   }

   if (g_staticParams.variableModParameters.bVarModSearch && status == true)
   {
      status = HandleMods();

      // Allocate array
         _ppbFirstRowModOccurrences = new UINT*[maxNumThreads];

         for (i=0; i < maxNumThreads; i++)
         {
            try
            {
               _ppbFirstRowModOccurrences[i] = new UINT[FIRSTROW_CHARS];
               memset(_ppbFirstRowModOccurrences[i], 0x0, sizeof(UINT) * FIRSTROW_CHARS);
            }
            catch (std::bad_alloc& ba)
            {
               char szErrorMsg[256];
               sprintf(szErrorMsg,  " Error - new(_ppbFirstRowModOccurrences[%d]). bad_alloc: %s.\n", FIRSTROW_CHARS, ba.what());
               sprintf(szErrorMsg+strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
               logerr(szErrorMsg);
               status = false;
            }
         }

         // Allocate array
         _ppbSecondRowModOccurrences = new UINT*[maxNumThreads];

         for (i=0; i < maxNumThreads; i++)
         {
            try
            {
               _ppbSecondRowModOccurrences[i] = new UINT[MAX_RANGE];
               memset(_ppbSecondRowModOccurrences[i], 0x0, sizeof(UINT) * MAX_RANGE);
            }
            catch (std::bad_alloc& ba)
            {
               char szErrorMsg[256];
               sprintf(szErrorMsg,  " Error - new(_ppbSecondRowModOccurrences[%d]). bad_alloc: %s.\n", MAX_RANGE, ba.what());
               sprintf(szErrorMsg+strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
               logerr(szErrorMsg);
               status = false;
            }
         }
   }

#endif /* PREPROCESS == 1 */

   return status;
}

static bool DeallocateMemory(int maxNumThreads)
{
   int i;

   /* Delete the memory allocated by preprocess stage */
   delete [] _pbSearchMemoryPools;
   delete [] _pbSignalPool;

   for (i = 0; i < maxNumThreads; i++)
   {
      delete[] _ppbDuplFragmentArray[i];
   }
   delete[] _ppbDuplFragmentArray;

   if (g_staticParams.variableModParameters.bVarModSearch)
   {
      for (i = 0; i < maxNumThreads; i++)
      {
         delete[] _ppbFirstRowModOccurrences[i];
      }
      delete[] _ppbFirstRowModOccurrences;
   }

   for (i = 0; i < maxNumThreads; i++)
   {
      delete[] _ppbFirstRowOccurrences[i];
   }
   delete[] _ppbFirstRowOccurrences;

   if (g_staticParams.variableModParameters.bVarModSearch)
   {
      for (i = 0; i < maxNumThreads; i++)
      {
         delete[] _ppbSecondRowModOccurrences[i];
      }
      delete[] _ppbSecondRowModOccurrences;
   }

   for (i=0; i<maxNumThreads; i++)
   {
      delete [] _ppbSecondRowOccurrences[i];
   }
   delete [] _ppbSecondRowOccurrences;

   return true;
}

static bool HandleMods()
{
   bool status = true;
   int entries = g_staticParams.iNumberOfMods;
   int allowed = (int) (ALLOWED_SEQUENCES);

   if (g_staticParams.iNumberOfMods != 0)
   {
      raw_mods = new DataRaw[((g_staticParams.iNumberOfMods) / (ALLOWED_SEQUENCES)) + 1];
      bwa_mods = new DataBWT[((g_staticParams.iNumberOfMods) / (ALLOWED_SEQUENCES)) + 1];

      varModEntries = new varModEntry*[((g_staticParams.iNumberOfMods) / (ALLOWED_SEQUENCES)) + 1];

      for (UINT num = 0; num < ((g_staticParams.iNumberOfMods) / (ALLOWED_SEQUENCES)) + 1; num++)
      {
         DataRaw *raw = (raw_mods + num);
         UINT size = 0;

         if (entries - allowed > 0)
         {
            size = (ALLOWED_SEQUENCES);
            entries -= (ALLOWED_SEQUENCES);
         }
         else
         {
            size = entries;
         }

         /* Allocate memory for modified peptides */
         try
         {
            raw->dataArr = new Peak[size * NUM_PEAKS];
            raw->upByte = new Peak[size * NUM_PEAKS];
            raw->num_entries = 0;
            raw->isFull = false;
         } catch (std::bad_alloc& ba)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg, " Error - new(raw_mods->dataArr[%d]). bad_alloc: %s.\n", size, ba.what());
            sprintf(szErrorMsg + strlen(szErrorMsg), "SLM BETA ran out of memory. Look into \"spectrum_batch_size\"\n");
            sprintf(szErrorMsg + strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            status = false;
         }

         try
         {
            varModEntries[num] = new varModEntry[size];
         }
         catch (std::bad_alloc& ba)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg, " Error - new(varModEntries[%d]). bad_alloc: %s.\n", size, ba.what());
            sprintf(szErrorMsg + strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
            logerr(szErrorMsg);
            status = false;
         }
         raw_mods_size++;
      }
   }

   else
   {
      status = false;
   }

   return status;
}

static bool InitiatePopulateStructures(int minNumThreads, int maxNumThreads, int iPercentStart, int iPercentEnd)
{
   bool status = true;
   sDBEntry *dbe;
   char szBuf[SIZE_BUF];
   FILE *fptr;
   int iTmpCh;
   long lEndPos;
   long lCurrPos;
   bool bTrimDescr;
   UINT dbEntryNum = 0;
   iPercentEnd = 100;
   // Create the thread pool containing g_staticParams.options.iNumThreads,
   // each hanging around and sleeping until asked to so a search.
   // NOTE: We don't want to read in ALL the database sequences at once or we
   // will run out of memory for large databases, so we specify a
   // "maxNumParamsToQueue" to indicate that, at most, we will only read in
   // and queue "maxNumParamsToQueue" additional parameters (1 in this case)
   ThreadPool<SearchThreadData *> *pSearchThreadPool = new ThreadPool<SearchThreadData *>(CometSearch::PreprocessThreadProc,
       minNumThreads, maxNumThreads, maxNumThreads /*maxNumParamsToQueue*/);

   /* Initialize the parameters */
   iNumSequences = 0;
   g_staticParams.databaseInfo.uliTotAACount = 0;
   g_staticParams.databaseInfo.iTotalNumProteins = 0;

   /* Read the database file */
   if ((fptr=fopen(g_staticParams.databaseInfo.szDatabase, "rb")) == NULL)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg, " Error - cannot read database file \"%s\".\n", g_staticParams.databaseInfo.szDatabase);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   if (!g_staticParams.options.bOutputSqtStream)
   {
      logout("     - Preprocess progress: ");
      fflush(stdout);
   }

   /* Get the end position of file and rewind the file pointer */
   fseek(fptr, 0, SEEK_END);
   lEndPos=ftell(fptr);
   rewind(fptr);

   // Load database entry header.
   lCurrPos = ftell(fptr);

   /* Count the number of peptide sequences
    * in the database file */
   while (!feof(fptr))
   {
      if(getc(fptr) == '>')
      {
         iNumSequences++;
      }
   }

   rewind(fptr);
   iTmpCh = getc(fptr);

   /* The number of sequences allowed are restricted to 2GB */
   if (iNumSequences > 0 && iNumSequences < (ALLOWED_SEQUENCES))
   {
      /* Allocate memory for protein idx */
      try
      {
         proteins = new sDBEntry[iNumSequences];
      }

      catch (std::bad_alloc& ba)
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg, " Error - new(proteins[%d]). bad_alloc: %s.\n", iNumSequences, ba.what());
         sprintf(szErrorMsg + strlen(szErrorMsg), "SLM BETA ran out of memory.\n");
         logerr(szErrorMsg);
         status = false;
      }

      status = AllocateMemory(g_staticParams.options.iNumThreads);

      if (!status)
      {
         cout << "\nFATAL: PreprocessDatabase: Allocate Memory Failed\n";
      }
   }
   else
   {
      status = false;
      cout << endl << "ABORT: The number of sequences in the database: "<< iNumSequences <<". Allowed limit: 1 - " << (ALLOWED_SEQUENCES) << endl;
   }

   if (status == true)
   {
      // Loop through entire database.
      while (!feof(fptr))
      {
         dbe = (proteins + dbEntryNum);
         dbe->strName = "";
         dbe->strSeq = "";
         dbe->iLenPep = 0;

         // Expect a '>' for sequence header line.
         if (iTmpCh != '>')
         {
            string strError = " Error - database file, expecting definition line here.\n";
            string strFormatError = strError + "\n";
            g_cometStatus.SetStatus(CometResult_Failed, strError);
            logerr(strFormatError.c_str());
            if (fgets(szBuf, SIZE_BUF, fptr) != NULL)
            {
               sprintf(szBuf, " %c%s", iTmpCh, szBuf);
               logerr(szBuf);
            }
            if (fgets(szBuf, SIZE_BUF, fptr) != NULL)
            {
               sprintf(szBuf, " %c%s", iTmpCh, szBuf);
               logerr(szBuf);
            }
            if (fgets(szBuf, SIZE_BUF, fptr) != NULL)
            {
               sprintf(szBuf, " %c%s", iTmpCh, szBuf);
               logerr(szBuf);
            }
            return false;
         }

         bTrimDescr = 0;

         if (!g_staticParams.options.bOutputSqtStream && !(iNumPeps % 200))
         {
            char szTmp[128];
            lCurrPos = ftell(fptr);
            // go from iPercentStart to iPercentEnd, scaled by lCurrPos/iEndPos
            sprintf(szTmp, "%3d%%",
                  (int) (((double) (iPercentStart
                        + (iPercentEnd - iPercentStart) * (double) lCurrPos / (double) lEndPos))));
            logout(szTmp);
            fflush(stdout);
            logout("\b\b\b\b");
         }

         while (((iTmpCh = getc(fptr)) != '\n') && (iTmpCh != '\r') && (iTmpCh != EOF))
         {
            // Don't bother storing description text past first blank.
            if (isspace(iTmpCh) || iscntrl(iTmpCh))
               bTrimDescr = 1;

            if (!bTrimDescr && dbe->strName.size() < (WIDTH_REFERENCE - 1))
               dbe->strName += iTmpCh;
         }

         /* HM: The point where the Protein starts */
         dbe->iSeqFilePosition = ftell(fptr);  // grab sequence file position here

         // Load sequence.
         while (((iTmpCh = getc(fptr)) != '>') && (iTmpCh != EOF))
         {
            if (33 <= iTmpCh && iTmpCh <= 126) // ASCII physical character range.
            {
               // Convert all sequences to upper case.
               dbe->strSeq += toupper(iTmpCh);
               /* Get the length of sequence */
               dbe->iLenPep++;
               g_staticParams.databaseInfo.uliTotAACount++;
            }
         }

         if (dbe->iLenPep <= MAX_PEPTIDE_LEN && dbe->iLenPep >= MIN_PEPTIDE_LEN)
         {
            double dCalcPepMass = g_staticParams.precalcMasses.dOH2ProtonCtermNterm;
            char *szProteinSeq = (char *)dbe->strSeq.c_str();

            for (int i = 0; i <= dbe->iLenPep - 1; i++)
            {
               dCalcPepMass += g_staticParams.massUtility.pdAAMassParent[(int) szProteinSeq[i]];
            }

            int inRange = CheckMassTolerance(dCalcPepMass);

            if (inRange != -1)
            {
               pepEntry *entry = (pepEntries + dbEntryNum);
               *entry = dCalcPepMass;
               iNumPeps++;
               dbEntryNum++;
               g_staticParams.databaseInfo.iTotalNumProteins++;
            }
         }
      }

      omp_set_num_threads(maxNumThreads);

      /* call Qsort  */
      /* The sorting is done in a parallel region */
      #pragma omp parallel
      {
        /* But we only want to sort the list once, so the first call
        * to Qsort is done only once thanks to the single parameter
        */
        #pragma omp single
          Qsort(0, (iNumPeps - 1));
      }

      iPercentStart = 0;
      iPercentEnd = 100;


      for (UINT i = 0; i < iNumPeps; i++)
      {
         sDBEntry *database = proteins + i;

         pSearchThreadPool->WaitForQueuedParams();

         SearchThreadData *pSearchThreadData = new SearchThreadData(*database);
         pSearchThreadData->sdbEntryNum = i;
         pSearchThreadPool->Launch(pSearchThreadData);

         if (!g_staticParams.options.bOutputSqtStream && !(i % (iNumPeps/100)))
         {
            char szTmp[128];
            lCurrPos = ftell(fptr);
            // go from iPercentStart to iPercentEnd, scaled by lCurrPos/iEndPos
            sprintf(szTmp, "%3d%%",
                  (int) (((double) (iPercentStart
                        + (iPercentEnd - iPercentStart) * (double) i / (double) iNumPeps))));
            logout(szTmp);
            fflush(stdout);
            logout("\b\b\b\b");
         }
         status = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();

         if (!status)
         {
            break;
         }
      }
   }

   // Check for errors one more time since there might have been an error
   // while we were waiting for the threads.
   if (status)
   {
      status = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();
   }

   // Wait for active search threads to complete processing.
   pSearchThreadPool->WaitForThreads();

   raw_data.num_entries = iNumPeps * NUM_PEAKS;

   delete pSearchThreadPool;
   pSearchThreadPool = NULL;

   fclose(fptr);

   if (!g_staticParams.options.bOutputSqtStream)
   {
      char szTmp[12];
      sprintf(szTmp, "%3d%%\n", iPercentEnd);
      logout(szTmp);
      fflush(stdout);
   }

   return status;
}

void CometSearch::PreprocessThreadProc(SearchThreadData *pSearchThreadData)
{
   // Grab available array from shared memory pool.
   int i;
   Threading::LockMutex(g_searchMemoryPoolMutex);
   for (i=0; i < g_staticParams.options.iNumThreads; i++)
   {
      if (!_pbSearchMemoryPools[i])
      {
         _pbSearchMemoryPools[i] = true;
         break;
      }
   }
   Threading::UnlockMutex(g_searchMemoryPoolMutex);

   // Fail-safe to stop if memory isn't available for the next thread.
   // Needs better capture and return?
   if (i == g_staticParams.options.iNumThreads)
   {
      printf("Error with memory pool.\n");
      exit(1);
   }

   // Give memory manager access to the thread.
   pSearchThreadData->pbSearchMemoryPool = &_pbSearchMemoryPools[i];

   CometSearch sqSearch;
   // DoSearch now returns true/false, but we already log errors and set
   // the global error variable before we get here, so no need to check
   // the return value here.

   sqSearch.sdbEntryNum = pSearchThreadData->sdbEntryNum;
   sqSearch.ith = i;
   sqSearch.PreprocessDatabase(pSearchThreadData->dbEntry, _ppbDuplFragmentArray[i]);

   delete pSearchThreadData;
   pSearchThreadData = NULL;
}

static bool Burrows_Wheeler_Transform()
{
   bool status = false;
   UINT *counter = new UINT[FIRSTROW_CHARS];
   UINT *occs = new UINT[MAX_RANGE];
   int threads = g_staticParams.options.iNumThreads;

   UINT *san = new UINT[raw_data.num_entries];

   if (raw_data.num_entries != 0)
   {
      UINT *SA = new UINT [raw_data.num_entries];

      status = divsufsort(raw_data.dataArr, (int *)SA, raw_data.num_entries);

      if (status == 0)
      {
         status = true;
         bwa_data.BWTposition = SA;
      }
      else
      {
         status = false;
         delete[] SA;
      }

      if (status == false)
      {
         cout << "FATAL: divsufsort failed\n";
         return status;
      }

      counter[0] = bwa_data.BWTfirstRow.occs[0];

      for (unsigned int j = 1; j < (FIRSTROW_CHARS); j++)
      {
         counter[j] = (UINT) bwa_data.BWTfirstRow.occs[j] + counter[j - 1];
      }

      cnt_data = new UINT[MAX_RANGE];

      cnt_data[0] = bwa_data.BWTSecondRow.occs[0];

      for (unsigned int j = 1; j < MAX_RANGE; j++)
      {
         cnt_data[j] = (UINT) bwa_data.BWTSecondRow.occs[j] + cnt_data[j - 1];
      }

#ifdef DEBUG
      for (unsigned int j = 1; j < (MAX_RANGE); j++)
      {
         if (cnt_data[j] < cnt_data[j-1] || cnt_data[j] > raw_data.num_entries)
         {
            cout << "LAKH DI LAANAT";
            cout << endl<< j;
            cout << "\t" <<cnt_data[j]<<endl;
         }
      }
#endif

      memset(occs, 0x0, MAX_RANGE * sizeof(UINT));

      UINT idd;
      UINT idx = 0;

      for (idd = 0; idd < FIRSTROW_CHARS; idd++)
      {
         for (idx = (idd == 0) ? 0 : counter[idd - 1]; idx < counter[idd]; idx++)
         {
#ifdef DEBUG
            if (raw_data.dataArr[bwa_data.BWTposition[idx]] != idd)
            {
               cout << "PAINCHOD!!!!" <<endl;
               exit(0);
            }
#endif
            UCHAR second = raw_data.upByte[bwa_data.BWTposition[idx]];
            USHORT val = second * 256 + idd;
            if (val == 0)
            {
               san[occs[val]] = bwa_data.BWTposition[idx] / NUM_PEAKS;
            }
            else
            {
               san[cnt_data[val - 1] + occs[val]] = bwa_data.BWTposition[idx] / NUM_PEAKS;
            }

            occs[val]++;
         }
      }
      /* Delete the raw arrays as no longer needed */
      delete[] raw_data.dataArr;
      delete[] raw_data.upByte;
      raw_data.isFull = true;
      raw_data.dataArr = NULL;
      raw_data.upByte = NULL;

      delete[] bwa_data.BWTposition;

      bwa_data.BWTposition = san;

   }

   if (g_staticParams.variableModParameters.bVarModSearch)
   {
      for (int n = 0; n < raw_mods_size; n++)
      {
         if (raw_mods[n].num_entries != 0)
         {
            UINT *san_mods = new UINT[raw_mods[n].num_entries];

            UINT *SA2 = new UINT[raw_mods[n].num_entries];

            status = divsufsort(raw_mods[n].dataArr, (int *) SA2, raw_mods[n].num_entries);

            if (status == 0)
            {
               status = true;
               bwa_mods[n].BWTposition = SA2;
            }
            else
            {
               status = false;
               delete[] SA2;
            }

            if (status == false)
            {
               cout << "FATAL: divsufsort - varMods failed\n";
               return status;
            }

            counter[0] = bwa_mods[n].BWTfirstRow.occs[0];

            for (unsigned int j = 1; j < (FIRSTROW_CHARS); j++)
            {
               counter[j] = (UINT) bwa_mods[n].BWTfirstRow.occs[j] + counter[j - 1];
            }

            UINT *cnt = new UINT[MAX_RANGE];

            cnt[0] = bwa_mods[n].BWTSecondRow.occs[0];

            for (unsigned int j = 1; j < (MAX_RANGE); j++)
            {
               cnt[j] = (UINT) bwa_mods[n].BWTSecondRow.occs[j] + cnt[j - 1];
            }

            memset(occs, 0x0, MAX_RANGE * sizeof(UINT));

            UINT idd;
            UINT idx = 0;
            for (idd = 0; idd < 256; idd++)
            {
               for (idx = (idd == 0) ? 0 : counter[idd - 1]; idx < counter[idd]; idx++)
               {
#ifdef DEBUG
                  if (raw_mods[n].dataArr[bwa_mods[n].BWTposition[idx]] != idd)
                  {
                     cout << "Please don't do this!!!" << endl;
                     cout << idd << "\t" << counter[idd - 1] << "\t" << counter[idd] << "\t\n";
                     exit(0);
                  }
#endif
                  UCHAR second = raw_mods[n].upByte[bwa_mods[n].BWTposition[idx]];
                  USHORT val = second * 256 + idd;
                  if (val == 0)
                  {
                     san_mods[occs[val]] = bwa_mods[n].BWTposition[idx] / NUM_PEAKS;
                  }
                  else
                  {
                     san_mods[cnt[val - 1] + occs[val]] = bwa_mods[n].BWTposition[idx] / NUM_PEAKS;
                  }

                  occs[val]++;
               }
            }
            /* Free the raw_mods memory and mark it full */
            delete[] raw_mods[n].dataArr;
            delete[] raw_mods[n].upByte;
            raw_mods[n].isFull = true;
            raw_mods[n].dataArr = NULL;
            raw_mods[n].upByte = NULL;
            delete[] bwa_mods[n].BWTposition;

            bwa_mods[n].BWTposition = san_mods;

            cnt_mods.push_back(cnt);

         }
      }
   }

   delete[] occs;
   delete[] counter;

   omp_set_num_threads(g_staticParams.options.iNumThreads);

   cout << " Sorting Data Buckets\n";

#pragma omp parallel for schedule(static,256)
   for (UINT t = 1; t < MAX_RANGE; t++)
   {

      if (cnt_data[t] <= raw_data.num_entries - 1)
      {
         if (cnt_data[t] - cnt_data[t - 1] > 1)
         {
#ifdef DEBUG
            cout << "t:\t" << t << endl;
            cout << "cnt_data[t-1]:\t" << cnt_data[t - 1] << endl;
            cout << "cnt_data[t]:\t" << cnt_data[t] << endl;
#endif
            int stt = (int)(cnt_data[t - 1]);

            if (stt <= raw_data.num_entries -1)
            {
               quicksort((int *) bwa_data.BWTposition, stt, (int) cnt_data[t] - 1);
            }
         }
      }
   }

   if (g_staticParams.variableModParameters.bVarModSearch)
   {
      for (int s = 0; s < raw_mods_size; s++)
      {
         if (raw_mods[s].num_entries > 0)
         {
#pragma omp parallel for schedule(static,256)
            for (UINT t = 1; t < MAX_RANGE; t++)
            {
               if (cnt_mods[s][t] <= raw_mods[s].num_entries - 1)
               {
                  if (cnt_mods[s][t] - cnt_mods[s][t - 1] > 2)
                  {
                     int stt = (int) (cnt_mods[s][t - 1]);

                     if (stt <= raw_mods[s].num_entries - 1)
                     {
                        quicksort((int *) bwa_mods[s].BWTposition, (int) cnt_mods[s][t - 1], (int) cnt_mods[s][t] - 1);
                     }
                  }
               }
            }
         }
	  }
   }

   return status;
}

static void PopulatePeptides(unsigned int* spectrum, UINT currEntries, UINT ii)
{
   unsigned idd = 0;
   int id = 0;

   merge_sort((int*)spectrum, PP_ARRAYSIZE/2);
   merge_sort((int*)(spectrum + PP_ARRAYSIZE/2), PP_ARRAYSIZE/2);

   for (id = 0; id < PP_ARRAYSIZE/2 && idd < NUM_PEAKS/2; id++)
   {
      if (spectrum[id] != 0)
      {
         Peak delta = (Peak) (spectrum[id] % FIRSTROW_CHARS);
         raw_data.dataArr[currEntries + idd] = delta;
         raw_data.upByte[currEntries + idd] = (Peak) (spectrum[id] / FIRSTROW_CHARS);
         idd++;

         /* Add the number of occurence into the occ array */
         _ppbFirstRowOccurrences[ii][delta]++;
         _ppbSecondRowOccurrences[ii][spectrum[id] & 0xFFFF]++;
      }
   }

   for (id = PP_ARRAYSIZE/2; id < PP_ARRAYSIZE && idd < NUM_PEAKS; id++)
   {
      if (spectrum[id] != 0)
      {
         Peak delta = (Peak) (spectrum[id] % FIRSTROW_CHARS);
         raw_data.dataArr[currEntries + idd] = delta;
         raw_data.upByte[currEntries + idd] = (Peak) (spectrum[id] / FIRSTROW_CHARS);
         idd++;

         /* Add the number of occurence into the occ array */
         _ppbFirstRowOccurrences[ii][delta]++;
         _ppbSecondRowOccurrences[ii][spectrum[id] & 0xFFFF]++;
      }
   }

   for (; idd < NUM_PEAKS; idd++)
   {
      raw_data.dataArr[currEntries + idd] = 0;
      raw_data.upByte[currEntries + idd] = 0;

      /* Add the number of occurence into the occ array */
      _ppbFirstRowOccurrences[ii][0]++;
      _ppbSecondRowOccurrences[ii][0]++;
   }

   return;
}

#ifdef VARMOD_SEARCH
static void PopulateModifiedPeptides(unsigned int* spectrum, UINT currEntries, UINT splt, UINT ii)
{
   unsigned idd = 0;
   int id = 0;

   merge_sort((int*)spectrum, PP_ARRAYSIZE/2);
   merge_sort((int*)(spectrum + PP_ARRAYSIZE/2), PP_ARRAYSIZE/2);

   for (id = 0; id < PP_ARRAYSIZE/2 && idd < NUM_PEAKS/2; id++)
   {
      if (spectrum[id] != 0)
      {
         Peak delta = (Peak) (spectrum[id] % FIRSTROW_CHARS);
         raw_mods[splt].dataArr[currEntries + idd] = delta;
         raw_mods[splt].upByte[currEntries + idd] = (Peak) (spectrum[id] / FIRSTROW_CHARS);
         idd++;

         /* Add the number of occurence into the occ array */
         _ppbFirstRowModOccurrences[ii][delta]++;
         _ppbSecondRowModOccurrences[ii][spectrum[id] & 0xFFFF]++;
      }
   }

   for (id = PP_ARRAYSIZE/2; id < PP_ARRAYSIZE && idd < NUM_PEAKS; id++)
   {
      if (spectrum[id] != 0)
      {
         Peak delta = (Peak) (spectrum[id] % FIRSTROW_CHARS);
         raw_mods[splt].dataArr[currEntries + idd] = delta;
         raw_mods[splt].upByte[currEntries + idd] = (Peak) (spectrum[id] / FIRSTROW_CHARS);
         idd++;

         /* Add the number of occurence into the occ array */
         _ppbFirstRowModOccurrences[ii][delta]++;
         _ppbSecondRowModOccurrences[ii][spectrum[id] & 0xFFFF]++;
      }
   }

   for (; idd < NUM_PEAKS; idd++)
   {
      raw_mods[splt].dataArr[currEntries + idd] = 0;
      raw_mods[splt].upByte[currEntries + idd] = 0;

      /* Add the number of occurence into the occ array */
      _ppbFirstRowModOccurrences[ii][0]++;
      _ppbSecondRowModOccurrences[ii][0]++;
   }


   return;
}
#endif

/* Allocate Memory for Peptides */
bool AllocateArraysForPeptides(UINT peptides)
{
   bool status = true;

   /* Allocate memory for regular peptides */
   try
   {
      raw_data.dataArr = new Peak[peptides * NUM_PEAKS];
      raw_data.upByte =  new Peak[peptides * NUM_PEAKS];
      raw_data.num_entries = 0;
      raw_data.isFull = false;
   }

   catch(std::bad_alloc& ba)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg,  " Error - new(iMSQPeptide[%d]). bad_alloc: %s.\n", peptides, ba.what());
      sprintf(szErrorMsg+strlen(szErrorMsg), "Comet ran out of memory. Look into \"spectrum_batch_size\"\n");
      sprintf(szErrorMsg+strlen(szErrorMsg), "parameters to address mitigate memory use.\n");
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      status = false;
   }

   return status;
}

// Break each protein into peptides and preprocess each peptide entry
bool CometSearch::PreprocessDatabase(sDBEntry dbe, bool *pbDuplFragment)
{

   // Standard protein database search.
   if (g_staticParams.options.iWhichReadingFrame == 0)
   {
      _proteinInfo.iProteinSeqLength = dbe.iLenPep;
      _proteinInfo.iSeqFilePosition = dbe.iSeqFilePosition;

      if (!GenerateAndPopulatePeptides((char *)dbe.strSeq.c_str(),
                             (char *)dbe.strName.c_str(),
                             0,
                             pbDuplFragment))
      {
         return false;
      }

      if (g_staticParams.options.bClipNtermMet && dbe.strSeq[0]=='M')
      {
         _proteinInfo.iProteinSeqLength -= 1;

         if (!GenerateAndPopulatePeptides((char *)dbe.strSeq.c_str()+1,
                                (char *)dbe.strName.c_str(),
                                1,
                                pbDuplFragment))
         {
            return false;
         }
      }
   }

   return true;
}

// Do the stuff around here
bool CometSearch::GenerateAndPopulatePeptides(char *szProteinSeq,
                                              char *szProteinName,
                                              bool bNtermPeptideOnly,
                                              bool *pbDuplFragment)
{
   int iStartPos = 0;
   int iEndPos = 0;
   int piVarModCounts[VMODS] = {0};
   int iWhichIonSeries;
   int ctIonSeries;
   int ctLen;
   int ctCharge;
   double dCalcPepMass = 0.0;


#ifndef VARMOD_SEARCH
   iMSQ_UNUSUED_PARAM(szProteinName);
#endif

   /* HM: TODO: Don't know the actual position of the peptide in
    *     a protein anymore so can't add any C/N Protein termini
    *     modifications
    */
   /*
   if (iStartPos == 0)
      dCalcPepMass += g_staticParams.staticModifications.dAddNterminusProtein;
   if (iEndPos == iProteinSeqLengthMinus1)
      dCalcPepMass += g_staticParams.staticModifications.dAddCterminusProtein;
   */

   /* HM: We are getting peptide */
   iEndPos = _proteinInfo.iProteinSeqLength - 1;

   //iLenPeptide = iEndPos-iStartPos+1;

   if (_proteinInfo.iProteinSeqLength > 0)
   {
      dCalcPepMass = g_staticParams.precalcMasses.dOH2ProtonCtermNterm;

      for (int i = iStartPos; i <= iEndPos; i++)
      {
         dCalcPepMass += g_staticParams.massUtility.pdAAMassParent[(int) szProteinSeq[i]];

         if (g_staticParams.variableModParameters.bVarModSearch)
         {
            CountVarMods(piVarModCounts, szProteinSeq[i], i);
         }
      }
   }

   UINT currEntries = sdbEntryNum * NUM_PEAKS;

#if (PREPROCESS == 1)

   /* TODO: This can be revisited in order to add more ion-series and charge states.
    * Private Linear Array to hold Database preprocess information */
   /* 3 (MAX_FRAGMENT_CHARGE) * 2 (NUMBER OF ION SERIES B|Y) * 64 (MAX PEPTIDE LENGTH) */
   unsigned int _uiBinnedIonsPreprocessing[PP_ARRAYSIZE] = { 0 };

   /* Generate Theoretical Spectrum here
    Calculate ion series.
    */
   int iLenMinus1 = iEndPos - iStartPos; // Equals iLenPeptide minus 1.
   int i;
   double dBion = g_staticParams.precalcMasses.dNtermProton;
   double dYion = g_staticParams.precalcMasses.dCtermOH2Proton;

   /* HM: TODO: Don't know the actual position of the peptide in
    *     a protein anymore so can't add any C/N Protein termini
    *     modifications

    if (iStartPos == 0)
    dBion += g_staticParams.staticModifications.dAddNterminusProtein;
    if (iEndPos == iProteinSeqLengthMinus1)
    dYion += g_staticParams.staticModifications.dAddCterminusProtein;*/

   for (i = iStartPos; i <= iEndPos; i++)
   {
      dBion += g_staticParams.massUtility.pdAAMassFragment[(int) szProteinSeq[i]];
      _pdAAforward[i] = dBion;

      dYion += g_staticParams.massUtility.pdAAMassFragment[(int) szProteinSeq[iEndPos - i]];
      _pdAAreverse[i] = dYion;
   }

   // Now get the set of binned fragment ions once to compare this peptide against all matching spectra.
   for (ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
   {
      iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];
      for (ctCharge = 1; ctCharge <= g_massRange.iMaxFragmentCharge; ctCharge++)
      {
         for (ctLen = 0; ctLen < iLenMinus1; ctLen++)
            pbDuplFragment[BIN(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse))] =
                  false;
      }
   }

   for (ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
   {
      iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];

      for (ctCharge = 1; ctCharge <= g_massRange.iMaxFragmentCharge; ctCharge++)
      {
         // As both _pdAAforward and _pdAAreverse are increasing, loop through
         // iLenPeptide-1 to complete set of internal fragment ions.
         for (ctLen = 0; ctLen < iLenMinus1; ctLen++)
         {
            int iVal = BIN(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse));

            if (pbDuplFragment[iVal] == false)
            {
               _uiBinnedIonsPreprocessing[ctIonSeries * 192 + (ctCharge - 1) * 64 + ctLen] = iVal;
               pbDuplFragment[iVal] = true;
            }
            else
            {
               _uiBinnedIonsPreprocessing[ctIonSeries * 192 + (ctCharge - 1) * 64 + ctLen] = 0;
            }
         }
      }
   }

   /* Populate the data into idx */
   (void) PopulatePeptides(_uiBinnedIonsPreprocessing, currEntries, ith);
#endif

   // If any variable mod mass is negative, consider adding to iEndPos as long
   // as peptide minus all possible negative mods is less than the dMaxMass????
   //
   // Otherwise, at this point, peptide mass is too big which means should be ok for varmod search.

   if (g_staticParams.variableModParameters.bVarModSearch)
   {
      if (HasVariableMod(piVarModCounts, iStartPos, iEndPos))
      {
         VarModSearch(szProteinSeq, szProteinName, piVarModCounts, iStartPos, iEndPos, pbDuplFragment);
      }
   }

   if (bNtermPeptideOnly)
   {
      return true;
   }

   return true;
}

/* HM: This API is analogous to the GenerateAndPopulatePeptides for modified peptides */
// FIX: 'false' is never returned by this function, why?
bool CometSearch::GenerateAndPopulateModPeptide(char *szProteinSeq,
                                                int iWhichQuery,
                                                bool *pbDuplFragment)
{
   char pcVarModSites[MAX_PEPTIDE_LEN_P2];
   int piVarModCharIdx[VMODS];
   int ctIonSeries;
   int ctLen;
   int ctCharge;
   int iWhichIonSeries;
   int i;
   int j;

   iMSQ_UNUSUED_PARAM(iWhichQuery);

   // at this point, need to compare current modified peptide
   // against all relevant entries

   // but first, calculate modified peptide mass as it could've changed
   // by terminating earlier than start/end positions defined in VarModSearch()
   double dCalcPepMass = g_staticParams.precalcMasses.dNtermProton + g_staticParams.precalcMasses.dCtermOH2Proton - PROTON_MASS;

   int iLenMinus1 = _varModInfo.iEndPos - _varModInfo.iStartPos;     // equals iLenPeptide-1
   int iLenPeptide = iLenMinus1+1;

   // contains positional coding of a variable mod at each idx which equals an AA residue
   memset(pcVarModSites, 0, _iSizepcVarModSites);
   memset(piVarModCharIdx, 0, sizeof(piVarModCharIdx));

   // deal with n-term mod
   for (j=0; j<VMODS; j++)
   {
      if ( strchr(g_staticParams.variableModParameters.varModList[j].szVarModChar, 'n')
            && !isEqual(g_staticParams.variableModParameters.varModList[j].dVarModMass, 0.0)
            && (_varModInfo.varModStatList[j].iMatchVarModCt > 0) )
      {
         if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
         {
            if (pcVarModSites[iLenPeptide] != 0)  // conflict in two variable mods on n-term
               return true;

            // store the modification number at modification position
            pcVarModSites[iLenPeptide] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
            dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;;
         }
         piVarModCharIdx[j]++;
      }
   }

   // deal with c-term mod
   for (j=0; j<VMODS; j++)
   {
      if ( strchr(g_staticParams.variableModParameters.varModList[j].szVarModChar, 'c')
            && !isEqual(g_staticParams.variableModParameters.varModList[j].dVarModMass, 0.0)
            && (_varModInfo.varModStatList[j].iMatchVarModCt > 0) )
      {
         if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
         {
            if (pcVarModSites[iLenPeptide+1] != 0)  // conflict in two variable mods on c-term
               return true;

            // store the modification number at modification position
            pcVarModSites[iLenPeptide+1] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
            dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;
         }
         piVarModCharIdx[j]++;
      }
   }



   // Generate pdAAforward for _pResults[0].szPeptide
   for (i=_varModInfo.iStartPos; i<=_varModInfo.iEndPos; i++)
   {
      int iPos = i - _varModInfo.iStartPos;

      dCalcPepMass += g_staticParams.massUtility.pdAAMassParent[(int)szProteinSeq[i]];

      // This loop is where all individual variable mods are combined
      for (j=0; j<VMODS; j++)
      {
         if (!isEqual(g_staticParams.variableModParameters.varModList[j].dVarModMass, 0.0)
               && (_varModInfo.varModStatList[j].iMatchVarModCt > 0)
               && strchr(g_staticParams.variableModParameters.varModList[j].szVarModChar, szProteinSeq[i]))
         {
            if (g_staticParams.variableModParameters.varModList[j].iVarModTermDistance == -1)
            {
               if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
               {
                  if (pcVarModSites[iPos] != 0)  // conflict in two variable mods on same residue
                     return true;

                  // store the modification number at modification position
                  pcVarModSites[iPos] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
                  dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;
               }
               piVarModCharIdx[j]++;
            }
            else  // terminal distance constraint specified
            {
               if (g_staticParams.variableModParameters.varModList[j].iWhichTerm == 0)      // protein N
               {
                  if (i <= g_staticParams.variableModParameters.varModList[j].iVarModTermDistance)
                  {
                     if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
                     {
                        if (pcVarModSites[iPos] != 0)  // conflict in two variable mods on same residue
                           return true;

                        // store the modification number at modification position
                        pcVarModSites[iPos] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
                        dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;
                     }
                     piVarModCharIdx[j]++;
                  }
               }
               else if (g_staticParams.variableModParameters.varModList[j].iWhichTerm == 1) // protein C
               {
                  if (i + g_staticParams.variableModParameters.varModList[j].iVarModTermDistance >= _proteinInfo.iProteinSeqLength-1)
                  {
                     if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
                     {
                        if (pcVarModSites[iPos] != 0)  // conflict in two variable mods on same residue
                           return true;

                        // store the modification number at modification position
                        pcVarModSites[iPos] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
                        dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;
                     }
                     piVarModCharIdx[j]++;
                  }
               }
               else if (g_staticParams.variableModParameters.varModList[j].iWhichTerm == 2) // peptide N
               {
                  if (iPos <= g_staticParams.variableModParameters.varModList[j].iVarModTermDistance)
                  {
                     if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
                     {
                        if (pcVarModSites[iPos] != 0)  // conflict in two variable mods on same residue
                           return true;

                        // store the modification number at modification position
                        pcVarModSites[iPos] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
                        dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;
                     }
                     piVarModCharIdx[j]++;
                  }
               }
               else if (g_staticParams.variableModParameters.varModList[j].iWhichTerm == 3) // peptide C
               {
                  if (iPos + g_staticParams.variableModParameters.varModList[j].iVarModTermDistance >= iLenMinus1)
                  {
                     if (_varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
                     {
                        if (pcVarModSites[iPos] != 0)  // conflict in two variable mods on same residue
                           return true;

                        // store the modification number at modification position
                        pcVarModSites[iPos] = _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]];
                        dCalcPepMass += g_staticParams.variableModParameters.varModList[j].dVarModMass;
                     }
                     piVarModCharIdx[j]++;
                  }
               }
            }
         }
      }
   }

   // Check to see if required mods are satisfied
   if (g_staticParams.variableModParameters.bRequireVarMod)
   {
      for (j=0; j<VMODS; j++)
      {
         if (g_staticParams.variableModParameters.varModList[j].bRequireThisMod
               && !isEqual(g_staticParams.variableModParameters.varModList[j].dVarModMass, 0.0))
         {
            bool bPresent = false;

            // if mod is required, see if it is present in the peptide
            for (i=_varModInfo.iStartPos; i<=_varModInfo.iEndPos; i++)
            {
               int iPos = i - _varModInfo.iStartPos;

               //FIX:  add in logic to check distance constraints
               if (pcVarModSites[iPos] == _varModInfo.varModStatList[j].iVarModSites[piVarModCharIdx[j]])
               {
                  bPresent = true;
                  break;
               }
            }

            if (!bPresent)
               return true;
         }
      }
   }

   if (dCalcPepMass < g_massRange.dMaxMass)
   {
      Threading::LockMutex(g_raw_modsMutex);

      if ((int) (split * (ALLOWED_SEQUENCES) + iNumModPeps) >= g_staticParams.iNumberOfMods)
      {
         cout << "\n\t" << split << "\t" <<  ALLOWED_SEQUENCES << "\t" << iNumModPeps << "\t" <<g_staticParams.iNumberOfMods << endl;
         cout << endl << " ERROR: The Number of varmods exceeded max_var_mods.\n"
               " Please increase the number of mods and run again\n"
               " Set the macro PREPROCESS = 0 and run to get the number of mods generated\n";
         exit(0);
      }

#if (PREPROCESS == 1)

      if (iNumModPeps == (ALLOWED_SEQUENCES)-1)
      {
         iNumModPeps = 0;
         split++;

         if (g_staticParams.variableModParameters.bVarModSearch)
         {
            for (int kk = 0; kk < g_staticParams.options.iNumThreads; kk++)
            {
               _pbSignalPool[kk] = true;
            }
         }
      }

      if (g_staticParams.variableModParameters.bVarModSearch)
      {
         if (_pbSignalPool[ith] == true)
         {
            UpdateTable(split - 1, ith);
            _pbSignalPool[ith] = false;
         }
      }

      /* Get the varModPeptide entry and populate it */
      varModEntry *entry = (varModEntries[split] + iNumModPeps);
      UINT splt = split;
      UINT currEntries = raw_mods[splt].num_entries;
      raw_mods[splt].num_entries += NUM_PEAKS;

#endif /* PREPROCESS == 1 */
      iNumModPeps++;

      Threading::UnlockMutex(g_raw_modsMutex);

#if (PREPROCESS == 1)
      /*
       * Private Linear Array to hold Database preprocess information
       * 3 (MAX_FRAGMENT_CHARGE) * 2 (NUMBER OF ION SERIES B|Y) * 64 (MAX PEPTIDE LENGTH)
       */
      unsigned int _uiBinnedIonsPreprocessing[PP_ARRAYSIZE] = {0};

      entry->indexProtein = this->sdbEntryNum;
      entry->iPepMass = dCalcPepMass;

      /*
       * Generate the MODSITES bitmask based on pcVarModSites
       */
      UINT mod_nums = 0;
      for (int s = 0; s < MAX_PEPTIDE_LEN_P2; s++)
      {
         /* Check if this bit needs to be set? */
         if (pcVarModSites[s] != 0)
         {
            /* Set the first 64 bits */
            if (s < MAX_PEPTIDE_LEN)
            {
               entry->sites.lower64 |= MODSITE(s);
            }
            /* Set the rest of 2 bits here */
            else
            {
               entry->sites.upper2 |= MODSITE(s-MAX_PEPTIDE_LEN);
            }

            if (mod_nums < 10)
            {
               entry->sites.upper2 |= (pcVarModSites[s] << (3 * mod_nums + 2));
               mod_nums++;
            }
            else
            {
               int temp = 1;
               cout << "SLM BETA: More than 10 mods in the peptide" << endl;
               while (temp);
            }
         }
      }

      double dBion = g_staticParams.precalcMasses.dNtermProton;
      double dYion = g_staticParams.precalcMasses.dCtermOH2Proton;

      if (_varModInfo.iStartPos == 0)
         dBion += g_staticParams.staticModifications.dAddNterminusProtein;
      if (_varModInfo.iEndPos == _proteinInfo.iProteinSeqLength - 1)
         dYion += g_staticParams.staticModifications.dAddCterminusProtein;

      // variable N-term
      if (pcVarModSites[iLenPeptide] > 0)
         dBion += g_staticParams.variableModParameters.varModList[pcVarModSites[iLenPeptide] - 1].dVarModMass;

      // variable C-term
      if (pcVarModSites[iLenPeptide + 1] > 0)
         dYion += g_staticParams.variableModParameters.varModList[pcVarModSites[iLenPeptide + 1] - 1].dVarModMass;

      // Generate pdAAforward for _pResults[0].szPeptide
      for (i = _varModInfo.iStartPos; i < _varModInfo.iEndPos; i++)
      {
         int iPos = i - _varModInfo.iStartPos;

         dBion += g_staticParams.massUtility.pdAAMassFragment[(int) szProteinSeq[i]];

         if (pcVarModSites[iPos] > 0)
            dBion += g_staticParams.variableModParameters.varModList[pcVarModSites[iPos] - 1].dVarModMass;

         _pdAAforward[iPos] = dBion;

         dYion += g_staticParams.massUtility.pdAAMassFragment[(int) szProteinSeq[_varModInfo.iEndPos - i
               + _varModInfo.iStartPos]];

         iPos = _varModInfo.iEndPos - i;
         if (pcVarModSites[iPos] > 0)
            dYion += g_staticParams.variableModParameters.varModList[pcVarModSites[iPos] - 1].dVarModMass;

         _pdAAreverse[i - _varModInfo.iStartPos] = dYion;

      }

      // now get the set of binned fragment ions once for all matching peptides

      // initialize pbDuplFragment here
      for (ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
      {
         iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];

         for (ctCharge = 1; ctCharge <= g_massRange.iMaxFragmentCharge; ctCharge++)
         {
            for (ctLen = 0; ctLen < iLenMinus1; ctLen++)
               pbDuplFragment[BIN(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse))] =
                     false;
         }
      }

      // set pbDuplFragment[bin] to true for each fragment ion bin
      for (ctIonSeries = 0; ctIonSeries < g_staticParams.ionInformation.iNumIonSeriesUsed; ctIonSeries++)
      {
         iWhichIonSeries = g_staticParams.ionInformation.piSelectedIonSeries[ctIonSeries];

         for (ctCharge = 1; ctCharge <= g_massRange.iMaxFragmentCharge; ctCharge++)
         {
            // as both _pdAAforward and _pdAAreverse are increasing, loop through
            // iLenPeptide-1 to complete set of internal fragment ions
            for (ctLen = 0; ctLen < iLenMinus1; ctLen++)
            {
               int iVal = BIN(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse));

               if (pbDuplFragment[iVal] == false)
               {
                  _uiBinnedIonsPreprocessing[ctIonSeries * 192 + (ctCharge - 1) * 64 + ctLen] = iVal;
                  pbDuplFragment[iVal] = true;
               }
               else
                  _uiBinnedIonsPreprocessing[ctIonSeries * 192 + (ctCharge - 1) * 64 + ctLen] = 0;
            }
         }
      }

      /* Call the function to sort the uIBinnedIonsPreprocessing and put the data into idx */
      PopulateModifiedPeptides(_uiBinnedIonsPreprocessing, currEntries, splt, ith);

      #endif

   }

   return true;
}

/*
 * Preprocess the database here
 */
bool CometSearchManager::PreprocessDatabase()
{
   bool status = true;

   /* Initialize the Semaphores for raw and modified peptides */
   Threading::CreateMutex(&g_raw_modsMutex);

   /* Adding temporary Percentages */
   status = InitiatePopulateStructures(g_staticParams.options.iNumThreads, g_staticParams.options.iNumThreads, 0, 100);

   int ii = 0;
   int thread = 0;

   omp_set_num_threads(g_staticParams.options.iNumThreads);

   if (status)
   {
      iNumModPeps += (split) * (ALLOWED_SEQUENCES);

      printf("     - Number of peps = %u\n", iNumPeps);
      printf("     - Number of modified peps = %u\n", iNumModPeps);
   }
   else
   {
      cout << "\nFATAL: PreprocessDatabase: Populate Structures Failed\n";
      exit(0);
   }

   for (thread = 0; thread < g_staticParams.options.iNumThreads; thread++)
   {
      int spll = split;

      if (g_staticParams.variableModParameters.bVarModSearch)
      {
         if (_pbSignalPool[thread] == true)
         {
            spll--;
            _pbSignalPool[thread] = false;
         }
      }
#pragma omp parallel for
      for (ii = 0; ii < FIRSTROW_CHARS; ii++)
      {
         bwa_data.BWTfirstRow.occs[ii] += _ppbFirstRowOccurrences[thread][ii];
         _ppbFirstRowOccurrences[thread][ii] = 0;

         if (g_staticParams.variableModParameters.bVarModSearch)
         {
            bwa_mods[spll].BWTfirstRow.occs[ii] += _ppbFirstRowModOccurrences[thread][ii];
            _ppbFirstRowModOccurrences[thread][ii] = 0;
         }
      }
#pragma omp parallel for
      for (ii = 0; ii < MAX_RANGE; ii++)
      {
         bwa_data.BWTSecondRow.occs[ii] += _ppbSecondRowOccurrences[thread][ii];
         _ppbSecondRowOccurrences[thread][ii] = 0;

         if (g_staticParams.variableModParameters.bVarModSearch)
         {
            bwa_mods[spll].BWTSecondRow.occs[ii] += _ppbSecondRowModOccurrences[thread][ii];
            _ppbSecondRowModOccurrences[thread][ii] = 0;
         }
      }
   }

#if (PREPROCESS == 1)
   if (status)
   {
      if (!g_staticParams.options.bOutputSqtStream)
      {
         logout(" Creating Suffix Arrays \n");
         fflush(stdout);
      }

      status = Burrows_Wheeler_Transform();
   }
#endif /* PREPROCESS == 1 */

   if (status)
   {
      status = DeallocateMemory(g_staticParams.options.iNumThreads);
   }
   else
   {
      cout << "\nFATAL: PreprocessDatabase: Burrows-Wheeler Transform Failed\n";
      exit(0);
   }

   if (!status)
   {
      cout << "\nFATAL: PreprocessDatabase: DeAllocateMemoryFailed\n";
      exit(0);
   }

   if (status)
   {
      if (!g_staticParams.options.bOutputSqtStream)
      {
         logout(" Allocating Memory for Results.. \n");
         fflush(stdout);
      }

      results = new short* [1/*g_staticParams.options.iNumThreads*/];

      UINT size = raw_data.num_entries;
      size /= NUM_PEAKS;

      for (int tt = 0; tt < 1/*g_staticParams.options.iNumThreads*/; tt++)
      {
         results[tt] = new short[size];
         memset(results[tt], 0x0, sizeof(short) * size);
      }

      if (raw_mods_size > 1)
      {
         size = (ALLOWED_SEQUENCES);
      }
      else if (raw_mods_size == 1)
      {
         size = raw_mods[0].num_entries / NUM_PEAKS;
      }
      else
      {
         size = 0;
      }

      if (size > 0)
      {
         results_mods = new short* [1/*g_staticParams.options.iNumThreads*/];
         for (int tt = 0; tt < 1/*g_staticParams.options.iNumThreads*/; tt++)
         {
            results_mods[tt] = new short[size];
            memset(results_mods[tt], 0x0, sizeof(short) * size);
         }
      }
   }

   Threading::DestroyMutex(g_raw_modsMutex);

#if (PREPROCESS == 0)

   exit(0);

#endif /* PREPROCESS == 1 */

   return status;
}

/**********************************************************************************************/
void merge_sort(int a[], int length)
{
   merge_sort(a, 0, length - 1);
}

void merge_sort(int a[], int low, int high)
{
   if (low >= high)                  //Base case: 1 value to sort->sorted
      return;                         //(0 possible only on initial call)
   else
   {
      int mid = (low + high) / 2;       //Approximate midpoint*
      merge_sort(a, low, mid);        //Sort low to mid part of array
      merge_sort(a, mid + 1, high);     //Sort mid+1 to high part of array
      merge(a, low, mid, mid + 1, high); //Merge sorted subparts of array
   }
}

void merge(int a[], int left_low, int left_high, int right_low, int right_high)
{
   int length = right_high - left_low + 1;
   int temp[length];
   int left = left_low;
   int right = right_low;
   for (int i = 0; i < length; ++i)
   {
      if (left > left_high)
         temp[i] = a[right++];
      else if (right > right_high)
         temp[i] = a[left++];
      else if (a[left] <= a[right])
         temp[i] = a[left++];
      else
         temp[i] = a[right++];
   }

   for (int i = 0; i < length; ++i)
      a[left_low++] = temp[i];
}

int CheckMassTolerance(double dCalcPepMass)
{
   unsigned int iPos = -1;

   if (dCalcPepMass >= g_massRange.dMinMass
       && dCalcPepMass <= g_massRange.dMaxMass)
   {
      iPos = 0;
   }

   return iPos;
}

/* Regular quiksort algorithm, with the only exception that
 * the recursive step is done in parallel with openmp tasks
 */
void Qsort(int first, int last) {
  float pivot, temp;
  sDBEntry pivote, tempe;

  int i_pivot, left, right;
  if (first >= last) return; // no need to sort
  // otherwise select a pivot
  i_pivot = (first + last) / 2;
  pivot = pepEntries[i_pivot];
  pivote = proteins[i_pivot];
  left = first;
  right = last;
  while (left <= right) {
    if (pepEntries[left] > pivot) { // swap left element with right element
       temp = pepEntries[left];
       tempe = proteins[left];
       pepEntries[left] = pepEntries[right];
       proteins[left] = proteins[right];

       pepEntries[right] = temp;
       proteins[right] = tempe;
       if (right == i_pivot) {
        i_pivot = left;
       }
       right--;
    }
    else {
      left++;
    }
  }
  // place the pivot in its place (i.e. swap with right element)
  temp = pepEntries[right];
  tempe = proteins[right];

  pepEntries[right] = pivot;
  proteins[right] = pivote;

  pepEntries[i_pivot] = temp;
  proteins[i_pivot] = tempe;
  // sort two sublists in parallel;

  /* The recursive steps in quicksort execution is implemented as separate tasks */
  #pragma omp task
    Qsort(first, (right - 1));
  #pragma omp task
    Qsort((right + 1), last);
}

void UpdateTable(int spl, int ith)
{
#if 0
   if (spl < 0 || ith < 0)
   {
      cout << "\nAKHIR KAAR BEGHAIRTI !!\n";
      exit(0);
   }
#endif

   int iJ = 0;
   for (iJ = 0; iJ < FIRSTROW_CHARS; iJ++)
   {
      bwa_mods[spl].BWTfirstRow.occs[iJ] += _ppbFirstRowModOccurrences[ith][iJ];
      _ppbFirstRowModOccurrences[ith][iJ] = 0;
   }
   for (iJ = 0; iJ < MAX_RANGE; iJ++)
   {
      bwa_mods[spl].BWTSecondRow.occs[iJ] += _ppbSecondRowModOccurrences[ith][iJ];
      _ppbSecondRowModOccurrences[ith][iJ] = 0;
   }
}

void quicksort(int *arr, int low, int high)
{
  int pivot, i, j, temp;
  if(low < high) {
    pivot = low; // select a pivot element
    i = low;
    j = high;
    while(i < j) {
      // increment i till you get a number greater than the pivot element
      while(arr[i] <= arr[pivot] && i <= high)
        i++;
      // decrement j till you get a number less than the pivot element
      while(arr[j] > arr[pivot] && j >= low)
        j--;
      // if i < j swap the elements in locations i and j
      if(i < j) {
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
      }
    }

    // when i >= j it means the j-th position is the correct position
    // of the pivot element, hence swap the pivot element with the
    // element in the j-th position
    temp = arr[j];
    arr[j] = arr[pivot];
    arr[pivot] = temp;
    // Repeat quicksort for the two sub-arrays, one to the left of j
    // and one to the right of j
    quicksort(arr, low, j-1);
    quicksort(arr, j+1, high);
  }
}
