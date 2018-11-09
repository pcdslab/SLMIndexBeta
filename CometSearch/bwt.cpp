// C program to find Burrows Wheeler transform of
// a given text
#include "bwt.h"
#include "SLMIndex.h"

unsigned int size;
Peak *text_ptr;

// Compares the rotations and
// sorts the rotations alphabetically
// Compares the rotations and
// sorts the rotations alphabetically
int cmpfunc(const void *x, const void *y)
{
   int result;
   struct rotation *rx = (struct rotation *) x;
   struct rotation *ry = (struct rotation *) y;
   Peak  *ptra = text_ptr + rx->index;
   Peak  *ptrb = text_ptr + ry->index;

   unsigned int idx1 = rx->index;
   unsigned int idx2 = ry->index;

   /* do something here */

   for (unsigned inc1 = 0, inc2 = 0; inc1 < size && inc2 < size; inc1++, inc2++)
   {
      if (inc1 == size - 1 && inc2 == size - 1)
      {
         result = 0;
         break;
      }

      if (*ptra == *ptrb)
      {
         if (idx1 < (size) && idx2 < (size))
         {
            ptra++;
            ptrb++;
            idx1++;
            idx2++;
         }

         else if (idx1 == (size))
         {
            ptra = text_ptr;
            ptrb++;
            idx2++;
            idx1 = 0;
         }

         else if (idx2 == (size))
         {
            ptrb = text_ptr;
            ptra++;
            idx1++;
            idx2 = 0;
         }

         else
         {
            bool error = true;
            while (error);
         }
      }

      else if (*ptra < *ptrb)
      {
         result = -1;
         break;
      }

      else if (*ptra > *ptrb)
      {
         result = 1;
         break;
      }

      else
      {
         bool error = true;
         while (error);
      }
   }

   return result;
}

// Takes text to be transformed and its length as
// arguments and returns the corresponding suffix array
unsigned int *computeSuffixArray(Peak *input_text, unsigned int len_text)
{
   // Array of structures to store rotations and
   // their indexes
   size = len_text;
   text_ptr = input_text;

   struct rotation *suff = new struct rotation[len_text];

   // Structure is needed to maintain old indexes of
   // rotations after sorting them
   for (unsigned int i = 0; i < len_text; i++)
   {
      suff[i].index = i;
      suff[i].suffix = (input_text + i);
   }

   cout <<"\n     - Sorting for BWT/SA now" << endl;

   // Sorts rotations using comparison function defined above.
   // Can do the merge sort here as well.
   BWTmergeSort(suff, 0, len_text);

   cout <<"     - Success: SA/BWT ready now" << endl << endl;

   // Stores the indexes of sorted rotations
   UINT *suffix_arr = new unsigned int[len_text];

   for (unsigned int i = 0; i < len_text; i++)
   {
      suffix_arr[i] = suff[i].index;
   }

   delete[] suff;
   // Returns the computed suffix array
   return suffix_arr;
}

// Takes suffix array and its size as arguments and returns
// the Burrows - Wheeler Transform of given text
unsigned char *findLastChar(Peak *input_text, unsigned int *suffix_arr, unsigned int n)
{
   // Iterates over the suffix array to find
   // the last char of each cyclic rotation
   Peak *bwt_arr = new unsigned char[n];
   unsigned int i;
   for (i = 0; i < n; i++)
   {
      // Computes the last char which is given by
      // input_text[(suffix_arr[i] + n - 1) % n]
      int j = suffix_arr[i] - 1;
      if (j < 0)
         j = j + n;

      bwt_arr[i] = input_text[j];
   }

   // Returns the computed Burrows - Wheeler Transform
   return bwt_arr;
}

unsigned char *findFirstChar(Peak *input_text, unsigned int *suffix_arr, unsigned int n)
{
   // Iterates over the suffix array to find
   // the last char of each cyclic rotation
   Peak *sorted_arr = new unsigned char[n];
   unsigned int i;
   for (i = 0; i < n; i++)
   {
      sorted_arr[i] = input_text[suffix_arr[i]];
   }

   // Returns the computed Burrows - Wheeler Transform
   return sorted_arr;
}

void BWTmergeSort(struct rotation arr[], int l, int r)
{
   if (l < r)
   {
       // Same as (l+r)/2, but avoids overflow for
       // large l and h
       int m = l+(r-l)/2;

       // Sort first and second halves
       BWTmergeSort(arr, l, m);
       BWTmergeSort(arr, m+1, r);

       BWTmerge(arr, l, m, r);
   }
}

void BWTmerge(struct rotation arr[], int l, int m, int r)
{
   int i, j, k;
   int n1 = m - l + 1;
   int n2 =  r - m;

   /* create temp arrays */
   struct rotation *L = new struct rotation[n1];
   struct rotation *R = new struct rotation[n2];

   /* Copy data to temp arrays L[] and R[] */
   for (i = 0; i < n1; i++)
       L[i] = arr[l + i];
   for (j = 0; j < n2; j++)
       R[j] = arr[m + 1+ j];

   /* Merge the temp arrays back into arr[l..r]*/
   i = 0; // Initial index of first subarray
   j = 0; // Initial index of second subarray
   k = l; // Initial index of merged subarray
   while (i < n1 && j < n2)
   {
       if (cmpfunc((void *)(L+i), (void *)(R+j)) == -1)
       {
           memcpy((void *)(arr+k), (void *)(L+i), sizeof(struct rotation));
           i++;
       }
       else
       {
           memcpy((void *)(arr+k), (void *)(R+j), sizeof(struct rotation));
           j++;
       }
       k++;
   }

   /* Copy the remaining elements of L[], if there
      are any */
   while (i < n1)
   {
       memcpy((void *)(arr+k), (void *)(L+i), sizeof(struct rotation));
       i++;
       k++;
   }

   /* Copy the remaining elements of R[], if there
      are any */
   while (j < n2)
   {
       memcpy((void *)(arr+k), (void *)(R+j), sizeof(struct rotation));
       j++;
       k++;
   }
   delete[] L;
   delete[] R;
}
