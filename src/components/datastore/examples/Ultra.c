#include "Vista.h"
#include "View.h"
#include "VHashTable.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

/**************************************************************************
 * Subroutine:  DumpUltra
 * Purpose   :  Dump tagged problem data to an ultra file 
 *************************************************************************/

void DumpUltra(View_t prob)
{
   int i ;
   FILE *fp ;
   char fname[100] ;
   char *tail ;
   VHashTraverse_t content ;

   View_t elem = prob.view("elem") ;

   strcpy(fname, prob.name()) ;

   /* Skip past the junk */
   for (tail=fname; isalpha(*tail); ++tail) ;

   sprintf(tail, "_%04d\0", prob.paramInt("cycle")) ;

   if ((fp = fopen(fname, "w")) == NULL)
   {
      printf("Could not open file %s. Aborting.\n", fname) ;
      exit (-1) ;
   }

   fprintf(fp, "# Ultra Plot File\n") ;
   fprintf(fp, "# Problem: %s\n", prob.name()) ;

   content = VHashTraverseInit() ;

   while (VHashTraverseGetNext(prob.params().itemTable, content) == true)
   {
      Param *currParam = (Param *) VHashTraverseEntry(content) ;
      switch(currParam->type())
      {
         case OTC_INT :
         {
            fprintf(fp, "# %s = %d\n", currParam->name(), currParam->Int()) ;
         }
         break ;

         case OTC_REAL :
         {
            fprintf(fp, "# %s = %f\n", currParam->name(), currParam->Real()) ;
         }
         break ;

         default:
            /* skip and warn */
            printf("Add another case in DumpUltra Scalar\n") ;
      }
   }

   VHashTraverseFree(content) ;

   content = VHashTraverseInit() ;

   while (VHashTraverseGetNext(elem.fields().itemTable, content) == true)
   {
      int plotlen = elem.length() ;
      Field *currField = (Field *) VHashTraverseEntry(content) ;
      switch(currField->type())
      {
         case OTC_INT :
         {
            int *data = currField->Int() ;
            fprintf(fp, "# %s\n", currField->name()) ;
            for (i=0; i<plotlen; ++i)
               fprintf(fp, "%f %f\n", (double) i, (double) data[i]) ;
            fprintf(fp, "\n") ;
         }
         break ;

         case OTC_REAL :
         {
            double *data = currField->Real() ;
            fprintf(fp, "# %s\n", currField->name()) ;
            for (i=0; i<plotlen; ++i)
               fprintf(fp, "%f %f\n", (double) i, data[i]) ;
            fprintf(fp, "\n") ;
         }
         break ;

         default:
            /* skip and warn */
            printf("Add another case in DumpUltra Scalar\n") ;
      }
   }

   VHashTraverseFree(content) ;

   fclose(fp) ;

   return ;
}

