/*I/O for FITS files with Healpix nested ordering
	
	Author:	Jason Wyse,
			School of Computer Science and Statistics,
			Lloyd Institute,
			Trinity College,
			Dublin 2,
			Ireland.
			mailto: wyseja@tcd.ie
*/
	
#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "healpix_utils.h"
#include "chealpix.h"
#include "defs.h"
	
static void astro_extract_nested_healpix_patch( char **filename, int *nside, int *patch, int *nrow, int *rows, int *ncol, int *cols, double *x, int *status, int *hdu_no, int *col_no );
	
static void astro_create_nested_healpix_file( char **filename, int *nside, int *nfields, int *status );
	
static void astro_write_nested_healpix_patch( char **filename, double *x, int *nside, int *nwrite, int *rows, int *cols, int *colnum, int *status, int  *patch );
	
static void astro_write_vector_to_healpix_file_0( char **filename, double *x, int *nside, int *colnum, int *status );	
	
	
static const R_CMethodDef cMethods[] = {
	  {"astro_extract_nested_healpix_patch", (DL_FUNC) &astro_extract_nested_healpix_patch, 11},
	  {"astro_create_nested_healpix_file", (DL_FUNC) &astro_create_nested_healpix_file, 4},
	  {"astro_write_nested_healpix_patch", (DL_FUNC) &astro_write_nested_healpix_patch, 9},
	  {"astro_write_vector_to_healpix_file_0", (DL_FUNC) &astro_write_vector_to_healpix_file_0, 5},
	  NULL
};

void R_init_HealR( DllInfo *info )
{
  R_registerRoutines( info, cMethods, NULL, NULL, NULL );
  R_useDynamicSymbols( info, FALSE );
}


static void astro_extract_nested_healpix_patch( char **filename, int *nside, int *patch, int *nrow, int *rows, int *ncol, int *cols, double *x, int *status, int *hdu_no, int *col_no )
{

	//x will be a vector of size length(rows)*length(cols) to be stacked column major

	int  	i, j, k, xi, yi, z, stat = 0, ord, 
			hdu_num = *hdu_no, hdu_type, any_null, col_num  = *col_no ;
			
	double *y, nulval = FITS_NULL_VAL ;
	
	char *pixtype, *comment ;
	
	fitsfile *fp;
	
	if( fits_open_file( &fp, *filename, READONLY, &stat ) ) { status[0] = 104; return; } 
	
	stat = 0 ;
	
	LONGLONG n = (LONGLONG) (*nside) * (*nside), elem = 1, row = ( *patch ) * n + 1; 
	
	y = calloc( (*nside)*(*nside) , sizeof(double) ) ;
	
	if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &stat ) ) { free(y); status[0] = 301; return; }
	
	stat = 0 ;
	
	if( fits_read_col( fp, TDOUBLE, col_num, row, elem, n , &nulval, y, &any_null, &stat ) ) { free(y); status[0] = 302; return;} ;
		
	//now pull out the relevant entries and place in to x (column major)
	
	int cc = 0, rc ;
	
	for( j=0; j<*ncol; j++ )
	{
		rc = 0;
		for( i=0; i<*nrow; i++ )
		{
			xi = rows[i]; yi = cols[j];
			
			z = xyf2nest2( *nside, xi, yi , 0);
			
			x[ (*nrow) * cc + rc ] = y[z] ; 
			
			rc++;
		}
		
		cc++;
	}
	
	fits_close_file( fp, &stat );
		
	free(y);
	
	status[0] = 0;
	
	return;

}

static void astro_create_nested_healpix_file( char **filename, int *nside, int *nfields, int *status )
{
	
	healpix_utils_create_nested_healpix_file( *filename, *nside, *nfields, status ) ;

}

static void astro_write_nested_healpix_patch( char **filename, double *x, int *nside, int *nwrite, int *rows, int *cols, int *colnum, int *status, int *patch )
{
	
	fitsfile *fp;
	int i, j;
	int a, k, l, stat = 0, any_null, hdunum=2, hdutype, count,
	   elem=1, nelem = ( *nside ) * ( *nside ), col_num = *colnum,
	   row = ( *patch ) * nelem + 1 , longnull=0;
	
	double *y = calloc( nelem, sizeof(double) ), nulval = FITS_NULL_VAL;
	
	if( fits_open_file( &fp, *filename , READWRITE, &stat ) ) { free(y); status[0] = 104; return; }
	if( fits_movabs_hdu( fp, hdunum, &hdutype, &stat ) ) { free(y); status[0] = 301; return; }
	
	if( fits_read_col( fp, TDOUBLE, col_num, row, elem, nelem , &nulval, y, &any_null, &stat ) ) { free(y); status[0] = 302; return;} ;
      
  // now overwrite the appropriate values in y
  for( l=0; l<*nwrite; l++ )
  {
    i = rows[l];
    j = cols[l];
    k = xyf2nest2( *nside, i, j, 0);
    y[k] = x[l];
  }

	if( fits_write_col( fp, TDOUBLE, col_num, row, elem, nelem, (void *)y, &stat ) ) { status[0] = 106; return;} ;
	
	fits_close_file( fp, &stat );
	
	free(y);
	
	*status = stat;
	
	return;
}

static void astro_write_vector_to_healpix_file_0( char **filename, double *x, int *nside, int *colnum, int *status )
{
	
	fitsfile *fp;
	int i, j;
	int a, k, l, stat = 0, anynul, hdunum=2, hdutype, count,
	   felem=1, nelem = 12 * ( *nside ) * ( *nside ), column = *colnum,
	   frow =1 , longnull=0;

	
	if( fits_open_file( &fp, *filename , READWRITE, &stat ) )
		{ *status = 104; return; }
	
	if( fits_movabs_hdu( fp, hdunum, &hdutype, &stat ) ) 
		{ *status = 301; return; }
	
	if( fits_write_col( fp, TDOUBLE, column, frow, felem, nelem, (void  *)x, &stat ) )
	  { *status = 302; return; }
	
	fits_close_file( fp, &stat );
	
	*status = stat;
	
	return;

}



