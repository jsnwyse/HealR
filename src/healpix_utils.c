#include "healpix_utils.h"

void healpix_utils_create_nested_healpix_file( char *filename, int nside, int nfields, int *status )
{
	
	fitsfile *fp;
	int stat=0, anynul, hdunum, hdutype;
	double nulval;
	int naxis2 = 12 * nside * nside  ;

	char **ttype;
	char **tform;
	char **tunit;
	char *extname = "--";

	int tfields = nfields;

	ttype = malloc( tfields * sizeof(char *) );
	tform = malloc( tfields * sizeof(char *) );
	tunit = malloc( tfields * sizeof(char *) );
	
	ttype[0] = "ttype0";
	ttype[1] = "ttype1";

	tform[0] = "E";
	tform[1] = "E";

	tunit[0] = "tunit0";
	tunit[1] = "tunit1";

	char *pixname = "PIXTYPE";
	char *orderingname = "ORDERING";
	char *pixtype = "HEALPIX";
	char *ordering = "NESTED";
	char *comment = NULL;
	char *nsidename = "NSIDE";
	int *firstpix;
	int *lastpix;
	int *Nside;
	int fpix = 0;
	int lpix = naxis2-1;
	int nsi = nside;
	
	Nside = &nsi;
	firstpix = &fpix;
	lastpix = &lpix;

	int k ;
	
	if( fits_create_file( &fp, filename, &stat) )
	{
		*status = 105; return;
	}
	
		
	if( fits_create_tbl( fp, BINARY_TBL, naxis2, tfields, ttype, tform, tunit, extname, &stat ) )
	{
		*status = 106; return;
	}
		
	fits_write_key( fp, TSTRING, pixname, pixtype, comment, &stat );
	fits_write_key( fp, TSTRING, orderingname, ordering, comment, &stat );
	fits_write_key( fp, TINT, nsidename, Nside, comment, &stat );
	fits_write_key( fp, TINT, "FIRSTPIX", firstpix, comment, &stat );
	fits_write_key( fp, TINT, "LASTPIX", lastpix, comment, &stat );
		
	fits_close_file( fp, &stat );
		
	return;

}

