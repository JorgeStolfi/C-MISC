#include <stdio.h>
#include "typedefs.h"
#include "constants.h"
#include "funcdefs.h"
#include "maindecl.h"

void main(void);

void main(void)
{
	int	line_y;
	int	pixel_x;
	t_3d	scrnx;
	t_3d	scrny;
	t_3d	firstray;
	t_3d	ray;
	t_color	color;
	double	dis;
	double	line[SCREENWIDTH][3];

	setup  ();
        viewingdefault();
    
	printf ("\nFile name = %s", outfilename);
	printf ("\nNObject = %d", nobject);
 	printf ("\nNLight = %d", nlight);
	printf ("\nNSurface = %d", nsurface);
	printf ("\nScreen = %d,%d", sizex, sizey);
	printf ("\nFov = %f,%f", vfov, hfov);
	printf ("\nMaxlevel = %d", maxlevel);
	printf ("\nEyep = %f,%f,%f", eyep.x, eyep.y, eyep.z);
	printf ("\nLookp = %f,%f,%f", lookp.x, lookp.y, lookp.z);
	printf ("\nUp = %f,%f,%f", up.x, up.y, up.z);
    
	viewing  (&scrnx, &scrny, &firstray);
	startpic (outfilename, sizey, sizex);

    	for (line_y = 0; line_y < sizey; line_y++)
	{
		for(pixel_x = 0; pixel_x < sizex; pixel_x++)
		{
	    		ray.x = firstray.x + pixel_x * scrnx.x - line_y * scrny.x;
	    		ray.y = firstray.y + pixel_x * scrnx.y - line_y * scrny.y;
	    		ray.z = firstray.z + pixel_x * scrnx.z - line_y * scrny.z;
	    
	    		normalize( &ray );

			/* actual ray trace */
	    		dis = intersect( object, nobject, -1, &eyep, &ray, &hit, &nrm);

			/* find color at point of intersection */
			shade( &hit, &ray, &nrm, &object[objhit], color );
	    
	    		if(dis > 0)	/* ray intersected object */
	    		{ 
				line[pixel_x][0] = color.r;
				line[pixel_x][1] = color.g;
				line[pixel_x][2] = color.b;
	    		} /* if */
	    		else
	    		{  		/* use background color */
				line[pixel_x][0] = background.r;
				line[pixel_x][1] = background.g;
				line[pixel_x][2] = background.b;
	    		} /* else */
		} /* for */

		linepic( line );  /* Output line of pixels */
		if (line_y % 10 == 0)
		{
	    		printf( "\nDone line %d", line_y );
	    		fflush( stdout );
		} /* if */
	} /* for */
	endpic();   /* done with picture */
} /* main */
/* */
/* */

