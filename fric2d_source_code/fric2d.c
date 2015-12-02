#define PROGRAM fric2d.c			/* Info for output/resume files		*/
#define VERSION 3.2.7				/*		"			*/
#define REVISE_DATE April 2014			/*		"			*/


/******************************* program: fric2d *****************************
*
* This program calculates two-dimensional fracture propagation paths in an
* elastic body.  Stress and displacement calculations are based on a
* boundary element program employing displacement discontinuities (TWODD)
* found in:
*
*						Crouch and Starfield, 1983
*				Boundary Element Methods in Solid Mechanics
*		With Applications in Rock Mechanics and Geological Engineering
*					George Allen & Unwin, London, 322 p.
*
*	
*****************************************************************************/

/**************************** REVISION HISTORY *******************************
*
* If you change this code in ANY way, you MUST describe and initial those
* changes below.
*
* DATE      CHANGE                                                      INTLS
* --------  ----------------------------------------------------------  -----
* Jan 93	Version 1.0 (NON-ANSI) released.			ALT
* Dec 93	Converted to ANSI					MLC
* Jan 94	Added elastic joint elements with Mohr Coulomb slip	MLC
* Jun 95	Added development of opening mode joints from faults	MLC
* Aug 99	Stop convergence after 100 iterations			MLC	
* Sep 02        Improved gravitational stress state with kratio         MLC
* Dec 11	Fixed unintialized variable problem in faultSlip	MLC
* Feb 13	Implemented more robust tension and slip tests		MLC 
*		Now faults will open when normal stresses > 0 
*		version 3.0
* Jul 13	Changed the Matrix Inversion to Numerical Recipes	MLC 
*		LU Decomposition. Added calculation of condition num 
*		Added dynamic friction with critical slip 
* Aug 13	3.1: Added inherent shear strength which drops to cohesion
*		upon slip, friction evolution				MLC
* Oct 13	3.1.2 Fixed issue with KI and Kii reporting on fault tips MLC
* Nov 13	3.2 Checked for oscillating convergence			MLC
* 		added fault tensile strength so faults open when normal
*		stresses exceed Tensile strength 			MLC
* Mar 14        3.2.4 Implemented reassignment of tensile strength	MLC 
*		after tensile failure 
* Apr 14	3.2.6 fixed reporting of opening failure		MLC
* Apr 14	3.2.7 Shear failure if the normal stress is greater 
* 		than shear strength divided by friction			MLC
*				 
* Known compiler warnings but code should still run ok: 
*getoption.c: In function ‘add_arg’:
*
*
* ALT = Andrew L. Thomas
* MLC = Michele Lynn Cooke
*
*****************************************************************************/
 
 
/******************************* includes/defines ***************************/
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "getoption.h"
#include "nrutil.h"
#if !defined(__APPLE__)				/*on macs malloc is superceded by stdlib */
   #include <malloc.h>
#endif
#define COMMENT_CHAR '*'			/* char starting infile comment lines	*/
#define DB "**db**"				/* string starting DPRINTF stmts	*/
#define DPRINTF if (debug) printf		/* debugging printf macro		*/
#define FALSE 0					/* false flag for boolean expressions	*/
#define ISECT -1				/* tip intersected flag			*/
#define ISECTTXT "intersected"			/* tip intersected flag in input file	*/
#define MAXELTS 5000				/* max # of BE's (boundaries & fracs)	*/
#define MAXFILE 50				/* max length of file names		*/
#define MAXFNAME 50				/* max length of fracture/fault name		*/
#define MAXLINE 120				/* max length of an input line		*/
#define MAXWORDS 20				/* max # of words on an input line	*/
#define MAXTITLE 80				/* max length of problem titles		*/
#define PROPYESTXT "yes"			/* allow tip prop flag in input file	*/
#define PROPYES 1				/* allow tip propagation flag		*/
#define PROPNO  0				/* disallow tip propagation flag	*/
#define PROPNOTXT "no"				/* disallow tip prop flag in input file	*/
#define TRUE 1					/* true flag for boolean expressions	*/
#define TENSILE 2				/* Flag for tension on friction elements*/
#define VFPRINTF if ((verbose) &&		/* verbose option fprintf macro		*/ \
	(ofp != NULL)) fprintf			/*     (continued)						*/

#define MAX_PDFITERS 30
 
						/* 	symmetry condition flags	*/
#define NOSYM 1					/* 	no symmetry			*/
#define XLSYM 2					/* 	symmetry about a line x = xsym	*/
#define YLSYM 3					/* 	symmetry about a line y = ysym	*/
#define XYLSYM 4				/* 	symmetry about two lines x & y	*/
#define XYASYM 5				/* 	anti-symmetry about a point x,y	*/
 
/***************************** structured types *****************************/
 
/*-----------------------------------------------------------------------
FRACTURE POINTERS: Because the boundary elements that make up a fracture
do not occupy contiguous positions in the BE arrays, a linked list of
pointers is necessary to keep track of which elements belong to which
fracture.  This makes it possible to generate an "resume" file at the end
of program execution which can used as an input file for continued
execution at some later time.  Below, a "segment" is a contiguous block
of elements in the BE arrays that belong to a single fracture.  A segment
is demarked by its starting and ending elements, as all elements in
between belong to the segment.
------------------------------------------------------------------------*/
 
struct frac {				/* linked list of ptrs to fracture	*/
					/* 	locations in arrays			*/
	char name[MAXFNAME];		/*		fracture name					*/
	int tip_prop[2];		/*		tip propagation flag array	*/
	struct frac_seg *first_seg;	/*		ptr to first segment		*/
	struct frac_seg *last_seg;	/*		ptr to last segment		*/
	struct frac *next;		/*		ptr to next fracture		*/
};
 
struct frac_seg {			/* linked list of ptrs to frac segments	*/
	int start_elt;			/*		array index of starting element	*/
	int end_elt;			/*		array index of ending element	*/
	struct frac_seg *next;		/*		pointer to next segment		*/
};
 
struct oline {				/* linked list of observation lines	*/
	double xbeg, ybeg;		/*		coords of starting point	*/
	double xend, yend;		/*		   "   "  ending point		*/
	int numpb;			/*		no. of pts. bet. start and end	*/
	struct oline *next;		/*		ptr to next obs. line		*/
};

struct fault{				/* linked list of ptrs to fault	*/
					/* 	locations in arrays			*/
	char name[MAXFNAME];		/*		fault name			*/
	int tails;			/*		flag for tail crack growth	*/
	int end_prop[2];		/*		tip propagation flag array	*/
	struct fault_seg *first_seg;	/*		ptr to first segment		*/
	struct fault_seg *last_seg;	/*		ptr to last segment		*/
	struct fault *next;		/*		ptr to next fracture		*/
};
 
struct fault_seg {			/* linked list of ptrs to frac segments	*/
	int start_elt;			/*		array index of starting element	*/
	int end_elt;			/*		array index of ending element	*/
	struct fault_seg *next;		/*		pointer to next segment		*/
};
	
/***************************** function declarations ************************/
/* 					(See function code for detailed comments) 				*/
 
void 	add_BEs(int num, double xbeg, double xend, double ybeg, double yend,
		int kode,double bvsh, double bvno, double bvshM, double bvnoM);
void 	addTail(FILE *ofp,int tipelt, int end, double pangle);
void 	adjustCoeff();
void	append_resumefile(FILE *);
void 	calc_coeff(double xi, double yi, double xj, double yj, double aj,
		   double cosbj, double sinbj);
void 	calc_BE_BCs(int i, double * usneg, double * unneg,
		    double * uspos, double * unpos, double * sigs, double * sign);
void	calc_inf_coeffs();
void	coeff(double x, double y, double cx, double cy, double a,
	      double cosb, double sinb, int msym);
int	disps_stresses_at_pt(double xp, double yp, double *ux, double *uy,
			     double *sigxx, double *sigyy, double *sigxy, int distcheck);
double	dist_to_closest_BE(double x, double y, int dont_check_BE, int *closest_BE);
void	findSlipped(FILE *ofp, int *prop_possible2);
double  fail(FILE *ofp, int tipelt, int top, double *pangle);
void  	faultSlip(FILE *ofp, int *contin);
void	fracs(FILE *ofp, FILE *gfp, int *niters, int pdfpos, int step,int nsteps);
void 	increment_BEs(int step, int nsteps);
extern 	int getwords(FILE *infp, char *line, int maxline, char *word[], int maxwords);
void	grow_fracs(FILE *ofp, int *prop_possible);
extern  void ludcmp(double **a, int n, int *indx, double *d);
extern  void lubksb(double **a, int n, int *indx, double b[]);
int	open_files(char *infile, FILE **ifp, char *outfile,
		   FILE **ofp, char *resumefile, FILE **rfp,
		   char *graphicsfile, FILE **gfp);
void	print_data(FILE *ofp, int iter, int bdata, int fracData, int faultData, int odata, int step);
void	print_graphics_start(FILE *gfp);
void	print_graphics_frame(FILE *gfp, int step, int increment);
void	print_graphics_end(FILE *gfp);
void	print_oline_data(FILE *ofp, int linenum, double xbeg, double ybeg,
			 double xend, double yend, int numpb, int titles);
void	print_BE_data(FILE *ofp, int first_elt, int last_elt, int titles);
void	print_BE_disps_stresses(FILE *ofp, int first_elt, int last_elt, int titles);
void 	print_fault_data(FILE *ofp, int first_elt, int last_elt, int titles);
void	print_prob_info(FILE *ofp);
double  propagation(FILE *ofp, int tipelt, double *pangle);
void 	quicksolve(int n, float *condnum);
int	readFault(FILE *ifp, char *line, char **word, int numwords);
int	readFracture(FILE *ifp, char *line, char **word, int nuwords);
void	reassign_BEs();
int	startup(FILE *ifp, FILE *rfp, int *niters, int *nsteps);

/****************************** global variables ****************************/
 
/* constants */
double	con;
double	cons;
int		dstrobe = 0;			/* write data every dstrobe iterations	*/
int		pdfiter[MAX_PDFITERS];		/* print data from iteration __ array	*/
int		bdata = TRUE;			/* write boundary BE data flag		*/
int		fracData = TRUE;			/* write fracture BE data flag		*/
int		faultData = TRUE;		/* write fault BE flag			*/
int		odata = TRUE;			/* write observation line data flag	*/
double	e;					/* young's modulus			*/
double	pr;					/* poisson's ratio			*/
double	k1c;					/* mode-I fracture toughness		*/
double tensileStrength;				/* tensile strength of rock 		*/
int  		ksym;				/* symmetry cond flag (see define stmts)*/
int		nblines = 0;			/* number of straight boundary lines	*/
int	   	nbBEs = 0;			/* number of boundary BEs		*/
int		nfracs = 0;			/* number of fractures			*/
int		numFaults = 0;			/* number of faults			*/
int		numFaultElements = 0;		/* number of fault elements		*/
int		nolines;			/* number of observation lines		*/
double	pi;					/* 3.14159.....				*/
double	pr1;					/* (1-2v) v- poissons ratio  		*/
double	pr2;					/* 2(1-v) v- poissons ratio  		*/
double	xsym;					/* x-coord of symmetry point or line	*/
double	ysym;					/* y-coord "     "       "   "    "	*/
int		gravity=0;			/* Flag to add gravitational stress	*/
double	kratio;				/* horizontal:vertical gravity stress  */
double	density;				/* Density of rock			*/
double  toler;					/* tolerance % for stopping interations */
int             tolCounter;                     /*a counter that stops interations      */

/* Used for iterative determination of influence coeffs */
double	sxxn;
double	sxxs;
double	sxyn;
double	sxys;
double	syyn;
double	syys;
double	uxn;
double	uxs;
double	uyn;
double	uys;
 
/* Boundary element arrays */
double	a[MAXELTS];				/* half-length of BE			*/
double	b[MAXELTS*2];				/* stress or displ BC magnitude on BE	*/
double	bvsS[MAXELTS];				/* STATIC shear stress or displ. BC mag.*/
double	bvnS[MAXELTS];				/* normal  "    "   "     "     "	*/
double	bvsM[MAXELTS];				/* MONO. shear stress or displ. BC mag.	*/
double	bvnM[MAXELTS];				/* normal  "    "   "     "     "	*/
double 	stepbvs[MAXELTS];			/* shear stress or displ. BC magnitude	*/
double	stepbvn[MAXELTS];			/* normal  "    "    "     "    "	*/
double  stiffS[MAXELTS];			/* Shear fault stiffness		*/
double  stiffN[MAXELTS];			/* Normal fault stiffness		*/
double  faultTen[MAXELTS];			/* Tensile strength of fault elements   */
double  strength[MAXELTS];			/* Inherent shear strength of material  */	
double	cohesion[MAXELTS];			/* cohesion of fault element		*/
double  friction[MAXELTS];                      /* coefficient of friction of fault el  */
double 	fric_dyn[MAXELTS];                      /* dynamic friction                     */
double 	fric_stat[MAXELTS];                     /* static friction                      */
double  crit_slip[MAXELTS];			/* distance over which friction drops   */
int		slipped[MAXELTS];		/* slip on fault flag			*/
double   *c[MAXELTS*2];				/* influence coeff matrix		*/
double	cosbet[MAXELTS];			/* cos(angle bet elt and x1-axis)	*/
double	d[MAXELTS*2];				/* pseudo stress or displ on BE		*/
int		kod[MAXELTS];			/* type of BC's (stress/displ) on BE	*/
double	sinbet[MAXELTS];			/* sin(angle bet elt and x1-axis)	*/
double	xm[MAXELTS];				/* x-coord of middle of BE		*/
double	ym[MAXELTS];				/* y-coord of middle of BE		*/
int		permSlip[MAXELTS];		/* whether element ever slipped 	*/
int		permOpen[MAXELTS];		/* whether an element ever open		*/
int		TCposs[MAXELTS][2];		/* whether element can have tail cracks	*/
 
/* Misc. */
int		debug = FALSE;			/* debugging flag			*/
int		verbose = FALSE;		/* verbose output option flag		*/
char	*global_argument;			/* global arg for getopt function	*/
int		numbe;				/* number of BE's			*/
char	title1[MAXTITLE];			/* char array for 1st line of title	*/
char	title2[MAXTITLE];			/* char array for 2nd line of title	*/
struct frac	*first_frac = NULL;		/* ptr to 1st frac in linked list	*/
struct frac	*last_frac = NULL;		/* ptr to last frac in linked list	*/
struct oline 	*first_oline = NULL;		/* ptr to 1st obs line in linked list	*/
struct fault	*first_fault = NULL;		/* ptr to 1st frac in linked list	*/
struct fault	*last_fault = NULL;		/* ptr to last frac in linked list	*/
	
/*************************** function: main *****************************/
int
main(int argc, char **argv)
{
	int		continue_ok = FALSE;	/* continue after iters completed flag	*/
	int		exit;			/* exit from loop flag			*/
	int		i;			/* loop control variable		*/
	FILE   *ifp = NULL;			/* input file pointer			*/
	char	infile[MAXFILE];		/* input file name			*/
	char	line[MAXLINE];			/* line of text				*/
	FILE   *ofp = NULL;			/* output file pointer			*/
	int		optarg;			/* cmd line option argument		*/
	char	outfile[MAXFILE];		/* ouput file name			*/
	FILE   *gfp = NULL;			/* graphics file pointer		*/
	char	graphicsfile[MAXFILE];		/* graphics file name			*/
	FILE   *rfp = NULL;			/* new input file pointer		*/
	char	resumefile[MAXFILE];		/* new input file name			*/
	int		nsteps;			/* # of loading steps to be completed	*/
	int		step;			/* current step being processes		*/
	int		niters;			/* # of iterations to be completed	*/
	int		additers;		/* additional # of iters to be completed*/
	static int	pdfpos;			/* vaule of next iteration to print	*/
	char   *word[MAXWORDS];			/* array of ptrs to wds. on line[]	*/
	
 	/*--------------------------
 	Calloc space for BE matrix 
 	---------------------------*/
 	for (i = 0; i < MAXELTS*2; i++)
 		c[i] = (double *) calloc((size_t) MAXELTS*2, sizeof(double));
 	
	/*----------------------------------------------------
	Set 1st char of filenames and problem titles to NULL
	----------------------------------------------------*/
	infile[0]     = '\0';
	outfile[0]    = '\0';
	graphicsfile[0]    = '\0';
	resumefile[0] = '\0';
	title1[0]     = '\0';
	title2[0]     = '\0';	
	pdfiter[0]    =0;
	
	/*------------------------------------------------
	Parse the command line arguments using getopt()
	------------------------------------------------*/
	exit = FALSE;
	while ((optarg = getoption("vdci:o:g:r:",argc,argv))
		!= NO_MORE_ARGS && !exit) {
		switch (optarg) {
			case NO_SUCH_ARG:
			case FILE_ARG:
				exit = TRUE;
				break;
			/* -v turns on verbose mode */
			case 'v':
				verbose = TRUE;
				break;
			/* -c prompts for addit. iters after requested iters completed */
			case 'c':
				continue_ok = TRUE;
				break;
			/* -d turn on debugging mode */
			case 'd':
				debug = TRUE;
				DPRINTF("%s IN MAIN: debugging turned on\n",DB);
				break;
			/* -i <filename> names the input file */
			case 'i':
				strncpy(infile,global_argument,MAXFILE);
				break;
			/* -o <filename> names the output file */
			case 'o':
				strncpy(outfile,global_argument,MAXFILE);
				break;
			/* -g <filename> names the graphics file */
			case 'g':
				strncpy(graphicsfile,global_argument,MAXFILE);
				break;
			/* -r <filename> names the resume file */
			case 'r':
				strncpy(resumefile,global_argument,MAXFILE);
				break;
		} /*switch*/
	} /*while*/
				
	DPRINTF("MADE IT THROUGH THE COMMAND LINE PARSING\n");			
 
	/* If invalid cmd line arg. given or no input file specified, quit */
	if (exit || infile[0] == '\0') {
		fprintf(stderr,"\n");
		fprintf(stderr,"Usage: fric2d -i <infile> [-o <outfile>] [-g <graphicsfile>]");
		fprintf(stderr,"[-r <resumefile>] [-c] [-d] [-v]\n");
		fprintf(stderr,"       (arguments may occur in any order)\n\n");
		return(-1);
	} /*if*/
 
 	DPRINTF("TRYING TO OPEN FILES\n");
 	
	/*----------------------------------------------
	Open the input and output files. Exit on error.
	-----------------------------------------------*/
	if (open_files(infile, &ifp, outfile, &ofp, resumefile, &rfp,
		     graphicsfile, &gfp) == -1)
		return(-1);
 
 	DPRINTF("CALLING STARTUP()\n");
 	
	/*-----------------------------------------------------------
	Start up the problem by reading the input file and setting up
	the initial BE arrays
	------------------------------------------------------------*/
	startup(ifp, rfp, &niters, &nsteps);
		
	DPRINTF("PRINTING PROBLEM INFO\n");
		
	/*-----------------------------------
	Print out the problem information
	------------------------------------*/
	if (ofp != NULL)
		print_prob_info(ofp);
 
	/*-----------------------------------
	Print starting stuff to gfp:
	------------------------------------*/
	if (gfp != NULL)
	  print_graphics_start(gfp);
 
	/*--------------
	Do the problem
	---------------*/
	for (step = 1; step <= nsteps; step++) 
	{

		fprintf(ofp,"\n\n");
		fprintf(ofp,"===============================\n");
		fprintf(ofp,"DATA FROM LOADING STEP #%d OF %d\n",step,nsteps);
		fprintf(ofp,"===============================\n\n");	

/*temp print statement*/
		printf("\nThe loading step is now %d out of %d.", step, nsteps);

		increment_BEs(step, nsteps);

		pdfpos = 0;		/* zero out the value of next iteration to print */
		
		for (;;)
		{
			DPRINTF("ENTERING FRACS()\n");
		
			fracs(ofp, gfp, &niters, pdfpos, step, nsteps);
		
			DPRINTF("OUT OF FRACS()\n");
 
			if (continue_ok) {
				printf("\nA total of %d iterations have been completed.\n",
				niters);
				printf("How many additional iterations (<return> to quit)? ");

				if (0 == getwords(stdin,line,MAXLINE,word,MAXWORDS))
					break;
				else
					additers = atoi(word[0]);
					niters += additers;
			} /*if*/
			else
				break;
		} /*for*/
	} /*for*/
 
 
	/*-----------------------------------------------------
	If resume file specified, append fracture info to file
	------------------------------------------------------*/
	if (rfp != NULL)
		append_resumefile(rfp);

	if (gfp != NULL)
	  print_graphics_end(gfp);
}
 
/*************************** function: fracs ********************************
*
*****************************************************************************/
void fracs(FILE *ofp, FILE *gfp, int *niters, int pdfpos, int step, int nsteps)
{
	int		i, contin, count;
	static	int oldniters = 1;
	int		prop_possible, prop_possible2;
	float	condnum;

	for (i=oldniters; i<=*niters; i++) {
		VFPRINTF(ofp,"\n\nStarting fracture growth increment #%d of %d\n",i,*niters);
		VFPRINTF(ofp,"----------------------------------\n");
/*temp*/	printf("\n   Starting fracture growth increment #%d of %d",i,*niters);
/*temp*/	printf("\n      Iteration: ");
		contin= TRUE;
		count = 0; 

 	while(contin == TRUE && count <= tolCounter){
		 
		++count;
/*temp*/	printf("%d ", count);	
/*temp*/	fflush(stdout);

		VFPRINTF(ofp,"\nStarting iteration #%d of frictional element convergence\n", count);

		/*-----------------------------------------------------
		Set-up the boundary condition arrays according to the step
		of loading
		-------------------------------------------------------*/
		reassign_BEs();

		/*------------------------------------------------------------
		Compute influence coefficients and set up system of algebraic
		equations
		-------------------------------------------------------------*/
		VFPRINTF(ofp,"  Calculating influence coefficients\n");
		calc_inf_coeffs();

		/*----------------------------------------
		Alter coefficients for fault elements
		---------------------------------------*/
		VFPRINTF(ofp,"  Altering coefficients for fault elements\n");
		adjustCoeff();
fflush(ofp);
		/*----------------------------------
		Solve system of algebraic equations
		-----------------------------------*/
		VFPRINTF(ofp,"  Solving system of linear equations\n");
		quicksolve(2*numbe, &condnum);
fflush(ofp); 
		/*----------------------------------------------
		Check fault elements for slip
		Change boundary stress if slipped
		----------------------------------------------*/	
		faultSlip(ofp, &contin);
	}/*while*/
	/*Convergence completed*/
	fprintf(ofp, "\nCONDITION NUMBER = %8.2e AFTER %d FRICTIONAL SLIP ITERATIONS:",condnum, count);
	fprintf(ofp, "\n=================================================================================");
	fprintf(ofp, "\n  If %d = %d convergence was incomplete and some interface elements", count, tolCounter);
	fprintf(ofp, "\n  do not match the solution. This may be resolved with a lower tolerance or");
	fprintf(ofp, "\n  higher tolCounter.  If some faults are slipping a lot (>0.5 the element length)");
	fprintf(ofp, "\n  the locally high strains will prevent convergence.  In this case try");
	fprintf(ofp, "\n  smaller boudary loads/displacements or simpler fault geometry.\n\n");
	printf("\n");
		/*----------------------------------------------------------
 		Print data after requested iterations
 		------------------------------------------------------------*/
 		if ((ofp != NULL) && (i == pdfiter[pdfpos]) ) {
 			print_data(ofp, i, bdata, fracData, faultData, odata, step);
			pdfpos++ ;
		} 		
		/*----------------------------------------------
		Check fracture tips for growth
		Add new elements at fracture tips if necessary
		----------------------------------------------*/
		grow_fracs(ofp, &prop_possible);

		/*-------------------------------------------------------
		Find the slipped elements and evaluate fracture growth 
		---------------------------------------------------------*/
		findSlipped(ofp, &prop_possible2);
	
		/* Output stuff to graphics file */
		if (gfp != NULL)
		  print_graphics_frame(gfp, step, i);

		/* Stop the iterations if all fractures are propagated or intersected*/
		if ((!prop_possible)&& (!prop_possible2) && (ofp != NULL)) {
			fprintf(ofp,"\n\n*****************************************\n");
			fprintf(ofp,"Increment #%d\n",i);
			fprintf(ofp,"ALL FRACTURE TIPS ARE EITHER INTERSECTED\n");
			fprintf(ofp,"OR ARE NOT ALLOWED TO PROPAGATE.\n");
			fprintf(ofp,"*****************************************\n");
			break;
		}/*if*/

	} /*for*/
	
}
 
/**************************** function: grow_fracs **************************
*
*****************************************************************************/
void grow_fracs(FILE *ofp, int *prop_possible)
{
	struct	frac		*current_frac;
	struct	frac_seg	*temp_seg;
 
	int		first_tip;
	double	dist;
	double	pangle;	
	int		i, j;
	int		closest_BE;
	int		tipelt;
	double	xbbeg, xbend;
	double	xbeg, xend;
	double	ybeg, yend;
	double	ybbeg, ybend;
	double	cosb, sinb;
	int		kode;
	double	bvshS, bvnoS, bvshM, bvnoM;
	double	failure;
 
	*prop_possible = FALSE;
	current_frac = first_frac;
	for (i=1; i<=nfracs; i++) {
 
		VFPRINTF(ofp,"  Checking fracture #%d of %d\n",i,nfracs);
 
		for (j=0; j<2; j++) {
 
			VFPRINTF(ofp,"\n    Checking fracture tip #%d\n",j+1);
 
			if (current_frac->tip_prop[j] == PROPNO) {
				VFPRINTF(ofp,"      Tip not allowed to propagate\n");
			}
		    	else if (current_frac->tip_prop[j] == ISECT) {
				VFPRINTF(ofp,"      Tip intersected\n");
			}
			else {
 
				if (first_tip == (j==0))
					tipelt = current_frac->first_seg->start_elt;
				else
					tipelt = current_frac->last_seg->end_elt;
	
				failure = propagation(ofp, tipelt, &pangle);
 
				if (failure >= 1.0) {
 					*prop_possible = TRUE;
					if (first_tip) {
						xbbeg = -a[tipelt] * (1 + 2*cos(pangle));
						ybbeg = -2 * a[tipelt] * sin(pangle);
						xbend = -a[tipelt];
						ybend = 0;
					}
					else {
						xbbeg = a[tipelt];
						ybbeg = 0;
						xbend = a[tipelt] * (1 + 2*cos(pangle));
						ybend = 2 * a[tipelt] * sin(pangle); 	
					}
 
					sinb = sinbet[tipelt];
					cosb = cosbet[tipelt];
 
					xbeg = xm[tipelt] + (xbbeg)*cosb - (ybbeg)*sinb;
					ybeg = ym[tipelt] + (xbbeg)*sinb + (ybbeg)*cosb;
					xend = xm[tipelt] + (xbend)*cosb - (ybend)*sinb;
					yend = ym[tipelt] + (xbend)*sinb + (ybend)*cosb;
 
					/* new crack tip BC's = former crack tip BC's */
					kode = kod[tipelt];
					bvshS = bvsS[tipelt];
					bvnoS = bvnS[tipelt];
					bvshM = bvsM[tipelt];
					bvnoM = bvnM[tipelt];
 
					VFPRINTF(ofp,"      Adding new tip element\n");
					VFPRINTF(ofp,"      (xbeg,ybeg) = %f,  %f\n",
						xbeg, ybeg);
					VFPRINTF(ofp,"      (xend,yend) = %f,  %f\n",
						xend, yend);
 
					add_BEs(1,xbeg,xend,ybeg,yend,kode,bvshS,bvnoS,bvshM,bvnoM);
 
					VFPRINTF(ofp,"      Adjusting fracture segment pointers\n");
 
					temp_seg = (struct frac_seg *)
						malloc(sizeof(struct frac_seg));
					if (first_tip) {
						temp_seg->start_elt = numbe;
						temp_seg->end_elt   = numbe;
						temp_seg->next = current_frac->first_seg;
						current_frac->first_seg = temp_seg;
					}
					else {
						current_frac->last_seg->next = temp_seg;
						current_frac->last_seg = current_frac->last_seg->next;
						current_frac->last_seg->start_elt = numbe;
						current_frac->last_seg->end_elt   = numbe;
						current_frac->last_seg->next      = NULL;
					}
 
					/* if end of new fracture tip is within one       */
					/* BE length of the middle of any other element,   */
					/* 	consider it intersected and set flag.		   */
					if (first_tip)
						dist = dist_to_closest_BE(xbeg,ybeg,numbe,&closest_BE);
					else
						dist = dist_to_closest_BE(xend,yend,numbe,&closest_BE);
					if ((dist < 2*a[tipelt]) || (dist < 2*a[closest_BE])) {
						current_frac->tip_prop[j] = ISECT;
						VFPRINTF(ofp,"      New tip intersected\n");
					} /*if*/
				}
			}
		}
	current_frac = current_frac->next;
	}
}

/**************************** function:propagation **************************
*
*************************************************************************/
/* This function determines the failure criterion for a given tip element */
double  propagation(FILE *ofp, int tipelt, double *pangle)
{
	double	unneg, usneg,unpos, uspos, sign, sigs;
	double	D1, D2, k1, k2;
	double  phi;
	double	cp, sp, failure;

	calc_BE_BCs(tipelt,&usneg,&unneg,&uspos,&unpos,&sigs,&sign);
	D1 = unneg - unpos;
	D2 = usneg - uspos;
	
		
	/* Determine Ki using Olsen's method of last DD */
	k1 = -(0.806*D1*e*sqrt(M_PI)) / (4*(1-pr*pr)*sqrt(2*a[tipelt]));
	k2 = -(0.806*D2*e*sqrt(M_PI)) / (4*(1-pr*pr)*sqrt(2*a[tipelt]));
	VFPRINTF(ofp,"\n      tip BE  = %d is a fault or fracture end\n",tipelt);


	/* Mode II */
	if (k1 <= 0 || k2 > 500*k1) 
                    *pangle = (k2 > 0) ? -70.5*M_PI/180 : 70.5*M_PI/180;

	/* Mixed Mode I+II */
               else {
                  	phi = atan(3*k2/k1);
                   	*pangle = asin((k2/k1)*cos(phi)) - phi;
               } /*else*/
 
	cp = cos(*pangle/2);
	sp = sin(*pangle);

	/* Failure Criterion */ 
	failure = cp * (k1/k1c*(cp*cp) - 1.5*k2/k1c*sp);
 
	VFPRINTF(ofp,"      k1      = %f\n",k1);
	VFPRINTF(ofp,"      k2      = %f\n",k2);
	VFPRINTF(ofp,"      pangle  = %f\n",(*pangle)*180/M_PI);				
	VFPRINTF(ofp,"      failure = %f   (propagation if > 1.0)\n", failure);

	return failure;
}

/**************************** function:failure **************************
*
*************************************************************************/
/* This function determines the failure criterion for a given tip element */
double  fail(FILE *ofp, int tipelt, int top, double *pangle)
{
	double	unneg, usneg,unpos, uspos, sign, sigs;
	double	leftU, rightU;
	double  gradient, tangent;
	double	diff, maxPrinc, failure;
	
	/* Find shear displacements to the left of the element in question */
	calc_BE_BCs(tipelt-1,&usneg,&unneg,&uspos,&unpos,&sigs,&sign);
	leftU = (top) ? uspos : usneg;
	
	/* Find shear displacements to the right of the element in question */
	calc_BE_BCs(tipelt+1,&usneg,&unneg,&uspos,&unpos,&sigs,&sign);
	rightU = (top) ? uspos : usneg;
	
	calc_BE_BCs(tipelt,&usneg,&unneg,&uspos,&unpos,&sigs,&sign);
		
	/* find slip gradient using center method */
	gradient=(rightU-leftU)/(a[tipelt+1]+2*a[tipelt]+a[tipelt-1]);
	

	/* The tangent is gradient component and influence of sign
		in this case sign alrady includes necessary gravity 
		effects.  They are added in calc_Bes */
	tangent = e/(1-pr*pr) * gradient + pr/(1-pr)*sign;

	/* Failure Criterion */ 
	diff = (tangent-sign);
	maxPrinc=(tangent+sign)/2+sqrt(diff*diff/4+sigs*sigs);
	failure = (maxPrinc > tensileStrength);
 
	/*Determine angle*/
	*pangle = 0.5 * atan2( 2*sigs, diff) + M_PI/2;
	
 	VFPRINTF(ofp,"    %f     %f    %7.1f  ", tangent, maxPrinc, *pangle*180/M_PI);
	
	return failure;
}

/**************************** function:findSlipped **************************
*
*****************************************************************************/
void findSlipped(FILE *ofp, int *prop_possible2)
{
	struct fault		*current_fault;
	struct fault_seg	*temp_seg;

	static int firstTime = 0;
	int n,i,j;
	int startFault,endFault;
	int TOP = 1, BOTTOM = 0;
	int LEFT =1, RIGHT =2;
	double  pangle;
	float  epsilon= -1e-8;
        double usneg,unneg,uspos,unpos, sigs, sign;

	*prop_possible2 = FALSE;
	current_fault = first_fault;

	/* Set default TC possible accordingly to user specs for first time only 
		End elements are at start_elt and end_elt. */
	if(firstTime==0) {
		for (i=1; i<=numFaults; i++) {
			startFault = current_fault->first_seg->start_elt - nbBEs;
			endFault = current_fault->last_seg->end_elt - nbBEs;
			
			if(current_fault->tails){
				for( n=startFault+1; n<=endFault; n++) {
					TCposs[n][0] =  TCposs[n][1] = TRUE;
				}/*for*/
			}else {
				for( n=startFault+1; n<=endFault; n++) {
					TCposs[n][0] =  TCposs[n][1] = FALSE;
				}/*for*/
			}/*else*/

			/* Tail Crack Possible for first and last elements of fault line ?*/
			TCposs[startFault][0]=(current_fault->end_prop[0]==PROPYES) ? TRUE:FALSE;
			TCposs[endFault][0]=(current_fault->end_prop[1]==PROPYES) ? TRUE:FALSE;

			current_fault = current_fault->next;
		}/*for*/
		firstTime++;
	}/*if*/

	current_fault = first_fault;
	for(i=1;i<=numFaults;i++) { 
		startFault = current_fault->first_seg->start_elt - nbBEs;
		endFault = current_fault->last_seg->end_elt - nbBEs;

		/* Determine failure at start of fault */
		if(TCposs[startFault][0] &&
			propagation(ofp, startFault+nbBEs, &pangle) >1.0) {
			*prop_possible2 = TRUE;
			addTail(ofp,startFault+nbBEs,LEFT,pangle);
			TCposs[startFault][0] = FALSE;
		}/*if*/
		
		/* Determine failure at end of fault */
		if( (TCposs[endFault][0]) && 
			(propagation(ofp, nbBEs+endFault, &pangle) >1.0) ) {
			*prop_possible2 = TRUE;
			addTail(ofp,nbBEs+endFault,RIGHT,pangle);
			TCposs[endFault][0] = FALSE;
		}

		VFPRINTF(ofp,"\n  BE    Tangential   maxPrinc(top)  angle(top)  Tangential  maxPrinc(bot)  angle(bot)");

		/* Determine if stress concentration is at top or bottom of
			interface element, find failure and grow crack if neccesary.
			Note: This doesn't assess the crack tip elements.*/
		for(j=startFault+1;j<=endFault-1; j++) {
			if( TCposs[j][0] || TCposs[j][1])
				VFPRINTF(ofp,"\n %3d", nbBEs+j);
			if( TCposs[j][0] && fail(ofp,(nbBEs+j),TOP,&pangle) ) {
				*prop_possible2 = TRUE;
				addTail(ofp,(nbBEs+j-1),RIGHT,pangle);
	printf("\n New crack from left top of element %d", nbBEs+j);
				TCposs[j][0] = FALSE;
			}/*if*/
			if( TCposs[j][1] && fail(ofp,(nbBEs+j),BOTTOM,&pangle) ) {
				*prop_possible2 = TRUE;
				addTail(ofp,(nbBEs+j),LEFT,pangle);
	printf("\n New crack from left bottom of element %d", nbBEs+j);
				TCposs[j][1] = FALSE;
			}/*if*/

		}/*for*/
		current_fault = current_fault->next;
	}/*for*/

	/* Reset permSlip to false for next loading increment */
	fprintf(ofp,"\n\nThe following frictional elements slipped in this loading step:\n");
	
	current_fault = first_fault;
	for(i=1;i<=numFaults;i++) { 
		startFault = current_fault->first_seg->start_elt - nbBEs;
		endFault = current_fault->last_seg->end_elt - nbBEs;

		fprintf(ofp,"%s:",current_fault->name);
		for( n=startFault; n<=endFault; n++){
			if(permSlip[n]) {
				fprintf(ofp, " %d",n + nbBEs );
				permSlip[n]=FALSE;
			}/*if*/
		} /*for*/
		fprintf(ofp, "\n");
		current_fault = current_fault->next;
	}/*for*/

        /* Print if the element failed in tension */ 
        fprintf(ofp,"\n\nThe following frictional elements opened in this loading step:\n");

        current_fault = first_fault;
        for(i=1;i<=numFaults;i++) {
                startFault = current_fault->first_seg->start_elt - nbBEs;
                endFault = current_fault->last_seg->end_elt - nbBEs;

                fprintf(ofp,"%s:",current_fault->name);
                for( n=startFault; n<=endFault; n++){

                        if(permOpen[n]) {
                                fprintf(ofp, " %d",n + nbBEs );
                        }/*if*/
                } /*for*/
                fprintf(ofp, "\n");
                current_fault = current_fault->next;
        }/*for*/

}

/*************************** function: addTail ********************************
*
*******************************************************************************/
void addTail(FILE *ofp, int tipelt, int end, double pangle)
{

	int static counter=0;
	int LEFT=1, RIGHT=2;
	double xbbeg, ybbeg, xbend,ybend;
	double	cosb, sinb;
	double xbeg,ybeg,xend,yend;
	double dist;
	int	closest_BE;

	counter++;

	if(end==LEFT) {
		xbbeg = -a[tipelt];
		ybbeg = 0;
		xbend = -a[tipelt] * (1 + 2*cos(pangle));
		ybend = -2 * a[tipelt] * sin(pangle);
	}else if (end==RIGHT) {
		xbbeg = a[tipelt];
		ybbeg = 0;
		xbend = a[tipelt] * (1 + 2*cos(pangle));
		ybend = 2 * a[tipelt] * sin(pangle);
	} else printf("\n ERROR in code \n");

	sinb = sinbet[tipelt];
	cosb = cosbet[tipelt];
 
	xbeg = xm[tipelt] + (xbbeg)*cosb - (ybbeg)*sinb;
	ybeg = ym[tipelt] + (xbbeg)*sinb + (ybbeg)*cosb;
	xend = xm[tipelt] + (xbend)*cosb - (ybend)*sinb;
	yend = ym[tipelt] + (xbend)*sinb + (ybend)*cosb;
 					
	VFPRINTF(ofp,"\n      Adding tail crack as new fracture \n");
	VFPRINTF(ofp,"      (xbeg,ybeg) = %f,  %f\n", xbeg, ybeg);
	VFPRINTF(ofp,"      (xend,yend) = %f,  %f\n", xend, yend);

	/* New tail crack is traction free */ 					
	add_BEs(1,xbeg,xend,ybeg,yend,1,0,0,0,0);
 
	nfracs++;

	/* If this is the first fracture allocate to both first and last fracs) */
	if (last_frac == NULL)
		first_frac = last_frac = (struct frac *) malloc(sizeof(struct frac));
	/* If not, allocate last frac and point from one before it*/
	else
		last_frac = last_frac->next = (struct frac *) malloc(sizeof(struct frac));

	last_frac->next = NULL;
	last_frac->first_seg = NULL;
	last_frac->last_seg  = NULL;

	/* Tail cracks can only grow at tip farthest from fault */
	sprintf(last_frac->name, "tailCrack%d",counter);
	last_frac->tip_prop[0] = PROPNO;
	last_frac->tip_prop[1] = PROPYES;

 
	/* Allocate frac segments first=last=only one */
	last_frac->first_seg = last_frac->last_seg =
		(struct frac_seg *) malloc(sizeof(struct frac_seg));

	last_frac->last_seg->start_elt = numbe;
	last_frac->last_seg->end_elt   = numbe;
	last_frac->last_seg->next      = NULL;


	/* If end of new tailcrack tip is within ONE-HALF      	*/
	/* BE length of the middle of any other element,  	*/ 
	/*	consider it itersected and set flag.		*/
	dist = dist_to_closest_BE(xend,yend,numbe,&closest_BE);
	if ((dist < a[tipelt]) || (dist < a[closest_BE])) {
		last_frac->tip_prop[1] = ISECT;
		VFPRINTF(ofp,"      New tip intersected\n");
	} /*if*/

}

/*************************** function: faultSlip ********************************
*
*****************************************************************************/
void faultSlip(FILE *ofp, int *contin)
{				
	int i, n;
	double usneg,unneg,uspos,unpos, sigs, sign;
	double grav, sheargrav=0, normgrav=0;
	double slip, shear_strength;
	double  epsilon= -1e-10, temp_fric, signLimit;

	*contin = FALSE;
	
	for(i=1; i<=numFaultElements; i++) {
		n = nbBEs + i;

		/* Find the stresses on the element */
		calc_BE_BCs(n,&usneg,&unneg,&uspos,&unpos,&sigs,&sign);
		/* Find the slip but remove the elastic give*/
		slip = (usneg - uspos) + sigs/stiffS[i];

		if(gravity ) {
                      grav = density * 9.81 * ym[n] / 1000000;
                      sheargrav = -(grav-kratio*grav)*sinbet[n]*cosbet[n] ;
                      normgrav = -grav*(kratio*sinbet[n]*sinbet[n] + cosbet[n]*cosbet[n]);
		} /*if*/

		/*If the fault has never slipped before then I will use the inherent
                 shear strength instead of cohesion*/
		if(permSlip[i] != TRUE){
			shear_strength = strength[i];
		}else{
			shear_strength = cohesion[i];
		}

                /*Figure out the evolution of friction*/
                if(crit_slip[i] == 0.0) friction[i] = fric_dyn[i];
                if( (crit_slip[i]-fabs(slip)) > epsilon ){
                        temp_fric = fric_stat[i] - ((fric_stat[i]-fric_dyn[i])/crit_slip[i])*fabs(slip);
                /*sometimes the elements unslip due to fractures and so we need to
                 make sure that the friction doesn't go back up*/
                       if(temp_fric < friction[i]) friction[i] = temp_fric;
                 }else{
                       friction[i] = fric_dyn[i];
                 }

		/* Set the upper limit of normal stress for the frictional slip to be valid*/	
		signLimit = shear_strength/friction[i];

		/* OPTION 1: if the normal stresses are tensile... */
		if ( sign >= epsilon + faultTen[i]) {
			/*printf("\n Fault element %d failed in tension ", n);*/
			faultTen[i] = 0.0;
			b[2*n-1] = sheargrav;
			b[2*n] = normgrav;
			slipped[i] = TENSILE;
	              	permOpen[i] = TRUE;
 
			/*Because once an element starts to open, it can influence neighboring
			  elements, we need to test convergence here as well.*/ 
			if ( toler < fabs(sign) + epsilon ) {
                                *contin = TRUE;
                        }
		
		/* OPTION 2: The shear stress is compared to the frictional resistance.  If the
		slip condition is met, the shear stress is constrained to the yield
		value. The continue condition is set to true to iterate process and
		the slip condition is set to true for altering influence coefficients. 
 	        this equation for slip is usually written
                         ( fabs( sigs ) >= cohesion[i] - sign * friction[i])
                but different machines sometimes treat == differently. So instead I will 
		formulate this using epsilon for the difference. */

		}else if(fabs( sigs ) - shear_strength + sign * friction[i] > epsilon ) {

			/*Since we want the element to go from strength to cohesion as soon as it
			  slips, we can use cohesion in the following equations.*/

			if(sigs <0)
				b[2*n-1] = - (cohesion[i] - sign * friction[i])+sheargrav;
			else 
				b[2*n-1] = cohesion[i] - sign * friction[i]+sheargrav;
			
			/*If the normal stress is grreater than the limit then tau = 0 */
			if(sign >= signLimit) b[2*n-1] = 0;

			b[2*n] = normgrav;
			slipped[i] = TRUE;
			permSlip[i] = TRUE;

			/* The condition for complete convergence is that the difference between 
				the new stress after slip is close to the old stress */
			if ( toler <  fabs(  b[2*n-1]-sheargrav - sigs) / (1+fabs(sigs)) ) {
				*contin = TRUE;
			}

		/* OPTION 3: If the slip condition is not met then the element could have 'over slipped'
                in this case the tolerance is not met and th element should slip again */
		} else if(slipped[i] == TRUE && toler <  fabs(  b[2*n-1]-sheargrav - sigs) / (1+fabs(sigs)) ) {

                        if(sigs <0)
                                b[2*n-1] = - (cohesion[i] - sign * friction[i])+sheargrav;
                        else
                                b[2*n-1] = cohesion[i] - sign * friction[i]+sheargrav;

                        b[2*n] = normgrav;
                        slipped[i] = TRUE;
			*contin = TRUE;

                /* OPTION 4: If the slip condition is not met, the element never slipped and the tolerance is fine,
		then the normal BC is set to zero and
                the shear is set to the total minus the elastic contribution.
                Since the elastic components are added into the adjust_Coeff sunroutine
                by way of the stiffnesses these are not specified a priori. */

		} else {
			b[2*n-1] = sigs + d[2*n-1] * stiffS[i]+sheargrav;
			b[2*n] = normgrav;
			slipped[i] = FALSE;
		} /*else*/

	} /*for*/
}


/*************************** function: increment_BEs ********************************
*
*****************************************************************************/
void increment_BEs(int step, int nsteps)
{				
	int i;
	double stepLoad;

	/* Determine the % of final load during incremental loading */
	stepLoad = (double) step / nsteps;

	/* Increment Conditions on Boundary Elements */
	for ( i = 1; i <= nbBEs; i++ ) {
		stepbvs[i] = stepLoad * bvsM[i];
		stepbvn[i] = stepLoad * bvnM[i];
	} /*for*/


	/* Increment Conditions on Fracture Elements */
	for ( i = nbBEs+numFaultElements+1; i <= numbe; i++ ) {
		stepbvs[i] = stepLoad * bvsM[i];
		stepbvn[i] = stepLoad * bvnM[i];
	} /*for*/
}


/*************************** function: reassign_BEs ***************************
*
*****************************************************************************/
void reassign_BEs()
{	
	int i;	

	/* This loop assigns the boundary values 
		Need to reassign boundary conditions each iteration 
		because quicksolve doesn't return the same values of b[] as given.
		b[] for fault elements are reassigned in the fault slip function 
		according to the slip history of the element. */	

	/* STATIC and MONOTONIC boundary conditions are added */

	/* Reset conditions on Boundary Elements */
	for ( i = 1; i <= nbBEs; i++ ){
		b[(2*i-1)] = bvsS[i] + stepbvs[i];
		b[(2*i)] = bvnS[i] + stepbvn[i];
	   } /*for*/

	/* Reset Conditions on Fracture Elements */
	for ( i = nbBEs+numFaultElements+1; i <= numbe; i++ ) {
		b[(2*i-1)] = bvsS[i] + stepbvs[i];
		b[(2*i)] = bvnS[i] + stepbvn[i];
	   } /*for*/

}
/**************************** function: calc_coeff **************************
*
****************************************************************************/
void calc_coeff(double xi, double yi, double xj, double yj, double aj,
		double cosbj, double sinbj)
{
 
	double xji, yji;				/* coords for image elements */
 
	coeff(xi,yi,xj,yj,aj,cosbj,sinbj,1);
	xji = 2.*xsym - xj;
	yji = 2.*ysym - yj;
	switch (ksym) {
		case NOSYM:
			break;
		case XLSYM:
			coeff(xi,yi,xji,yj,aj,cosbj,-sinbj,-1);
			break;
		case YLSYM:
			coeff(xi,yi,xj,yji,aj,-cosbj,sinbj,-1);
			break;
		case XYLSYM:
			coeff(xi,yi,xji,yj,aj,cosbj,-sinbj,-1);
			coeff(xi,yi,xj,yji,aj,-cosbj,sinbj,-1);
			coeff(xi,yi,xji,yji,aj,-cosbj,-sinbj,1);
			break;
		case XYASYM:
			coeff(xi,yi,xji,yji,aj,-cosbj,-sinbj,1);
			break;
	} /*switch*/
}
 
 
/***************************** function: coeff ******************************
*
****************************************************************************/
void coeff(double x, double y, double cx, double cy, double a,
	   double cosb, double sinb, int msym)
{
	/****    global variables    ****/
	/* uses "common blocks" s1 & s2 */
 
	/**** automatic variables ****/
	double cos2b, sin2b, cosb2, sinb2;
	double xb, yb, r1s, r2s, fl1, fl2;
	double fb2, fb3, fb4, fb5, fb6, fb7;
	double uxds, uxdn, uyds, uydn;
	double sxxds, sxxdn, syyds, syydn, sxyds, sxydn;
 
	/**** function body ****/
	cos2b = cosb*cosb - sinb*sinb;
	sin2b = 2.*sinb*cosb;
	cosb2 = cosb*cosb;
	sinb2 = sinb*sinb;
 
	xb = (x-cx)*cosb + (y-cy)*sinb;
	yb = -(x-cx)*sinb + (y-cy)*cosb;
 
	r1s = (xb-a)*(xb-a) + yb*yb;
	r2s = (xb+a)*(xb+a) + yb*yb;
	fl1 = 0.5 * log(r1s);
	fl2 = 0.5 * log(r2s);
	fb2 = con*(fl1-fl2);
 
	if (yb == 0) {
		fb3 = 0;
		if (fabs(xb) < a)
			fb3 = con * M_PI;
	} /*if*/
	else
		fb3 = -con * (atan((xb+a)/yb) - atan((xb-a)/yb));
 
	fb4 = con * (yb/r1s - yb/r2s);
	fb5 = con * ((xb-a)/r1s - (xb+a)/r2s);
	fb6 = con * (((xb-a)*(xb-a)-yb*yb)/(r1s*r1s)
				- ((xb+a)*(xb+a)-yb*yb)/(r2s*r2s));
	fb7 = 2.*con*yb*((xb-a)/(r1s*r1s) - (xb+a)/(r2s*r2s));
 
	uxds = -pr1*sinb*fb2 + pr2*cosb*fb3 + yb*(sinb*fb4-cosb*fb5);
	uxdn = -pr1*cosb*fb2 - pr2*sinb*fb3 - yb*(cosb*fb4+sinb*fb5);
	uyds =  pr1*cosb*fb2 + pr2*sinb*fb3 - yb*(cosb*fb4+sinb*fb5);
	uydn = -pr1*sinb*fb2 + pr2*cosb*fb3 - yb*(sinb*fb4-cosb*fb5);
 
	sxxds = cons * (2.*cosb2*fb4 + sin2b*fb5 + yb*(cos2b*fb6 - sin2b*fb7));
	sxxdn = cons * (-fb5 + yb*(sin2b*fb6 + cos2b*fb7));
	syyds = cons * (2.*sinb2*fb4 - sin2b*fb5 - yb*(cos2b*fb6 - sin2b*fb7));
	syydn = cons * (-fb5 - yb*(sin2b*fb6 + cos2b*fb7));
	sxyds = cons * (sin2b*fb4 - cos2b*fb5 + yb*(sin2b*fb6+cos2b*fb7));
	sxydn = cons * (-yb * (cos2b*fb6 - sin2b*fb7));
 
	uxs += msym*uxds;
	uxn += uxdn;
	uys += msym*uyds;
	uyn += uydn;
 
	sxxs += msym*sxxds;
	sxxn += sxxdn;
	syys += msym*syyds;
	syyn += syydn;
	sxys += msym*sxyds;
	sxyn += sxydn;
 
}

	
/*************************** function: adjustCoeff **************************
*
****************************************************************************/
void adjustCoeff()
{
	/**** automatic variables ****/
	int i, n;

	for (i=1; i<=numFaultElements; i++)
		{n = nbBEs + i;

		/* The normal stresses are related to stiffness if non-tensile */ 
		
		if(slipped[i] != TENSILE) /* PUT bacK IN */
			c[2*n][2*n] = c[2*n][2*n] + stiffN[i];
		
		/* The shear stresses are related to stiffness when there is no slip */
		if(slipped[i] == FALSE)
			c[2*n-1][2*n-1] = c[2*n-1][2*n-1] + stiffS[i];
	
} /*for*/
} 


 
/************************** function calc_inf_coeffs ************************
*
*****************************************************************************/
void calc_inf_coeffs()
{
	/**** automatic variables ****/
	int i, in, is;
	double cosbi, sinbi;
	int kode;
	int j, jn, js;
 
	for (i=1; i<=numbe; i++) {
		in = 2*i;
		is = in-1;
		cosbi = cosbet[i];
		sinbi = sinbet[i];
		kode = kod[i];
		for (j=1; j<=numbe; j++) {
			jn = 2*j;
			js = jn-1;
			sxxs = sxxn = syys = syyn = sxys = sxyn = 0;
			uxs = uxn = uys = uyn = 0;
			calc_coeff(xm[i],ym[i],xm[j],ym[j],a[j],cosbet[j],
					sinbet[j]);
 
			switch (kode) {
				case 1:
					c[is][js] = (syys-sxxs)*sinbi*cosbi
								+ sxys*(cosbi*cosbi-sinbi*sinbi);
					c[is][jn] = (syyn-sxxn)*sinbi*cosbi
								+ sxyn*(cosbi*cosbi-sinbi*sinbi);
					c[in][js] = sxxs*sinbi*sinbi - 2.*sxys*sinbi*cosbi
								+ syys*cosbi*cosbi;
					c[in][jn] = sxxn*sinbi*sinbi - 2.*sxyn*sinbi*cosbi
								+ syyn*cosbi*cosbi;
					break;
				case 2:
					c[is][js] = uxs*cosbi + uys*sinbi;
					c[is][jn] = uxn*cosbi + uyn*sinbi;
					c[in][js] = -uxs*sinbi + uys*cosbi;
					c[in][jn] = -uxn*sinbi + uyn*cosbi;
					break;
				case 3:
					c[is][js] = uxs*cosbi + uys*sinbi;
					c[is][jn] = uxn*cosbi + uyn*sinbi;
					c[in][js] = sxxs*sinbi*sinbi - 2.*sxys*sinbi*cosbi
								+ syys*cosbi*cosbi;
					c[in][jn] = sxxn*sinbi*sinbi - 2.*sxyn*sinbi*cosbi
								+ syyn*cosbi*cosbi;
					break;
				case 4:
					c[is][js] = (syys-sxxs)*sinbi*cosbi
								+ sxys*(cosbi*cosbi-sinbi*sinbi);
					c[is][jn] = (syyn-sxxn)*sinbi*cosbi
								+ sxyn*(cosbi*cosbi-sinbi*sinbi);
					c[in][js] = -uxs*sinbi + uys*cosbi;
					c[in][jn] = -uxn*sinbi + uyn*cosbi;
					break;
			} /*switch*/
		} /*for*/
	} /*for*/
}
 
 
/************************** function: calc_BE_BCs *************************
*
****************************************************************************/
void calc_BE_BCs(int i, double * usneg, double * unneg,
		 double * uspos, double * unpos, double * sigs, double * sign)
{
	/**** automatic variables ****/
	int in, is;
	double cosbi, sinbi;
	double uxneg, uyneg, sigxx, sigyy, sigxy, grav;
 
	in = 2*i;
	is = in-1;
	cosbi = cosbet[i];
	sinbi = sinbet[i];

	disps_stresses_at_pt(xm[i],ym[i],&uxneg,&uyneg,&sigxx,&sigyy,&sigxy,
		FALSE);

	if(gravity) {
		grav = density * 9.81 * ym[i] / 1000000;
		sigyy = sigyy + grav;
		sigxx = sigxx + kratio * grav;
	}
 
	*usneg = uxneg*cosbi + uyneg*sinbi;		/* these x and y are global coordinate system */
	*unneg = -uxneg*sinbi + uyneg*cosbi;
	*uspos = *usneg - d[is];
	*unpos = *unneg - d[in];
	*sigs = (sigyy-sigxx)*sinbi*cosbi + sigxy*(cosbi*cosbi-sinbi*sinbi);
	*sign = sigxx*sinbi*sinbi - 2.*sigxy*sinbi*cosbi + sigyy*cosbi*cosbi;
}
 
 
/********************** function: disps_stresses_at_pt **********************
*
***************************************************************************/
int disps_stresses_at_pt(double xp, double yp, double *ux, double *uy,
			 double *sigxx, double *sigyy, double *sigxy, int distcheck)
{
	/**** automatic variables ****/
	int j, jn, js;
	int closest_BE;
 
 
	*ux = *uy = 0;
	*sigxx = 0.0;
	*sigyy = 0.0;
	*sigxy = 0.0;
	
	if (distcheck) {
		if ((dist_to_closest_BE(xp,yp,0,&closest_BE)) < a[closest_BE])
			return -1;
	} /*if*/
	for (j=1; j<=numbe; j++) {
		jn = 2*j;
		js = jn-1;
		sxxs = sxxn = syys = syyn = sxys = sxyn = 0;
		uxs = uxn = uys = uyn = 0;
		calc_coeff(xp,yp,xm[j],ym[j],a[j],cosbet[j],sinbet[j]);
 
		*ux += uxs*d[js] + uxn*d[jn];
		*uy += uys*d[js] + uyn*d[jn];
		*sigxx += sxxs*d[js] + sxxn*d[jn];
		*sigyy += syys*d[js] + syyn*d[jn];
		*sigxy += sxys*d[js] + sxyn*d[jn];
	} /*for*/

	return 1;
}
 
/************************ function: dist_to_closest_BE ***********************
*
*****************************************************************************/
double dist_to_closest_BE(double x, double y, int dont_check_BE, int *closest_BE)
 
#define DIST(L1,L2) sqrt((L1)*(L1) + (L2)*(L2))
#define MIN(A,B) (A < B) ? A:B
 
{
	double xi, yi;
	double mindist = -1;
	double oldmindist;
	int i;
 
	xi = 2.*xsym - x;
	yi = 2.*ysym - y;
	for (i=1; i<=numbe; i++) {
		oldmindist = mindist;
		if (i != dont_check_BE) {
			if (mindist == -1)
				mindist = DIST(x-xm[i], y-ym[i]);
			else
				mindist = MIN(mindist, DIST(x-xm[i], y-ym[i]));
 
			switch (ksym) {
				case NOSYM:
					break;
				case XLSYM:
					mindist = MIN(mindist, DIST(xi-xm[i], y-ym[i]));
					break;
				case YLSYM:
					mindist = MIN(mindist, DIST(x-xm[i], yi-ym[i]));
					break;
				case XYLSYM:
					mindist = MIN(mindist, DIST(xi-xm[i], y-ym[i]));
					mindist = MIN(mindist, DIST(x-xm[i], yi-ym[i]));
					mindist = MIN(mindist, DIST(xi-xm[i], yi-ym[i]));
					break;
				case XYASYM:
					mindist = MIN(mindist, DIST(xi-xm[i], yi-ym[i]));
					break;
			} /*switch*/
		} /*if*/
		if (oldmindist != mindist)
			*closest_BE = i;
	} /*for*/
	DPRINTF("%s IN DIST_TO_CLOSEST_BE: distance = %f\n",DB,mindist);
	return (mindist);
}
 
/**************************** function: add_BEs ****************************
*
* Define locations, sizes, orientations, and boundary conditions of
* boundary elements
*
****************************************************************************/
void add_BEs(int num, double xbeg, double xend, double ybeg, double yend,
	     int kode,double bvshS, double bvnoS,double bvshM, double bvnoM)
{
	int ne, m;
	double xd, yd, sw;
	double grav;
 	static int counter= 0;		/* counts times entered this function*/

	++counter;
	xd = (xend-xbeg)/num;
	yd = (yend-ybeg)/num;
	sw = sqrt(xd*xd + yd*yd);
	for (ne=1; ne<=num; ne++) {
		m = ++numbe;
		xm[m] = xbeg + 0.5*(2.*ne-1)*xd;
		ym[m] = ybeg + 0.5*(2.*ne-1)*yd;
		a[m] = 0.5*sw; 			/*half length of element*/
		sinbet[m] = yd/sw;		/* a, sinbet and cosbet remain */
		cosbet[m] = xd/sw;		/* constant for one straight length */
		
		kod[m] = kode;
	

 		/*If there is gravity, then along elements 
		  the negative of gravitational stresses
                  should be entered as the STATIC boundary conditions */

                  if(gravity ) {
                      grav = density * 9.81 * ym[m] / 1000000;

                      switch (kod[m]) {
                          case 1:
                                bvsS[m] = bvshS-(grav-kratio*grav)*sinbet[m]*cosbet[m] ;
                                bvnS[m] = bvnoS-grav*(kratio*sinbet[m]*sinbet[m] + cosbet[m]*cosbet[m]);
				break;
                          case 2:
                                break;
                          case 3:
                                bvnS[m] = bvnoS-grav*(kratio*sinbet[m]*sinbet[m] + cosbet[m]*cosbet[m]);
                                break;
                          case 4:
                                bvsS[m] = bvshS-(grav-kratio*grav)*sinbet[m]*cosbet[m] ;
                                break;
             		} /*switch*/
			/* If the element is a fault we need to adjust the BEs for gravity
                		before the influence coefficient matrix is inverted.
				Even though the b[] are assigned for boundaries and
				fractures in reassign_BEs.  We do all the elements here
				since fault cannot just be distinguished at this point. */
				b[(2*m-1)] = bvsS[m] ;
               			b[(2*m)] = bvnS[m] ;

		}else{

		   bvsS[m] = bvshS;	/* store STATIC boundary conditions*/
		   bvnS[m] = bvnoS; 
  		   }/*else*/

		bvsM[m] = bvshM;	/* store final boundary conditions*/
		bvnM[m] = bvnoM; 	/*	for MONOTONIC loading	  */	
			
	}  /*for*/
}
 
 
/************************* function: quicksolve ****************************
*
****************************************************************************/
 
/* This function solves the system of algebraic equations for the vector d[].
	The vector d[] contains the conditions along the boundary elements 
	July 2013. Updating this to LUDCMP from numerical recipes*/ 

void quicksolve(int n, float *condnum)
{
	int i,imax,j,k;
	int indx[n+1];
	double ri, sum=0.0;
        double normOriginal = 0.0, normInverse = 0.0;
	double col[n+1];
	double** original = dmatrix(1, n, 1, n);
	double** inverse = dmatrix(1, n, 1, n);

	/*printf("\n");*/
	/*Store the original C  matrix in a local matrix and find norm*/
	for (i = 1; i<=n; i++) {
                for (j=1; j<=n; j++) {
                        original[i][j] = c[i][j];
                        /*printf("%f ",c[i][j]);*/
			sum = sum + fabs(c[i][j]);
                }
		/*printf("\n");*/
                if(sum > normOriginal){
                        normOriginal = sum;
                }
                sum=0.0;
        }
        /*printf("norm = %f \n", normOriginal);*/

	/*The matrix c and number n are provided to ludcmp, which returns the
	index and the row_interchange, ri
	how is ri used?*/
	ludcmp(c, n, indx, &ri);
	/*Now C is factorized by LU decomposition */
	
	/*Now run through LUBKSB in order to solve for d*/
	lubksb(c, n, indx, &b[0]);

	/*Now the vector b contains the solution for d*/
	for(i=1; i<=n; i++) {
		d[i] = b[i];
	}
	
	/*Now we can use LU decomposed C matrix to find the inverse by solving
		for the matrix that would when multiplied give the indentity matrix*/ 
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++) col[i]=0.0;
		col[j]=1.0;
		lubksb(c,n,indx,col);
		for(i=1;i<=n;i++) inverse[i][j] = col[i];
	} 

	sum=0.0;
        for (i = 1; i<=n; i++) {
                for (j=1; j<=n; j++) {
                    	sum = sum + fabs(inverse[i][j]);
                }
                if(sum > normInverse){
                        normInverse = sum;
                }
                sum=0.0;
        }
        /*printf("norm Inverse = %e \n", normInverse);*/

	*condnum = normOriginal * normInverse;

	free_matrix(original, 1, n, 1, n);
	free_matrix(inverse, 1, n, 1, n);
}


 /*********************** function: print_fault_data ****************************
*
****************************************************************************/
void print_fault_data(FILE *ofp, int first_elt, int last_elt, int titles)
{
	int i, n;
	double xbeg, ybeg, xend, yend;
	double angle;
 
	if (titles) {
		fprintf(ofp," Elt     X-beg     Y-beg     X-end ");
		fprintf(ofp,"    Y-end     Length   angle     Ks  ");
		fprintf(ofp,"     Kn       Ten     So    cohes     fric   mu-stat mu-dyn  slip-dis\n");
		fprintf(ofp,"---- --------- --------- --------- ");
		fprintf(ofp,"--------- ----------  -------  ------- ");
		fprintf(ofp," -------  ------  ------  ------   ------   -----   -----  --------\n");
	}
 
	for (i=first_elt; i<=last_elt; i++) {
		xbeg = xm[i] - a[i]*cosbet[i];
		ybeg = ym[i] - a[i]*sinbet[i];
		xend = xm[i] + a[i]*cosbet[i];
		yend = ym[i] + a[i]*sinbet[i];
		angle = atan2(sinbet[i],cosbet[i]) * 180/M_PI;
		fprintf(ofp,"%4d %9.5f %9.5f %9.5f %9.5f %10.6f %8.3f ",
			i, xbeg, ybeg, xend, yend, 2*a[i], angle);
		n = i- nbBEs;
		fprintf(ofp,"%8.1e %8.1e %7.3f %7.3f %7.3f ",
			stiffS[n], stiffN[n], faultTen[n], strength[n], cohesion[n]);
                fprintf(ofp,"%8.4f %7.3f %7.3f %9.2e\n",
                        friction[n], fric_stat[n], fric_dyn[n], crit_slip[n]);

	} /*for*/
}
 
/*********************** function: print_BE_data ****************************
*
****************************************************************************/
void print_BE_data(FILE *ofp, int first_elt, int last_elt, int titles)
{
	int i;
	double xbeg, ybeg, xend, yend;
	double angle;
 
	if (titles) {
		fprintf(ofp,"\t\t\t\t\t\t\t\t\t\t\t\t    STATIC \t\t MONOTONIC\n");
		fprintf(ofp," Elt        X-beg        Y-beg        X-end ");
		fprintf(ofp,"       Y-end       Length    Angle kode ");
		fprintf(ofp," US or Sig-S  UN or Sig-N");
		fprintf(ofp," US or Sig-S  UN or Sig-N\n");
		fprintf(ofp,"---- ------------ ------------ ------------ ");
		fprintf(ofp,"------------ ------------ -------- ---- ");
		fprintf(ofp,"------------ ------------ ------------ ------------\n");
	}
 
	for (i=first_elt; i<=last_elt; i++) {
		xbeg = xm[i] - a[i]*cosbet[i];
		ybeg = ym[i] - a[i]*sinbet[i];
		xend = xm[i] + a[i]*cosbet[i];
		yend = ym[i] + a[i]*sinbet[i];
		angle = atan2(sinbet[i],cosbet[i]) * 180/M_PI;
		fprintf(ofp,"%4d %12.5f %12.5f %12.5f %12.5f %12.5f %8.3f ",
			i, xbeg, ybeg, xend, yend, 2*a[i], angle);
		fprintf(ofp,"%4d %12.5e %12.5e %12.5e %12.5e\n",
			kod[i], bvsS[i], bvnS[i], stepbvs[i], stepbvn[i]);
	} /*for*/
}
 
/******************* function: print_BE_disps_stresses ********************
*
* Prints the displacements and stresses at the midpoints of the boundary
* elements first_elt through last_elt.
*
****************************************************************************/
void print_BE_disps_stresses(FILE *ofp, int first_elt, int last_elt, int titles)
{
	double usneg, unneg, uspos, unpos;
	double sigs, sign ;
	int i,in,is;
	
	if (titles) {
		fprintf(ofp," Elt           DS        US(-)        US(+) ");
		fprintf(ofp,"          DN        UN(-)        UN(+)      Sigma-S ");
		fprintf(ofp,"     Sigma-N\n");
		fprintf(ofp,"---- ------------ ------------ ------------ ");
		fprintf(ofp,"------------ ------------ ------------ ------------ ");
		fprintf(ofp,"------------\n");
	} /*if*/
 
	for (i=first_elt; i<=last_elt; i++) {
		in = 2*i;
		is = in-1;
		calc_BE_BCs(i,&usneg,&unneg,&uspos,&unpos,&sigs,&sign);
 		
		fprintf(ofp,"%4d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e ",
			i, d[is], usneg, uspos, d[in], unneg, unpos);
		fprintf(ofp,"%12.5e %12.5e\n",sigs, sign);
	} /*for*/
}
 
 
/************************* function: print_prob_info ***********************
*
****************************************************************************/
void
print_prob_info(FILE *ofp)
{
	int i;


	fprintf(ofp,"This file was created with Fric2D version 3.2.7 - revision date April 2014 \n");	
	fprintf(ofp,"\nTITLE1: %s\n", title1);
	fprintf(ofp,"TITLE2: %s\n\n", title2);
 
	fprintf(ofp,"SYMMETRY: ");
	switch (ksym) {
		case NOSYM:
			fprintf(ofp,"none");
			break;
		case XLSYM:
			fprintf(ofp,"x = %f is a line of symmetry", xsym);
			break;
		case YLSYM:
			fprintf(ofp,"y = %f is a line of symmetry", ysym);
			break;
		case XYLSYM:
			fprintf(ofp,"x = %f & y = %f are lines of symmetry", xsym, ysym);
			break;
		case XYASYM:
			fprintf(ofp,"(x,y) = (%f,%f) is a point of anti-symmetry",
					xsym, ysym);
	} /*switch*/
	fprintf(ofp,"\n\n");

	fprintf(ofp,"Data will be written from the following iterations:\n");
	if (pdfiter[0] == 0)
		fprintf(ofp,"    (never)\n\n");
	else {
		fprintf(ofp,"    ");
		for (i=0; pdfiter[i] != 0; i++)
			fprintf(ofp,"%d  ",pdfiter[i]);
		fprintf(ofp,"\n\n");
	}		
	
	if (pdfiter[0] != 0) {
		fprintf(ofp,"     Boundary Data:..........%s\n",(bdata) ? "yes":"no");
		fprintf(ofp,"     Fracture Data:..........%s\n",(fracData) ? "yes":"no");
		fprintf(ofp,"     Fault Data:.............%s\n",(faultData) ? "yes":"no");
		fprintf(ofp,"     Observation Line Data:..%s\n\n",(odata) ? "yes":"no");
	} /*if*/

	fprintf(ofp,"     Superpose Gravitational Stresses ....%s\n",(gravity) ? "yes":"no");
	fprintf(ofp,"    %12f = Ratio of horizontal:vertical gravitational stresses\n",kratio);
	fprintf(ofp,"    %12f = Rock Density\n\n",density);

	fprintf(ofp,"%12d = number of straight-line boundary segments\n",nblines);
	fprintf(ofp,"%12d = number of fractures\n",nfracs);
	fprintf(ofp,"%12d = number of faults\n",numFaults);
	fprintf(ofp,"%12d = number of observation lines\n\n",nolines);
	fprintf(ofp,"%12f = Poisson's ratio\n",pr);
	fprintf(ofp,"%12f = Young's modulus\n",e);
	fprintf(ofp,"%12f = Fracture Toughness\n",k1c);
	fprintf(ofp,"%12f = Tensile Strength\n", tensileStrength);
	fprintf(ofp,"%12f = Convergence Tolerance\n",toler);
	fprintf(ofp,"%12d = Iteration Limit\n",tolCounter);
}
 
 
 
/*************************** function: startup ***************************
*
* Reads the problem info from the input file.  Sets up initial boundary/
* fracture BE array.  Sets fracture pointers.
*
**************************************************************************/
int startup(FILE *ifp, FILE *rfp, int *niters, int *nsteps)
{
	int i, numwords;		/* numwords - # of words gotten by getwords()	*/
	char line[MAXLINE], *word[MAXWORDS];
	int num, kode;
	double xbeg, ybeg, xend, yend, bvshS, bvnoS, bvshM, bvnoM;
	struct oline *current_oline = NULL;
 
 	DPRINTF("INSIDE STARTUP...STARTING VARIABLE PROCESSING LOOP\n");
 	

	/* read variable list */
	for (;;) {
 
 		DPRINTF("AT TOP OF PROCESSING LOOP...READING A LINE\n");
 		
		/* get line, exit function on EOF */
		if ((numwords = getwords(ifp,line,MAXLINE,word,MAXWORDS)) == EOF) {
			fprintf(stderr,"\nerror: Unexpected EOF in input file\n");
			return (-1);
		} /*if*/
		
		DPRINTF("LINE = %s\n",line);
 
		/* write line to resume file, if requested */
		if (rfp != NULL) {
			if ((fputs(line, rfp)) == EOF) {
				fprintf(stderr,"\nerror: Error writing to resume file\n");
				return (-1);
			} /*if*/
		} /*if*/
 
		/* skip blank and comment lines */
		if (numwords == 0 || word[0][0] == COMMENT_CHAR)
			;
 
		/* exit loop when end of variable list reached */
		else if (!strcmp(word[0],"end"))
			break;
 
		/* parse variables */
		else if (numwords >= 3 && !strcmp(word[0],"title1")) {
			strncpy(title1,word[2],MAXTITLE-1);
			title1[MAXTITLE-1] = '\0';
		}
		else if (numwords >= 3 && !strcmp(word[0],"title2")) {
			strncpy(title2,word[2],MAXTITLE-1);
			title2[MAXTITLE-1] = '\0';
		}
		else if (numwords >= 3 && !strcmp(word[0],"ksym"))
			ksym = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"xsym"))
			xsym = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"ysym"))
			ysym = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"pr"))
			pr = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"e"))
			e = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"tensileStrength"))
			tensileStrength = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"gravity"))
			gravity = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"kratio"))
			kratio = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"density"))
			density = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"nolines"))
			nolines = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"nblines"))
			nblines = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"k1c"))
			k1c = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"bdata"))
			bdata = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"fracData"))
			fracData = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"faultData"))
			faultData = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"odata"))
			odata = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"nsteps"))
			*nsteps = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"niters"))
			*niters = atoi(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"tolerance"))
			toler = atof(word[2]);
                else if (numwords >= 3 && !strcmp(word[0],"tolCounter"))
                        tolCounter = atof(word[2]);
		else if (numwords >= 3 && !strcmp(word[0],"pdfiter")) {
			for (i=3; i<=numwords ;i++)
				pdfiter[i-3] = atoi(word[i-1]);
			pdfiter[numwords-2] = 0;
		} /*if*/	
	} /*for*/

 
	/*---------------------------------
	 Calculate some needed constants
	---------------------------------*/
	con = 1./(4.*M_PI*(1.-pr));
	cons = e/(1.+pr);
	pr1 = 1.-2.*pr;
	pr2 = 2.*(1.-pr);
 
	/* create a linked list of observation lines */
	i = 1;
	while (i <= nolines) {
 
		/* get line, exit function on EOF */
		if ((numwords = getwords(ifp,line,MAXLINE,word,MAXWORDS)) == EOF) {
			fprintf(stderr,"\nerror: Unexpected EOF in input file\n");
			return (-1);
		} /*if*/
 
		/* write line to resume file, if requested */
		if (rfp != NULL) {
			if ((fputs(line, rfp)) == EOF) {
				fprintf(stderr,"\nnerror: Error writing to resume file\n");
				return (-1);
			} /*if*/
		} /*if*/
 
		/* skip blank and comment lines */
		if (numwords == 0 || word[0][0] == COMMENT_CHAR)
			;
 
		/* create the new node */
		else {
			/* if first obs line, create first node in linked list */
			if (current_oline == NULL) {
				current_oline = (struct oline *) malloc(sizeof(struct oline));
				first_oline = current_oline;
			}
			else {
				current_oline->next = (struct oline *)
					malloc(sizeof(struct oline));
				current_oline = current_oline->next;
			}
			current_oline->xbeg  = atof(word[0]);
			current_oline->ybeg  = atof(word[1]);
			current_oline->xend  = atof(word[2]);
			current_oline->yend  = atof(word[3]);
			current_oline->numpb = atoi(word[4]);
			current_oline->next  = NULL;
			i++;
		}
	}

					
	/* add the boundary lines to the BE arrays */
	i = 1;
	while (i <= nblines) {
 
		/* get line, exit function on EOF */
		if ((numwords = getwords(ifp,line,MAXLINE,word,MAXWORDS)) == EOF) {
			fprintf(stderr,"\nerror: Unexpected end of file\n");
			return (-1);
		} /*if*/
 
		/* write line to resume file, if requested */
		if (rfp != NULL) {
			if ((fputs(line, rfp)) == EOF) {
				fprintf(stderr,"\nnerror: Error writing to resume file\n");
				return (-1);
			} /*if*/
		} /*if*/
 
		/* skip blank and comment lines */
		if (numwords == 0 || word[0][0] == COMMENT_CHAR)
			;
 
		else {
			num  = atoi(word[0]);
			xbeg = atof(word[1]);
			ybeg = atof(word[2]);
			xend = atof(word[3]);
			yend = atof(word[4]);
			kode = atoi(word[5]);
			bvshS = atof(word[6]);
			bvnoS = atof(word[7]);
			bvshM = atof(word[8]);
			bvnoM = atof(word[9]);
			add_BEs(num,xbeg,xend,ybeg,yend,kode,bvshS,bvnoS,bvshM,bvnoM);
			nbBEs += num;
			i++;
		}
	}
 
	/* Read fractures and faults */
	for (;;) {
		int grabLine = 1;

		/*
		 * readFracture (see below) always reads one too many lines because
		 * it doesn't know how many fracture segments there are in the input
		 * file.  So, we only grab a line if we didn't just read a fracture:
		 */
		if (grabLine)
			numwords = getwords(ifp, line, MAXLINE, word, MAXWORDS);

		/* Exit loop if we get end of file */
		if (numwords == EOF) break;

		/* Skip blank lines and comment lines */
		if (numwords == 0 || word[0][0] == COMMENT_CHAR) continue;

		/* If the first word on the line is fracture... */
		if (strcmp(word[0], "fracture") == 0) {
			int readError = readFracture(ifp, line, word, numwords);
			if (readError == -1) break;
			grabLine = 0;
		}
		/* If the first word on the line is fault... */
		else if (strcmp(word[0], "fault") == 0) {
			int readError = readFault(ifp,line, word, numwords);
			if (readError == -1) break;
			grabLine = 1;
		}
	}
}

/*
 * This routine is called after the first line of a fracture (with the "fracture"
 * keyword) is read.  The input file, a line buffer used to read into, an array of
 * words on that line, and the number of words read are returned.  This routine
 * makes sure the right stuff was read on the first line and then grabs the next
 * line of input into line/word/numwords, and allocates a new fracture with the right
 * info.  This routine returns 0 if everything is read OK, -1 if there was a read error.
 * The loop for reading fracture segments is stopped when it encounters a blank line.
 * This insures that data isn't lost by being read before the break command. 
 */
int
readFracture(FILE *ifp, char *line, char **word, int numwords)
{
	int num, kode, i;
	double xbeg, ybeg, xend, yend, bvshS, bvnoS, bvshM, bvnoM;

	nfracs++;
	if (last_frac == NULL)
		first_frac = last_frac = (struct frac *)
			malloc(sizeof(struct frac));
	else
		last_frac = last_frac->next = (struct frac *)
			malloc(sizeof(struct frac));
	last_frac->next = NULL;
	last_frac->first_seg = NULL;
	last_frac->last_seg  = NULL;		
 
	if (numwords != 4) {
		fprintf(stderr,"\nerror: Missing fracture name or propagation ");
		fprintf(stderr,"flags in input file\n");
		return (-1);
	} /*if*/
 
	strncpy(last_frac->name,word[1],MAXFNAME);
	for (i=0; i<2; i++) {
		if (!strcmp(word[2+i],PROPYESTXT))
			last_frac->tip_prop[i] = PROPYES;
		else if (!strcmp(word[2+i],ISECTTXT))
			last_frac->tip_prop[i] = ISECT;
		else
			last_frac->tip_prop[i] = PROPNO;
	} /*for*/
 
	for (;;) {
		/* get line, exit loop on EOF */
		if ((numwords = getwords(ifp,line,MAXLINE,word,MAXWORDS)) == EOF)
			return(0);
 
		/* skip comment lines */
		else if (word[0][0] == COMMENT_CHAR)
			continue;  /* continues for loop */

 		/* break with blank lines */
		else if (numwords == 0) break;

		else {
			if (numwords != 10) {
				fprintf(stderr, "Error: fracture %s, incomplete segment\n",
					last_frac->name);
				return -1;
			}
			num  = atoi(word[0]);
			xbeg = atof(word[1]);
			ybeg = atof(word[2]);
			xend = atof(word[3]);
			yend = atof(word[4]);
			kode = atoi(word[5]);
			bvshS = atof(word[6]);
			bvnoS = atof(word[7]);
			bvshM = atof(word[8]);
			bvnoM = atof(word[9]);
			add_BEs(num,xbeg,xend,ybeg,yend,kode,bvshS,bvnoS,bvshM,bvnoM);
 
			if (last_frac->last_seg == NULL)
				last_frac->first_seg = last_frac->last_seg =
					(struct frac_seg *)
					malloc(sizeof(struct frac_seg));
			else
				last_frac->last_seg = last_frac->last_seg->next =
					(struct frac_seg *)
					malloc(sizeof(struct frac_seg));
			last_frac->last_seg->start_elt = numbe - num + 1;
			last_frac->last_seg->end_elt   = numbe;
			last_frac->last_seg->next      = NULL;
		} /*else*/
	} /*for*/

	return 0;  /* everything read OK */
}
 
/*
 * This routine is called after the first line of a fault (with the "fault"
 * keyword) is read.  The input file, a line buffer used to read into, an array of
 * words on that line, and the number of words read are returned.  This routine
 * makes sure the right stuff was read on the line.
 * This routine returns 0 if everything is read OK, -1 if there was a read error.\
 */
int
readFault(FILE *ifp, char *line, char **word, int numwords)
{

	int num, i;
	double xbeg, ybeg, xend, yend;
 
	numFaults++;

	/* Allocate memory and define structures for another fault line */
	if (last_fault == NULL)
		first_fault = last_fault = (struct fault *)
			malloc(sizeof(struct fault));
	else
		last_fault = last_fault->next = (struct fault *)
			malloc(sizeof(struct fault));
	last_fault->next = NULL;
	last_fault->first_seg = NULL;
	last_fault->last_seg  = NULL;
	
	/* Check if enough parameters in first line of fault data */
	if (numwords != 5) {
		fprintf(stderr,"\nerror: Missing fault name or propagation ");
		fprintf(stderr,"flags in input file\n");
		return (-1);
	} /*if*/

	strncpy(last_fault->name,word[1],MAXFNAME);
	if (!strcmp(word[2],PROPYESTXT))
			last_fault->tails= TRUE;
		else
			last_fault->tails = FALSE;
	for (i=0; i<2; i++) {
		if (!strcmp(word[3+i],PROPYESTXT))
			last_fault->end_prop[i] = PROPYES;
		else
			last_fault->end_prop[i] = PROPNO;
	} /*for*/

	for (;;) {
		/* get line, exit loop on EOF */
		if ((numwords = getwords(ifp,line,MAXLINE,word,MAXWORDS)) == EOF)
			return(0);
 
		/* skip comment lines */
		else if ( word[0][0] == COMMENT_CHAR)
			continue;  /* continues for loop */
 
		/* stop reading if gets to a blank line */
		else if (numwords == 0 ) break ;
		
		else {
			if (numwords != 13) {
				fprintf(stderr, "Error: Fault %s, incomplete segment\n",
					last_fault->name);
				return -1;
			}
			num  = atoi(word[0]);
			xbeg = atof(word[1]);
			ybeg = atof(word[2]);
			xend = atof(word[3]);
			yend = atof(word[4]);
		
			for(i=numFaultElements + 1; i<= numFaultElements + num +1; i++){
				stiffS[i] = atof(word[5]);
				stiffN[i] = atof(word[6]);
				faultTen[i] = atof(word[7]);
				strength[i] = atof(word[8]);
				cohesion[i] = atof(word[9]);
				fric_stat[i] = atof(word[10]);
				friction[i] = atof(word[10]);
				fric_dyn[i] = atof(word[11]);
				crit_slip[i] = atof(word[12]);
				slipped[i] = FALSE;
			}/*for*/

			add_BEs(num,xbeg,xend,ybeg,yend,1,0,0,0,0);
			
			numFaultElements += num;
 
			if (last_fault->last_seg == NULL)
				last_fault->first_seg = last_fault->last_seg =
					(struct fault_seg *)
					malloc(sizeof(struct fault_seg));
			else
				last_fault->last_seg = last_fault->last_seg->next =
					(struct fault_seg *)
					malloc(sizeof(struct fault_seg));
			last_fault->last_seg->start_elt = numbe - num + 1;
			last_fault->last_seg->end_elt   = numbe;
			last_fault->last_seg->next      = NULL;
		} /*else*/
	} /*for*/

	return 0;  /* everything read OK */
}
 
/************************** function: open_files ***************************
*
****************************************************************************/
int
open_files(char *infile, FILE **ifp, char *outfile,
	   FILE **ofp, char *resumefile, FILE **rfp, char *graphicsfile, FILE **gfp)
{
	/* Open the input file */
	if ((*ifp = fopen(infile,"r")) == NULL) {
		fprintf(stderr,"\nerror: cannot open the file %s\n",infile);
		return(-1);
	} /*if*/
 
	/* If output file specified, open the file.  If output file */
	/* requested is stdout, set ofp = stdout.					*/
	if (outfile[0] != '\0')  {
		if (!strcmp(outfile,"stdout"))
			*ofp = stdout;
		else if ((*ofp = fopen(outfile,"w")) == NULL) {
			fprintf(stderr,"\nerror: cannot open the file %s\n",outfile);
			return(-1);
		} /*else if*/
	} /*if*/
 
	/* If resume file specified, open the file */
	if (resumefile[0] != '\0') {
		if ((*rfp = fopen(resumefile,"w")) == NULL) {
			fprintf(stderr,"\nerror: cannot open the file %s\n",resumefile);
			return(-1);
		} /*if*/
	} /*if*/
	/* If graphics file specified, open the file */
	if (graphicsfile[0] != '\0') {
		if ((*gfp = fopen(graphicsfile,"w")) == NULL) {
			fprintf(stderr,"\nerror: cannot open the file %s\n",graphicsfile);
			return(-1);
		} /*if*/
	} /*if*/
	return 0;
}
 
/********************** function: append_resumefile ***********************
*
***************************************************************************/
void append_resumefile(FILE *rfp)
 
{
	struct frac *current_frac;
	struct frac_seg *current_seg;
	int i, j;
	double xbeg, xend, ybeg, yend;
 
	current_frac = first_frac;
	for (i=1; i<=nfracs; i++) {
 
		current_seg = current_frac->first_seg;
		fprintf(rfp,"\nfracture %s    ",current_frac->name);
 
		for (j=0; j<2; j++) {
			if (current_frac->tip_prop[j] == PROPYES)
				fprintf(rfp,"%s   ",PROPYESTXT);
			else if (current_frac->tip_prop[j] == ISECT)
				fprintf(rfp,"%s   ",ISECTTXT);
			else
				fprintf(rfp,"%s   ",PROPNOTXT);
		} /*for*/
		fprintf(rfp,"\n");
 
		for (;;) {
			for (j = current_seg->start_elt; j <= current_seg->end_elt; j++) {
				xbeg = xm[j] - a[j]*cosbet[j];
				ybeg = ym[j] - a[j]*sinbet[j];
				xend = xm[j] + a[j]*cosbet[j];
				yend = ym[j] + a[j]*sinbet[j];
				fprintf(rfp,"1  %12f  %12f  %12f  %12f  %1d  %12f  %12f\n",
					xbeg, ybeg, xend, yend, kod[j], stepbvs[j], stepbvn[j]);
			} /*for*/
 
			if ((current_seg->next) != NULL)
				current_seg = current_seg->next;
			else
				break;
		} /*for*/
 
		current_frac = current_frac->next;
	} /*for*/
}
	
			
/************************ function print_data *****************************
*
***************************************************************************/
void print_data(FILE *ofp, int iter, int bdata, int fracData, int faultData, int odata, int step)
{
	struct frac *current_frac;
	struct frac_seg *current_seg;
	struct fault *current_fault;
	struct fault_seg *current_fault_seg;
	struct oline *current_oline;
	int i, j;
	double xbeg, ybeg, xend, yend;
	int numpb;
	int titles;

	fprintf(ofp,"DATA FROM CRACK GROWTH ITERATION #%d\n",iter);
	fprintf(ofp,"======================================\n\n");
 
	if (bdata) {
		fprintf(ofp,"\n--------------\n");
		fprintf(ofp,"BOUNDARY LINES\n");
		fprintf(ofp,"--------------\n");
		if (nbBEs > 0) {
			if (step == 1) print_BE_data(ofp, 1, nbBEs, TRUE);
			fprintf(ofp,"\n");
			print_BE_disps_stresses(ofp, 1, nbBEs, TRUE);
			fprintf(ofp,"\n");
		}
		else
			fprintf(ofp,"There are no boundary lines in this problem.\n\n");
	} /*if*/
 
	if (fracData) {
		fprintf(ofp,"\n--------------\n");
		fprintf(ofp,"FRACTURES     \n");
		fprintf(ofp,"--------------\n\n");
		current_frac = first_frac;
		if (nfracs > 0) {
			for (i=1; i<=nfracs; i++) {
				fprintf(ofp,"FRACTURE: %s\n\n",current_frac->name);
				for (j=1; j<=2; j++) {
					titles = TRUE;
					current_seg = current_frac->first_seg;
					for (;;) {
	
						if (j==1)
							print_BE_data(ofp, current_seg->start_elt,
								current_seg->end_elt, titles);
						if (j==2)
							print_BE_disps_stresses(ofp, current_seg->start_elt,
								current_seg->end_elt, titles);
	
						titles = FALSE;
	
						if (current_seg->next != NULL)
							current_seg = current_seg->next;
						else
							break;
					} /*for*/
					fprintf(ofp,"\n");
				} /*for*/
 
				current_frac = current_frac->next;
			} /*for*/
		} /*if*/
		else
			fprintf(ofp,"There are no fractures in this problem.\n\n");
	} /*if*/
 
	if (faultData) {
		fprintf(ofp,"\n--------------\n");
		fprintf(ofp,"FAULTS     \n");
		fprintf(ofp,"--------------\n\n");
		current_fault = first_fault;
		if(numFaults >0) {
			for (i=1; i<=numFaults; i++) {
				fprintf(ofp,"FAULT: %s\n\n",current_fault->name);
				for (j=1; j<=2; j++) {
					titles = TRUE;
					current_fault_seg = current_fault->first_seg;
					for (;;) {
	
						if (j==1 && step==1)
							print_fault_data(ofp, current_fault_seg->start_elt,
								current_fault_seg->end_elt, titles);
						if (j==2)
							print_BE_disps_stresses(ofp, current_fault_seg->start_elt,
								current_fault_seg->end_elt, titles);
	
						titles = FALSE;
	
						if (current_fault_seg->next != NULL)
							current_fault_seg = current_fault_seg->next;
						else
							break;
					} /*for*/
					fprintf(ofp,"\n");
				} /*for*/
 
				current_fault = current_fault->next;
			} /*for*/
		} /*if*/
		else
			fprintf(ofp,"There are no faults in this problem.\n\n");
	} /*if*/

				
	if (odata) {
		fprintf(ofp,"\n-----------------\n");
		fprintf(ofp,"OBSERVATION LINES\n");
		fprintf(ofp,"-----------------\n\n");
 
		titles = TRUE;
		current_oline = first_oline;
		if (nolines > 0) {
			for (i=1; i<=nolines; i++) {
				xbeg  = current_oline->xbeg;
				ybeg  = current_oline->ybeg;
				xend  = current_oline->xend;
				yend  = current_oline->yend;
				numpb = current_oline->numpb;
				print_oline_data(ofp, i, xbeg, ybeg, xend, yend, numpb, titles);
				titles = FALSE;
 
				current_oline = current_oline->next;
			} /*for*/
			fprintf(ofp,"\n");
		} /*if*/
		else
			fprintf(ofp,"There are no observation lines in this problem.\n\n");
	} /*if*/
}
			
 
/************************ function: print_oline_data *************************
*
* Compute displacements and stresses at specified points in body
*
*****************************************************************************/
void print_oline_data(FILE *ofp, int linenum, double xbeg, double ybeg,
		      double xend, double yend, int numpb, int titles)
{
	int i;
	int nump;
	double delx, dely;
	double xp, yp, lithv, lithh;
	int exit;
	double ux, uy, sigxx, sigyy, sigxy;
 
/* I have temporarily removed the line numbers and the row of dashes
so that Spyglass can more easily read the data in. */

	if (titles) {
		fprintf(ofp,"DISPLACEMENTS AND STRESSES AT SPECIFIED POINTS ");
		fprintf(ofp,"IN THE BODY\n\n");
		fprintf(ofp,"   X-coord      Y-coord         UX          UY");
		fprintf(ofp,"        Sigma-XX     Sigma-YY      Sigma-XY");
		
		if (gravity) {
                        fprintf(ofp,"    Sigma-XX(w/ grav)   Sigma-YY(w/ grav)\n");
                }else{
                        fprintf(ofp,"\n");
                }/*else if*/
	} /*if*/
	
	nump = numpb+1;
	delx = (xend-xbeg)/nump;
	dely = (yend-ybeg)/nump;
	if (numpb > 0)
		nump++;
	if ((delx*delx + dely*dely) == 0)
		nump = 1;
	for (i=1; i<=nump; i++) {
		xp = xbeg + (i-1)*delx;
		yp = ybeg + (i-1)*dely;
		exit = disps_stresses_at_pt(xp,yp,&ux,&uy,&sigxx,&sigyy,&sigxy,
			TRUE);
 

		if (exit != -1) {

                        fprintf(ofp,"%12.5f %12.5f %12.5e %12.5e %12.5e ",
                                xp, yp, ux, uy, sigxx);
                        fprintf(ofp,"%12.5e %12.5e ",
                                sigyy , sigxy );

			if (gravity) {
                                lithv = density * 9.8 * yp/1000000;
				lithh = kratio * lithv;
			
			fprintf(ofp," %16.5e  %16.5e\n", 
				 sigxx + lithh, sigyy + lithv );
                        }else{
				fprintf(ofp,"\n");
			}/*else if*/
		
		} /*if*/
	} /*for*/
}

/**********************************************************************/
/* This function initiates the graphics file */
void
print_graphics_start(FILE *gfp)
{
  fprintf(gfp, "#Inventor V2.0 ascii\n");
  fprintf(gfp, "Separator { \n");
  fprintf(gfp, " Font { size 16 name Helvetica }\n");
  fprintf(gfp, " Blinker { speed 0.2\n");

  print_graphics_frame(gfp, 0, 0);
}

/* This function prints the graphics to the file */
void
print_graphics_frame(FILE *gfp, int step, int increment)
{
  int i;
  fprintf(gfp, "  Separator {\n");
  fprintf(gfp, "   Text2 { string [ \"Loading Step %d\", \"Fracture Increment %d\" ] \n", 
		step, increment);
  fprintf(gfp, "           justification CENTER }\n");

  /* Boundary elements: */
  fprintf(gfp, "   BaseColor { rgb 0 0 1 }\n");
  fprintf(gfp, "   Coordinate3 { point [\n");
  for (i = 1; i <= nbBEs; i++) {
    float x1, y1, x2, y2;
    x1 = xm[i]-cosbet[i]*a[i];
    y1 = ym[i]-sinbet[i]*a[i];
    x2 = xm[i]+cosbet[i]*a[i];
    y2 = ym[i]+sinbet[i]*a[i];
    fprintf(gfp, "%f %f 0.0, %f %f 0.0,\n", x1, y1, x2, y2);
  }
  fprintf(gfp, "   ] }\n");

  fprintf(gfp, "   LineSet { numVertices [ ");
  for (i = 1; i <= nbBEs; i++) fprintf(gfp, "2,");
  fprintf(gfp, "]\n");
  fprintf(gfp, "   }\n");

  /* Faults */
  fprintf(gfp, "   BaseColor { rgb 0 1 0 }\n");
  fprintf(gfp, "   Coordinate3 { point [\n");
  for (i = nbBEs+1; i <= nbBEs+numFaultElements; i++) {
    float x1, y1, x2, y2;
    x1 = xm[i]-cosbet[i]*a[i];
    y1 = ym[i]-sinbet[i]*a[i];
    x2 = xm[i]+cosbet[i]*a[i];
    y2 = ym[i]+sinbet[i]*a[i];
    fprintf(gfp, "%f %f 0.0, %f %f 0.0,\n", x1, y1, x2, y2);
  }
  fprintf(gfp, "   ] }\n");

  fprintf(gfp, "   LineSet { numVertices [ ");
  for (i = nbBEs+1; i <=nbBEs+numFaultElements; i++) fprintf(gfp, "2,");
  fprintf(gfp, "]\n");
  fprintf(gfp, "   }\n");

 /* Fractures */
  fprintf(gfp, "   BaseColor { rgb 1 0 0 }\n");
  fprintf(gfp, "   Coordinate3 { point [\n");
  for (i = nbBEs+ numFaultElements+1; i <= numbe; i++) {
    float x1, y1, x2, y2;
    x1 = xm[i]-cosbet[i]*a[i];
    y1 = ym[i]-sinbet[i]*a[i];
    x2 = xm[i]+cosbet[i]*a[i];
    y2 = ym[i]+sinbet[i]*a[i];
    fprintf(gfp, "%f %f 0.0, %f %f 0.0,\n", x1, y1, x2, y2);
  }
  fprintf(gfp, "   ] }\n");

  fprintf(gfp, "   LineSet { numVertices [ ");
  for (i = nbBEs+ numFaultElements+1; i <= numbe; i++) fprintf(gfp, "2,");
  fprintf(gfp, "]\n");
  fprintf(gfp, "   }\n");


  fprintf(gfp, "  }\n");
}

void
print_graphics_end(FILE *gfp)
{
  fprintf(gfp, " }\n");
  fprintf(gfp, "}\n");
}
