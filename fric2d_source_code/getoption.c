#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "getoption.h"
#if !defined(__APPLE__)                         /*on macs malloc is superceded by stdlib */
   #include <malloc.h>
#endif

struct arg_rec *read_args(int argc, char *argv[], char *arg_string);
static struct arg_rec *add_arg(	struct arg_rec	**top,
				int				i,
				char			arg_type,
				char			*argv[]);

/************************** function: getoption ******************************
*
* Creates a linked list of arguments to the program based on the arg_string
* passed in.  If an argument has a ':' after it the next argv argument belongs
* with that flag.  It tags the structures as 4 types: FILE_ARG, meaning that
* there is no flag that goes with it; FLAG, a -x argument, with no file that
* goes with it; a FLAG_AND_FILE, which has a -x argument followed by a file
* argument; and finally an error argument NO_SUCH_ARG.  The structure of an
* argument list is defined in getoption.h as follows:
* 
*		struct arg_rec{
*			char flag_name;
*			char *file_arg;
*			char arg_type;
*			struct arg_rec *next_rec;
*		};
*
* By: Dan Kolkowitz, Nov. 1990.
*     Written as part of a solution set for CS040, Stanford University.
* 
****************************************************************************/
#ifdef ANSI
int getoption(char *arg_string, int argc, char *argv[])
#else
int getoption(arg_string, argc, argv)
char *arg_string;
int argc;
char *argv[];
#endif

{
	/**** automatic variables ****/
	static int args_parsed=0;
	static struct arg_rec *next_arg, *current_arg;

	/**** function body ****/

	/*****************************************************************
	* The first time getoption is called, call read_args() to set up
	* the linked list of arguments. Also set the next_arg pointer.
	*****************************************************************/
	if (!args_parsed)
		if ((next_arg = read_args(argc,argv,arg_string)) 
			== (struct arg_rec *) NO_SUCH_ARG)
			return (NO_SUCH_ARG);
	args_parsed++;

	if (!next_arg)
		return (NO_MORE_ARGS); 

	current_arg = next_arg;
	next_arg = next_arg->next_rec;

	switch (current_arg->arg_type) {
		case FILE_ARG:
			global_argument = current_arg->file_arg;
			return (FILE_ARG);
		case FLAG:
			global_argument = NULL;
			return (current_arg->flag_name);
		case FLAG_AND_FILE:
			global_argument = current_arg->file_arg;
			return (current_arg->flag_name);
		case NO_SUCH_ARG:
			return (NO_SUCH_ARG);
	} /*switch*/
}



/************************** function: read_args ************************
* 
* Parses each argument and determines what type to add.
* 	Calls add_arg() to allocate the actual structure and set the
* 	actual field types.
* 
************************************************************************/
#ifdef ANSI
struct arg_rec *read_args(int argc, char *argv[], char *arg_string)
#else
struct arg_rec *read_args(argc, argv, arg_string)
int argc;
char *argv[];
char *arg_string;
#endif

{
	/**** automatic variables ****/
	int i=1;
	char *c_ptr, newarg;
	static struct arg_rec *top;
	
	/**** function body ****/
	while (i < argc) {
		if (argv[i][0] == '-') {
			newarg = argv[i][1];
			c_ptr = (char *) strchr(arg_string, newarg);
			/* c_ptr = (char *) index(arg_string, newarg);  */
			if (!c_ptr)
				return ((struct arg_rec *) NO_SUCH_ARG);
			if (*(c_ptr+1) == ':') {
				add_arg(&top, i, FLAG_AND_FILE, argv);
				i+=2;
			} else add_arg(&top, i++, FLAG, argv);
		} else add_arg(&top, i++, FILE_ARG, argv);
	}
	return (top);
}

/**************************** function: add_arg **************************
* 
* Allocates an actual node.
* 	Mallocs the space for the structure.
* 	Sets a pointer to the argument if it is required.
* 
**************************************************************************/
#ifdef ANSI
static struct arg_rec *add_arg(	struct arg_rec	**top,
								int				i,
								char			arg_type,
								char			*argv[])
#else
static struct arg_rec *add_arg(top, i, arg_type, argv)
struct arg_rec **top;
int i;
char arg_type;
char *argv[];
#endif

{
	/**** automatic varibales ****/
	static struct arg_rec *last_arg;
	struct arg_rec *test_var;
	int *test_var2;

	/**** function body ****/
	
	/****************************************************************
	* The first time add_arg is called, top points nowhere, and the
	* the first node in the linked list of arguments must be created.
	* 
	* Allocate the memory needed for the new argument record and
	* reset the last_arg pointer.
	****************************************************************/
	if (!*top) {
		*top = (struct arg_rec *) malloc(sizeof(struct arg_rec));
		test_var = (struct arg_rec *) malloc(sizeof(struct arg_rec));
		test_var2 = (int *) malloc(sizeof(int));
		last_arg = *top;
	} /*if*/
	else {
		last_arg->next_rec = (struct arg_rec *) malloc(sizeof(struct arg_rec));
		last_arg = last_arg->next_rec;
	} /*else*/

	/**************************************************************
	* Put the appropriate info into the new (last) argument record
	**************************************************************/
	last_arg->next_rec = NULL;
	last_arg->arg_type = arg_type;

	switch (arg_type) {
		case FILE_ARG:
			last_arg->file_arg = argv[i];
			break;
		case FLAG_AND_FILE:
			last_arg->file_arg = argv[i+1];
		case FLAG:
			last_arg->flag_name = argv[i][1];
			break;
	} /*switch*/
	return last_arg;
}

