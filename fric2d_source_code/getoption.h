/******************************* getopt.h *********************************
*
* Data declarations for use with getopt.c
* NOTE: the first define of getopt_h_ guarantees that things will only
*       be declared once.
*
***************************************************************************/

#ifndef getopt_h_
#define getopt_h_

#define FILE_ARG 0
#define FLAG 1
#define FLAG_AND_FILE 2
#define NO_MORE_ARGS -1
#define NO_SUCH_ARG -2

struct arg_rec {
	char flag_name;
	char *file_arg;
	char arg_type;
	struct arg_rec *next_rec;
};

#ifdef ANSI
int getoption(char *arg_string, int argc, char *argv[]);
#endif

extern char *global_argument;
#endif
