#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <ctype.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
#define MAXWORDS 20 

/**************************** function: getwords ***********************
*
* Reads a line of input into the character array line from the file infp
* using fgets.  Copies line to the static array wline.  Then sets the 
* character pointers word[0]..word[n-1] to point to the first letters
* of the 1st..nth words in line.  Changes the space following each word 
* in line to '\0', so that each word[i] will be a null terminated, 
* one-word string.
*
* Returns the number of words on the line or EOF on end of file.
*
* By: Andrew L. Thomas (February, 1991)
*
************************************************************************/
#ifdef ANSI
int getwords(	FILE	*infp,
				char	*line,
				int		maxline,
				char	*word[],
				int		maxwords)
#else
int getwords(infp, line, maxline, word, maxwords)
FILE *infp;
char *line, *word[];
int maxline, maxwords;
#endif

{
	int i, j;
	int inword = FALSE;
	int inquotes = FALSE;
	extern int debug;
	static char *wline = NULL;
	static int oldmaxline;


	j = 0;										/* j = index for word[]		*/

	/*------------------------------------------
	free/allocate memory for wline if necessary
	-------------------------------------------*/
	if ((wline == NULL) || (maxline > oldmaxline)) {
		if (wline != NULL)
			free(wline);
		wline = (char *) malloc((size_t) maxline);
		oldmaxline = maxline;
	}

	/*------------------------------------
	read line from file and copy to wline
	-------------------------------------*/
	if (fgets(line, maxline, infp) == NULL) 
		return( EOF );
	strcpy(wline,line);

	/*--------------------------------------
	Step through line and set word pointers
	---------------------------------------*/
	for (i = 0; wline[i] != '\0'; i++) {
		if (wline[i]=='"') {
			if (!inquotes) {
				if (inword) {
					word[j] = &wline[i+1];
					inword = FALSE;
				} /*if*/
				else
					word[j++] = &wline[i+1];
				inquotes = TRUE;
			} /*if*/
			else {
				wline[i] = '\0';
				inquotes = FALSE;
			} /*else*/
		} /*if*/
		else if (!inquotes && !inword && !isspace(wline[i])) {
			if(j< MAXWORDS){
				word[j++] = &wline[i];
				inword = TRUE;
			}else{
				fprintf(stderr,"\nToo many words within a line of the input file\n");
				return (-1);
			}
		} /*if*/
		else if (inword && isspace(wline[i])) {
			wline[i] = '\0';
			inword = FALSE;	
		} /*else if*/
	} /*for*/

	return j;
}
