/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   testvit.c
**      Purpose: driver for testing the Viterbi code.
**      Organization: University of Maryland
**
**	Update:
**	Author:	Tapas Kanungo
**	Purpose: run both viterbi with probabilities and
**		viterbi with log, change output etc.
**      $Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <glib.h> //added
#include "nrutil.h"
#include "hmm.h"
static char rcsid[] = "$Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $";

int main (int argc, char **argv)
{
	int 	t, T,T_new,repeat_len,prefix_idx;
	HMM  	hmm;
	int	*O;	/* observation sequence O[1..T] */
	int	*q;	/* state sequence q[1..T] */
	double **delta;
	int	**psi;
	int **traceback_dir; //added, for Viterbi modification
	double 	proba, logproba;
	FILE	*fp, *traceback_del, *traceback, *delta_out;
	int i; //added
	int n; //added
	int offset; //added
	char *S; //added, input sequence as chars
	//GList *list = NULL; //added, input to Viterbi modification
	GList *list_log = NULL;
	GList *likelihoods_ints = NULL;
	GList *likelihoods_dec = NULL;



	if (argc != 6) { //changed to 4 so we can include T and O separately, might need to allow this to be longer if it reads the array as separate arguments, changed to 5 to include repeat length, changed to 8 to include debug (remove before running through python script) 

		printf("%d",argc);
		printf("Usage error \n");
		printf("Usage: testvit <model.hmm> <obs.seq length> <obs.seq> \n");
		exit (1);
	}
	//open HMM file
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[1]);
		exit (1);
	}
	//read HMM file according to ReadHMM
	ReadHMM(fp, &hmm);
	fclose(fp);

	//read seq file THIS IS THE PART WE WANT TO EDIT TO ALLOW RUNNING ACROSS MORE THAN ONE SEQUENCE/READ
//	fp = fopen(argv[2], "r");
//	if (fp == NULL) {
//		fprintf(stderr, "Error: File %s not found\n", argv[2]);
//		exit (1);
//	}

	//ReadSequence(fp, &T, &O); //I have to change this method to change how the program reads the sequences
	//fclose(fp); //I don't think I'll need this after I change it

	T = atoi(argv[2]); //read T from input
	S = argv[3];
	//printf("%s",S);
	repeat_len =  atoi(argv[4]); // for calculation of RDn case
	prefix_idx = atoi(argv[5]); // for calculation of RDn case
	//printf("%s",S);
	//printf("%d\n",repeat_len);
	//scanf(fp, "T= %d\n", pT);
	O = ivector(1,T);
	//sscanf(S,"%d", &i);
	//printf("%d",i);
	i=1;
	while (sscanf(S, "%d%n", &n, &offset)==1){ //added this to read string of sequence
		O[i] = n;
		S += offset;
		//printf("%d\n",O[i]);
		i++;
	}




	q = ivector(1,T);

	traceback_dir = imatrix(1, T, 1, hmm.N); //added
	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);

	//printf("------------------------------------\n");
	//printf("Viterbi using direct probabilities\n");
	//Viterbi(&hmm, T, O, delta, psi, traceback_dir,&list, &T_new,q, &proba);
	//fprintf(stdout, "Viterbi  MLE log prob = %E\n", log(proba));
	//fprintf(stdout, "Optimal state sequence:\n");
	//PrintSequence(stdout, T, q);
	//printf("%d",T_new);
	//printf("%d\n", g_list_length(list));
	//PrintPath(stdout, T_new, list); //function I am writing to read back the linked list
	//printf("------------------------------------\n");
	//printf("Viterbi using log probabilities\n");
	/* note: ViterbiLog() returns back with log(A[i][j]) instead
	** of leaving the A matrix alone. If you need the original A,
	** you can make a copy of hmm by calling CopyHMM */
	ViterbiLog(&hmm, T, O, delta, psi, traceback_dir,&list_log,&likelihoods_ints,&likelihoods_dec,&T_new, q, &logproba, repeat_len, prefix_idx);

	fprintf(stdout, "Viterbi  MLE log prob = %E\n", logproba);
	fprintf(stdout, "Optimal state sequence:\n");
	//PrintSequence(stdout, T, q);
	PrintPath(stdout, T_new, list_log);
	PrintLikelihoods(stdout, T_new, likelihoods_ints);
	PrintLikelihoods(stdout, T_new, likelihoods_dec);

	//for debuggings: *traceback_del, *traceback, *delta_out
	// traceback_del = fopen(argv[5], "w");
	// traceback = fopen(argv[6], "w");
	// delta_out = fopen(argv[7], "w");
	// PrintTraceBack(traceback_del,traceback_dir,&hmm,T); //print traceback matrix for deletions
	// PrintTraceBack(traceback,psi,&hmm,T); // print out total traceback matrices
	// PrintDelta(delta_out,delta,&hmm,T);
	//PrintAllOuts(stdout, T_new,list_log,delta,hmm.N, T);
	printf("------------------------------------\n");
	printf("The two log probabilites and optimal state sequences\n");
	printf("should identical (within numerical precision). \n");

	free_ivector(q, 1, T);
	free_ivector(O, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_imatrix(traceback_dir,1,T,1,hmm.N); //added
	g_list_free(list_log);
	g_list_free(likelihoods_ints);
	g_list_free(likelihoods_dec);
	free_dmatrix(delta, 1, T, 1, hmm.N);
	FreeHMM(&hmm);
}
