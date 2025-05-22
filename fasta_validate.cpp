/* 
 * Validate a fasta file
 *
 * Rob Edwards January 2019.
 *
 * We exit with:
 *	0: this is a valid fasta file
 *	1: the first line does not start with a >
 *	2: the ids are not unique
 *	4: lines in the sequence (that do not start >) contain characters that do not match the perl regexp /[A-Z][a-z] /
 *
 *	Over 200:
 *	internal errors, eg. unable to allocate memory, etc.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <search.h>
#include <zlib.h>

#include "fasta_validate.h"


int contains_non_word_characters(char *seq, int verbose) {
	/*
	 * we consider word characters to be ASCII:
	 *      64 < char < 91   (A=65; Z=90)
	 *      96 < char < 123  (a=97; z=122)
	 */

	
	if (seq == NULL) {
		fprintf(stderr, "Empty line receieved. Empty string?\n");
		return 8;
	}

	for (int i=0; i<strlen(seq); i++) {
		if ((int) seq[i] < 65) {
			if (((int) seq[i] != 10) && ((int) seq[i] != 13))
				return 1;
		}
		else if (((int) seq[i] > 90) && ((int) seq[i] < 97))
			return 2;
		else if ((int) seq[i] > 122)
			return 3;
	}
	return 0;
}

void help(char *nm) {
	// printf("%s version %d\n\nCheck and validate a fasta file\n", nm, VERSION);
	printf("Rob Edwards.\n\n");
	printf("%s checks your fasta file and exits with:\n", nm);
	printf("\t0: this is a valid fasta file\n");
	printf("\t1: the first line does not start with a >\n");
	printf("\t2: the ids are not unique\n");
	printf("\t4: lines in the sequence (that do not start >) contain characters that do not match the perl regexp /[A-Z][a-z]/\n");
	printf("\t8: There is a sequence with zero length in it\n\n");
	printf("\tOver 200:\n");
	printf("\t\tinternal errors, eg. unable to allocate memory, etc.\n\n");
	printf("Usage: %s [options] [fasta file]\n\t-v: verbose output\n\t-V: current version\n\t-h: this help\n", nm);
}


int fasta_validator_main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("%s [-h] [-v] [-V] [fasta file]\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	int verbose = 0;

	if (strcmp(filename, "-h") == 0) {
		help(argv[0]);
		return 0;
	}

	// if (strcmp(filename, "-V") == 0) {
	// 	printf("%s version %d.\nFor more information please see https://github.com/linsalrob/fasta_validator\n", argv[0], VERSION);
	// 	return 0;
	// }


	if (strcmp(filename, "-v") == 0) {
		verbose = 1;
		filename = argv[2];
	}

	return validate_fasta_core(filename, verbose);
}


int validate_fasta_file(const char *filename, int verbose) {
	if (strlen(filename) >= 3 && strcmp(filename + strlen(filename) - 3, ".gz") == 0) {
    	return validate_fasta_core_gzip(filename, verbose);
	}
    return validate_fasta_core(filename, verbose);
}


int validate_fasta_core(const char *filename, int verbose) {
    FILE *fp;
    char line[MAXLINELEN];
    int firstline = 1;
    int seqcount = 0;

    int hc = hcreate(NUMSEQS);
    if (hc == 0) {
        fprintf(stderr, "Unable to create the hash table\n");
        return -1;
    }

    if ((fp = fopen(filename, "r")) == NULL) {
        if (verbose)
            fprintf(stderr, "Can't open file %s\n", filename);
        return 1;
    }

    while ((fgets(line, MAXLINELEN, fp)) != NULL) {
        if (line[0] == '>') {
            if (!firstline && seqcount == 0) {
                if (verbose)
                    fprintf(stderr, "ERROR: We have an empty sequence\n");
                return 8;
            }
            firstline = 0;
            seqcount = 0;

            char *p = strchr(line, ' ');
            if (p) *p = '\0';

            ENTRY item;
            item.key = strdup(line);
            ENTRY *found_item;
            if ((found_item = hsearch(item, FIND)) != NULL) {
                if (verbose) {
                    fprintf(stderr, "ERROR: Found a duplicate id: |%s|\n", line);
                    fprintf(stderr, "ERROR: Found a duplicate id: |%s|\n", found_item->key);
                }
                return 2;
            }
            (void) hsearch(item, ENTER);
        } else {
            if (firstline > 0) {
                if (verbose)
                    fprintf(stderr, "ERROR: The first line should start with a >\n");
                return 1;
            }
            if (contains_non_word_characters(line, verbose)) {
                if (verbose)
                    fprintf(stderr, "ERROR: We have a non word character!\n");
                return 4;
            }
            seqcount += strlen(line);
        }
    }

    if (seqcount == 0) {
        if (verbose)
            fprintf(stderr, "ERROR: At end: We have an empty sequence\n");
        return 8;
    }

    if (fp != NULL) {
        fclose(fp);
    }

    hdestroy();
    return 0;
}

int validate_fasta_core_gzip(const char *filename, int verbose) {
    gzFile fp;
    char line[MAXLINELEN];
    int firstline = 1;
    int seqcount = 0;

    if (hcreate(NUMSEQS) == 0) {
        fprintf(stderr, "Unable to create the hash table\n");
        return -1;
    }

    fp = gzopen(filename, "rb");
    if (fp == NULL) {
        if (verbose)
            fprintf(stderr, "Can't open file %s\n", filename);
        return 1;
    }

    while (gzgets(fp, line, MAXLINELEN) != NULL) {
        if (line[0] == '>') {
            if (!firstline && seqcount == 0) {
                if (verbose)
                    fprintf(stderr, "ERROR: We have an empty sequence\n");
                gzclose(fp);
                return 8;
            }
            firstline = 0;
            seqcount = 0;

            char *p = strchr(line, ' ');
            if (p) *p = '\0';

            ENTRY item;
            item.key = strdup(line);  // must be strdup for hsearch
            ENTRY *found_item;
            if ((found_item = hsearch(item, FIND)) != NULL) {
                if (verbose) {
                    fprintf(stderr, "ERROR: Found a duplicate id: |%s|\n", line);
                    fprintf(stderr, "ERROR: Found a duplicate id: |%s|\n", found_item->key);
                }
                gzclose(fp);
                return 2;
            }
            (void) hsearch(item, ENTER);
        } else {
            if (firstline > 0) {
                if (verbose)
                    fprintf(stderr, "ERROR: The first line should start with a >\n");
                gzclose(fp);
                return 1;
            }
            if (contains_non_word_characters(line, verbose)) {
                if (verbose)
                    fprintf(stderr, "ERROR: We have a non word character!\n");
                gzclose(fp);
                return 4;
            }
            seqcount += strlen(line);
        }
    }

    gzclose(fp);
    hdestroy();

    if (seqcount == 0) {
        if (verbose)
            fprintf(stderr, "ERROR: At end: We have an empty sequence\n");
        return 8;
    }
    

    return 0;
}