#include <errno.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERSION "0.1.2"
//#define DEBUG 1
void usage(const char* cmd) {
  	fprintf(stderr, "https://github.com/yttria-aniseia/fastq-lengths (v%s)\n\
usage: %s [subcommand [stopafter]] FASTQ_FILE\n\
\tsubcommand: lengths | median | count | summary\n\
\tstopafter:  only judge this many entries\n\
\n\
fastq-lengths lengths 100 example.fq\n\
32  1\n\
150 99\n", VERSION, cmd);
}
enum subcommand { LENGTHS, MEDIAN, SUMMARY, COUNT };

long long med_seq = 0;
long median = 0;
long num_nodes = 0;

typedef struct length_el {
  long length;
  long count;
} length_el;

int node_compare(const void *node1, const void *node2) {
  long difference = ((length_el *)node2)->length - ((length_el *)node1)->length;
  return (difference > 0) - (difference < 0); // sgn (avoid long to int cast)
}
void print_node(const void *ptr, VISIT order, __attribute__((unused)) int level) {
  if (order == postorder || order == leaf) {
	const length_el *el = *(const length_el **)ptr;
	num_nodes++;
	printf("%li\t%li\n", el->length, el->count);
  }
}
void find_median(const void *ptr, VISIT order, __attribute__((unused)) int level) {
  if (!median && (order == postorder || order == leaf)) {
	const length_el *el = *(const length_el **)ptr;
	num_nodes++;
	if ((med_seq -= el->count) <= 0)
	  median = el->length;
  }
}

int main(int argc, char** argv) {
  long long stop_after = -1;
  enum subcommand cmd = LENGTHS;
  switch (argc) {
  case 4: // subcommand firstN fastqfile
	if (sscanf(argv[2], "%lli", &stop_after) != 1) {
	  usage(argv[0]); exit(EINVAL);
	}
	__attribute__((fallthrough));
  case 3: // subcommand fastqfile
	if (strncmp("median", argv[1], 3) == 0)
	  cmd = MEDIAN;
	else if (strncmp("count", argv[1], 3) == 0)
	  cmd = COUNT;
	else if (strncmp("summary", argv[1], 3) == 0)
	  cmd = SUMMARY;
	else // omitting 'lengths' when specifying <stopafter> is not supported.
	  cmd = LENGTHS;
	break;
  case 2: // fastqfile
	break;
  default:
	usage(argv[0]); exit(EINVAL);
  }
  FILE *fq = fopen(argv[argc - 1], "r");
  if (fq == NULL) {
	fprintf(stderr, "%s does not exist or cannot be opened.\n", argv[argc - 1]);
	usage(argv[0]);
	exit(errno);
  }

  void *seq_lengths = NULL;
  length_el *element_ptr;
  void *node = NULL;

  long long total_seqs = 0;
  long seq_len = 0;
  long last_seq_len = -1;
  long newlines = 0;
  char c;
  do {
	seq_len = 0; newlines = 0;
	// @<seqname>
	if (fgetc(fq) != '@') {
	  if (feof(fq)) break;
	  fprintf(stderr, "error: expected @<seqname> line (pos %li): %s\n",
			  ftell(fq), argv[1]);
	  exit(1);
	}
	while ((c = (char)fgetc(fq)) != EOF && c != '\n');
	// <seq>
	while ((c = (char)fgetc(fq)) != EOF && c != '+')
	  if (c == '\n') newlines++; else seq_len++;

	// Update Counts
	total_seqs++;
	if (seq_len != last_seq_len) {
	  element_ptr = malloc(sizeof(length_el));
	  element_ptr->length = seq_len;
	  element_ptr->count = 1;
	  node = tsearch((void *)element_ptr, &seq_lengths, node_compare);
	  if (node == NULL) {
		fprintf(stderr, "out of memory while collecting sequence lengths\n");
		exit(ENOMEM);
	  } else if (*(length_el **)node != element_ptr) {
		(*(length_el **)node)->count++;
		free(element_ptr);
	  }
	  last_seq_len = seq_len;
	} else {
	  (*(length_el **)node)->count++;
	}

	// +[<seqname>] (+ already consumed)
	while ((c = (char)fgetc(fq)) != EOF && c != '\n');
	// <qual>
  } while (!feof(fq) && (total_seqs != stop_after) &&
		   fseek(fq, seq_len + newlines, SEEK_CUR) == 0);

  if (ferror(fq)) {
	fprintf(stderr, "error %i occurred while parsing %s", ferror(fq), argv[1]);
	exit(1);
  }
  fclose(fq);

  switch (cmd) {
  case LENGTHS:
	twalk(seq_lengths, print_node);
	break;
  case MEDIAN:
	med_seq = total_seqs / 2;
	twalk(seq_lengths, find_median);
	printf("%li\n", median);
	break;
  case COUNT: // COUNT
	printf("%lli\n", total_seqs);
	break;
  case SUMMARY:
	med_seq = total_seqs / 2;
	twalk(seq_lengths, find_median);
	printf("%li\t%li\t%lli\n", median, num_nodes, total_seqs);
	break;
  }
  return 0;
}
