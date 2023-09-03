#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

//#define DEBUG 1

long long med_seq = 0;
long median = 0;

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
	printf("%li\t%li\n", el->length, el->count);
  }
}

void find_median(const void *ptr, VISIT order, __attribute__((unused)) int level) {
  if (!median && (order == postorder || order == leaf)) {
	const length_el *el = *(const length_el **)ptr;
	if ((med_seq -= el->count) <= 0)
	  median = el->length;
  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
	fprintf(stderr, "usage: %s FASTQ_FILE\n", argv[0]);
	exit(EINVAL);
  }
  FILE *fq = fopen(argv[1], "r");
  if (fq == NULL)
	exit(errno);

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

#ifdef DEBUG
	printf("sequence length: %li\n", seq_len);
#endif
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
  } while (!feof(fq) && fseek(fq, seq_len + newlines, SEEK_CUR) == 0);

  if (ferror(fq)) {
	fprintf(stderr, "error %i occurred while parsing %s", ferror(fq), argv[1]);
	exit(1);
  }
  fclose(fq);

 //twalk(seq_lengths, print_node);
  med_seq = total_seqs / 2;
  twalk(seq_lengths, find_median);
  printf("%li", median);
  return 0;
}