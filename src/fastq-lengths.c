#include <errno.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERSION "0.2.8"
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


long long page_ofs = -1;
inline long long _pos(char *buf, const size_t bufsize, const char *nextc) {
  if (page_ofs == -1)
	return 0;
  return page_ofs * (long long)bufsize + (nextc - buf);
}
inline size_t readmore(char *buf, char **nextc, char **end, const size_t bufsize,
					FILE *fp) {
    size_t count = fread(buf, 1, bufsize, fp);
    if (count <= 0)
      return count;
    page_ofs++;
    *end = buf + count;
    *nextc = buf;
	return count;
}
inline int nextchar(char *buf, char **nextc, char **end, const size_t
                    bufsize, FILE *fp) {
  if (*nextc == *end) {
    if (readmore(buf, nextc, end, bufsize, fp) <= 0)
	  return EOF;
  }
  char c = **nextc;
  (*nextc)++;
  return c;
}
inline int skip(long skipn, char *restrict buf, char **nextc, char **end, const size_t
                bufsize, FILE *restrict fp) {
  (*nextc) += skipn - 1;
  if (*nextc >= *end) {
    skipn = *nextc - *end;
    if (readmore(buf, nextc, end, bufsize, fp) <= 0)
	  return EOF;
    return skip(skipn + 1, buf, nextc, end, bufsize, fp);
  }
  char c = **nextc;
  (*nextc)++;
  return c;
}
inline long long scanuntil(char delim, char *restrict buf, char **nextc, char **end, const size_t
						bufsize, FILE *restrict fp) {
  long long start = _pos(buf, bufsize, *nextc);
  do {
	if (*nextc == *end) {
	  if (readmore(buf, nextc, end, bufsize, fp) <= 0)
		break;
	}
  } while (*((*nextc)++) != delim);
  return _pos(buf, bufsize, *nextc) - start;
}

int main(int argc, char** argv) {
  const size_t BUFSIZE = 16 * 1024;
  char *const buf = malloc(BUFSIZE);
  char *end = buf; char *nextc = buf;
  if (buf == NULL) {
    fprintf(stderr, "could not allocate buffer!");
    exit(errno);
  }

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
  FILE *fq;
  if (strncmp("-", argv[argc - 1], 2) == 0)
    fq = stdin;
  else
    fq = fopen(argv[argc - 1], "r");
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
	if ((c = nextchar(buf, &nextc, &end, BUFSIZE, fq)) != '@') {
	  if (c == EOF) break;
	  fprintf(stderr, "error: expected @<seqname> line (pos %lli): %s\n",
			  _pos(buf, BUFSIZE, nextc), argv[argc - 1]);
	  exit(1);
	}
	(void)scanuntil('\n', buf, &nextc, &end, BUFSIZE, fq);
	// <seq>
	do {
	  seq_len += scanuntil('\n', buf, &nextc, &end, BUFSIZE, fq) - 1;
	  newlines++;
	  c = nextchar(buf, &nextc, &end, BUFSIZE, fq);
	  if (c == '+')
		break;
	  if (c == '\n')
		newlines++;
	  else
		seq_len++;
	} while (1);

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
	(void)scanuntil('\n', buf, &nextc, &end, BUFSIZE, fq);

	// <qual>
	c = skip(seq_len + newlines, buf, &nextc, &end, BUFSIZE, fq);
  } while (c != EOF && (total_seqs != stop_after));

  if (ferror(fq) || errno) {
	fprintf(stderr, "error: %i (errno %i) occurred while parsing %s",
			ferror(fq), errno, argv[1]);
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
	//MEDIAN, UNIQUES, TOTAL SEQS, BYTES READ 
	printf("%li\t%li\t%lli\t%lli\n", median, num_nodes, total_seqs, _pos(buf, BUFSIZE, nextc));
	break;
  }
  return 0;
}
