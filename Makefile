SRC=./src
BIN=./bin

CFLAGS=-O2 -Wextra -Wall -Wconversion -mtune=generic

.PHONY: all clean
all: fastq-lengths fastq-median

$(SRC)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $(SRC)/$@ $^

fastq-lengths: $(SRC)/fastq-lengths.o
	mkdir -p bin
	$(CC) $(CFLAGS) -o$(BIN)/$@ $^

fastq-median: $(SRC)/fastq-median.o
	mkdir -p bin
	$(CC) $(CFLAGS) -o$(BIN)/$@ $^


clean:
	rm -f $(SRC)/*.o
	rm -f $(BIN)/fastq-lengths
	rm -f $(BIN)/fastq-median