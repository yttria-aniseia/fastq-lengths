prefix=/usr/local
exec_prefix=$(prefix)
bindir=$(exec_prefix)/bin

SRC=./src
BIN=./bin
BIN_STATIC=$(BIN)/x86_64
NAME=fastq-lengths

CFLAGS=-O3 -Wextra -Wall -Wconversion -std=c2x -mtune=native
CFLAGS_STATIC=$(CFLAGS) -static -fdata-sections -ffunction-sections -Wl,--gc-sections -s -march=x86-64 -mtune=generic
INSTALL=install
INSTALL_PROGRAM=$(INSTALL)

.PHONY: all clean
all: fastq-lengths
static: fastq-lengths-static

$(SRC)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $(SRC)/$@ $^

fastq-lengths: $(SRC)/fastq-lengths.o
	mkdir -p $(BIN)
	$(CC) $(CFLAGS) -o$(BIN)/$@ $^
fastq-lengths-static: $(SRC)/fastq-lengths.o
	$(MAKE) BIN='$(BIN_STATIC)' CFLAGS='$(CFLAGS_STATIC)'\
		fastq-lengths

install: $(BIN)/$(NAME) installdirs
	$(INSTALL_PROGRAM) $< $(DESTDIR)$(bindir)
install-strip:
	$(MAKE) INSTALL_PROGRAM='$(INSTALL_PROGRAM) -s' \
		install
installdirs:
	mkdir -p $(DESTDIR)$(bindir)

uninstall:
	rm -f $(DESTDIR)$(bindir)/$(NAME)

clean:
	rm -f $(SRC)/*.o
	rm -rf $(BIN)/*
