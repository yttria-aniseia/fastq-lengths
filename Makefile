prefix=/usr/local
exec_prefix=$(prefix)
bindir=$(exec_prefix)/bin

SRC=./src
BIN=./bin
NAME=fastq-lengths

CFLAGS=-O3 -Wextra -Wall -Wconversion -std=c2x -mtune=generic

INSTALL=install
INSTALL_PROGRAM=$(INSTALL)

.PHONY: all clean
all: fastq-lengths

$(SRC)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $(SRC)/$@ $^

fastq-lengths: $(SRC)/fastq-lengths.o
	mkdir -p $(BIN)
	$(CC) $(CFLAGS) -o$(BIN)/$@ $^

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
	rm -f $(BIN)/*
