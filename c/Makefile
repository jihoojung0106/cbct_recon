CC = gcc
CFLAGS = -Wall -Wextra
LIBS = -lfftw3 -lm

main: main.c
	$(CC) $(CFLAGS) -o main main.c $(LIBS)

run: main
	./main