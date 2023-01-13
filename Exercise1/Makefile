all: interpolationspolynom


interpolationspolynom: interpolationspolynom.c
	gcc -Wall -O1 --std=c99 -o interpolationspolynom interpolationspolynom.c -lm

## unter MinGW evtl.:
#interpolationspolynom: interpolationspolynom.c
#	gcc -Wall -O1 --std=gnu99 -o interpolationspolynom interpolationspolynom.c -lm

clean:
	rm interpolationspolynom
