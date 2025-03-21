CC=g++ -g -flto -O3
CFLAGS=-c -I. -std=c++11

all: clean main

main: .obj/main.o .obj/Graph.o
	${CC} .obj/main.o .obj/Graph.o -o main
	rm .obj/*.o

.obj/main.o: main.cpp
	${CC} ${CFLAGS} -o .obj/main.o main.cpp

.obj/Graph.o: Graph.cpp
	${CC} ${CFLAGS} -o .obj/Graph.o Graph.cpp

clean:
	rm -rf *o .obj/
	mkdir .obj
