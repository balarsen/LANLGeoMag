all: Velocity


Velocity: Velocity.c
	gcc Velocity.c -Wall `pkg-config --cflags --libs lgm` -o Velocity

clean:
	rm -rf *~ *.o Velocity *.dSYM
