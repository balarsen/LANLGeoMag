all: ParticleDrift


ParticleDrift: ParticleDrift.c Makefile
	gcc -Wall ParticleDrift.c `pkg-config --cflags --libs lgm` -o ParticleDrift 



clean:
	rm -f *~ *.o ParticleDrift
