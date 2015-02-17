.SUFFIXES: .cu .cuh .h

INC = -I. -I./include
LIBDIR =
LIB =

CC = icc -openmp
CPP = h5pcc
MPICC = h5pcc

CFLAGS = $(INC) -std=c99
CPPFLAGS = $(INC)

srcdir = src
gpusrc = kernel_l1.c kernel_l2.c kernel_l3.c
cpusrc = host_main.c host_func.c host_launcher.c
exec = ising

default: $(exec) 


%.o: %.c
	$(CC) -O3 $(CFLAGS) -c $<

%.o: %.cpp
	$(CPP) -O3 $(CPPFLAGS) -c $<

mpiprocess.o: mpiprocess.c
	$(MPICC) -c $<


ising: host_main.o host_func.o host_launcher.o kernel_l1.o kernel_l2.o kernel_l3.o
	$(CC) -O3 -o $@ $^

mpi_ising: host_main.o mpiprocess.o host_func.o host_launcher.o kernel_l1.o kernel_l2.o kernel_l3.o
	$(MPICC) -O3 $(LIB) $(LIBDIR) -o $@ $^

prof: $(gpusrc) $(cpusrc)
	$(CPP) -O0 $(CPPFLAGS) $(PROFFLAG) $^

profclean: $(gpusrc) $(cpusrc)
	$(CPP) -O0 $(CPPFLAGS) $(PROFFLAG) -clean $^


clean:
	rm -r *.o $(exec)

cleanoutput:
	rm -r output_*

