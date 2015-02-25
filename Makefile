.SUFFIXES: .cu .cuh .h

INC = -I. -I./include
LIBDIR =
LIB =

#CC = gcc -O3
#CC = icc -openmp -O3
CC = clang -O3

#pollycc = clang -Xclang -load -Xclang /home/yfang11/compiler/llvm/3.4.2+polly/llvm_obj/lib/LLVMPolly.so
#CC = $(pollycc) -O3 -mllvm -polly -mllvm -polly-vectorizer=polly



CFLAGS = $(INC) -std=c99

srcdir = src
gpusrc = kernel_l1.c kernel_l2.c kernel_l3.c
cpusrc = host_main.c host_func.c host_launcher.c
exec = ising

default: $(exec) 


%.o: %.c
	$(CC) $(CFLAGS) -c $<

ising: host_main.o host_func.o host_launcher.o kernel_l1.o kernel_l2.o kernel_l3.o
	$(CC) -o $@ $^ -lm



clean:
	rm -r *.o $(exec)

cleanoutput:
	rm -r output_*

