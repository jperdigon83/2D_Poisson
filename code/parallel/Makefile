#Makefile for 2DPoisson code

CXX = g++
TARGET = 2DPoisson_parallel
Np = 4
N = 128 # Number of points in one direction of space
Lx = 1 # Space domain in x: [0,Lx]
Ly = 1 # Space domain in y: [0,Ly]

TARGET: 2DPoisson_parallel.cpp
	mpic++ -o $(TARGET) 2DPoisson_parallel.cpp

run:
	mpirun -n $(Np) ./$(TARGET) $(N) $(Lx) $(Ly)
clean:
	@echo "Deleting the executable"
	rm -f $(TARGET)
