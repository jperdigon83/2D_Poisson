#Makefile for 2DPoisson code

TARGET = 2DPoisson
N = 256 # Number of points in one direction of space
Lx = 1 # Space domain in x: [0,Lx]
Ly = 1 # Space domain in y: [0,Ly]

TARGET: 2DPoisson.cpp
	g++ -o $(TARGET) 2DPoisson.cpp

run:
	./$(TARGET) $(N) $(Lx) $(Ly)
clean:
	@echo "Deleting the executable"
	rm -f $(TARGET)
