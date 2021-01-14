#Makefile for 2DPoisson code


CXX = g++
TARGET = 2DPoisson
N = 128 # Number of points in one direction of space

TARGET: 2DPoisson.cpp
	g++ -o $(TARGET) 2DPoisson.cpp

run:
	./$(TARGET) $(N)
clean:
	@echo "Deleting the executable"
	rm -f $(TARGET)



