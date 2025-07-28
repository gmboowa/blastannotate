# Makefile for blastannotate

CXX = g++
CXXFLAGS = -std=c++17 -O3 -lpthread

TARGET = blastannotate
SRC = blastannotate.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
