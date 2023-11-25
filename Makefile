OBJ = raster3d.o raster_tri.o
EXE = raster3d
CXXFLAGS = -MMD -O3 -std=c++11 -g -fopenmp
CXX = g++

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE) -lm

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf $(EXE) $(OBJ)
