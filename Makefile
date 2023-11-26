OBJ = raster3d.o raster_tri.o
EXE = raster3d
CXXFLAGS = -MMD -O3 -std=c++11 -g
CXX = g++

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE) -lm -lSDL2

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf $(EXE) $(OBJ)
