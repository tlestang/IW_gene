CXXFLAGS=-O3
EXECNAME=run_IW

all: LatticeSite.o Walls.o topo.o LatticeBoltzmann.o write_vtk.o main.o
	g++ -o $(EXECNAME) LatticeSite.o Walls.o topo.o LatticeBoltzmann.o write_vtk.o main.o

LatticeSite.o: LatticeSite.cpp
	g++ -o LatticeSite.o -c LatticeSite.cpp $(CXXFLAGS)
Walls.o: Walls.cpp
	g++ -o Walls.o -c Walls.cpp $(CXXFLAGS)
topo.o: topo.cpp
	g++ -o topo.o -c topo.cpp $(CXXFLAGS)
LatticeBoltzmann.o: LatticeBoltzmann.cpp
	g++ -o LatticeBoltzmann.o -c LatticeBoltzmann.cpp $(CXXFLAGS)
write_vtk.o: write_vtk.cpp
	g++ -o write_vtk.o -c write_vtk.cpp
main.o: main.cpp
	g++ -o main.o -c main.cpp $(CXXFLAGS)

clean:
	rm -rf *.o
squeakyclean: clean
	rm -f $(EXECNAME)
