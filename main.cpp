#include <iostream>
#include <cstdlib>
#include "SimplexAbstract.hpp"
#include <math.h>
#include <list>
#include <time.h>
#include "SimplexAlpha.hpp"
//#include "RCC.hpp"
#include "Covering_CORE.hpp"
#include "OptimizedCovering.hpp"
#include "Lex.hpp"
#include <mpi.h>
#include <functional>

/*compile with*/
/*g++ main.cpp SimplexAbstract.cpp RCC.hpp SimplexAlpha.cpp LocalSearch.cpp Matrix.cpp Covering.hpp OptimizedCovering.hpp -o motion */
/*ffmpeg -r 36 -f image2 -s 500x500 -i topology%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4*/

using namespace std;

Complex komplex;
VectorInt v0, v1, v2, v3, v4;
MatrixVectorInt matV;
SimplexAlpha simplex, simplex0;
MatrixSimplexAlpha com;

Simplex s0, s1, s2, s3, s4;
SimplexProd subComplex;

SimplicialMap map0, copyMap, map1, map2;
ComplexProduct KxK;
SubComplexJ L, JsubL, JsubL1, JsubL2;
LocalSearch suchen, suchenCore1, suchen2; 
MatrixIntList ll;
MatrixInt Test;
SimplexAlpha testSimplex;

//RCC rcc, addFacets;
Covering c;
OptimizedCovering O;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char **argv) {
	
	srand(time(NULL));
	
	int M = 12000;
        double r = 0.1;
	int Ns = 10;

	int maxInt = 3;
	int numOfMaxSimplex = 4;
	komplex.initComplex(numOfMaxSimplex, maxInt + 1);

	komplex.K.A[0][0].initSimplex(3, 0, 1, 2);
	komplex.K.A[0][1].initSimplex(3, 0, 1, 3);
	komplex.K.A[0][2].initSimplex(3, 0, 2, 3);
	komplex.K.A[0][3].initSimplex(3, 1, 2, 3);

	//komplex.K.A[0][3].initSimplex(2, 1, 2);
	KxK.initComplexProduct(komplex);
	KxK.escMaximalSimplices();
	map1.projection1(KxK);
	map0.projection2(KxK);

	L.initSubComplexJ(96);
	int counting = 0;

	for (int i = 0; i < KxK.listOfFacets.m; i++) {
                for (int j = 0; j < KxK.listOfFacets.rowLength.getA(0, i); j++) {
                        L.initA(counting, KxK.listOfFacets.A[i][j]);
                        counting += 1;
                }
        }

	 JsubL.initSubComplexJ(9);
	 JsubL.initA(0, KxK.listOfFacets.A[3][1]);
	 JsubL.initA(1, KxK.listOfFacets.A[8][1]);
	 JsubL.initA(2, KxK.listOfFacets.A[7][1]);
	 JsubL.initA(3, KxK.listOfFacets.A[2][1]);
	 JsubL.initA(4, KxK.listOfFacets.A[7][0]);
	 JsubL.initA(5, KxK.listOfFacets.A[4][0]);
	 JsubL.initA(6, KxK.listOfFacets.A[3][0]);
	 JsubL.initA(7, KxK.listOfFacets.A[2][0]);
	 JsubL.initA(8, KxK.listOfFacets.A[0][0]);

	 //L.initZero_Skeleton();
	 JsubL.initZero_Skeleton();

	 komplex.initAdjMat();
	 komplex.graph.addWeight(0, 1, 1);
	 komplex.graph.addWeight(0, 2, 1);
	 komplex.graph.addWeight(1, 2, 1);

	//suchen.initLocalSearch(JsubL, komplex, map1, map0, M, r);	

	int ccount = 0;
	
	//c.initCovering(L, komplex, map1, map0, M, r, argc, argv);
        //c.runCovering(L, komplex, map1, map0, M, r, argc, argv);
        //c.endCovering();
	//suchen.updateLocalSearch(JsubL, komplex, map1, map0, M, r);
	
	//rcc.initRCC_CORE(L, komplex, map1, map0, M, r, argc, argv);
	//rcc.runRCC_CORE(L, komplex, map1, map0, M, r, argc, argv);
	//rcc.escRCCResults();

	//rcc.endRCC_CORE();
	
	O.initOptimizedCovering_CORE(L, komplex, map1, map0, M, r,  argc, argv);
	O.runOptimizedCovering_CORE(L, komplex, map1, map0, M, Ns, r, argc, argv);

	cout << "\n\n\n";
	
	return 0;
}
