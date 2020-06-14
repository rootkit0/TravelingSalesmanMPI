/* ******************************************************************** */
/*               Algoritmo Branch-And-Bound Secuencial                  */
/* ******************************************************************** */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include "libtsp.h"
 
using namespace std;
 
unsigned int NCIUDADES;
 
main (int argc, char **argv) {
        MPI::Init(argc,argv);
	switch (argc) {
		case 3:		NCIUDADES = atoi(argv[1]);
					break;
		default:	cerr << "La sintaxis es: bbseq <tamaño> <archivo>" << endl;
					exit(1);
					break;
	}
 
	int** tsp0 = reservarMatrizCuadrada(NCIUDADES);
	tNodo	nodo,         // nodo a explorar
			lnodo,        // hijo izquierdo
			rnodo,        // hijo derecho
			solucion;     // mejor solucion
	bool activo,        // condicion de fin
		nueva_U;       // hay nuevo valor de c.s.
	int  U;             // valor de c.s.
	tPila pila;         // pila de nodos a explorar
 
	U = INFINITO;                  // inicializa cota superior
	InicNodo (&nodo);              // inicializa estructura nodo
	PilaInic (&pila);              // inicializa pila
	LeerMatriz (argv[2], tsp0);    // lee matriz de fichero
	double t=MPI::Wtime();
	
	activo = !Inconsistente(tsp0);
	
	while (activo) {       // ciclo del Branch&Bound
		Ramifica (&nodo, &lnodo, &rnodo, tsp0);
		nueva_U = false;
		if (Solucion(&rnodo)) {
			if (rnodo.ci < U) {    // se ha encontrado una solucion mejor
				U = rnodo.ci;
				nueva_U = true;
				CopiaNodo (&rnodo, &solucion);
			}
		}
		else {                    //  no es un nodo solucion
			if (rnodo.ci < U) {     //  cota inferior menor que cota superior
				if (!PilaPush (&pila, &rnodo)) {
					printf ("Error: pila agotada\n");
					liberarMatriz(tsp0);
					exit (1);
				}
			}
		}
		if (Solucion(&lnodo)) {
			if (lnodo.ci < U) {    // se ha encontrado una solucion mejor
				U = lnodo.ci;
				nueva_U = true;
				CopiaNodo (&lnodo,&solucion);
			}
		}
		else {                     // no es nodo solucion
			if (lnodo.ci < U) {      // cota inferior menor que cota superior
				if (!PilaPush (&pila, &lnodo)) {
					printf ("Error: pila agotada\n");
					liberarMatriz(tsp0);
					exit (1);
				}
			}
		}
		if (nueva_U) PilaAcotar (&pila, U);
		activo = PilaPop (&pila, &nodo);
	}
 	t=MPI::Wtime()-t;
 	
    	MPI::Finalize();
	EscribeSolucion(&solucion, t);
	liberarMatriz(tsp0);
}