#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string.h>
#include <mpi.h>
#include "libtsp.h"
 
using namespace std;
 
unsigned int NCIUDADES;
int size, rango, tareaPeticion, optimoU;
int tagPeticion = 50;
int tagRespuesta = 60;
int tagOptimo = 70;
int tareasFinalizadas = 0;
bool fin = false;
MPI_Datatype nodoDatatype;
MPI_Request request;
MPI_Request requestOptimo;
MPI_Status status;

struct tNodoPar {
	int ci;
	int incl[100];
	int orig_excl;
	int dest_excl[100];
	int task_fin;
} nodoPar;

void asignarTrabajo(tPila *pila, tNodo *nodo) {
	int flag;
	//Comprobamos si tenemos alguna peticion de trabajo pendiente
	MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
	while(flag) {
		if(PilaPop(pila, nodo)) {
			//Adaptar el formato del struct para enviarlo
			nodoPar.ci = nodo->ci;
			for(int i=0; i<NCIUDADES; ++i) {
				nodoPar.incl[i] = nodo->incl[i];
			}
			nodoPar.orig_excl = nodo->orig_excl;
			for(int j=0; j<NCIUDADES; ++j) {
				nodoPar.dest_excl[j] = nodo->dest_excl[j];
			}
			nodoPar.task_fin = 0;
			//Enviar el struct a la tarea recibida
			MPI_Send(&nodoPar, 1, nodoDatatype, tareaPeticion, tagRespuesta, MPI_COMM_WORLD);
		}
		else {
			if(tareasFinalizadas < size-1) {
				++tareasFinalizadas;
				nodoPar.task_fin = 1;
				MPI_Send(&nodoPar, 1, nodoDatatype, tareaPeticion, tagRespuesta, MPI_COMM_WORLD);
			}
		}
		MPI_Irecv(&tareaPeticion, 1, MPI_INT, MPI_ANY_SOURCE, tagPeticion, MPI_COMM_WORLD, &request);
		MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
	}
}

void obtenerTrabajo(tPila *pila, tNodo *nodo) {
	//Enviar mensaje con el numero de tarea
	MPI_Send(&rango, 1, MPI_INT, 0, tagPeticion, MPI_COMM_WORLD);
	//Esperar recepcion del nodo
	if(MPI_Recv(&nodoPar, 1, nodoDatatype, 0, tagRespuesta, MPI_COMM_WORLD, &status) == MPI_SUCCESS) {
		//Adaptar el formato del struct y guardarlo
		nodo->ci = nodoPar.ci;
		for(int i=0; i<NCIUDADES; ++i) {
			nodo->incl[i] = nodoPar.incl[i];
		}
		nodo->orig_excl = nodoPar.orig_excl;
		for(int j=0; j<NCIUDADES; ++j) {
			nodo->dest_excl[j] = nodoPar.dest_excl[j];
		}
		nodo->task_fin = nodoPar.task_fin;
		if(nodo->task_fin == 1) {
			fin = true;
		}
		//Push a la pila con el nodo obtenido
		PilaPush(pila, nodo);
	}
	else {
		exit(2);
	}
}

void compartirOptimo(int U, bool nueva_U, tPila *pila) {
	int flag;
	//Si el proceso ha encontrado un optimo lo comparte con los otros procesos
	if(nueva_U) {
		for(int i=0; i<size; ++i) {
			if(i != rango) {
				MPI_Send(&U, 1, MPI_INT, i, tagOptimo, MPI_COMM_WORLD);
			}
		}
	}
	//Comprovar los optimos de los otros procesos
	MPI_Test(&requestOptimo, &flag, MPI_STATUS_IGNORE);
	while(flag) {
		if(optimoU < U) {
			//DEBUG
			printf("Proceso: %d, Nuevo optimo recibido: %d\n", rango, optimoU);
			//Actualizar el valor
			U = optimoU;
			//Acotar la pila con el nuevo optimo obtenido
			PilaAcotar (pila, U);
		}
		MPI_Irecv(&optimoU, 1, MPI_INT, MPI_ANY_SOURCE, tagOptimo, MPI_COMM_WORLD, &requestOptimo);
		MPI_Test(&requestOptimo, &flag, MPI_STATUS_IGNORE);
	}
}

int main (int argc, char **argv) {
    MPI::Init(argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rango);

	switch (argc) {
		case 3:
			NCIUDADES = atoi(argv[1]);
			break;
		default:
			cerr << "La sintaxis es: bbseq <tamaño> <archivo>" << endl;
			exit(1);
			break;
	}
 
	int** tsp0 = reservarMatrizCuadrada(NCIUDADES);
	tNodo	nodo,		// nodo a explorar
			lnodo,		// hijo izquierdo
			rnodo,		// hijo derecho
			solucion;
	bool activo,		// condicion de fin
		nueva_U;		// hay nuevo valor de c.s.
	int  U;				// valor de c.s.
	tPila pila;			// pila de nodos a explorar
	U = INFINITO;		// inicializa cota superior

	//Estructura de datos MPI
	int blocklenghts[5] = {1, 100, 1, 100, 1};
	MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Aint index;
	MPI_Type_extent(MPI_INT, &index);
	MPI_Aint offsets[5];
	offsets[0] = (MPI_Aint) 0;
	offsets[1] = index;
	offsets[2] = index + index*100;
	offsets[3] = index + index*100 + index;
	offsets[4] = index + index*100 + index + index*100;
	MPI_Type_struct(5, blocklenghts, offsets, types, &nodoDatatype);
	MPI_Type_commit(&nodoDatatype);
	InicNodo (&nodo);	// inicializa estructura nodo
	PilaInic (&pila);	// inicializa pila
	double t=MPI::Wtime();

	//El proceso master lee la matriz de fichero
	if(rango == 0) {
		LeerMatriz (argv[2], tsp0);
	}
	//Compartimos la matriz con todos los procesos
	MPI_Bcast(&tsp0[0][0], NCIUDADES*NCIUDADES, MPI_INT, 0, MPI_COMM_WORLD);
	activo = !Inconsistente(tsp0);

	if(rango == 0) {
		int tasks = size;
		if(size > 4) {
			tasks *= 5;
		}
		//El proceso master llena la pila hasta que el tamano sea igual al n. procesos
		while(PilaTamanio(&pila) < tasks) {
			Ramifica (&nodo, &lnodo, &rnodo, tsp0);
			PilaPush(&pila, &lnodo);
			PilaPush(&pila, &rnodo);
			PilaPop(&pila, &nodo);
		}
	}
	else {
		//Los procesos piden trabajo al master
		obtenerTrabajo(&pila, &nodo);
	}

	//Iniciamos las recepciones no bloqueantes de las peticiones y optimos
	MPI_Irecv(&tareaPeticion, 1, MPI_INT, MPI_ANY_SOURCE, tagPeticion, MPI_COMM_WORLD, &request);
	MPI_Irecv(&optimoU, 1, MPI_INT, MPI_ANY_SOURCE, tagOptimo, MPI_COMM_WORLD, &requestOptimo);

	PilaPop(&pila, &nodo);
	while (!fin) {
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

		//Compartir y obtener los optimos de los distintos procesos
		compartirOptimo(U, nueva_U, &pila);

		if(rango == 0) {
			//El proceso master atiende las peticiones de trabajo
			asignarTrabajo(&pila, &nodo);
		}
		else {
			if(PilaVacia(&pila)) {
				//Si la pila esta vacia pedimos trabajo al proceso master
				obtenerTrabajo(&pila, &nodo);
			}
		}

		//Sacamos el ultimo elemento de la pila
		PilaPop(&pila, &nodo);

		if(rango == 0 && tareasFinalizadas == size-1){
			fin = true;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
 	t=MPI::Wtime()-t;
    //El proceso master escribe la solucion
	if(rango == 0) {
		printf("Escribiendo la solucion, coste total(ci): %d\n", solucion.ci);
		EscribeSolucion(&solucion, t);
	}
	liberarMatriz(tsp0);

	MPI::Finalize();
	return 0;
}