#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

void Solver::Init(unsigned N, float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
}
/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/
void Solver::FreeData(void)
{
//TODO: Libera los buffers de memoria.
	if (u != NULL)
	free(u);
	if (v != NULL)
	free(v);
	if (dens != NULL)
	free(dens);
	if (u_prev != NULL)
	free(u_prev);
	if (v_prev != NULL)
	free(v_prev);
	if (dens_prev != NULL)
	free(dens_prev);
}

void Solver::ClearData(void)
{
//TODO: Borra todo el contenido de los buffers

	for (unsigned int i = 0; i < (N + 2)*(N + 2); i++)
	{
		u[i] = 0;
		v[i] = 0;
		dens[i] = 0;
		u_prev[i] = 0;
		v_prev[i] = 0;
		dens_prev[i] = 0;
	}

}

bool Solver::AllocateData(void)
{
//TODO:
//Reservamos memoria, en caso de fallo devlvemos false.
//Antes de devolver true, hay que limpiar la memoria reservada con un ClearData().



	this->u = (float *)malloc(sizeof(float)* (N+2) * (N + 2));
	this->v = (float *)malloc(sizeof(float)* (N + 2) * (N + 2));
	this->dens = (float *)malloc(sizeof(float)* (N + 2) * (N + 2));
	this->u_prev = (float *)malloc(sizeof(float)* (N + 2) * (N + 2));
	this->v_prev = (float *)malloc(sizeof(float)* (N + 2) * (N + 2));
	this->dens_prev = (float *)malloc(sizeof(float)* (N + 2) * (N + 2));

	if ((u_prev == NULL) || (v_prev == NULL) || (dens_prev == NULL) || (u == NULL) || (v == NULL) || (dens == NULL)) exit(1);
	ClearData();

	return true;
}

void Solver::ClearPrevData() 
{
//TODO: Borra el contenido de los buffers _prev
	for (unsigned int i = 0; i < (N + 2)*(N + 2); i++)
	{
		u_prev[i] = 0;
		v_prev[i] = 0;
		dens_prev[i] = 0;
	}
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
//TODO: Añade el valor de source al array de densidades. Sería interesante usar la macro: XY_TO_ARRAY
	dens_prev[XY_TO_ARRAY(x, y)] += source;
}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
//TODO: Añade el valor de fuerza a sus respectivos arrays. Sería interesante usar la macro: XY_TO_ARRAY
	u_prev[XY_TO_ARRAY(x, y)] += forceX;
	v_prev[XY_TO_ARRAY(x, y)] += forceY;
}

void Solver::Solve()
{
	VelStep();
	DensStep();
}

void Solver::DensStep()
{
	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
}

void Solver::VelStep()
{
	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev, v)
	Diffuse(1, u, u_prev);  
	Diffuse(2, v, v_prev); 
	Project(u, v, u_prev, v_prev);		//Mass conserving.
	SWAP (u_prev,u)			
	SWAP (v_prev,v)
	Advect(1, u, u_prev, u_prev, v_prev);
	Advect(2, v, v_prev, u_prev, v_prev);
	Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
//TODO: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
	
	for(unsigned int i=0; i < (N + 2)*(N + 2); i++)
		base[i] += source[i];

}


void Solver::SetBounds(int b, float * x)
{
/*TODO:
Input b: 0, 1 or 2.
	0: borders = same value than the inner value.
	1: x axis borders inverted, y axis equal.
	2: y axis borders inverted, x axis equal.
	Corner values allways are mean value between associated edges.
*/


	if (b == 1) {
		for (unsigned i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = x[XY_TO_ARRAY(i, 1)];//y normal x opuesto
			x[XY_TO_ARRAY(i, N + 1)] = x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(N + 1, i)] = -x[XY_TO_ARRAY(N, i)];
			x[XY_TO_ARRAY(0, i)] = -x[XY_TO_ARRAY(1, i)];
		}
		x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
		x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
		x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
		x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;
	}
	else if (b == 2) {
		for (unsigned i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = -x[XY_TO_ARRAY(i, 1)];//y opuesto x normal
			x[XY_TO_ARRAY(i, N + 1)] = -x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(N + 1, i)] = x[XY_TO_ARRAY(N, i)];
			x[XY_TO_ARRAY(0, i)] = x[XY_TO_ARRAY(1, i)];
		}
		x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
		x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
		x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
		x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;
	}
	else if (b == 0) {
		for (int i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = x[XY_TO_ARRAY(i, 1)];
			x[XY_TO_ARRAY(i, N + 1)] = x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(N + 1, i)] = x[XY_TO_ARRAY(N, i)];
			x[XY_TO_ARRAY(0, i)] = x[XY_TO_ARRAY(1, i)];
		}
		x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
		x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
		x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
		x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;
	}
}

/*
https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel <- Solución de valores independientes.
Despreciando posibles valores de x no contiguos, se simplifica mucho. Mirar diapositivas y la solución de Gauss Seidel de términos independientes.
Gauss Seidel -> Matrix x and x0
*/
void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii)
{
	for (int i = 0; i <= N + 1; i++) {
		for (int j = 0; j <= N + 1; j++)
		{
			x[(N + 2)*i + j] = x0[(N + 2)*i + j];
		}
	}
	for (int r = 0; r < 20; r++) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY(i, j)] = (1 / aii)*(aij*(x[XY_TO_ARRAY(i, j - 1)] + x[XY_TO_ARRAY(i - 1, j)] 
				+ x[XY_TO_ARRAY(i + 1, j )] + x[XY_TO_ARRAY(i, j + 1)]) + x0[XY_TO_ARRAY(i, j)]);
		END_FOR

			SetBounds(b, x);
	}
}

/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0)
{
//Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
	float aij = diff * dt * (N * N);
	float aii =   1 + 4 * diff * dt * (N * N) ;
	LinSolve(b, x, x0, aij, aii);

}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
void Solver::Advect(int b, float * d, float * d0, float * u, float * v)
{
//Se aplica el campo vectorial realizando una interploación lineal entre las 4 casillas más cercanas donde caiga el nuevo valor.

	float x, y, difa, difc;
	int a, c;

	//int posa, posc;
	//float factora, factorc, x1, y1;

	for (int i = 1; i < N + 1; i++) 
	{
		for (int j = 1; j < N + 1; j++)
		{
			x = (float)i - u[XY_TO_ARRAY(i, j)] * N *dt;
			y = (float)j - v[XY_TO_ARRAY(i, j)] * N *dt;

			if (x < 0.5f) x = 0.5f;
			if (y < 0.5f) y = 0.5f;
			if (x > N + 0.5f) x = N + 0.5f;
			if (y > N + 0.5f) y = N + 0.5f;

			a = (int)x;
			c = (int)y;
		
			difa = x - a;
			difc = y - c;

			//if (difa > 0.5)
			//	a = a + 1;

			//if (difc > 0.5)
			//	c = c + 1;

			//if (difa == 0.5)
			//	factora = 0.5;
			//else
			//factora = 2 * (0.5 - difa);

			//if (difc == 0.5)
			//	factorc = 0.5;
			//else
			//factorc = 2 * (0.5 - difc);


			//	if (difa <= 0.5)
			//	{
			//		posa = a + 1;
			//	}
			//	else
			//	{
			//		posa = a - 1;
			//	}

			//	if (difc <= 0.5)
			//	{
			//		posc= c + 1;
			//	}
			//	else
			//	{
			//		posc = c - 1;
			//	}

			d[XY_TO_ARRAY(i, j)] = ((1 - fabs(difa))* (1 - fabs(difc)))*d0[XY_TO_ARRAY(a, c)]
				+ ((1 - fabs(difa))* fabs(difc))*d0[XY_TO_ARRAY(a, c + 1)]
				+ (fabs(difa)* (1 - fabs(difc)))*d0[XY_TO_ARRAY(a + 1, c)]
				+ (fabs(difa)* fabs(difc))*d0[XY_TO_ARRAY(a + 1, c + 1)];



			//d[XY_TO_ARRAY(i, j)] = ((1 - fabs(factora))* (1 - fabs(factorc)))*d0[XY_TO_ARRAY(posa, posc)]
			//	+ ((1 - fabs(factora))* fabs(factorc))*d0[XY_TO_ARRAY(posa, c)]
			//	+ (fabs(factora)* (1 - fabs(factorc)))*d0[XY_TO_ARRAY(a, posc)]
			//	+ (fabs(factora)* fabs(factorc))*d0[XY_TO_ARRAY(a, c)];
		}
	}

	SetBounds(b, d);
}

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
		p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	SetBounds(0, div);
	SetBounds(0, p);

	LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
		v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	SetBounds(1, u);
	SetBounds(2, v);
}