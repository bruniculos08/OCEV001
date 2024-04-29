#include <bits/stdc++.h>

using namespace std;

#define GENERATIONS 500000
#define MUTATION 0.01

typedef struct Population
{
    void **matrix;
    int type_size, ind_dim, pop_size;
} population;

typedef struct Optimization
{
    double (*objective_function) (int, int);
    bool (*restriction_function) (int, int, int);
    int x_a, x_b, y_a, y_b;
    double restriction_value;
} optimization;


// Funções auxiliares:
int clampi(int x, int low, int high);
void **allocateGenericMatrix(int width, int height, void *value, int type_size);
void printIntMatrix(int **matrix, int width, int height);
void printBoolMatrix(bool **matrix, int width, int height);
void printDoubleMatrix(double **matrix, int width, int height);
void permute(int *list, int size);

// Funções para gerar população:
int **generateIntPop(int ind_dim, int pop_size, int low, int high);
int **generatePermPop(int ind_dim, int pop_size, int low, int high);
bool **generateBoolPop(int ind_dim, int pop_size);
double **generateDoublePop(int ind_dim, int pop_size, double low, double high);

/* --------------------- Funções para modelagem do problema de função real --------------------- */

// Função fitness: para avaliação de um indivíduo em relação a uma função f
double fitness(bool *number, int number_size, optimization &o);
// Torneio: selecionar n aleatórios e retornar o melhor.
int tournament(optimization &o, int sample_size, population &p);
// Mutação: 
void mutation(bool *number, int number_size);
// Cruzamento: misturar 2 individuos (50% de chance para cada alelo).
void *breed(int ind1, int ind2, population &p); 
// 
int evolve(optimization &o, population &p);
// 
void changeGeneration(vector<bool *> &children, population &p);
// 
double profit(int x, int y);
// 
bool workers(int x, int y, int restriction_value);
