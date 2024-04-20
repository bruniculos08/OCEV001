#include <bits/stdc++.h>

using namespace std;

#define GENERATIONS 20000
#define POP_SIZE 20
double C_min = 0;

typedef struct Population
{
    void **matrix;
    int type_size, ind_dim, pop_size;
} population;

typedef struct Function_info
{
    double (*f) (double);
    double a, b, C;
    int number_size, decimal_houses;
} function_info;


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
double fitness(bool *number, function_info &info);
// Torneio: selecionar n aleatórios e retornar o melhor.
int tournament(int sample_size, population &p, function_info &info);
// Cruzamento: misturar 2 individuos (50% de chance para cada alelo).
void *breed(int ind1, int ind2, population &p);
// Mutação: alterar randomicamente alelos de um indivíduo.
void mutation(bool *dna, population &p);
// Trocar a velha geração pela nova geração:
void changeGeneration(vector<bool *> &children, population &p);
// Densenvolver gerações:
int evolve(function_info &info, population &p);
// Normalizar o valor do(s) número(s) representados por cada indivíduo:
double normalize(function_info &info, double x);
// Função real a ser otmizada:
double f(double x);
// Encontrar C_min para função fitness:
double findCMin(function_info &info, int sample_size);