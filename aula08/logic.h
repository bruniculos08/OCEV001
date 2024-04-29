#include <bits/stdc++.h>

using namespace std;

#define GENERATIONS 20000
#define POP_SIZE 20
#define MUTATION_PROBABILITY 0.00001
double C_min = 0;

union allele
{
    int integer;
    bool boolean;
    double real;
};

typedef struct Population
{
    vector<vector<allele>> matrix;
    int ind_dim, pop_size;
} population;

typedef struct Function_info
{
    double (*f) (double);
    double a, b, C;
    int number_size, decimal_houses;
} function_info;


// Funções auxiliares:
int clampi(int x, int low, int high);
vector<vector<allele>> allocateGenericMatrix(int width, int height, allele &value);
void printIntMatrix(vector<vector<allele>> &matrix);
void printBoolMatrix(vector<vector<allele>> &matrix);
void printDoubleMatrix(vector<vector<allele>> &matrix);
void permute(vector<allele> &v);

// Funções para gerar população:
void generateIntPop(int width, int height, vector<vector<allele>> &v, allele low, allele high);
void generatePermPop(int width, int height, vector<vector<allele>> &v, allele low, allele high);
void generateBoolPop(int width, int height, vector<vector<allele>> &v);
void generateDoublePop(int width, int height, vector<vector<allele>> &v, allele low, allele high);

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