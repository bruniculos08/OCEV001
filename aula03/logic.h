#include <bits/stdc++.h>

using namespace std;
double c = 0;

typedef struct Population
{
    void **matrix;
    int type_size, ind_dim, pop_size;
} population;

typedef struct Clause
{
    int variables[3];
} clause;

typedef struct Formula
{
    int formumla_size, number_of_variables;
    clause *expr;
} formula;

// Funções auxiliares:
int clampi(int x, int low, int high);
void **allocateGenericMatrix(int width, int height, void *value, int type_size);
void printIntMatrix(int **matrix, int width, int height);
void printBoolMatrix(bool **matrix, int width, int height);
void printDoubleMatrix(double **matrix, int width, int height);
void permute(int *list, int size);
int removeRandomId(vector<int> &v);

// Funções para gerar população:
int **generateIntPop(int ind_dim, int pop_size, int low, int high);
int **generatePermPop(int ind_dim, int pop_size, int low, int high);
bool **generateBoolPop(int ind_dim, int pop_size);
double **generateDoublePop(int ind_dim, int pop_size, double low, double high);

/* --------------------- Funções para modelagem do problema 3-SAT --------------------- */

// Gerar problema 3SAT de n clausulas (note que "variables_num" deve ser igual à "pop_size"):
void generate3SAT(formula &sat, int formula_size, int variables_num);
// Gerar problema 3SAT a partir de arquivo:
void fromFileGenerate3SAT(formula &sat, string file_path);
// Densenvolver gerações:
int evolve(formula &sat, population &p);
// Trocar a velha geração pela nova geração:
void changeGeneration(vector<bool *> &children, population &p);
// Printar fórmula 3SAT:
void print3SAT(formula &sat);
// Função Fitness: para SAT (avalia o número de clausulas satisfeitas por um indivíduo).
int fitness(formula &sat, int ind_id, population &p);
// Torneio: selecionar n aleatórios e retornar o melhor.
int tournament(formula &sat, int sample_num, population &p, int exclude_id);
// Cruzamento: misturar 2 individuos (50% de chance para cada alelo).
void *crossover(int ind1, int ind2, population &p);
// Mutação: alterar randomicamente alelos de um indivíduo.
void mutation(bool *dna, population &p);
