#include <bits/stdc++.h>
#include <gnuplot-iostream.h>

using namespace std;

#define POPULATION_SIZE 30
#define MUTATION_PROBABILITY 0.20
#define CUT_POINTS_CROSSOVER 10
#define TOURNAMENT_SAMPL 0
#define GENERATIONS 2000

union allele
{
    int integer;
    bool boolean;
    double real;
};

typedef struct Population
{
    vector<allele> *matrix;
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

vector<pair<double, double>> data_plot;

// Funções auxiliares:
int clampi(int x, int low, int high);
vector<allele> *allocateGenericMatrix(int width, int height, allele value);
void printIntMatrix(vector<allele> *matrix, int width, int height);
void printBoolMatrix(vector<allele> *matrix, int width, int height);
void printDoubleMatrix(vector<allele> *matrix, int width, int height);
void permute(vector<int> &dna);
template <typename T> void eraseFast(vector<T> &v, int index);
template <typename T> T getRandom(vector<T> &v);
void generateGraph();

// Funções para gerar população:
vector<allele> *generateIntPop(int ind_dim, int pop_size, allele low, allele high);
vector<allele> *generatePermPop(int ind_dim, int pop_size, allele low, allele high);
vector<allele> *generateBoolPop(int ind_dim, int pop_size);
vector<allele> *generateDoublePop(int ind_dim, int pop_size, allele low, allele high);

/* --------------------- Funções para modelagem do problema 3-SAT --------------------- */

// Gerar problema 3SAT de n clausulas (note que "variables_num" deve ser igual à "pop_size"):
void generate3SAT(formula &sat, int formula_size, int variables_num);
// Gerar problema 3SAT a partir de arquivo:
void fromFileGenerate3SAT(formula &sat, string file_path);
// Densenvolver gerações:
int evolve(formula &sat, population &p);
// Trocar a velha geração pela nova geração:
void changeGeneration(vector<allele> *children, population &p);
// Printar fórmula 3SAT:
void print3SAT(formula &sat);
// Função Fitness: para SAT (avalia o número de clausulas satisfeitas por um indivíduo).
int fitness(formula &sat, int ind_id, population &p);
// Torneio: selecionar n aleatórios e retornar o melhor.
int tournament(formula &sat, int samples_num, population &p, vector<int> &single_ind);
// Roleta:
int roulette(vector<int> &single_ind, formula &sat, population &p);
// Realiza crossover entre todos:
void crossoverAll(formula &sat, int samples_num, population &p, int cuts_num);
// Cruzamento: misturar 2 individuos (50% de chance para cada alelo).
void crossover(int ind1, int ind2, population &p, int cuts_num);
// Mutação: alterar randomicamente alelos de um indivíduo.
void mutation(vector<allele> &dna, population &p);
// Pegar o melhor indivíduo da geração:
int getElite(formula &sat, population &p, vector<allele> &elite_place);
// Pega o id do pior indivíduo:
int getWorst(formula &sat, population &p);
