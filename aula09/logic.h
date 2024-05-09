#include <bits/stdc++.h>
// #include <gnuplot-iostream.h>

using namespace std;

#define POPULATION_SIZE 30
#define MUTATION_PROBABILITY 0.001
//  Se o seguinte valor for 0 ou maior que a dimensão dos individuos - 2, 
//  faz-se crossover uniforme:
#define CUT_POINTS_CROSSOVER 2 
//  Se este valor for 0 ou maior que POPULATION_SIZE, 
//  todos os individuos participam do torneio:
#define TOURNAMENT_SAMPL 0 
#define GENERATIONS 1000
//  Se o seguinte valor for true será feita seleção por torneio, 
//  e se for false será feita seleção por roleta viciada:
#define SELECTION_OPT false 

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

typedef struct NQueens
{
    int **board;
    int dim;
} nqueens;

vector<pair<double, double>> data_plot;

//  Funções auxiliares:
int clampi(int x, int low, int high);
void allocateGenericMatrix(int width, int height, allele value, vector<vector<allele>> &v);
void printIntMatrix(vector<vector<allele>> &matrix);
void printBoolMatrix(vector<vector<allele>> &matrix, int width, int height);
void printDoubleMatrix(vector<vector<allele>> &matrix, int width, int height);
void permute(vector<int> &dna);
template <typename T> void eraseFast(vector<T> &v, int index);
template <typename T> T getRandom(vector<T> &v);
void generateGraph();
double randomDouble(double max, double min);

//  Funções para gerar população:
void generateIntPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);
void generatePermPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);
void generatePairPermPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);
void generateBoolPop(int ind_dim, int pop_size, vector<vector<allele>> &v);
void generateDoublePop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);

//  Funções para modelagem do problema N-Rainhas:
//  - Gerar problema da N-Rainhas a partir de arquivo:
void fromFileGenerateNQueens(nqueens &b, string file_path);
//  - Gerar tabuleiro:
int **createBoard(int width, int height);
//  - Densenvolver gerações:
int evolve(nqueens &b, population &p);
//  - Trocar a velha geração pela nova geração:
void changeGeneration(vector<allele> *children, population &p);
//  - Printar tabuleiro:
void printBoard(nqueens &b);
//  - Função Fitness: para SAT (avalia o número de clausulas satisfeitas por um indivíduo).
int fitness(nqueens &b, int ind_id, population &p);
//  - Torneio: selecionar n aleatórios e retornar o melhor.
int tournament(nqueens &b, int samples_num, population &p, vector<int> &single_ind);
//  - Roleta viciada:
int roulette(vector<int> &single_ind, nqueens &sat, population &p);
//  - Realiza crossover entre todos:
void crossoverAll(nqueens &b, int samples_num, population &p, bool genetic_operator);
//  - Cruzamento: misturar 2 individuos gerando 2 filhos usando método PMX.
void crossoverPMX(int ind1, int ind2, population &p, int cuts_num);
//  - Cruzamento: misturar 2 individuos gerando 2 filhos usando método CX.
void crossoverCX(int ind1, int ind2, population &p, int cuts_num);
//  - Mutação: alterar randomicamente alelos de um indivíduo.
void mutation(vector<allele> &dna, population &p);
//  - Pegar o melhor indivíduo da geração:
int getElite(nqueens &b, population &p, vector<allele> &elite_place);
//  - Pega o id do pior indivíduo:
int getWorst(nqueens &b, population &p);
//  - Calcula a posição equivalente de um número no tabuleiro
pair<int, int> getPosition(nqueens &b, int id);
//  - Verifica o número de conflitos nas diagonais de uma rainha:
int lookDiag(vector<allele> &ind, nqueens &b, allele actual_queen);
//  - Verifica o número de conflitos na colunas e na linha de uma rainha:
int lookRowAndCol(vector<allele> &ind, nqueens &b, allele actual_queen);