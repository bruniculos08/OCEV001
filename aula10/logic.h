#include <bits/stdc++.h>
#include <gnuplot-iostream.h>
using namespace std;


#define POPULATION_SIZE 30
#define MUTATION_PROBABILITY 0.01
//  Se o seguinte valor for 0 ou maior que a dimensão dos individuos - 2, 
//  faz-se crossover uniforme:
#define CUT_POINTS_CROSSOVER 2
//  Se este valor for 0 ou maior que POPULATION_SIZE, 
//  todos os individuos participam do torneio:
#define TOURNAMENT_SAMPL 5
#define GENERATIONS 1000
//  Se o seguinte valor for "true" será feita seleção por torneio, 
//  e se for "false" será feita seleção por roleta viciada:
#define SELECTION_OPT false
//  Se o seguinte valor for true será feito crossover PMX, e se
//  for "false" será feita crossover CX:
#define CROSSOVER_OPT false
//  Coeficiente de relação entre o fitness escalonado máximo e o
//  fitness escalonado médio:
#define SCALING_CONSTANT 1.8
// 6 de 20 (150 gerações, 4 pontos de corte (crossover PMX), torneio com 5 amostras, taxa de mutação = 0.01, população = 100 e C = 1.5)
// 2 de 5 (150 gerações, 4 pontos de corte (crossover PMX), torneio com 5 amostras, taxa de mutação = 0.02, população = 100 e C = 1.8)
// 2 de 5 (150 gerações, 4 pontos de corte (crossover PMX), torneio com 5 amostras, taxa de mutação = 0.02, população = 100 e C = 1.2)
// ? de ? (150 gerações, 2 pontos de corte (crossover PMX), roleta, taxa de mutação = 0.02, população = 100 e C = 1.2)
// 2 de 5 (150 gerações, crossover CX, roleta, taxa de mutação = 0.00001, população = 100 e C = 1.2)
// 0 de 2 (150 gerações, crossover CX, torneio com 15 amostras, taxa de mutação = 0.00001, população = 100 e C = 1.2)
// 8 de 20 (150 gerações, crossover CX, torneio com 15 amostras, taxa de mutação = 0.001, população = 100 e C = 1.2)
// 0 de 0 (150 gerações, crossover CX, torneio com 15 amostras, taxa de mutação = 0.001, população = 100 e C = 1.2)
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
    double alpha, beta;
} population;

typedef struct NQueens
{
    int **board;
    int dim;
} nqueens;

// Variáveis globais:
vector<pair<double, double>> data_plot;
vector<mutex> control_crossover(POPULATION_SIZE);

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
int linearSearch(vector<allele> &v, allele e, int start, int end);

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
void printBoard(nqueens &b, vector<allele> &ind);
//  - Função Fitness: para SAT (avalia o número de clausulas satisfeitas por um indivíduo).
int fitness(nqueens &b, int ind_id, population &p);
//  - Torneio: selecionar n aleatórios e retornar o melhor.
int tournament(nqueens &b, int samples_num, population &p, vector<int> &single_ind);
//  - Roleta viciada:
int roulette(vector<int> &single_ind, nqueens &b, population &p);
//  - Realiza crossover entre todos:
void crossoverAll(nqueens &b, int samples_num, population &p, int cuts_num);
//  - Cruzamento: misturar 2 individuos gerando 2 filhos usando método PMX.
void crossoverPMX(int ind1, int ind2, population *p, int cuts_num);
//  - Cruzamento: misturar 2 individuos gerando 2 filhos usando método CX.
void crossoverCX(int ind1, int ind2, population *p, int cuts_num);
//  - Mutação: alterar randomicamente alelos de um indivíduo.
void mutation(vector<allele> &dna, population &p);
//  - Mutação: alterar randomicamente alelos de um indivíduo.
void mutateAllele(vector<allele> &dna, population &p, int allele_index, vector<int> &poll);
//  - Pegar o melhor indivíduo da geração:
int getElite(nqueens &b, population &p, vector<allele> &elite_place);
//  - Pega o id do pior indivíduo:
int getWorst(nqueens &b, population &p);
//  - Calcula a posição equivalente de um número no tabuleiro
pair<int, int> getPosition(nqueens &b, allele queen, vector<allele> &v);
//  - Verifica o número de conflitos nas diagonais de uma rainha:
int lookDiag(vector<allele> &ind, nqueens &b, allele actual_queen);
//  - Verifica o número de conflitos na colunas e na linha de uma rainha:
int lookRowAndCol(vector<allele> &ind, nqueens &b, allele actual_queen);

// Funções para implementação de escalonamento linear:
//  - Média da função fitness:
double fitnessAvg(nqueens &b, population &p);
//  - Valor máximo da função fitness:
double fitnessMin(nqueens &b, population &p);
//  - Valor mínimo da função fitness:
double fitnessMax(nqueens &b, population &p);
//  - Calculo dos valores alpha e beta:
void calcEsc(nqueens &b, population &p, double c);
//  - Função fitness com escalonamento:
double fitnessEsc(nqueens &b, int ind_id, population &p);
