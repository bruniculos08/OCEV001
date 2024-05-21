#include <bits/stdc++.h>
#include <dirent.h>
#include <regex>
#include <gnuplot-iostream.h>

using namespace std;

enum crossover_mode {pmx, cx, generic};

#define POPULATION_SIZE 50
// #define MUTATION_PROBABILITY 0.005
// #define MUTATION_PROBABILITY 0.002
#define MUTATION_PROBABILITY 0.05
//  Se o seguinte valor for 0 ou maior que a dimensão dos individuos - 2, 
//  faz-se crossover uniforme:
#define CUT_POINTS_CROSSOVER 0
//  Se este valor for 0 ou maior que POPULATION_SIZE, 
//  todos os individuos participam do torneio:
#define TOURNAMENT_SAMPL 5
// #define GENERATIONS 160000
#define GENERATIONS 1200
//  Se o seguinte valor for "true" será feita seleção por torneio, 
//  e se for "false" será feita seleção por roleta viciada:
#define SELECTION_OPT false
#define CROSSOVER_OPT generic
//  Coeficiente de relação entre o fitness escalonado máximo e o
//  fitness escalonado médio:
#define SCALING_CONSTANT 1.2
//  Coeficientes de penalidade:
#define R_PENALTY -5.0
#define C_PENALTY 1.0
//  Flag especifiica para problema de programação liner (se o valor for
//  "true" os valores das variáveis serão considerados como inteiros, e se
//  for "false" os valores serão considerados como reais):
#define LINEAR_PROG_OPT true
//  Intervalo de gerações em que será printado os valores do melhor e do
//  pior (ou aleatório) indivíduo:
#define PRINT_INTERVAL 500
//  Flag para definir se o elitismo do melhor indivíduo irá repor o pior
//  indivíduo da nova geração (true) ou um indivíduo aleatório (false):
#define ELITISM_OPT true
//  Flag para ligar ou desligar todos os prints da função evolve:
#define PRINTS_OPT false

// #define INTERVAL_Y_AXIS 8
// #define INTERVAL_Y_AXIS 10
#define INTERVAL_Y_AXIS 0.1

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

typedef struct LinearProg
{
    //  Coeficientes da função objetivo:
    vector<double> obj_coeffs;
    //  Coeficientes das restrições (o último é 
    //  o valor do lado direito da restrição "<="):
    vector<vector<double>> leq_restr;
    //  Intervalos de cada variável (na forma [a, b]):
    vector<pair<double, double>> intervals;
    
    //  Número de casas inteiras de cada número:
    int integers_places, decimal_places;
} linearProg;

typedef struct Clause
{
    int variables[3] = {0, 0, 0};
} clause;

typedef struct Formula
{
    int formumla_size, number_of_variables;
    clause *expr;
} formula;

//  Vetor para salvar o valor do melhor indivíduo de cada geração:
vector<pair<double, double>> best_plot;
//  Vetor para salvar os valor médios dos indivíduos de cada geração:
vector<pair<double, double>> average_plot;
//  Vetor para salvar o valor do melhor pior de cada geração:
vector<pair<double, double>> worst_plot;
//  Vetor para salvar vetores de pontos de gráfico:
vector<vector<pair<double, double>>> bests_buffer;
vector<vector<pair<double, double>>> averages_buffer;
vector<vector<pair<double, double>>> worsts_buffer;

double worst_value = 0;
double best_value = 0;

//  Funções auxiliares:
int clampi(int x, int low, int high);
void allocateGenericMatrix(int width, int height, allele value, vector<vector<allele>> &v);
void printIntMatrix(vector<vector<allele>> &matrix);
void printBoolMatrix(vector<vector<allele>> &matrix, int width, int height);
void printDoubleMatrix(vector<vector<allele>> &matrix, int width, int height);
void permute(vector<int> &dna);
template <typename T> void eraseFast(vector<T> &v, int index);
template <typename T> T getRandom(vector<T> &v);
void generateGraph(vector<vector<pair<double, double>>> data_plot, string name, int opt);
double randomDouble(double max, double min);
int linearSearch(vector<allele> &v, allele e, int start, int end);

//  Funções para gerar população:
void generateIntPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);
void generatePermPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);
void generatePairPermPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);
void generateBoolPop(int ind_dim, int pop_size, vector<vector<allele>> &v);
void generateDoublePop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v);

//  Funções para modelagem do problema N-Rainhas:
//  - Gerar tabuleiro (considerando o problema das N-Rainhas):
int **createBoard(int width, int height);
//  - Printar solução encontrada para problema de programação linear:
void printLinearProg(linearProg &b, vector<allele> &ind);
//  - Densenvolver gerações (função genérica para qualquer problema desde que haja overload
//  de funções):
template <typename ProblemType>
int evolve(ProblemType &b, population &p);
//  - Trocar a velha geração pela nova geração (função genérica para qualquer 
//  problema desde que haja overload de funções):
void changeGeneration(vector<vector<allele>> &children, population &p);
//  - Printar tabuleiro (considerando o problema das N-Rainhas):
void printBoard(nqueens &b, vector<allele> &ind);
//  - Função Fitness: para SAT (avalia o número de clausulas satisfeitas por um indivíduo).
int fitness(formula &b, int ind_id, population &p);
//  - Função Fitness: para N-Rainhas (avalia o número de conflitos no tabuleiro).
int fitness(nqueens &b, int ind_id, population &p);
//  - Função Fitness: para programação linear (avalia o valor da função objetivo).
double fitness(linearProg &b, int ind_id, population &p);
//  - Torneio: selecionar n aleatórios e retornar o melhor (função genérica para qualquer 
//  problema desde que haja overload de funções):
template <typename ProblemType>
int tournament(ProblemType &b, int samples_num, population &p, vector<int> &single_ind);
//  - Roleta viciada (função genérica para qualquer problema desde que haja overload de funções):
template <typename ProblemType>
int roulette(vector<int> &single_ind, ProblemType &b, population &p);
//  - Realiza crossover entre todos (função genérica para qualquer problema desde que haja overload de funções):
template <typename ProblemType>
void crossoverAll(ProblemType &b, int samples_num, population &p, int cuts_num, vector<vector<allele>> &children);
//  - Cruzamento: misturar 2 individuos gerando 2 filhos usando método PMX (para problemas envolvendo índividuos permutados).
void crossoverPMX(int ind1, int ind2, population &p, int cuts_num, vector<vector<allele>> &children);
//  - Cruzamento: misturar 2 individuos gerando 2 filhos usando método CX (para problemas envolvendo índividuos permutados).
void crossoverCX(int ind1, int ind2, population &p, int cuts_num, vector<vector<allele>> &children);
//  - Cruzamento: misturar 2 indivíduos gerando 2 filhos usando o método genérico:
void crossoverGeneric(int ind1, int ind2, population &p, int cuts_num, vector<vector<allele>> &children);
//  - Mutação: alterar randomicamente indivíduo binário:
void mutationBin(vector<allele> &dna, population &p);
//  - Mutação: alterar randomicamente alelos de um indivíduo permutado:
void mutationPerm(vector<allele> &dna, population &p);
//  - Mutação: alterar randomicamente alelos de um indivíduo permutado:
void mutateAllelePerm(vector<allele> &dna, population &p, int allele_index, vector<int> &poll);
//  - Pegar o melhor indivíduo da geração (função genérica para qualquer problema desde que haja overload de funções):
template <typename ProblemType> 
int getElite(ProblemType &b, population &p, vector<allele> &elite_place);
//  - Pega o id do pior indivíduo:
template <typename ProblemType> 
int getWorst(ProblemType &b, population &p);
//  - Calcula a posição equivalente de um número no tabuleiro (considerando o problema das N-Rainhas):
pair<int, int> getPosition(nqueens &b, allele queen, vector<allele> &v);
//  - Verifica o número de conflitos nas diagonais de uma rainha (considerando o problema das N-Rainhas):
int lookDiag(vector<allele> &ind, nqueens &b, allele actual_queen);
//  - Verifica o número de conflitos na colunas e na linha de uma rainha (considerando o problema das N-Rainhas):
int lookRowAndCol(vector<allele> &ind, nqueens &b, allele actual_queen);

//  Funções para implementação de escalonamento linear (funções genéricas para qualquer problema desde que haja
//  overload de funções):
//  - Média da função fitness:
template <typename ProblemType> 
double fitnessAvg(ProblemType &b, population &p);
//  - Valor máximo da função fitness:
template <typename ProblemType> 
double fitnessMin(ProblemType &b, population &p);
//  - Valor mínimo da função fitness:
template <typename ProblemType> 
double fitnessMax(nqueens &b, population &p);
//  - Calculo dos valores alpha e beta:
template <typename ProblemType> 
void calcEsc(ProblemType &b, population &p, double c);
//  - Função fitness com escalonamento:
template <typename ProblemType> 
double fitnessEsc(ProblemType &b, int ind_id, population &p);

//  Funções para o problema 3-SAT
//  - Gerar problema 3SAT a partir de arquivo:
void fromFileGenerate3SAT(formula &sat, string file_path);
