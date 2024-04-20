#include <bits/stdc++.h>

using namespace std;


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