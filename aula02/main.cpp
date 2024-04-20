#include "logic.h"

int main(void)
{
    srand(time(0));

    double **matrix;
    matrix = generateDoublePop(5, 2, 1.0, 5.0);
    printDoubleMatrix(matrix, 5, 2);

    return 0;
}

int clampi(int k, int low, int high)
{
    assert(low < high);
    if(k < low) k = low;
    if(k > high) k = high;
    return k;
}

void **allocateGenericMatrix(int width, int height, void *value, int type_size)
{
    type_size /= sizeof(char);

    void **matrix;
    matrix = (void **) malloc(sizeof(void *) * height);
    char **aux;
    aux = (char **) matrix;

    for (int y = 0; y < height; y++)
    {
        *(matrix + y) = malloc(type_size * width);
        for(int x = 0; x < width; x++)
        {
            if (value != NULL) memcpy((void *) (*(aux + y) + x * type_size), value, type_size);
        }
        // cout << "\n";
    }
    return matrix;
}

void printIntMatrix(int **matrix, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            cout << matrix[y][x] << " " ;
        }
        cout << endl;
    }
}

void printBoolMatrix(bool **matrix, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            cout << matrix[y][x] << " " ;
        }
        cout << endl;
    }
}

void printDoubleMatrix(double **matrix, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            cout << matrix[y][x] << " " ;
        }
        cout << endl;
    }
}

void permute(int *list, int size)
{
    for (int i = 0; i < size; i++)
    {
        int j = rand() % size;
        int aux = list[i];
        list[i] = list[j];        
        list[j] = aux;        
    }
}

int **generateIntPop(int ind_dim, int pop_size, int low, int high)
{
    assert(low < high);
    srand(time(0));

    int **matrix;
    matrix = (int **) allocateGenericMatrix(ind_dim, pop_size, (void *) &low, sizeof(int));

    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x++)
        {
            matrix[y][x] = (rand() % (high - low)) + low;
        }
    }

    return matrix;
}

int **generatePermPop(int ind_dim, int pop_size, int low, int high)
{
    assert(low < high);
    assert(high-low >= ind_dim);
    int **matrix;
    matrix = (int **) allocateGenericMatrix(ind_dim, pop_size, NULL, sizeof(int));
    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x ++)
        {
            try_again:
            int num = (rand() % (high - low)) + low;
            for (int i = 0; i < x; i++)
            {
                if (matrix[y][i] == num) goto try_again;
            }
            matrix[y][x] = num;
        }
    }
    return matrix;
}

bool **generateBoolPop(int ind_dim, int pop_size)
{
    // srand(time(0));
    bool temp = true;

    bool **matrix;
    matrix = (bool **) allocateGenericMatrix(ind_dim, pop_size, (void *) &temp, sizeof(bool));


    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x++)
        {
            matrix[y][x] = (rand() % 2 == 0);
        }
    }

    return matrix;
}

double **generateDoublePop(int ind_dim, int pop_size, double low, double high)
{
    assert(low < high);

    double **matrix;
    matrix = (double **) allocateGenericMatrix(ind_dim, pop_size, (void *) NULL, sizeof(double));

    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x++)
        {
            matrix[y][x] = ((double) rand() / RAND_MAX) * (high - low) + low;
        }
    }

    return matrix;
}