#include "logic.h"

int main(void)
{
    srandom(time(0));

    population p;
    p.pop_size = 10;
    p.ind_dim = 5 + 5;
    p.type_size = sizeof(bool);

    p.matrix = (void **) generateBoolPop(p.ind_dim, p.pop_size);
    printBoolMatrix((bool **) p.matrix, p.ind_dim, p.pop_size);

    optimization o;
    o.objective_function = profit;
    o.restriction_function = workers;
    o.restriction_value = 40;

    int xmen = evolve(o, p);
    cout << "final best = " << xmen << endl;
    cout << "final best fitness function = " << fitness((bool *)p.matrix[xmen], p.ind_dim, o) * 1360.0 << endl;
    printBoolMatrix((bool **) p.matrix, p.ind_dim, p.pop_size);

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
            if (value != NULL) 
                memcpy((void *) (*(aux + y) + x * type_size), value, type_size);
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
        int j = random() % size;
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
            matrix[y][x] = (random() % (high - low)) + low;
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
            int num = (random() % (high - low)) + low;
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
            matrix[y][x] = (bool) (random() % 2 == 0);
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
            matrix[y][x] = ((double) random() / RAND_MAX) * (high - low) + low;
        }
    }

    return matrix;
}

/* --------------------- Funções para modelagem do problema dos rádios --------------------- */

double profit(int x, int y)
{
    // Substituindo x por 24 (x_max) e y por 16 (y_max) tem-se 1360 como resultado...
    // ... por isso para se normalizar o resultado entre [0,1] se divide por 1360:
    return (double) (30 * x + 40 * y) / 1360.0;
}

bool workers(int x, int y, int restriction_value)
{
    // Não é necessário normalizar essa função pois ela serve apenas para verificar...
    // ... se os valores de x e y atendem à restrição:
    return x + 2 * y <= restriction_value;
}

int evolve(optimization &o, population &p)
{ 
    int current_best_id = 0, counter = 0;
    vector<bool *> children;
    do
    {
        int ind1 = tournament(o, p.pop_size/2, p), ind2 = tournament(o, p.pop_size/2, p);
        bool *son;
        son = (bool *) breed(ind1, ind2, p);
        mutation(son, p.ind_dim);
        children.push_back(son);

        if(children.size() == p.pop_size)
        {
            changeGeneration(children, p);
            children.clear();
        }

        for (int i = 0; i < p.pop_size; i++)
        {
            if (fitness((bool *) p.matrix[i], p.ind_dim, o) >= fitness((bool *) p.matrix[current_best_id], p.ind_dim, o))
            {
                current_best_id = i;
            }
        }
        counter++;
    } while (counter != GENERATIONS);
    return current_best_id;
}

void changeGeneration(vector<bool *> &children, population &p)
{
    int current_id = 0;
    for (bool *dna : children)
    {
        for (int i = 0; i < p.ind_dim; i++)
        {
            memcpy((char *) *(p.matrix + current_id), dna, p.ind_dim);
        }
        current_id++;
    }
}

double fitness(bool *number, int number_size, optimization &o)
{
    int x = 0, y = 0;
    double penalty = 0.0;
    for (int i = number_size/2 - 1; i >= 0; i--)
    {
        if(number[i + number_size/2]) x += (int) pow(2.0, (double) i);
        if(number[i]) y += (int) pow(2.0, (double) i);
    }

    // Colocar os valores das variáveis dentro do intervalo de restrição [0, 24] e [0, 16], seguindo...
    // as fórmulas de normalização:
    // x = x_min + (x_max - x_min)/((2^L) - 1) * d
    // y = y_min + (y_max - y_min)/((2^L) - 1) * d
    // Onde L é o número de bits de cada número (nesse caso ambos possuem 5 bits).
    x = ((24.0/(pow(2.0, (double) (number_size/2)) - 1.0)) * (double) x);
    y = ((16.0/(pow(2.0, (double) (number_size/2)) - 1.0)) * (double) y);
    // Obs.: o tamanho de cada número na verdade é number_size/2.

    // A função de penalidade também deve ser normalizada para o intervalo [0,1]:
    penalty = max(0.0, ((double) x + 2.0* (double) y - 40.0)/16.0);
    // Note que substituindo-se x por 24 (x_max) e y por 16 (y_max), temos x + 2*y - 40 = 16, por isso...
    // ... se divide por 16 para normalizar.
    return (o.restriction_function(x, y, o.restriction_value)) ? o.objective_function(x, y) - penalty : - penalty;
}

int tournament(optimization &o, int sample_size, population &p)
{
    int choosen_ind_id = -1, rand_id; 
    double best_fitness = 0.0;
    for(int i = 0; i < sample_size; i++)
    {
        rand_id = (random() % p.pop_size);
        if(i == 0 || fitness((bool *) p.matrix[rand_id], p.ind_dim, o) > best_fitness) 
            choosen_ind_id = rand_id;
    }
    return choosen_ind_id;
}

void mutation(bool *number, int number_size)
{
    for (int i = 0; i < number_size; i++)
    {
        bool *b;
        b = number + i;
        // if(1.0 / (double) number_size <= (double) random()/ (double) RAND_MAX)
        if((double) random()/ (double) RAND_MAX <= MUTATION)
            *b = !(*b);
    }
}

void *breed(int ind1, int ind2, population &p)
{
    void *son;
    son = malloc(p.type_size * p.ind_dim);
    for (int i = 0; i < p.ind_dim; i++)
    {
        ((random() % 2) == 0) 
            ? memcpy((char *) son + i * p.type_size, ((char *) *(p.matrix + ind1) + i * p.type_size), p.type_size)
            : memcpy((char *) son + i * p.type_size, ((char *) *(p.matrix + ind2) + i * p.type_size), p.type_size);  
    }
    return son;
}