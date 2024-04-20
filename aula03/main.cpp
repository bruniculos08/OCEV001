#include "logic.h"

int main(void)
{
    srandom(time(0));

    population p;
    p.pop_size = 2;
    p.ind_dim = 100;
    p.type_size = sizeof(bool);

    p.matrix = (void **) generateBoolPop(p.ind_dim, p.pop_size);
    printBoolMatrix((bool **) p.matrix, p.ind_dim, p.pop_size);

    formula sat;
    generate3SAT(sat, 430, p.ind_dim);
    print3SAT(sat);

    int xmen = evolve(sat, p);
    cout << "final best = " << xmen << endl;
    cout << "final best fitness function = " << fitness(sat, xmen, p) << endl;
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

/* --------------------- Funções para modelagem do problema 3-SAT --------------------- */

void generate3SAT(formula &sat, int formula_size, int variables_num)
{
    sat.formumla_size = formula_size;
    sat.expr = (clause *) malloc(sizeof(clause) * formula_size);
    for(int i = 0; i < formula_size; i++)
    {
        for(int j = 0; j < 3; j++) sat.expr[i].variables[j] = random() % variables_num; 
    }
}

int evolve(formula &sat, population &p)
{
    int best_fitness = 0, current_best_id;
    vector<bool *> children;
    do
    {
        int ind1 = tournament(sat, 5, p), ind2 = tournament(sat, 2, p);
        bool *son;
        son = (bool *) breed(ind1, ind2, p);
        mutation(son, p);
        children.push_back(son);

        if(children.size() == p.pop_size)
        {
            changeGeneration(children, p);
            children.clear();
        }

        for (int i = 0; i < p.pop_size; i++)
        {
            if (fitness(sat, i, p) > best_fitness)
            {
                best_fitness = fitness(sat, i, p);
                current_best_id = i;
            }
        }

        cout << "best fitness = " << best_fitness << endl;

    } while (best_fitness != sat.formumla_size);
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

void print3SAT(formula &sat)
{
    for (int i = 0; i < sat.formumla_size; i++)
    {
        cout << "(x" << sat.expr[i].variables[0] << " \\/ x" << sat.expr[i].variables[1]
            << " \\/ x" << sat.expr[i].variables[2] << ")";   
        if(i != sat.formumla_size - 1) cout << " /\\ ";  
    }
    cout << endl;
}

int fitness(formula &sat, int ind_id, population &p)
{
    int count = 0;
    bool **matrix = (bool **) p.matrix;
    for(int i = 0; i < sat.formumla_size; i++)
    {
        clause phi = sat.expr[i];
        bool value = matrix[ind_id][phi.variables[0]] || matrix[ind_id][phi.variables[1]] || matrix[ind_id][phi.variables[2]];
        if(value) count++; 
    }
    return count;
}

int tournament(formula &sat, int samples_num, population &p)
{
    int choosen_ind_id = -1, rand_id, best_fitness = -1;
    for(int i = 0; i < samples_num; i++)
    {
        rand_id = (random() % p.pop_size);
        if(best_fitness == -1 || fitness(sat, rand_id, p) > best_fitness) 
            choosen_ind_id = rand_id;
    }
    return choosen_ind_id;
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

void mutation(bool *dna, population &p)
{
    for (int i = 0; i < p.ind_dim; i++)
    {
        bool *b;
        b = dna + i;
        if(1.0 / (double) p.ind_dim <= (double) random()/ (double) RAND_MAX)
            *b = !(*b);
    }
}