#include "logic.h"

int main(void)
{
    srandom(time(0));

    function_info info;
    info.f = f;
    info.a = -2.0;
    info.b = 2.0;
    info.number_size = 11;
    info.decimal_houses = 5;
    C_min = findCMin(info, 100);

    // É necessário garantir que o tamanho do intervalo seja representável com o número de bits escolhido:
    assert(pow(2.0, (double) info.number_size - 1) / pow(2.0, (double) info.decimal_houses - 1) >= max(fabs(info.a), fabs(info.b)));

    population p;
    p.pop_size = POP_SIZE;
    p.ind_dim = info.number_size;
    p.type_size = sizeof(bool);
    p.matrix = (void **) generateBoolPop(p.ind_dim, p.pop_size);

    int lisan_al_gaib = evolve(info, p);
    cout << "Best fitness (individual " << lisan_al_gaib << ") = " 
    << fitness((bool *) p.matrix[lisan_al_gaib], info) - C_min << endl;

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

/* --------------------- Funções para modelagem do problema de otimização de função real --------------------- */

double f(double x)
{
    return cosf64(20 * x) - (fabs(x)/2.0) + ((x * x * x)/4.0);
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

double fitness(bool *number, function_info &info)
{
    double x = 0, exp;
    bool signal = !number[0];
    int integer_houses = (info.number_size - 1) - info.decimal_houses; 
    for (int i = 0; i < info.number_size - 1; i++)
    {
        (number[info.number_size - (i + 1)]) ? exp = (double) i : exp = 0.0; 
        (signal) ? x += pow(2.0, exp) : x -= pow(2.0, exp); 
    }
    x /= pow(2.0, (double) info.decimal_houses);

    // cout << "x = " << x << endl;
    x = normalize(info, x);
    // cout << "x_norm = " << x << endl;
    return info.f(x) + C_min;
}

double normalize(function_info &info, double x)
{
    double max_value, min_value, x_min = info.a, x_max = info.b;
    
    max_value = pow(2.0, (double) (info.number_size) - 1) - 1;
    max_value /= pow(2.0, (double) info.decimal_houses);
    min_value = - max_value;

    return (x - min_value)/(max_value - min_value) * (info.b - info.a) + info.a;
}

double findCMin(function_info &info, int sample_size)
{
    double x, fx, fx_min = info.f(0);
    for (int i = 0; i < sample_size; i++)
    {
        x = ((double) random()/ (double) RAND_MAX) * (info.b - info.a) + info.a;
        fx = info.f(x);
        if(fx < fx_min) fx_min = fx;
    }
    return (fx_min < 0) ? fabs(fx_min) : 0;
}

int tournament(int sample_size, population &p, function_info &info)
{
    int choosen_ind_id = -1, rand_id, best_fitness = -1;
    for(int i = 0; i < sample_size; i++)
    {
        rand_id = (random() % p.pop_size);
        if(choosen_ind_id == -1 || fitness((bool *) p.matrix[rand_id], info) > best_fitness) 
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

int evolve(function_info &info, population &p)
{
    int counter = 0, current_best_id = 0;
    vector<bool *> children;
    do
    {
        int ind1 = tournament(2, p, info), ind2 = tournament(2, p, info);
        bool *son;
        son = (bool *) breed(ind1, ind2, p);
        mutation(son, p);
        children.push_back(son);
        
        if(children.size() == p.pop_size)
        {
            changeGeneration(children, p);
            children.clear();
            counter++;
            cout << "Best fitness (individual " << current_best_id << ") = " 
            << fitness((bool *) p.matrix[current_best_id], info) - C_min << endl;
        }

        for (int i = 0; i < p.pop_size; i++)
        {
            if (fitness((bool *) p.matrix[i], info) > fitness((bool *) p.matrix[current_best_id], info))
                current_best_id = i;
        }
        
    } while (counter != GENERATIONS);
    return current_best_id;
}