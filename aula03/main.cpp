#include "logic.h"

int main(void)
{
    srandom(time(0));

    formula sat;
    fromFileGenerate3SAT(sat, "data.txt");
    print3SAT(sat);

    population p;
    p.pop_size = 40;
    p.ind_dim = sat.number_of_variables;
    p.type_size = sizeof(bool);

    p.matrix = (void **) generateBoolPop(p.ind_dim, p.pop_size);
    printBoolMatrix((bool **) p.matrix, p.ind_dim, p.pop_size);

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
    srand(time(0));
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
        for(int j = 0; j < 3; j++) sat.expr[i].variables[j] = (random() % variables_num) + 1; 
    }
}

void fromFileGenerate3SAT(formula &sat, string file_path)
{
    ifstream my_file (file_path); 
    string line;
    char *cstr = new char[1024];

    assert(my_file.is_open());

    getline(my_file, line);
    strcpy(cstr, line.c_str());

    int number_of_variables, number_of_clauses;
    sscanf(cstr, "%*s %*s %i %i\n", &number_of_variables, &number_of_clauses);

    sat.formumla_size = number_of_clauses;
    sat.number_of_variables = number_of_variables;
    sat.expr = (clause *) malloc(sizeof(clause) * number_of_clauses);

    for(int i = 0; i < number_of_clauses && getline(my_file, line); i++)
    {  
        strcpy(cstr, line.c_str());
        sscanf(cstr, "%d %d %d %*s", &sat.expr[i].variables[0], &sat.expr[i].variables[1], &sat.expr[i].variables[2]);
    }
}

int evolve(formula &sat, population &p)
{
    int best_fitness = 0, current_best_id;
    vector<bool *> children;
    do
    {
        int ind1 = tournament(sat, 10, p, -1), ind2 = tournament(sat, 10, p, ind1);
        bool *son;
        son = (bool *) crossover(ind1, ind2, p);
        mutation(son, p);
        children.push_back(son);

        if(children.size() == p.pop_size)
        {
            changeGeneration(children, p);
            cout << "last generation best fitness = " << fitness(sat, current_best_id, p) << endl;
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
        cout << "(";
        (sat.expr[i].variables[0] > 0) ? cout << "x" << sat.expr[i].variables[0] : cout << "¬x" << abs(sat.expr[i].variables[0]);
        (sat.expr[i].variables[1] > 0) ? cout << " \\/ x" << sat.expr[i].variables[1] : cout << " \\/ ¬x" << abs(sat.expr[i].variables[1]);
        (sat.expr[i].variables[2] > 0) ? cout << " \\/ x" << sat.expr[i].variables[2] : cout << " \\/ ¬x" << abs(sat.expr[i].variables[2]);
        cout << ")";   
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
        bool value = false;
        value = value || (phi.variables[0] >= 0) ? matrix[ind_id][phi.variables[0] - 1] : !(matrix[ind_id][abs(phi.variables[0]) - 1]);
        value = value || (phi.variables[1] >= 0) ? matrix[ind_id][phi.variables[1] - 1] : !(matrix[ind_id][abs(phi.variables[1]) - 1]);
        value = value || (phi.variables[2] >= 0) ? matrix[ind_id][phi.variables[2] - 1] : !(matrix[ind_id][abs(phi.variables[2]) - 1]);
        if(value) count++; 
    }
    return count;
}

int tournament(formula &sat, int samples_num, population &p, int exclude_id)
{
    int choosen_ind_id = -1, rand_id, best_fitness = -1;
    for(int i = 0; i < samples_num; i++)
    {
        rand_id = (random() % p.pop_size);
        if(rand_id == exclude_id)
            if(exclude_id == p.pop_size - 1) rand_id--;
            else rand_id++; 

        if(best_fitness == -1 || fitness(sat, rand_id, p) > best_fitness) 
            choosen_ind_id = rand_id;        
    }
    return choosen_ind_id;
}

void *crossover(int ind1, int ind2, population &p)
{
    void *son;
    son = malloc(p.type_size * p.ind_dim);
    for (int i = 0; i < p.ind_dim; i++)
    {
        (random() % 2 == 0)
            ? memcpy((char *) son + i * p.type_size, ((char *) *(p.matrix + ind1) + i * p.type_size), p.type_size)
            : memcpy((char *) son + i * p.type_size, ((char *) *(p.matrix + ind2) + i * p.type_size), p.type_size);  
    }
    return son;
}

void mutation(bool *dna, population &p)
{
    for (int i = 0; i < p.ind_dim; i++)
    {
        if(1.0 / (double) p.ind_dim <= (double) random()/ (double) RAND_MAX)
            dna[i] = !dna[i];
    }
}