#include "logic.h"

int main(void)
{
    srandom(time(0));

    formula sat;
    fromFileGenerate3SAT(sat, "simple.txt");
    // print3SAT(sat);

    population p;
    p.pop_size = POPULATION_SIZE;
    p.ind_dim = sat.number_of_variables;

    p.matrix = generateBoolPop(p.ind_dim, p.pop_size);
    // printBoolMatrix(p.matrix, p.ind_dim, p.pop_size);

    int xmen = evolve(sat, p);
    cout << "final best = " << xmen << endl;
    cout << "final best fitness function = " << fitness(sat, xmen, p) << endl;
    printBoolMatrix(p.matrix, p.ind_dim, p.pop_size);

    generateGraph();

    return 0;
}

int clampi(int k, int low, int high)
{
    assert(low < high);
    if(k < low) k = low;
    if(k > high) k = high;
    return k;
}

vector<allele> *allocateGenericMatrix(int width, int height, allele value)
{
    vector<allele> *ptr;
    ptr = (vector<allele> *) malloc(sizeof(vector<allele>) * height);
    return ptr;
}

void printIntMatrix(vector<allele> *matrix, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (allele field : matrix[y])
        {
            cout << field.integer << " " ;
        }
        cout << endl;
    }
}

void printBoolMatrix(vector<allele> *matrix, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (allele field : matrix[y])
        {
            cout << field.boolean << " " ;
        }
        cout << endl;
    }
}

void printDoubleMatrix(vector<allele> *matrix, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (allele field : matrix[y])
        {
            cout << field.real << " " ;
        }
        cout << endl;
    }
}

void permute(vector<allele> &dna)
{
    for (int i = 0; i < (int) dna.size(); i++)
    {
        int j = random() % dna.size();
        int aux = dna[i].integer;
        dna[i].integer = dna[j].integer;        
        dna[j].integer = aux;        
    }
}

vector<allele> *generateIntPop(int ind_dim, int pop_size, allele low, allele high)
{
    vector<allele> *matrix;
    matrix = new vector<allele>[pop_size];

    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x++)
        {
            matrix[y][x].integer = (random() % (high.integer - low.integer)) + low.integer;
        }
    }

    return matrix;
}

vector<allele> *generatePermPop(int ind_dim, int pop_size, allele low, allele high)
{
    assert(low.integer < high.integer && (high.integer - low.integer) >= ind_dim);

    vector<allele> *matrix;
    matrix = new vector<allele>[pop_size];
    
    allele temp;
    temp.integer = 0;

    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x ++)
        {
            matrix[y].push_back(temp);
            try_again:
            int num = (random() % (high.integer - low.integer)) + low.integer;
            for (int i = 0; i < x; i++)
            {
                if (matrix[y][i].integer == num) goto try_again;
            }
            matrix[y][x].integer = num;
        }
    }
    return matrix;
}

vector<allele> *generateBoolPop(int ind_dim, int pop_size)
{
    allele temp;
    temp.boolean = true;

    vector<allele> *matrix;
    matrix = new vector<allele>[pop_size];

    for (int y = 0; y < pop_size; y++)
    {
        vector<allele> v;
        matrix[y] = v;

        for (int x = 0; x < ind_dim; x++)
        {
            temp.boolean = (bool) (random() % 2 == 0);
            matrix[y].push_back(temp);
        }
    }
    return matrix;
}

vector<allele> *generateDoublePop(int ind_dim, int pop_size, allele low, allele high)
{
    assert(low.real < high.real);

    vector<allele> *matrix;
    matrix = new vector<allele>[pop_size];

    allele temp;

    for (int y = 0; y < pop_size; y++)
    {
        for (int x = 0; x < ind_dim; x++)
        {
            temp.real = ((double) random() / RAND_MAX) * (high.real - low.real) + low.real;
            matrix[y].push_back(temp);
        }
    }
    return matrix;
}

template <typename T> void eraseFast(vector<T> &v, int index)
{
    int last_index = v.size() - 1;
    v[index] = v[last_index];
    v.pop_back(); 
}

template <typename T> T getRandom(vector<T> &v)
{
    int random_index = random() % v.size();
    T item = v[random_index];
    eraseFast(v, random_index);
    return item;
}

void generateGraph()
{
    Gnuplot gp("gnuplot");
    gp << "set title 'Gráfico de convergência com " << GENERATIONS  << " gerações'\n";
    gp << "set terminal 'svg'\n";
    // gp << "set output \"GraphicFiles/convergencia\"\n";
    gp << "set output \"convergencia.svg\"\n";

    gp << "set xlabel \"geração\"\n";
    gp << "set ylabel \"fitness do melhor individuo\"\n";
    gp << "set grid\n";

    gp << "plot '-' with lines title 'points'\n";
    gp.send1d(data_plot);
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
    int best_fitness = 0, best_id = 0, worst_id = 0, counter = GENERATIONS;
    vector<allele> best_individual;
    best_id = getElite(sat, p, best_individual);
    do
    {
        crossoverAll(sat, p.pop_size, p, CUT_POINTS_CROSSOVER);
        worst_id = getWorst(sat, p);
        p.matrix[worst_id].swap(best_individual);

        if(GENERATIONS > 0) counter--;
        
        best_id = getElite(sat, p, best_individual);
        cout << fitness(sat, best_id, p) << " <- current best result" << endl;
        data_plot.push_back(make_pair((double) (GENERATIONS - counter), (double) fitness(sat, best_id, p)));

    } while (best_fitness < sat.formumla_size && counter != 0);
    return best_id;
}

void changeGeneration(vector<allele> *children, population &p)
{
    for(int i = 0; i < p.pop_size; i++)
        p.matrix[i].clear();
    free(p.matrix);
    p.matrix = children;
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
    vector<allele> *matrix;
    matrix = p.matrix;

    for(int i = 0; i < sat.formumla_size; i++)
    {
        clause phi = sat.expr[i];
        bool value = false;
        value = value || (phi.variables[0] >= 0) ? matrix[ind_id][phi.variables[0] - 1].boolean : !(matrix[ind_id][abs(phi.variables[0]) - 1]).boolean;
        value = value || (phi.variables[1] >= 0) ? matrix[ind_id][phi.variables[1] - 1].boolean : !(matrix[ind_id][abs(phi.variables[1]) - 1]).boolean;
        value = value || (phi.variables[2] >= 0) ? matrix[ind_id][phi.variables[2] - 1].boolean : !(matrix[ind_id][abs(phi.variables[2]) - 1]).boolean;
        if(value) count++; 
    }
    return count;
}

int tournament(formula &sat, int samples_num, population &p, vector<int> &single_ind)
{
    int choosen_ind_id = -1, rand_id, best_fitness = -1;
    for(int i = 0; i < min(samples_num, (int) single_ind.size()); i++)

    if(samples_num <= 0) samples_num = single_ind.size();
    
    {
        rand_id = getRandom<int>(single_ind);
        if(best_fitness == -1 || fitness(sat, rand_id, p) > best_fitness) 
            choosen_ind_id = rand_id;        
        else single_ind.push_back(rand_id);
    }
    return choosen_ind_id;
}

int roulette(vector<int> &single_ind, formula &sat, population &p)
{
    int summ = 0;
    for(int index : single_ind)
        summ += fitness(sat, index, p);
    
    double prob; 
    vector<double> roulette_chances;
    for(int index : single_ind)
    {
        prob = 100.0 * (double) fitness(sat, index, p) / (double) summ;
        roulette_chances.push_back(prob);
    }

    double winner = 100.0 * (double) ((double) random() / (double) RAND_MAX);
    double acc = 0.0;
    int result = 0;
    for(int i = 0; i < (int) roulette_chances.size(); i++)
    {
        if(acc <= winner && winner <= acc + roulette_chances[i])
        {
            result = single_ind[i];
            eraseFast(single_ind, i);
            // cout << "Individuo selecionado pela roleta = " 
            //      << single_ind[i] << " (fitness = " << fitness(sat, result, p) << ")" << endl;
            return result;
        }
        acc += roulette_chances[i];
    }
    return 0;
}

void crossoverAll(formula &sat, int samples_num, population &p, int cuts_num)
{
    vector<int> single_ind;
    for(int i = 0; i < p.pop_size; i++) single_ind.push_back(i);

    int parent1_index, parent2_index;
    while(single_ind.size() >= 2)
    {
        // parent1_index = tournament(sat, samples_num, p, single_ind);
        // parent2_index = tournament(sat, samples_num, p, single_ind);
        parent1_index = roulette(single_ind, sat, p);
        parent2_index = roulette(single_ind, sat, p);
        crossover(parent1_index, parent2_index, p, cuts_num);
    }
}

void crossover(int ind1, int ind2, population &p, int cuts_num)
{
    vector<allele> son1;
    vector<allele> son2;

    vector<allele> *parent1;
    vector<allele> *parent2;
    parent1 = p.matrix + ind1;
    parent2 = p.matrix + ind2;

    bool flag = false;

    if(cuts_num <= 0 || cuts_num >= p.ind_dim - 1)
        for (int i = 0; i < p.ind_dim; i++)
        {
            if(random() % 2 == 0)
            {
                son1.push_back((*parent1)[i]);
                son2.push_back((*parent2)[i]);
            }
            else
            {
                son1.push_back((*parent2)[i]);
                son2.push_back((*parent1)[i]);
            }      
        }
    else
        for (int i = 0, cut_point = 0, old_cut_point = 0; i <= cuts_num; i++)
        {
            (i == cuts_num)
            ? cut_point = p.ind_dim
            : cut_point += random() % ((int) ceil((float) p.ind_dim / (float) cuts_num)) + 1;
            flag = !flag;
            for(int j = old_cut_point; j < cut_point; j++)
            {
                if(flag) 
                {
                    son1.push_back((*parent1)[j]);
                    son2.push_back((*parent2)[j]);
                }
                else
                {
                    son1.push_back((*parent2)[j]);
                    son2.push_back((*parent1)[j]);
                }
            }
            old_cut_point = cut_point;
        }
    mutation(son1, p);
    mutation(son2, p);
    (*parent1).swap(son1); 
    (*parent2).swap(son2); 

}

void mutation(vector<allele> &dna, population &p)
{
    for (int i = 0; i < p.ind_dim; i++)
    {
        // if(1.0 / (double) p.ind_dim <= (double) random()/ (double) RAND_MAX)
        if((double) random()/ (double) RAND_MAX <= MUTATION_PROBABILITY)
            dna[i].boolean = !dna[i].boolean;
    }
}

int getElite(formula &sat, population &p, vector<allele> &elite_place)
{
    int best_fitness = fitness(sat, 0, p), current_best_id = 0;
    for (int i = 0; i < p.pop_size; i++)
        {
            if (fitness(sat, i, p) > best_fitness)
            {
                best_fitness = fitness(sat, i, p);
                current_best_id = i;
            }
        }
    elite_place = p.matrix[current_best_id];
    return current_best_id;
}

int getWorst(formula &sat, population &p)
{
    int worst_fitness = fitness(sat, 0, p), current_worst_id = 0;
    for (int i = 0; i < p.pop_size; i++)
        {
            if (fitness(sat, i, p) < worst_fitness)
            {
                worst_fitness = fitness(sat, i, p);
                current_worst_id = i;
            }
        }
    return current_worst_id;
}
