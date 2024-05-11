#include "logic.h"

int main(void)
{
    srandom(time(0));

    //  Criar o tabuleiro:
    nqueens b;
    b.dim = 64;
    b.board = createBoard(b.dim, b.dim);

    //  Gerar a população:
    allele low, high;
    low.integer = 0;
    high.integer = b.dim;
    population p;
    p.ind_dim = b.dim;
    p.pop_size = POPULATION_SIZE;
    generatePermPop(b.dim, p.pop_size, low, high, p.matrix);

    int result = evolve(b, p);
    cout << "best ind = " << result << endl;
    cout << "best ind fitness = " << fitness(b, result, p) << endl;
    generateGraph();

    // printIntMatrix(p.matrix);
    cout << "Resultado: " << endl;
    // printBoard(b, p.matrix[result]);
    return 0;
}

void allocateGenericMatrix(int width, int height, allele value, vector<vector<allele>> &v)
{
    vector<allele> line(width, value);
    vector<vector<allele>> temp(height, line);
    v.swap(temp);
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

void generatePermPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v)
{
    assert(low.integer < high.integer && (high.integer - low.integer) >= ind_dim);
    vector<allele> line(ind_dim, low);
    vector<vector<allele>> matrix(pop_size, line);
    
    // std::random_device rd;
    // std::mt19937 gen(rd());
    
    vector<int> poll(high.integer - low.integer), temp_poll;
    for(int i = 0; i < (int) poll.size(); i++) poll[i] = low.integer + i;
    
    int rand_index;
    for (int y = 0; y < pop_size; y++)
    {
        temp_poll = poll;
        for (int x = 0; x < ind_dim; x ++)
        {
            rand_index = random() % temp_poll.size();
            matrix[y][x].integer = temp_poll[rand_index];
            eraseFast(temp_poll, rand_index);
        }
        // shuffle(begin(matrix[y]), end(matrix[y]),  gen);
        permute(matrix[y]);
    }
    v.swap(matrix);
}

void generateIntPop(int ind_dim, int pop_size, allele low, allele high, vector<vector<allele>> &v)
{
    vector<allele> line(ind_dim, low);
    vector<vector<allele>> matrix(pop_size, line);
    
    for(vector<allele> &row : matrix)
        for(allele &num : row)
            num.integer = (low.integer + (random() % (high.integer - low.integer)));
}

void printIntMatrix(vector<vector<allele>> &matrix){
    for(vector<allele> &row : matrix){
        for(allele &num : row){
            cout << num.integer << " ";
        }
        cout << endl;
    }
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

double randomDouble(double max, double min)
{
    return (((double) random()) / ((double) RAND_MAX)) * (max - min) + min;
}

int linearSearch(vector<allele> &v, allele e, int start, int end)
{
    assert((size_t) end < v.size());
    for (int i = start; i <= end; i++) if(v[i].integer == e.integer) return i;
    return -1;
}

/* ----------------- Funções para modelagem do problema N-Rainhas --------------- */

int **createBoard(int width, int height)
{
    int **matrix;
    matrix = new int*[height];
    for(int i = 0; i < height; i++)
    {
        matrix[i] = new int[width];
        for(int j = 0; j < width; j++)
            matrix[i][j] = 0;
    }
    return matrix;
}

int evolve(nqueens &b, population &p)
{
    calcEsc(b, p, SCALING_CONSTANT);
    // double best_fitness = 0.0;
    int best_id = 0, worst_id = 0, counter = GENERATIONS;
    vector<allele> best_individual;
    best_id = getElite(b, p, best_individual);
    cout << "Started evolving" << endl;
    do
    {
        calcEsc(b, p, SCALING_CONSTANT);
        crossoverAll(b, p.pop_size, p, CUT_POINTS_CROSSOVER);
        worst_id = getWorst(b, p);
        // cout << fitnessEsc(b, worst_id, p) << " <- old (worst) result" << endl;
        if(counter % 200 == 0) 
            cout << fitness(b, worst_id, p) << " <- old (worst) result" << endl;
        p.matrix[worst_id] = best_individual;


        best_id = getElite(b, p, best_individual);
        // best_fitness = fitnessEsc(b, best_id, p); 

        // cout << best_fitness << " <- current best result" << endl;
        if(counter % 200 == 0)
            cout << fitness(b, best_id, p) << " <- current best result" << endl;
        
        if(GENERATIONS > 0) counter--;
        data_plot.push_back(make_pair((double) (GENERATIONS - counter), (double) fitness(b, best_id, p)));

    } while (fitness(b, best_id, p) != (b.dim * (b.dim - 1)) && counter != 0);
    return best_id;
}

void changeGeneration(vector<vector<allele>> &children, population &p)
{
    for(int i = 0; i < p.pop_size; i++)
        p.matrix[i] = children[i];
    children.clear();
}

void printBoard(nqueens &b, vector<allele> &ind)
{
    b.board = createBoard(b.dim, b.dim);
    for(allele &queen : ind)
    {
        pair<int, int> pos = getPosition(b, queen, ind);
        b.board[pos.first][pos.second] = 1;
    }
    for(int i = 0; i < b.dim; i++)
    { 
        for(int j = 0; j < b.dim; j++) cout << b.board[i][j] << " ";
        cout << endl;
    }
}

int fitness(nqueens &b, int ind_id, population &p)
{
    vector<allele> pieces;
    pieces = p.matrix[ind_id];

    int conflicts = 0;
    for(allele &queen : pieces){
        // conflicts += lookRowAndCol(pieces, b, queen);
        conflicts += lookDiag(pieces, b, queen);
    }
    return (b.dim * (b.dim - 1)) - conflicts;
}

int tournament(nqueens &b, int samples_num, population &p, vector<int> &single_ind)
{
    vector<int> poll_rest;
    int choosen_ind_id = -1, rand_id;
    double best_fitness = -1.0;
    if(samples_num <= 0 || samples_num >= (int) single_ind.size()) 
        samples_num = single_ind.size();

    for(int i = 0; i < samples_num; i++)
    {
        rand_id = getRandom<int>(single_ind);
        if(choosen_ind_id == -1 || fitnessEsc(b, rand_id, p) > best_fitness)
        {
            if(choosen_ind_id != -1) poll_rest.push_back(choosen_ind_id);
            choosen_ind_id = rand_id;        
            best_fitness = fitnessEsc(b, choosen_ind_id, p);
        }
        else poll_rest.push_back(rand_id);
    }

    for(int id : poll_rest) single_ind.push_back(id);
    poll_rest.clear(); 
    return choosen_ind_id;
}

int roulette(vector<int> &single_ind, nqueens &b, population &p)
{
    double summ = 0.0;
    for(int index : single_ind){
        // pthread_mutex_lock(&mtx_crossover[index]);
        summ += (double) fitnessEsc(b, index, p);
        // pthread_mutex_unlock(&mtx_crossover[index]);
    }

    if(summ == 0.0){
        int rand_id = random() % single_ind.size();
        int result = single_ind[rand_id];
        eraseFast(single_ind, rand_id);
        return rand_id;
    }

    double prob; 
    vector<double> roulette_chances;
    for(int index : single_ind)
    {
        // pthread_mutex_lock(&mtx_crossover[index]);
        prob = 100.0 * fitnessEsc(b, index, p) / summ;
        // pthread_mutex_unlock(&mtx_crossover[index]);
        roulette_chances.push_back(prob);
    }

    double winner = 100.0 * (double) ((double) random() / (double) RAND_MAX), acc = 0.0;
    int result = 0;
    for(int i = 0; i < (int) roulette_chances.size(); i++)
    {
        if(acc <= winner && winner <= acc + roulette_chances[i])
        {
            result = single_ind[i];
            eraseFast(single_ind, i);
            return result;
        }
        acc += roulette_chances[i];
    }
    return 0;
}

void crossoverAll(nqueens &b, int samples_num, population &p, int cuts_num)
{
    vector<int> single_ind;
    for(int i = 0; i < p.pop_size; i++) single_ind.push_back(i);
    int thread_counter = 0;

    int parent1_index, parent2_index;
    while(single_ind.size() >= 2)
    { 
        args *t_arg;
        t_arg = new args;
        if(SELECTION_OPT == true)
        {
            parent1_index = tournament(b, samples_num, p, single_ind);
            parent2_index = tournament(b, samples_num, p, single_ind);
        }
        else
        {
            parent1_index = roulette(single_ind, b, p);
            parent2_index = roulette(single_ind, b, p);
        }
        t_arg->ind1 = parent1_index;
        t_arg->ind2 = parent2_index;
        t_arg->cuts_num = cuts_num;
        t_arg->p = &p;
        if(CROSSOVER_OPT)
            pthread_create(&threads[thread_counter], NULL, threadCrossoverPMX, (void *) t_arg);
            // crossoverPMX(parent1_index, parent2_index, p, cuts_num);
        else
            pthread_create(&threads[thread_counter], NULL, threadCrossoverCX, (void *) t_arg);
            // crossoverCX(parent1_index, parent2_index, p, cuts_num);
        thread_counter++;
    }
    for(int i = 0; i < thread_counter; i++)
    {
        pthread_join(threads[i], NULL);
    }
}

void crossoverPMX(int ind1, int ind2, population &p, int cut_size)
{
    allele none;
    none.integer = -1;
    
    vector<allele> f1 = p.matrix[ind1];
    vector<allele> f2 = p.matrix[ind2];
    vector<allele> s1(p.ind_dim, none);
    vector<allele> s2(p.ind_dim, none);

    int cut1 = random() % p.ind_dim, cut2 = random() % p.ind_dim;
    if(cut_size > 0) cut2 = min(cut1 + cut_size - 1, p.ind_dim - 1);

    // Colocando segmentos nos filhos:
    for(int i = cut1; i <= cut2; i++)
    {
        // (1) Copiar segmento do primeiro pai para o primeiro filho:
        s1[i] = f1[i];
        // (2) Copiar segmento do segundo pai para o segundo filho:
        s2[i] = f2[i];
    }
    // Colocando alelos por mapeamento:
    for(int i = cut1; i <= cut2; i++)
    {
        // (1) Mapeamento do segundo pai para o primeiro filho:
        if(linearSearch(s1, f2[i], cut1, cut2) == -1)
        {
            allele n = s1[i];
            allele m = f2[i];
            int n_pos_f2 = linearSearch(f2, n, 0, f2.size() - 1);
            if(s1[n_pos_f2].integer == -1)
                s1[n_pos_f2] = m;
            else{
                do{
                    allele k = s1[n_pos_f2];
                    int k_pos_f2 = linearSearch(f2, k, 0, f2.size() - 1);
                    if(s1[k_pos_f2].integer == -1)
                    {
                        s1[k_pos_f2] = m;
                        break;
                    }
                    n_pos_f2 = k_pos_f2;
                } while (true);
            }
        }
        // (2) Mapeamento do primeiro pai para o segundo filho:
        if(linearSearch(s2, f1[i], cut1, cut2) == -1)
        {
            allele n = s2[i];
            allele m = f1[i];
            int n_pos_f1 = linearSearch(f1, n, 0, f1.size() - 1);
            if(s2[n_pos_f1].integer == -1)
                s2[n_pos_f1] = m;
            else{
                do{
                    allele k = s2[n_pos_f1];
                    int k_pos_f1 = linearSearch(f1, k, 0, f1.size() - 1);
                    if(s2[k_pos_f1].integer == -1)
                    {
                        s2[k_pos_f1] = m;
                        break;
                    }
                    n_pos_f1 = k_pos_f1;
                } while (true);
            }
        }
    } 

    // (3) Completar os vetores:
    int index_s1 = 0, index_s2 = 0;
    for (int i = 0; i < p.ind_dim; i++)
    {
        if(linearSearch(s1, f2[i], 0, p.ind_dim - 1) == -1)
        {
            index_s1 = linearSearch(s1, none, 0, p.ind_dim - 1);
            s1[index_s1] = f2[i];
        }
        if(linearSearch(s2, f1[i], 0, p.ind_dim - 1) == -1)
        {
            index_s2 = linearSearch(s2, none, 0, p.ind_dim - 1);
            s2[index_s2] = f1[i];
        }
    }

    mutation(s1, p);
    mutation(s2, p);

    // Ponto crítico (para dados do individuo 1):
    // pthread_mutex_lock(&mtx_crossover[ind1]);
    p.matrix[ind1].swap(s1);
    // pthread_mutex_unlock(&mtx_crossover[ind1]);

    // Ponto crítico (para dados do individuo 2):
    // pthread_mutex_lock(&mtx_crossover[ind2]);
    p.matrix[ind2].swap(s2);
    // pthread_mutex_unlock(&mtx_crossover[ind2]);
}

void *threadCrossoverPMX(void *data)
{
    args *A;
    A = (args *) data;
    crossoverPMX(A->ind1, A->ind2, *(A->p), A->cuts_num);
    pthread_exit(NULL);
}

void crossoverCX(int ind1, int ind2, population &p, int cuts_num)
{
    allele none;
    none.integer = -1;
    allele start_value;
    int temp_pos;
    vector<allele> s1(p.ind_dim, none);
    vector<allele> s2(p.ind_dim, none);
    vector<allele> p1 = p.matrix[ind1];
    vector<allele> p2 = p.matrix[ind2];
    vector<allele> aux;
    for(int i = 0; i < p.ind_dim; i++)
    {
        // (1) Ciclo sobre o primeiro filho:
        if(s1[i].integer == -1)
        {
            temp_pos = i;
            start_value = p1[temp_pos];
            do
            {
                s1[temp_pos] = p1[temp_pos];
                temp_pos = linearSearch(p1, p2[temp_pos], 0, p.ind_dim - 1);
            } while (p1[temp_pos].integer != start_value.integer);      
        }
        // (2) Ciclo sobre o segundo filho:
        if(s2[i].integer == -1)
        {
            temp_pos = i;
            start_value = p2[temp_pos];
            do
            {
                s2[temp_pos] = p2[temp_pos];
                temp_pos = linearSearch(p2, p1[temp_pos], 0, p.ind_dim - 1);
            } while (p2[temp_pos].integer != start_value.integer);      
        }
        aux = p1;
        p1 = p2;
        p2 = aux;        
    }

    mutation(s1, p);
    mutation(s2, p);

    // pthread_mutex_lock(&mtx_crossover[ind1]);
    p.matrix[ind1] = s1;
    // pthread_mutex_unlock(&mtx_crossover[ind1]);

    // pthread_mutex_lock(&mtx_crossover[ind2]);
    p.matrix[ind2] = s2;
    // pthread_mutex_unlock(&mtx_crossover[ind2]);
}

void *threadCrossoverCX(void *data)
{
    args *A;
    A = (args *) data;
    crossoverCX(A->ind1, A->ind2, *(A->p), A->cuts_num);
    pthread_exit(NULL);
}

void mutation(vector<allele> &dna, population &p){
    vector<int> poll(p.ind_dim - 1);
    for(size_t i = 0; i < poll.size(); i++) poll[i] = (int) i;
    for (int i = 0; i < p.ind_dim; i++)
    {
        if(randomDouble(0.0, 1.0) <= MUTATION_PROBABILITY)
            mutateAllele(dna, p, i, poll);
    }
}

void mutateAllele(vector<allele> &dna, population &p, int allele_index, vector<int> &poll)
{
    int random_index = random() % poll.size();
    if(random_index >= allele_index)
        random_index += 1;
    allele aux = dna[allele_index];
    dna[allele_index] = dna[random_index];
    dna[random_index] = aux;
}

int getElite(nqueens &b, population &p, vector<allele> &elite_place)
{
    double best_fitness = fitnessEsc(b, 0, p);
    int current_best_id = 0;
    for (int i = 0; i < p.pop_size; i++)
        {
            if (fitnessEsc(b, i, p) > best_fitness)
            {
                best_fitness = fitnessEsc(b, i, p);
                current_best_id = i;
            }
        }
    elite_place = p.matrix[current_best_id];
    return current_best_id;
}

int getWorst(nqueens &b, population &p)
{
    double worst_fitness = fitnessEsc(b, 0, p), current_worst_id = 0;
    for (int i = 0; i < p.pop_size; i++)
        {
            if (fitnessEsc(b, i, p) < worst_fitness)
            {
                worst_fitness = fitnessEsc(b, i, p);
                current_worst_id = i;
            }
        }
    return current_worst_id;
}

pair<int, int> getPosition(nqueens &b, allele queen, vector<allele> &v)
{
    int col = queen.integer;
    int row = linearSearch(v, queen, 0, v.size() - 1);
    return make_pair(row, col);
}

int lookDiag(vector<allele> &ind, nqueens &b, allele actual_queen)
{
    pair<int, int> coord_ind, coord_temp;
    coord_ind = getPosition(b, actual_queen, ind);
    //  Obs.: chess_piece é o valor de inteiro do alelo.
    
    //  As diagonais que passam pela rainha que está sendo avaliada são as
    //  duas retas que satisfazem as equações de reta (com inclinação de 45º) 
    //  montadas a partir das coordenadas da rainhas:
    //  B = y_r - x_r, com A = 1 (Ax + B = y)
    int B_growing = coord_ind.first - coord_ind.second;
    //  B = y_r + x_r, com A = -1 (Ax + B = y)
    int B_decreasing = coord_ind.first + coord_ind.second;

    int conflicts = 0;
    for(allele &other_queen : ind)
    {   
        if(other_queen.integer == actual_queen.integer)
            continue;
        coord_temp = getPosition(b, other_queen, ind);
        if(coord_temp.first - coord_temp.second == B_growing || 
                coord_temp.first + coord_temp.second == B_decreasing)
                conflicts++;
    }
    return conflicts;
}

int lookRowAndCol(vector<allele> &ind, nqueens &b, allele actual_queen)
{
    pair<int, int> coord_ind, coord_temp;
    coord_ind = getPosition(b, actual_queen, ind);
    int conflicts = 0;
    for(allele &other_queen : ind)
    {   
        if(other_queen.integer == actual_queen.integer)
            continue;
        coord_temp = getPosition(b, other_queen, ind);
        if(coord_temp.first == coord_ind.first || coord_temp.second == coord_ind.second)
            conflicts++;
    }
    return conflicts;
}

// Funções para implementação de escalonamento linear:
//  - Média da função fitness:
double fitnessAvg(nqueens &b, population &p)
{
    double sum = 0.0;
    for(int i = 0; i < p.pop_size; i++)
        sum += (double) fitness(b, i, p);
    sum /= (double) p.pop_size;
    return sum;
}
//  - Valor máximo da função fitness:
double fitnessMin(nqueens &b, population &p)
{
    double m = (double) fitness(b, 0, p);
    for(int i = 0; i < p.pop_size; i++)
        m = min(m, (double) fitness(b, i, p)); 
    return m;
}
//  - Valor mínimo da função fitness:
double fitnessMax(nqueens &b, population &p)
{
    double m = (double) fitness(b, 0, p);
    for(int i = 0; i < p.pop_size; i++)
        m = max(m, (double) fitness(b, i, p)); 
    return m;
}
//  - Calculo dos valores alpha e beta:
void calcEsc(nqueens &b, population &p, double c)
{
    double f_min = fitnessMin(b, p),
    f_max = fitnessMax(b, p), f_avg = fitnessAvg(b, p), d;

    if(f_min > (c * f_avg - f_max)/(c - 1.0))
    {
        d = f_max - f_avg;
        p.alpha = f_avg * (c - 1) / d;
        p.beta = f_avg * (f_max - c * f_avg) / d;
    }
    else
    {
        d = f_max - f_min;
        p.alpha = f_avg / d;
        p.beta = (-f_min) * f_avg / d;
    }

    if(abs(d) <= 0.00001)
    {
        p.alpha = 1.0;
        p.beta = 0.0;
    }
}
//  - Função fitness com escalonamento:
double fitnessEsc(nqueens &b, int ind_id, population &p)
{
    return max(0.0, p.alpha * ((double) fitness(b, ind_id, p)) + p.beta);
}