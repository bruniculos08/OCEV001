#include "logic.h"

int main(void)
{
    srandom(time(0));

    //  Problema das N-Rainhas:
    // nqueens b;
    // b.dim = 128;
    // for(int i = 0; i < 10; i++)
    // {
    //     //  Criar o tabuleiro:
    //     b.board = createBoard(b.dim, b.dim);
    //     allele low, high;
    //     low.integer = 0;
    //     high.integer = b.dim;
    //     //  Gerar a população:
    //     population p;
    //     p.ind_dim = b.dim;
    //     p.pop_size = POPULATION_SIZE;
    //     generatePermPop(b.dim, p.pop_size, low, high, p.matrix);
    //     if(i == 0)
    //         worst_value = (double) fitnessMin(b, p);
    //     // printIntMatrix(p.matrix);
    //     int result = evolve(b, p);
    //     cout << "Resultado: " << fitness(b, result, p) << endl;
    //     // printBoard(b, p.matrix[result]);
    //     bests_buffer.push_back(best_plot);
    //     best_plot.clear();
    //     averages_buffer.push_back(average_plot);
    //     average_plot.clear();
    //     worsts_buffer.push_back(worst_plot);
    //     worst_plot.clear();
    // }
    // generateGraph(bests_buffer, "nqueens" + to_string(b.dim) + "-bests", 0);
    // generateGraph(averages_buffer, "nqueens" + to_string(b.dim) + "-averages", 1);
    // generateGraph(worsts_buffer, "nqueens" + to_string(b.dim) + "-worsts", 2);

    // Problema de programação linear:
    linearProg t;
    t.decimal_places =  0;
    t.integers_places = 6;
    t.obj_coeffs = {30.0, 40.0};
    t.intervals = {make_pair(0.0, 24.0), make_pair(0.0, 16.0)};
    t.leq_restr = {{1.0, 2.0, 40.0}};
    for(int i = 0; i < 10; i++)
    {
        population p;
        p.ind_dim = (t.integers_places + t.decimal_places) * (int) t.intervals.size();
        p.pop_size = POPULATION_SIZE;
        generateBoolPop(p.ind_dim, p.pop_size, p.matrix);
        if(i == 0)
            worst_value = (double) fitnessMin(t, p);
        int result = evolve(t, p);
        cout << "best ind = " << result << endl;
        cout << "best ind fitness = " << fitness(t, result, p) << endl;
        printLinearProg(t, p.matrix[result]);
        bests_buffer.push_back(best_plot);
        best_plot.clear();
        averages_buffer.push_back(average_plot);
        average_plot.clear();
        worsts_buffer.push_back(worst_plot);
        worst_plot.clear();
    }
    generateGraph(bests_buffer, "linearProg-bests", 0);
    generateGraph(averages_buffer, "linearProg-averages", 1);
    generateGraph(worsts_buffer, "linearProg-worsts", 2);

    // Problema 3-SAT:
    // formula sat;
    // fromFileGenerate3SAT(sat, "sat.txt");
    // for(int i = 0; i < 10; i++)
    // {
    //     population p;
    //     p.pop_size = POPULATION_SIZE;
    //     p.ind_dim = sat.number_of_variables;
    //     generateBoolPop(p.ind_dim, p.pop_size, p.matrix);
    //     if(i == 0)
    //         worst_value = (double) fitnessMin(sat, p);
    //     int result = evolve(sat, p);
    //     cout << "best ind = " << result << endl;
    //     cout << "best ind fitness = " << fitness(sat, result, p) << endl;
    //     bests_buffer.push_back(best_plot);
    //     best_plot.clear();
    //     averages_buffer.push_back(average_plot);
    //     average_plot.clear();
    //     worsts_buffer.push_back(worst_plot);
    //     worst_plot.clear();
    // }
    // generateGraph(bests_buffer, "3SAT-bests", 0);
    // generateGraph(averages_buffer, "3SAT-averages", 1);
    // generateGraph(worsts_buffer, "3SAT-worsts", 2);

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

void generateBoolPop(int ind_dim, int pop_size, vector<vector<allele>> &v)
{
    vector<allele> line(ind_dim);
    vector<vector<allele>> matrix(pop_size, line);
    for(int y = 0; y < pop_size; y++)
        for(int x = 0; x < ind_dim; x++)
            matrix[y][x].boolean = (bool) (random() % 2);
    v.swap(matrix);
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

void generateGraph(vector<vector<pair<double, double>>> data_plot, string name, int opt)
{
    regex file_pattern{"convergencia_" + name + "[0-9]+.svg"};
    DIR *dir_ptr;
    struct dirent *entry;
    dir_ptr = opendir("Imagens");
    int file_counter = 1;
    while ((entry = readdir(dir_ptr)) != NULL)
    {
        if(entry->d_type == DT_REG && regex_match(entry->d_name, file_pattern))
            file_counter++;
    }
    closedir(dir_ptr);

    Gnuplot gp("gnuplot");
    gp << fixed;
    gp << setprecision(4);
    gp << "set title \"Gráfico de convergência com " << GENERATIONS 
        << " gerações, população de tamanho " << POPULATION_SIZE
        << ",\\n coeficiente de fitness escalonado de " << setprecision(1) << SCALING_CONSTANT
        << ", coeficiente de penalidade (r) de " << setprecision(4) << R_PENALTY
        << ",\\n taxa de mutação de " << setprecision(4) << 100 * MUTATION_PROBABILITY << "\%"
        << ",\\n elitismo substituindo " << ((ELITISM_OPT) ? "o pior indivíduo" : "um indivíduo aleatório")
        << ", crossover " << ((CROSSOVER_OPT == cx) ? "CX" : ((CROSSOVER_OPT == pmx) ? "PMX" : "genérico"))
        << "\\n usando " << ((SELECTION_OPT) ? ("torneio de " +  to_string(TOURNAMENT_SAMPL) + "indivíduos")
                : ("roleta viciada"))
        << " e " << ((CROSSOVER_OPT == cx) ? "sem pontos de cortes no crossover (usando método CX)" 
                : ((CROSSOVER_OPT == pmx && CUT_POINTS_CROSSOVER > 0) ? ("crossover com corte de tamanho máx." + to_string(CUT_POINTS_CROSSOVER))
                : ((CROSSOVER_OPT == pmx) ? "crossover com corte de tamanho máx. aleatório" 
                : ((CUT_POINTS_CROSSOVER <= 0) ? "crossover uniforme" : "crossover com " + to_string(CUT_POINTS_CROSSOVER) + " pontos de corte" ))))
        << "\"\n";
    gp << "set terminal 'svg'\n";
    gp << "set output \"Imagens/convergencia_" << name << 
    ((file_counter < 10) ? "0" : "") << file_counter << ".svg\"\n";
    gp << "set xlabel \"geração\"\n";
    // gp << "set ytics " << INTERVAL_Y_AXIS << "\n";
    // gp << "set yrange [" << ((worst_value >= 0.0) ? 0 : (((int) worst_value) - INTERVAL_Y_AXIS)) << ":" << (((int) best_value) + INTERVAL_Y_AXIS) << "]\n";
    gp << "set yrange [" << ((worst_value == 0.0) ? 0 : (worst_value - INTERVAL_Y_AXIS)) << ":" << (best_value + INTERVAL_Y_AXIS) << "]\n";
    // cout << "(((int) best_value) + INTERVAL_Y_AXIS) = " << (best_value + INTERVAL_Y_AXIS) << endl;
    if(opt == 0)
        gp << "set ylabel \"fitness do melhor indivíduo\"\n";
    else if(opt == 1)
        gp << "set ylabel \"fitness médio da geração\"\n";
    else
        gp << "set ylabel \"fitness do pior indivíduo\"\n";
    gp << "set grid\n";

    string cores[10] = {
        "#FF0000",  // Vermelho
        "#00FF00",  // Verde
        "#0000FF",  // Azul
        "#FFFF00",  // Amarelo
        "#FF00FF",  // Magenta
        "#00FFFF",  // Ciano
        "#800000",  // Marrom
        "#008000",  // Verde escuro
        "#000080",  // Azul marinho
        "#FFA500"   // Laranja
    };

    gp << "plot ";
    for(int i = 0; i < data_plot.size(); i++)
        gp << "'-' lt rgb \"" << cores[i] << "\" with lines smooth sbezier notitle,";
        // gp << "'-' lt rgb \"" << cores[i] << "\" with lines smooth csplines notitle,";
    gp << "\n";
    for(int i = 0; i < data_plot.size(); i++)
        gp.send(data_plot[i]);
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

/* ----------------- Funções para modelagem do problema --------------- */

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

template <typename ProblemType> 
int evolve(ProblemType &b, population &p)
{
    calcEsc(b, p, SCALING_CONSTANT);
    vector<allele> best_individual;
    vector<vector<allele>> children;
    int best_id = getElite(b, p, best_individual), worst_id = 0, counter = GENERATIONS;
    best_individual = p.matrix[best_id];
    int worst_or_rand_id = 0;
    worst_value = min((double) fitness(b, worst_or_rand_id, p), worst_value);
    best_value = max((double) fitness(b, best_id, p), best_value);
    cout << "Started evolving" << endl;
    do
    {
        calcEsc(b, p, SCALING_CONSTANT);
        crossoverAll(b, p.pop_size, p, CUT_POINTS_CROSSOVER, children);
        
        changeGeneration(children, p);
        
        if(ELITISM_OPT) worst_or_rand_id = getWorst(b, p);
        else worst_or_rand_id = random() % p.pop_size;
        if(PRINTS_OPT && counter % PRINT_INTERVAL == 0){
            cout << fitness(b, worst_or_rand_id, p) << " <- old (worst or random) result" << endl;
            cout << worst_or_rand_id << " <- worst or random id" << endl;
        }
        p.matrix[worst_or_rand_id] = best_individual;
        // best_individual.clear();

        // if(counter % 200 == 0){
        //     cout << fitness(b, best_id, p) << " <- old best result" << endl;
        //     cout << best_id << " <- old best id" << endl;
        // }

        best_id = getElite(b, p, best_individual);
        best_individual = p.matrix[best_id];

        if(PRINTS_OPT && counter % PRINT_INTERVAL == 0){
            cout << fitness(b, best_id, p) << " <- current best result" << endl;
            cout << best_id << " <- current best id" << endl;
        }

        if(GENERATIONS > 0) counter--;
        best_plot.push_back(make_pair((double) (GENERATIONS - counter), (double) fitness(b, best_id, p)));
        average_plot.push_back(make_pair((double) (GENERATIONS - counter), (double) fitnessAvg(b, p)));
        worst_plot.push_back(make_pair((double) (GENERATIONS - counter), (double) fitnessMin(b, p)));
        worst_value = min((double) fitnessMin(b,p), worst_value);
        best_value = max((double) fitness(b, best_id, p), best_value);
    // } while (fitness(b, best_id, p) != (b.dim * (b.dim - 1)) && counter != 0);
    } while (counter != 0);
    return best_id;
}

void printLinearProg(linearProg &b, vector<allele> &ind)
{
    double obj_value = 0.0, expoent;
    int number_size = b.integers_places + b.decimal_places,
        number_variables = (int) b.obj_coeffs.size(),
        start;
    vector<double> values(number_variables, 0.0);

    for(int i = 0; i < number_variables; i++)
    {
        start = i * number_size;
        for(int j = start + number_size - 1; j >= start; j--)
        {
            expoent = (double) (start + number_size - 1 - j);
            if(ind[j].boolean == true) values[i] += pow(2.0, expoent);
        }
    }

    for(int i = 0; i < number_variables; i++)
    {
        double &xi_min = b.intervals[i].first;
        double &xi_max = b.intervals[i].second;
        values[i] = xi_min + 
            ((xi_max - xi_min) / 
            ((pow(2.0, (double) number_size) - 1) / pow(2.0, (double) b.decimal_places)))
            * values[i];
        if(LINEAR_PROG_OPT) values[i] = round(values[i]);
        cout << "normalized x_" << i + 1 << " = " << values[i] << endl;
    }

    for(int i = 0; i < number_variables; i++)
    {
        obj_value += b.obj_coeffs[i] * values[i];
    }
    cout << "obj_value = " << obj_value << endl;
}

void changeGeneration(vector<vector<allele>> &children, population &p)
{
    // for(int i = 0; i < p.pop_size; i++)
    //     p.matrix[i] = children[i];
    p.matrix.swap(children);
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

int fitness(formula &sat, int ind_id, population &p)
{
    int count = 0;
    vector<vector<allele>> &matrix = p.matrix;
    // matrix = p.matrix;

    for(int i = 0; i < sat.formumla_size; i++)
    {
        clause &phi = sat.expr[i];
        bool value = false;
        value = value || ((phi.variables[0] >= 0) ? matrix[ind_id][phi.variables[0] - 1].boolean : !((matrix[ind_id][abs(phi.variables[0]) - 1]).boolean));
        value = value || ((phi.variables[1] >= 0) ? matrix[ind_id][phi.variables[1] - 1].boolean : !((matrix[ind_id][abs(phi.variables[1]) - 1]).boolean));
        value = value || ((phi.variables[2] >= 0) ? matrix[ind_id][phi.variables[2] - 1].boolean : !((matrix[ind_id][abs(phi.variables[2]) - 1]).boolean));
        if(value == true) count++; 
    }
    return count;
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

double fitness(linearProg &b, int ind_id, population &p){
    double penalty = 0.0, expoent, summ, summ_max, obj_value, obj_value_max;
    int number_size = b.integers_places + b.decimal_places,
        number_variables = (int) b.obj_coeffs.size(),
        start, broken_restr = 0;
    bool constrain_flag = false;
    //  Verificar se todos os intervalos estão definidos:
    assert((int) b.intervals.size() == number_variables);
    //  Valores das variáveis:
    vector<double> values(number_variables, 0.0);

    for(int i = 0; i < number_variables; i++)
    {
        start = i * number_size;
        for(int j = start + number_size - 1; j >= start; j--)
        {
            expoent = (double) (start + number_size - 1 - j);
            if(p.matrix[ind_id][j].boolean == true) values[i] += pow(2.0, expoent);
        }
    }

    //  Normalizando variáveis para os seus respectivos intervalos:
    for(int i = 0; i < number_variables; i++)
    {
        double &xi_min = b.intervals[i].first;
        double &xi_max = b.intervals[i].second;
        values[i] = xi_min + 
            ((xi_max - xi_min) / 
            ((pow(2.0, (double) number_size) - 1) / pow(2.0, (double) b.decimal_places)))
            * values[i];
        if(LINEAR_PROG_OPT) values[i] = round(values[i]);
    }

    //  Calculando o valor da função objetivo com as variáveis normalizadas:
    obj_value = 0.0, obj_value_max = 0.0;
    for(int i = 0; i < (int) b.obj_coeffs.size(); i++)
    {
        obj_value += values[i] * b.obj_coeffs[i];
        obj_value_max += b.intervals[i].second * b.obj_coeffs[i];
    }
    //  Normalizando o valor da função objetivo:
    obj_value /= obj_value_max;

    //  Calculando o valor da penalidade normalizada
    for(int i = 0; i < (int) b.leq_restr.size(); i++)
    {
        summ = 0.0, summ_max = 0.0;
        for(int j = 0; j < (int) b.leq_restr[i].size() - 1; j++)
        {
            summ += b.leq_restr[i][j] * values[j];
            summ_max += b.leq_restr[i][j] * b.intervals[j].second;
        }
        if(summ > b.leq_restr[i].back()) {
            constrain_flag = true;
            broken_restr++;
            // O problema era estar somando a penalidade mesmo quando não era violada (isto
            // é, quando seu valor era negativo e portanto isto subtraia os valores positivos
            // dados pelas violações):
            penalty += (double) (summ - b.leq_restr[i].back()) / (double) (summ_max - b.leq_restr[i].back());
        } 
    }
    penalty = ((double) penalty) / ((double) broken_restr);
    
    return (constrain_flag == true) ? obj_value + (R_PENALTY * max(0.0, penalty)) : obj_value;
}

template <typename ProblemType>
int tournament(ProblemType &b, int samples_num, population &p, vector<int> &single_ind)
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

template <typename ProblemType>
int roulette(vector<int> &single_ind, ProblemType &b, population &p)
{
    double summ = 0.0;
    for(int index : single_ind){
        summ += (double) fitnessEsc(b, index, p);
    }

    if(summ == 0){
        int rand_id = random() % single_ind.size();
        int result = single_ind[rand_id];
        eraseFast(single_ind, rand_id);
        return result;
    }

    double prob; 
    vector<double> roulette_chances;
    for(int index : single_ind)
    {
        prob = 100.0 * fitnessEsc(b, index, p) / summ;
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

template <typename ProblemType>
void crossoverAll(ProblemType &b, int samples_num, population &p, int cuts_num, vector<vector<allele>> &children)
{
    vector<int> single_ind;
    for(int i = 0; i < p.pop_size; i++) single_ind.push_back(i); 

    int parent1_index, parent2_index, child_counter = 0;
    while(child_counter < p.pop_size)
    { 
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
        if(CROSSOVER_OPT == pmx)
            crossoverPMX(parent1_index, parent2_index, p, cuts_num, children);
        else if(CROSSOVER_OPT == cx)
            crossoverCX(parent1_index, parent2_index, p, cuts_num, children);
        else if(CROSSOVER_OPT == generic)
            crossoverGeneric(parent1_index, parent2_index, p, cuts_num, children);
        single_ind.push_back(parent1_index);
        single_ind.push_back(parent2_index);
        child_counter += 2;
    }

}

void crossoverPMX(int ind1, int ind2, population &p, int cut_size, vector<vector<allele>> &children)
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

    mutationPerm(s1, p);
    mutationPerm(s2, p);
    children.push_back(s1);
    children.push_back(s2);
}

void crossoverCX(int ind1, int ind2, population &p, int cuts_num, vector<vector<allele>> &children)
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
    mutationPerm(s1, p);
    mutationPerm(s2, p);
    children.push_back(s1);
    children.push_back(s2);
}

void crossoverGeneric(int ind1, int ind2, population &p, int cuts_num, vector<vector<allele>> &children)
{
    vector<allele> son1;
    vector<allele> son2;

    vector<allele> parent1;
    vector<allele> parent2;
    parent1 = p.matrix[ind1];
    parent2 = p.matrix[ind2];

    bool flag = false;
    if(cuts_num <= 0 || cuts_num >= p.ind_dim - 1)
        for (int i = 0; i < p.ind_dim; i++)
        {
            if(random() % 2 == 0)
            {
                son1.push_back(parent1[i]);
                son2.push_back(parent2[i]);
            }
            else
            {
                son1.push_back(parent2[i]);
                son2.push_back(parent1[i]);
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
                    son1.push_back(parent1[j]);
                    son2.push_back(parent2[j]);
                }
                else
                {
                    son1.push_back(parent1[j]);
                    son2.push_back(parent2[j]);
                }
            }
            old_cut_point = cut_point;
        }
    mutationBin(son1, p);
    mutationBin(son2, p);
    children.push_back(son1);
    children.push_back(son2);
}

void mutationPerm(vector<allele> &dna, population &p){
    vector<int> poll(p.ind_dim - 1);
    for(size_t i = 0; i < poll.size(); i++) poll[i] = (int) i;
    for (int i = 0; i < p.ind_dim; i++)
    {
        if(randomDouble(0.0, 1.0) <= MUTATION_PROBABILITY)
            mutateAllelePerm(dna, p, i, poll);
    }
}

void mutateAllelePerm(vector<allele> &dna, population &p, int allele_index, vector<int> &poll)
{
    int random_index = random() % poll.size();
    if(random_index >= allele_index)
        random_index += 1;
    allele aux = dna[allele_index];
    dna[allele_index] = dna[random_index];
    dna[random_index] = aux;
}

void mutationBin(vector<allele> &dna, population &p)
{
    double rand_num;
    for(int i = 0; i < (int) dna.size(); i++)
    {
        rand_num = randomDouble(0.0, 1.0);
        if(rand_num <= MUTATION_PROBABILITY)
            dna[i].boolean = !dna[i].boolean;
    }
}

template <typename ProblemType> 
int getElite(ProblemType &b, population &p, vector<allele> &elite_place)
{
    // double best_fitness = fitnessEsc(b, 0, p);
    double best_fitness = (double) fitness(b, 0, p);
    int current_best_id = 0;
    for (int i = 0; i < p.pop_size; i++)
        {
            // if (fitnessEsc(b, i, p) > best_fitness)
            if (((double) fitness(b, i, p)) > best_fitness)
            {
                // best_fitness = fitnessEsc(b, i, p);
                best_fitness = fitness(b, i, p);
                current_best_id = i;
            }
        }
    elite_place = p.matrix[current_best_id];
    return current_best_id;
}

template <typename ProblemType>
int getWorst(ProblemType &b, population &p)
{
    // double worst_fitness = fitnessEsc(b, 0, p), current_worst_id = 0;
    double worst_fitness = (double) fitness(b, 0, p);
    int current_worst_id = 0;
    for (int i = 0; i < p.pop_size; i++)
        {
            // if (fitnessEsc(b, i, p) < worst_fitness)
            if (fitness(b, i, p) < worst_fitness)
            {
                // worst_fitness = fitnessEsc(b, i, p);
                worst_fitness = fitness(b, i, p);
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

template <typename ProblemType> 
double fitnessAvg(ProblemType &b, population &p)
{
    double sum = 0.0;
    for(int i = 0; i < p.pop_size; i++)
        sum += (double) fitness(b, i, p);
    sum /= (double) p.pop_size;
    return sum;
}

template <typename ProblemType> 
double fitnessMin(ProblemType &b, population &p)
{
    double m = (double) fitness(b, 0, p);
    for(int i = 0; i < p.pop_size; i++)
        m = min(m, (double) fitness(b, i, p)); 
    return m;
}

template <typename ProblemType> 
double fitnessMax(ProblemType &b, population &p)
{
    double m = (double) fitness(b, 0, p);
    for(int i = 0; i < p.pop_size; i++)
        m = max(m, (double) fitness(b, i, p)); 
    return m;
}

template <typename ProblemType> 
void calcEsc(ProblemType &b, population &p, double c)
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

template <typename ProblemType>
double fitnessEsc(ProblemType &b, int ind_id, population &p)
{
    return max(0.0, p.alpha * ((double) fitness(b, ind_id, p)) + p.beta);
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
    sat.expr = new clause[number_of_clauses];

    for(int i = 0; i < number_of_clauses && getline(my_file, line); i++)
    {  
        strcpy(cstr, line.c_str());
        sscanf(cstr, "%d %d %d %*s", &sat.expr[i].variables[0], &sat.expr[i].variables[1], &sat.expr[i].variables[2]);
        // cout << sat.expr[i].variables[0] << ", " << sat.expr[i].variables[1] << ", " << sat.expr[i].variables[2] << endl;
    }
}