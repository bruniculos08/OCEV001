#include "logic.h"

int main(void)
{
    // srandom(time(0));

    //  Criar o tabuleiro:
    nqueens b;
    b.dim = 8;
    b.board = createBoard(8, 8);

    //  Gerar a população:
    allele low, high;
    low.integer = 0;
    high.integer = b.dim * b.dim;
    population p;
    generatePermPop(b.dim, 10, low, high, p.matrix);

    cout << "Fitness do individuo " << 1 << " = " << fitness(b, 1, p) << endl;
    for(allele &a : p.matrix[1]){
        pair<int,int> u = getPosition(b, a.integer);
        cout << "("<< u.first << ", " << u.second << ")" << " ";
    }
    cout << endl;

    printIntMatrix(p.matrix);
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
    
    std::random_device rd;
    std::mt19937 gen(rd());
    
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
        shuffle(begin(matrix[y]), end(matrix[y]),  gen);
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

double randomDouble(double max, double min)
{
    return (((double) random()) / ((double) RAND_MAX)) * (max - min) + min;
}

int fitness(nqueens &b, int ind_id, population &p)
{
    vector<allele> pieces;
    pieces = p.matrix[ind_id];

    int conflicts = 0;
    for(allele &queen : pieces)
        conflicts += lookRowAndCol(pieces, b, queen) + lookDiag(pieces, b, queen);

    return (b.dim * (b.dim - 1)) - conflicts;
}

int **createBoard(int width, int height)
{
    int **matrix;
    matrix = new int*[height];
    for(int j = 0; j < height; j++)
    {
        matrix[j] = new int[width];
    }
    return matrix;
}

pair<int, int> getPosition(nqueens &b, int id)
{
    int row = (int) floor((double) id / (double) b.dim);
    int col = id % b.dim;
    return make_pair(row, col);
}

int lookDiag(vector<allele> &ind, nqueens &b, allele actual_queen)
{
    pair<int, int> coord_ind, coord_temp;
    coord_ind = getPosition(b, actual_queen.integer);
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
        coord_temp = getPosition(b, other_queen.integer);
        if(coord_temp.first - coord_temp.second == B_growing || 
                coord_temp.first + coord_temp.second == B_decreasing)
                conflicts++;
    }
    return conflicts;
}

int lookRowAndCol(vector<allele> &ind, nqueens &b, allele actual_queen)
{
    pair<int, int> coord_ind, coord_temp;
    int conflicts = 0;
    for(allele &other_queen : ind)
    {   
        if(other_queen.integer == actual_queen.integer)
            continue;
        coord_temp = getPosition(b, other_queen.integer);
        if(coord_temp.first == coord_ind.first || coord_temp.second == coord_ind.second)
            conflicts++;
    }
    return conflicts;
}
