import numpy as np

tamanho_pop = 10
dim_pop = 15

cod_pop = np.zeros((tamanho_pop, dim_pop), dtype=int)

for i in range(0, cod_pop.shape[0]):
    for j in range(0, cod_pop.shape[1]):
        cod_pop[i][j] = (np.random.randint(0, 2))

print("código genético da população: ")
print(cod_pop)

cod_pop_int = np.zeros((tamanho_pop, dim_pop), dtype=int)

for i in range(0, cod_pop_int.shape[0]):
    for j in range(0, cod_pop_int.shape[1]):
        cod_pop_int[i][j] = (np.random.randint(-5,11))

print("código genético da população: ")
print(cod_pop_int)

# Inteiro permutado é inteiro sem repetição

# Anotações:
# Genótipo são as variáveis (um vetor v) e Fenótipo é a função f(v), isto é, o indivíduo montado...
# ... através dos genes.
# Um gene é uma posição de um vetor (que tem um determinado tipo) enquanto o alelo é o valor que...
# ... que se encontra nessa posição.
# O número de gerações é o número de "Survivor Selection".
# Deve-se verificar em cada cromossomo se ele leva a um indivíduo válido (uma possível solução para...
# ... o problema).
# As restrições as vezes tem de ser tratadas na função objetivo (aparentemente na maioria das vezes).

