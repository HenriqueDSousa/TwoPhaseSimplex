import numpy as np
import sys
from simplex_parsing import obj_func_parser, get_constraints



def zero_aprox(value):
    # Valores muito baixo sao aproximados para 0

    if (abs(value) < 1e-4):
        value = 0.0

    return value


def standart_pl(c,A,b, b_signs):
    """
    args= 
        c: funcao objetivo
        A: matrix de coeficientes para as constraints
        b: valores à direita das constraints

    return=
        c,A e b em um problema de programacao linear com variaveis basicas e nao basicas
    """
       

    for i in range(len(b_signs)):

        if b_signs[i] == ">=":
            A[i,:] = A[i,:] * (-1)
            b[i] = b[i] * (-1)
        
        if b_signs[i] == "==" :
            
            A = np.vstack((A, -A[i, :]))
            b.append(-b[i])

    m, n = A.shape

    A = np.hstack((A, np.eye(m)))  # Adicionar matrix identidade
    c = np.hstack((c, np.zeros(m)))  # Padding na funcao objetivo com 0s
   
    
    c = [x * (-1) for x in c]

    return c, A, b


class Tableau:

    class Result:
        EMPTY = -3
        UNBOUND = -2
        INVIABLE = -1
        OPTIMAL = 1
    


    def __init__(self):

        self.n = 0
        self.m = 0
        self.c = np.array([])
        self.A = np.array([])
        self.b = np.array([])
        self.status = self.Result.EMPTY
        self.solution = []
        self.certificate = []
        self.max = 1

    def set_tableau(self, input:str) :

        """
            A partir de uma string de entrada, obtem o
            valor das matrizes c, A e b e outras informacoes sobre o 
            LPP
        """

        c, self.max = obj_func_parser(input)
        A, b , b_signs = get_constraints(input)

        self.m, self.n = A.shape
        if self.max == 0:
            
            c = [(-1) * x for x in c]

        
        self.c, self.A, self.b = standart_pl(c, A, b, b_signs)
        
        self.m = (self.A.shape)[0]
        
        self.certificate = np.zeros(self.m)
        

    

    def create_aux_tableau(self):

        """
            Cria um tableau para o problema auxiliar
        """

        tb = np.zeros((self.m+1, self.n+ self.m+ self.m+1))
        
        tb[1:self.m+1 , :(self.n + self.m) ] = self.A
        
        for i in range(1, self.n+1):
            factor = 1
            if(self.b[i-1] < 0):
                factor  *= -1
            tb[i, :] = tb[i, :] * factor

        tb[1:, -1] = np.abs(self.b)
        
        # Adicionando variaveis auxiliares
        tb[1:, (self.n+self.m):-1] = np.eye(self.m)
        tb[0, (self.n+self.m):-1] = np.ones(self.m)

        for line in tb[1:,:]:
            tb[0,:] -= line        

        return tb
    
    
    def get_pivot_column_index(self, vecC):
        """
            Retorna a coluna do pivo, seguindo a regra de Bland
            Se retornar -1, nao ha como pivotear
        """

        for i in range(len(vecC)):
            if round(vecC[i], 5) < 0:
                return i
        
        return -1

    def get_pivot_row_index(self, column, b):

        """
            Funcao para descobrir a linha do pivo, seguindo a regra de Bland
            Se retornar -1, nao ha como pivotear
        """

        min_idx = -1
        min_val = -1
        for i in range(self.m):

            if round(column[i], 5) > 0:

                factor = round(b[i]/column[i], 5)
                if min_idx == -1 or factor < min_val:
                    min_idx = i
                    min_val = factor
    
        return min_idx


    def escalonar(self, pcol, prow, tb):
    
        tb[prow, :] = tb[prow, :] * (1/tb[prow, pcol])

        for i in range(self.m+1):
            if i == prow:
                continue
            factor = -tb[i, pcol]
            tb[i, :] += (factor*tb[prow,:])
    
    def solve_tableau(self, tb):

        vfunc = np.vectorize(zero_aprox)
        # pegar o indice da coluna que será pivoteada usando a regra de Bland
        pcol = self.get_pivot_column_index(tb[0, :-1])
        
        while pcol != -1:  

        # escolher qual a linha ideal para o pivot segundo Bland
            prow = self.get_pivot_row_index(tb[1:, pcol], tb[1:, -1])
            
            if(prow == -1):
                self.status = self.Result.UNBOUND
                #certificado de ilimitado
                self.certificate = tb[0:, pcol]
                return tb
            
            prow += 1
            
            # escalonar a matriz
            self.escalonar(pcol, prow, tb)
            
            
            tb = vfunc(tb)
            
            # pegar o novo pivot
            pcol = self.get_pivot_column_index(tb[0, :-1])
            
        
        return tb
    
    
    def index_of_unique_one(self, col):
        zeros_in_col = np.round(col, 5) == 0.0
        ones_in_col = np.round(col, 5) == 1.0
        
        if(np.count_nonzero(zeros_in_col) == (self.m-1) # count the quantity of zeros
        and np.count_nonzero(ones_in_col) == 1): # count the quantity of zeros
            return  np.where(ones_in_col)[0][0]
        
        return -1

    def canonical_form(self, tb):

        """
            Um tableau estara na forma canonica sse todas posicoes
            nas colunas basicas forem 0, assim como suas entradas 
            correspondentes na funcao objetivo, exceto pelo pivo, que devera ser 1.

        """

        matA = tb[1:, :self.n].T
        vecC = tb[0, :self.n]
        basis = set()

        for i in range(self.n):
            one_idx = self.index_of_unique_one(matA[i])
        if one_idx >= 0 and one_idx not in basis:
            basis.add(one_idx)

            tb[0,:] += ((-vecC[i])*tb[one_idx+1,:])


    def basic_cols(self, tb):
        """
            Funcao para encontrar as variaveis basicas
        """
        base = []
        for i in range(self.n + self.m):
            
            line = self.index_of_unique_one(tb[1:, i])
            if (line != -1) and i not in base:
                base.append((i, line))

        return base

    def get_solution(self, tb, base):

        RHS = tb[1:, -1]
        objective_vals = np.arange(self.n)
        solution = np.zeros(self.n)
 
        for col, line in base: 
            if col < self.n:
                solution[col] = RHS[line] 

        return solution

    def solve(self):
        
        # Fase 1
        tb = self.create_aux_tableau()
        
        

        tb = self.solve_tableau(tb)        

        if round(tb[0][-1], 5) == 0:
        # Fase 2

            # remover a parte auxiliar da matrix 
            tb = np.delete(tb, np.s_[(self.n+self.m):-1], axis=1)
            
            tb[0, :-1] = self.c
            
            self.canonical_form(tb)
            
            # inicio das iteracoes
            tb = self.solve_tableau(tb) 

            if (self.status == self.Result.EMPTY):
                
                self.status = self.Result.OPTIMAL
                self.certificate = tb[1, :-1]

            self.c =  -tb[0, :self.n]
            self.optimal = tb[0, -1]

            if self.max == 0:
                self.optimal *= (-1)



        else:
    
            self.certificate = tb[1, :-1]
            self.status = self.Result.INVIABLE
            
        # solucao
        base = self.basic_cols(tb)
        self.solution = self.get_solution(tb, base)

#----------------------------------------------------------   EXECUCAO   -----------------------------------------------------------------------------------------------#


# arquivo de entrada
if '-i' in sys.argv:

    filename = ""

    if sys.argv.index('-i') + 1 < len(sys.argv):
        filename = sys.argv[sys.argv.index('-i') + 1]
    
    else:
        print("Erro: arquivo de entrada não especificado")
        sys.exit(1)

    table = Tableau()   


    file = open(filename, "r")
    input = file.read()

    table.set_tableau(input)

    table.solve()

    # arquivo de output
    if '-o' in sys.argv:

        filename = filename = sys.argv[sys.argv.index('-o') + 1]


    outfile=  open(filename, "w")

    #UNBOUND
    if (table.status == -2):
        print("Status: ilimitado", file= outfile)
        print("Certificado:", file= outfile)
        print(table.certificate, file= outfile)

    elif (table.status == -1):

        print("Status: inviavel", file= outfile)
        print("Certificado:", file= outfile)
        print(table.certificate, file= outfile)

    elif (table.status == 1):

        print("Status: otimo", file= outfile)
        print("Objetivo: ", table.optimal, file= outfile)
        print("Solucao:", file= outfile)
        print(table.solution, file= outfile)
        print("Certificado:", file= outfile)
        print(table.certificate, file= outfile)
