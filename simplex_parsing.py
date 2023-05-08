import numpy as np


def insert_number_at_index(lst, num, index):
    """
    Insere um numero em um especifico index na lista, com 0s para posicoes nao inseridas anteriormente.
    
    Args:
        lst (list): Lista
        num (int ou float): numero a ser inserido
        index (int): Index a inserir o numero
        
    Returns:
        list: lista com a nova insercao
    """
    #padding
    if index >= len(lst):
        lst.extend([0] * (index - len(lst) + 1))

    # insere
    lst[index] = num

    return lst



def obj_func_parser(input):

    """
    Parser para a funcao objetivo
    
    Args:
        input(str): string com a função objetivo
        
    Returns:
        list: funcao objetivo em formato matricial 
        max: booleano para verificar se é tratado o caso de min ou max -> 1 = max e 0 = min              
    """
    input = input.split(sep="\n")[0]
    
    max = 1
    if (input.split()[0] == "MIN"):
        max = 0

    c = []

    next_negative = False
    for word in input.split()[1:]:

        if word[0] == "+":
            continue
        
        elif word[0] == "-":
            next_negative = True
            continue

        elif word[0] == "−":
            next_negative = True
            continue

        #caso de xN ser multiplicado por 1
        if word[0].isdigit() == False:


            if(next_negative == False):
                insert_number_at_index(c, 1, int(word[1:])-1)

            else:
                #caso seja negativo
                insert_number_at_index(c, -1, int(word[1:])-1) 
                next_negative = False
        

        #caso tenha um coeficiente
        elif word[0].isdigit() == True:
            
            coef = str("")
            div = str("")
            i = 0
            for caract in word:

                #coletando o coeficiente ate *
                if caract != "*" and caract != "/":
                    coef = coef + caract
                    i += 1

                elif caract == "/":
                    div = coef
                    coef = ""
                    i += 1

                else:

                    num = 0.0
                    if div != "":
                        num = float(div)/float(coef)
                    else:
                        num = float(coef)
                    i+=2
                    if next_negative == False:
                        insert_number_at_index(c, num , int(word[i:])-1)
                    else:
                        insert_number_at_index(c, (-1) * num, int(word[i:])-1)
                        next_negative = False


    return c, max



def add_constraint(input):

    new_A = []
    new_b = 0
    b_sign = ""
    next_negative = False

    input = input.split()

    for word in input[:-1]:
        
        if word == ">=":
            b_sign = ">="
            continue

        elif word == "<=":
            b_sign = "<="
            continue

        elif word == "=":
            b_sign = "="
            continue
        
        if word[0] == "+":
            continue
        
        if word[0] == "-":
            next_negative = True
            continue
        
        if word[0] == "−":
            next_negative = True
            continue
        
        if word[0].isdigit() == False:

            if(next_negative == False):
    
                insert_number_at_index(new_A, 1, int(word[1:])-1)

            else:
                #caso seja negativo
                insert_number_at_index(new_A, -1, int(word[1:])-1) 
                next_negative = False
        
        elif word[0].isdigit() == True:
        
            coef = str("")
            div = str("")
            i = 0
            for caract in word:

                #coletando o coeficiente ate *
                if caract != "*" and caract != "/":
                    coef = coef + caract
                    i += 1

                elif caract == "/":
                    div = coef
                    coef = ""
                    i += 1

                else:
                    num = 0.0
                    if div != "":
                        num = float(div)/float(coef)
                    else:
                        num = float(coef)
                    i+=2
                    if next_negative == False:
                        insert_number_at_index(new_A, num , int(word[i:])-1)
                    else:
                        insert_number_at_index(new_A, (-1) * num, int(word[i:])-1)
                        next_negative = False

        new_b =  float(input[-1])

    return new_A, new_b, b_sign

def get_constraints(input:str):

    A = []
    b = []
    b_signs = []

    for line in input.split(sep="\n")[1:]:
        new_A, new_b, b_sign = add_constraint(line)
        A.append(new_A)
        b.append(new_b)
        b_signs.append(b_sign)

    max_len = 0
    for sublist in A:
        
        if len(sublist) > max_len:
            max_len = len(sublist) 

    matrix = []
    # Padding para tornar todas listas com o mesmo tamanho
    for sublist in A :
        padded_sublist = sublist + [0] * (max_len - len(sublist))
        matrix.append(padded_sublist)

    A = np.array(matrix)
    
    return A, b, b_signs


