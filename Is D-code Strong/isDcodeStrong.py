from itertools import combinations
from sage.combinat import combination
import math



# проверка D-кода на удовлетворение условиям теоремы 3 (стойкие коды)
def is_theorem_3(d_code, m1, m2):

    i = 0
    while i < len(d_code):
        if d_code[i][0] >= m1 / 2 and d_code[i][1] >= m2 / 2:
            return True
        if i < len(d_code) - 1:
            j = i + 1
            while j < len(d_code):
                if d_code[i][0] + d_code[j][0] >= m1 and d_code[i][1] + d_code[j][1] >= m2:
                    return True
                j += 1
        i += 1
    return False

# проверка D-кода на удовлетворение условиям теоремы 4 (стойкие коды)
def is_theorem_4(d_code, m1, m2):
    if len(d_code) < 2:
        return False

    if d_code[0][0] < m1 / 2 and d_code[0][1] >= m2 / 2:
        j = 1
        while j < len(d_code):
            if d_code[j][0] < m1 / 2 or d_code[j][1] >= m2 / 2:
                return False
            j += 1
    else:
        return False

    i = 0
    while i < len(d_code) - 1:
        j = i + 1
        while j < len(d_code):
            if d_code[i][0] + d_code[j][0] < m1 or d_code[i][1] + d_code[j][1] >= m2:
                return False
            j += 1
        i += 1

    return True

# проверка D-кода на удовлетворение условиям теоремы 5 (стойкие коды)
def is_theorem_5(d_code, m1, m2):
    if len(d_code) < 2:
        return False

    if d_code[-1][0] >= m1 / 2 and d_code[-1][1] < m2 / 2:
        j = 0
        while j < len(d_code) - 1:
            if d_code[j][0] >= m1 / 2 or d_code[j][1] < m2 / 2:
                return False
            j += 1
    else:
        return False

    i = 0
    while i < len(d_code) - 1:
        j = i + 1
        while j < len(d_code):
            if d_code[i][0] + d_code[j][0] >= m1 or d_code[i][1] + d_code[j][1] < m2:
                return False
            j += 1
        i += 1

    return True
    

#Вычисляет степени заданного D-кода и проверяет делимость размерности, возвращает True, если одна из степеней раскладывается в прямую сумму
def dcode_degrees(dcode, m1, m2, print_mode=True):
    
    dim = dim_d_code(dcode, m1, m2)
    if print_mode:
        print(f"{1} degree: {dcode}, dim: {dim_d_code(dcode, m1, m2)}, div 2^m_1: {float(dim/(2**m1))}, div 2^m_2: {float(dim/(2**m2))}")
    
    prev_dcode = set(dcode)
    next_dcode = set()
    flag = True
    degree = 2
    while flag:
        
        #вычисляем степень
        for dcode_i in prev_dcode:
            for dcode_j in dcode:
                next_dcode.add((min(dcode_i[0]+dcode_j[0], m1), min(dcode_i[1]+dcode_j[1], m2)))
                if min(dcode_i[0]+dcode_j[0], m1) == m1 and min(dcode_i[1]+dcode_j[1], m2) == m2:
                    flag = False
                    if float(dim/(2**m1))%1 != 0 and float(dim/(2**m2))%1 != 0:
                        result_flg = False
                    else:
                        result_flg = True
    
        #убираем вложенные слагаемые
        result_dcode = set()
        sub_flg = False
        for dcode_i in next_dcode:
            for dcode_j in next_dcode:
                if dcode_i[0] <= dcode_j[0] and dcode_i[1] < dcode_j[1] or dcode_i[0] < dcode_j[0] and dcode_i[1] <= dcode_j[1]:
                    sub_flg = True
            if not sub_flg:
                result_dcode.add(dcode_i)
            sub_flg = False
            
        prev_dcode = result_dcode
        next_dcode = set()
        result_dcode = sorted(list(result_dcode), key = lambda element : element[0])
        
        dim_new = dim_d_code(result_dcode, m1, m2)
        if dim_new == dim:
            flag = False
            if float(dim/(2**m1))%1 != 0 and float(dim/(2**m2))%1 != 0:
                result_flg = False
            else:
                result_flg = True
            
        dim = dim_new
        if print_mode:
            print(f"{degree} degree: {result_dcode}, dim: {dim_d_code(result_dcode, m1, m2)}, div 2^m_1: {float(dim/(2**m1))}, div 2^m_2: {float(dim/(2**m2))}")
        degree += 1 
        
    return result_flg


def isDcodeStrong(dcode, m1, m2):
    
    theorem_3_flg = is_theorem_3(dcode, m1, m2)
    theorem_4_flg = is_theorem_4(dcode, m1, m2)
    theorem_5_flg = is_theorem_5(dcode, m1, m2)
    
    if (theorem_3_flg or theorem_4_flg or theorem_5_flg or not dcode_degrees(dcode, m1, m2, False)):
        print(f'D-code {dcode} (m1 = {m1}, m2 = {m2}) is strong')
    else:
        print(f'D-code {dcode} (m1 = {m1}, m2 = {m2}) is weak')
        
        
        
if __name__ == '__main__':
    
    dcode = [(0, 3), (3, 2), (5, 1)]
    m1 = 8
    m2 = 8
    
    isDcodeStrong(dcode, m1, m2)
    

