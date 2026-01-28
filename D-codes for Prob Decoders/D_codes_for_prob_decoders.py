from itertools import combinations
from sage.combinat import combination
from multiprocessing import Process
from multiprocessing import Manager
import random
import math
import time


M_1, M_2 = 5, 8

# проверка на то, что D-код составлен из подходящей последовательности
def is_correct_dots(d_code):
    i = 1
    while i < len(d_code):
        if d_code[i][0] <= d_code[i - 1][0] or d_code[i][1] >= d_code[i - 1][1]:
            return False
        i += 1
    return True


# вычисление размерности D-кода
def dim_d_code(d_code, m1=M_1, m2=M_2):
    i = 0
    d_b = set()
    d_b.add((0, d_code[0][1] + 1))
    while i < len(d_code) - 1:
        d_b.add((d_code[i][0] + 1, d_code[i+1][1] + 1))
        i += 1
    d_b.add((d_code[-1][0] + 1, 0))

    d = set()
    for (k, l) in d_b:
        i = k
        while i <= m1 + 1:
            j = l
            while j <= m2 + 1:
                d.add((i, j))
                j += 1
            i += 1
    i = 0
    dim = 0
    for (i, j) in d:
        dim += binomial(m1, i) * binomial(m2, j)

    return (2 ** m1) * (2 ** m2) - dim


# вычисление кодового растояния D-кода
def dist_d_code(d_code, m1=M_1, m2=M_2):
    i = 1
    min_dist = (2 ** (m1 - d_code[0][0])) * (2 ** (m2 - d_code[0][1]))
    while i < len(d_code):
        min_dist_i = (2 ** (m1 - d_code[i][0])) * (2 ** (m2 - d_code[i][1]))
        min_dist = min_dist_i if min_dist_i < min_dist else min_dist
        i += 1
    return min_dist


# вычисление вероятности успеха атаки ISD
def get_sec_lvl(n, dim, dist):
    return round(-math.log2(float(binomial(n - round((dist - 1) / 2), dim) / binomial(n, dim))), 2)


# вычисление вероятности P1 выбора хороших блоков
def get_correct_prob(n, blocks_count, count_bad):
    return float(binomial((n - count_bad), blocks_count)) / binomial(n, blocks_count)


# проверка D-кода на удовлетворение условиям теоремы 2 (нестойкие коды)
def is_theorem_2(d_code):
    ind = 0 if len(d_code)  < 2 else -1
    if (d_code[0][0] >= M_1 / 2 and d_code[0][1] < M_1 / 2 
        or  d_code[ind][0] < M_1 / 2 and d_code[ind][1] >= M_1 / 2):
        return true
    else:
        return false  


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


# перебор разных D-кодов, выбор наиболее стойких к обычной ISD среди стойких/нестойких к модифицированной атаке
# 
def get_d_codes(r1, r2, m1, m2, dots_count=0, strength=True, sec_lvl_flg=False, print_mode=True):
    
        
    n = (2 ** m1) * (2 ** m2)

    all_dots_ords = []
    i = 0
    while i <= m1:
        j = 0
        while j <= m2:
            if (i<=r1 and j<=r2):
                all_dots_ords.append((i, j))
            j += 1
        i += 1

    d_codes_count = 0
    choosen_d_codes_count = 0

    if dots_count == 0:
        i = 1
        finish = min(m1, m2, max(m1, m2))
    else:
        i = dots_count
        finish = dots_count

    result_dcodes=[]
    sec_d_code = []
    max_sec_lvl = 1

    while i <= finish:
        all_subsets = (subset for subset in combinations(all_dots_ords, i) if is_correct_dots(subset))
        for cur_dots_ords in all_subsets:
            theorem_2_flg = is_theorem_2(cur_dots_ords)
            theorem_3_flg = is_theorem_3(cur_dots_ords, m1, m2)
            theorem_4_flg = is_theorem_4(cur_dots_ords, m1, m2)
            theorem_5_flg = is_theorem_5(cur_dots_ords, m1, m2)
            if strength and (theorem_3_flg or theorem_4_flg or theorem_5_flg or not dcode_degrees(cur_dots_ords, m1, m2, False)):
                
                #if cur_dots_ords[0][0]==3 and cur_dots_ords[0][1]==3:
                #print(dcode_degrees(cur_dots_ords, m1, m2, True))
                #print(theorem_3_flg)
                #print(theorem_4_flg)
                #print(theorem_5_flg)
                
                result_dcodes.append(list(cur_dots_ords))
                dim = dim_d_code(cur_dots_ords, m1, m2)
                dist = dist_d_code(cur_dots_ords)
                if sec_lvl_flg:
                    sec_lvl = get_sec_lvl(n, dim, dist)
                    sec_d_code = list(cur_dots_ords) if sec_lvl < max_sec_lvl else sec_d_code
                    max_sec_lvl = sec_lvl if sec_lvl < max_sec_lvl else max_sec_lvl
                    sec_lvl_str = f"; sec_lvl = {sec_lvl}"
                else:
                    sec_lvl_str = ''
                
                strength_str = (' ' + str(3) if theorem_3_flg else '') + (' ' + str(4) if theorem_4_flg else '') \
                               + (' ' + str(5) if theorem_5_flg else '')
                if print_mode:
                    print(f"D-code: {cur_dots_ords}; strength:{strength_str}; dim = {dim}; "
                        f"dist = {dist}" + sec_lvl_str)
                choosen_d_codes_count += 1
            elif not strength and (dcode_degrees(cur_dots_ords, m1, m2, False) or theorem_2_flg):
                result_dcodes.append(list(cur_dots_ords))
                dim = dim_d_code(cur_dots_ords, m1, m2)
                dist = dist_d_code(cur_dots_ords)
                if sec_lvl_flg:
                    sec_lvl = get_sec_lvl(n, dim, dist)
                    sec_d_code = list(cur_dots_ords) if sec_lvl < max_sec_lvl else sec_d_code
                    max_sec_lvl = sec_lvl if sec_lvl < max_sec_lvl else max_sec_lvl
                    sec_lvl_str = f"; sec_lvl = {sec_lvl}"
                else:
                    sec_lvl_str = ''
                if print_mode:
                    print(f"D-code: {cur_dots_ords}; dim = {dim}; dist = {dist}" + sec_lvl_str)
                    
                
                choosen_d_codes_count += 1
            d_codes_count += 1
        i += 1
        
    if print_mode:
        print(f"N = {n}")
        print(f"D-codes count: {d_codes_count}, {'strength' if strength else 'weak'} D-codes count: {choosen_d_codes_count}")
    
    if sec_lvl_flg and choosen_d_codes_count > 0 and print_mode:
        print(f"Most strength: {sec_d_code}; dim = {dim_d_code(sec_d_code, m1, m2)}; dist = {dist_d_code(sec_d_code)}; "
              f"sec_lvl = {max_sec_lvl}")

    return result_dcodes


# вычисление размерности кода Рида-Маллера 
def dim_RM(r, m):
    i = 0
    dim = 0
    while i <= r:
        dim += binomial(m, i)
        i += 1
    return dim


# нахождение наиболее стойкого к атаке ISD кода Рида-Маллера с заданным параметром m
def sec_RM(m):
    i = 0
    max_sec = 1
    max_r = 0
    n = 2 ** m
    while i <= m:
        sec_lvl = get_sec_lvl(n, dim_RM(i, m), 2 ** (m - i))
        max_r = i if max_sec > sec_lvl else max_r
        max_sec = sec_lvl if max_sec > sec_lvl else max_sec
        i += 1
    return [max_r, max_sec]

# поиск порождающей матрицы кода Рида-Маллера
def get_rm_gen_matr(r, m):
    gen_matr = []
    gen_matr_1 = []
    for i in range(0, r + 1):
        if i == 0:
            gen_matr_i = [[1 for j in range(0, 2 ** m)]]
        elif i == 1:
            gen_matr_1 = [[int(bin(j)[2:].zfill(m)[i]) for j in range(0, 2 ** m)] for i in range(0, m)]
            gen_matr_i = gen_matr_1
        else:
            gen_matr_i = []
            for comb in combinations(gen_matr_1, i):
                v_i = [1 for i in range(0, 2 ** m)]
                for v in comb:
                    for dig in range(0, 2 ** m):
                        v_i[dig] *= v[dig]
                gen_matr_i.append(v_i)
        gen_matr.extend(gen_matr_i)

    return matrix(GF(2),gen_matr)


# вычисление количества подматриц (из count_iter подматриц), составленных из K случайно выбранных столбцов и имеющих ранг k
def count_full_rank_RM(matr, K, k, count_iter):
    
    counter = 0
    
    for i in range(0, count_iter):
        columns = set()
        while sum(1 for i in columns) < K:
            columns.add(random.randint(0, len(matr[0]) - 1))
        
        if matr.matrix_from_columns(columns).rank() == k:
            counter += 1
        #print(f"-- {i} finish, {counter} submatrix found --")
            
    return counter
    
    
# вычисление количества подматриц (из count_iter подматриц), 
# составленных из blocks_count случайных блоков длины n2 и имеющих ранг k
def count_full_rank_d_code(matr, dim, count_iter, all_blocks_count, cur_blocks_count, block_length, desc_flg=False, reverse_flg=False, print_mode=True):
    
    counter = 0
    
    
    
    for i in range(1, count_iter):
        blocks = set()
        while sum(1 for i in blocks) < cur_blocks_count:
            blocks.add(random.randint(0, all_blocks_count - 1))
        
        columns = []
        for j in blocks:
            for l in range(0, block_length):
                if not reverse_flg:
                    columns.append((block_length) * j + l)
                else:
                    columns.append((block_length) * l + j)
        if print_mode:
            print('-- columns found --')
            print(f"Columns num = {len(columns)}")
            print(f'Cur matr rank = {matr.matrix_from_columns(columns).rank()}')
            
        if matr.matrix_from_columns(columns).rank() == dim:
            counter += 1
            if desc_flg:
                return counter
        if print_mode:
            print(f"-- {i} finish, {counter} submatrix found --")
            
    return counter


# вычисление для матрицы matr минимального необходимого количества блоков длины n2, 
# при котором подматрица, составленная из такого количества блоков, будет иметь ранг k с непренебрежимой вероятностью
# для каждого количества блоков вычисляется вероятность P1
def find_blocks_count(matr, dim_full, dim_small, all_blocks_count, bad_blocks_count, block_length, reverse_flg, print_mode=True):
    
    min_blocks_count = min(ceil(dim_full / dim_small) + 1, all_blocks_count)
    
    good_block_count = all_blocks_count - bad_blocks_count
    
    
    for i in range(all_blocks_count, min_blocks_count, -1):
        if print_mode:
            print(f"--------FINDING MATR FROM {i} BLOCKS ---------")
            print(f"P1 = {'{0:.3E}'.format(get_correct_prob(all_blocks_count, i, bad_blocks_count))}")
        p2_cnt = count_full_rank_d_code(matr, dim_full, 100, all_blocks_count, i, block_length, true, reverse_flg, print_mode)
        if p2_cnt == 0:
            if i == all_blocks_count:
                if print_mode:
                    print('NOT FOUND')
                return -1
            else:
                if print_mode:
                    print('DONE')
                return i+1
    return -1


# Вычисление "порождающей" матрицы D-кода (неполного ранга)
def get_d_code_gen_matr(d_code, m1, m2, print_mode=True):
    
    if print_mode:
        print('-- start matr --')
    gen_matr = matrix(GF(2), 1, 2 ** (m1 + m2))
    for (r_1, r_2) in d_code:
        gen_matr_i = get_rm_gen_matr(r_1, m1).tensor_product(get_rm_gen_matr(r_2, m2))
        #print('-- tens prod --')
        gen_matr = block_matrix([[gen_matr],[gen_matr_i]])
        #print('-- concat --')
        #for row in gen_matr_i.rows():
            #gen_matr = matrix(gen_matr.rows() + [row])
    
    #print('-- finish concat --')
    
    #V = VectorSpace(GF(2), 2 ** (M_1 + M_2))
    #gen_matr = V.subspace(gen_matr.rows())
    #gen_matr = matrix(gen_matr.basis())
    
    if print_mode:
        print('-- finish matr --')
    return gen_matr


# вычисление ранга всех блоков длины n2
def d_code_rank_blocks():
    
    matr = get_d_code_gen_matr([[4, 3], [5, 2]])
    
    for j in range(2 ** M_1):
    
        columns = []
        for i in range(0, 2 ** M_2):
            columns.append((2 ** M_2) * j + i)
    
        print(columns)
        print(matr.matrix_from_columns(columns).rank())


# перебор тензорных произведений кодов Рида-Маллера,
# вычисление для них вероятностей p1, p2, p_attack,
# вычисление количества кодов, для которых Ng < k1 (таблица 1,2 статьи KoLe2022)
def print_tens_codes():
    
    counter_all = 0
    counter_bad = 0
    
    codes = []
    
    for m in range(7,9):
        
        r2 = ceil(m/2) - 1
        
        for r1 in range(ceil(m/2), m):
            
            code = []
            code.append(r1)
            code.append(m)
            code.append(r2)
            code.append(m)
            codes.append(code)
    
    
    codes.append([4,8,3,7])
    codes.append([4,8,2,8])
    print(codes)
    
 
    for code in codes:
        
        r1 = code[0]
        m1 = code[1]
        r2 = code[2]
        m2 = code[3]


        k1 = dim_RM(r1, m1)
        k2 = dim_RM(r2, m2)
        k = k1 * k2
        n = 2 ** (m1 + m2)
        d = 2 ** ((m1 - r1) + (m2 - r2))
        
        Ng_min = 2 ** m1 - floor((2 ** (m1 - r1) * 2 ** (m2 - r2) - 1)/(2 ** (m2 - r2) + 1))
        Ng_avg = 2 ** m1 - calc_avg_bad_blocks_cnt_prob(2 ** m1, 2 ** m2, d, 2 ** (m2 - r2), 8)
        
        p_isd = get_sec_lvl(n, k, d)
        
        if Ng_min < k1:
            counter_all += 1
            counter_bad += 1

        prev_p = 0

        for K in range(k1, 2 ** m1):
            counter_all += 1

            p1_min = float(binomial(Ng_min, K)/binomial(2 ** m1, K))
            p1_avg = float(binomial(Ng_avg, K)/binomial(2 ** m1, K))
            
            p2 = float(count_full_rank_RM(get_rm_gen_matr(r1,m1), K, k1, 10000) / 10000)
            p_min = p1_min * p2
            p_avg = p1_avg * p2
            
            #print(f"p1_min = {p1_min}, p2 = {p2} , p_min = {p_min}")

            print(f"$\\mathrm{{RM}}({r1},{m1}) \\otimes \\mathrm{{RM}}({r2},{m2})$ & {p2} & {K} & {'{0:.3E}'.format(p_isd)} & {'{0:.3E}'.format(p_min)} & {'{0:.3E}'.format(p_avg)}\\\\")

            if p_min < prev_p or p_min == 0:
                print('\\Xhline{1pt}')
                break
            else:
                print('\\hline')

            prev_p = p_min
            
    print(f"count_all = {counter_all}, count_bad = {counter_bad}")

    for code in codes:
        
        r1 = code[0]
        m1 = code[1]
        r2 = code[2]
        m2 = code[3]
                
        k1 = dim_RM(r1, m1)
        k2 = dim_RM(r2, m2)
        k = k1 * k2
        n = 2 ** (m1 + m2)
        d = 2 ** ((m1 - r1) + (m2 - r2))
        
        t = floor((d - 1)/2)
        t2 = floor((2 ** (m2 - r2) - 1)/2)
        
        Nb_max = floor((2 ** (m1 - r1) * 2 ** (m2 - r2) - 1)/(2 ** (m2 - r2) + 1))
        Nb_avg = calc_avg_bad_blocks_cnt_prob(2 ** m1, 2 ** m2, d, 2 ** (m2 - r2), 8)
        
        result = dict()
        rs = []
        rs.append(Nb_max)
        rs.append(Nb_avg)
        calc_avg_bad_blocks_cnt_prob_worker(t, t2, 2 ** m1, 2 ** m2, rs, result)
        all_comb = binomial(2 ** (m1 + m2), t)
        
        p_ng_min = result[Nb_max]*1.0/all_comb
        p_ng_avg = float(result[Nb_avg]*1.0/all_comb)
        
                
        print(f"$\\mathrm{{RM}}({r1},{m1}) \\otimes \\mathrm{{RM}}({r2},{m2})$ & {2 ** m1} & {2 ** m1 - Nb_max} & {2 ** m1 - Nb_avg} & {'{0:.3E}'.format(p_ng_min)} & {'{0:.3E}'.format(p_ng_avg)} \\\\")
        print('\\hline')


# вычисляет среднее количество плохих блоков (веса больше T2) длины N2 в коде длины N с количеством ошибок веса T за count_iter итераций
def find_avg_bad_blocks_cnt(N, N2, dist, dist_small, count_iter):
    
    T = floor((dist - 1)/2)
    T2 = floor((dist_small - 1)/2)
    
    error = []
    for i in range(0, T):
        error.append(1)
        
    for i in range(0, N-T):
        error.append(0)
        
    avg=0
    for i in range(0,count_iter):
        random.shuffle(error)
        #print(error)
        
        count_bad_blocks = 0
        block_weight = 0
        for j in range(0, N):
            
            if error[j] == 1:
                block_weight+=1
                
            if (j+1) % N2 == 0:
                #print(j)
                if block_weight > T2:
                    #print('bad_block_found')
                    count_bad_blocks += 1
                block_weight = 0
               
        avg+=count_bad_blocks
        #print(f"{i}-iter: bad_blocks_count = {count_bad_blocks}, avg = {avg*1.0/(i+1)}") 
        
    print(f"avg_bad_blocks_count = {avg*1.0/count_iter}")
    return avg*1.0/count_iter


# (многопоточное) воркер для find_avg_bad_blocks_cnt_mp
def calc_bad_blocks_count_worker(error, N, N2, T2, count_iter, result):
    
    total_bab_blocks_cnt = 0

    for i in range(0, count_iter):
        
        random.shuffle(error)
        count_bad_blocks = 0
        block_weight = 0
        
        for j in range(0, N):
            
            if error[j] == 1:
                block_weight+=1
                
            if (j+1) % N2 == 0:
                if block_weight > T2:
                    count_bad_blocks += 1
                    
                block_weight = 0
        
        total_bab_blocks_cnt+=count_bad_blocks
    
    result.append(total_bab_blocks_cnt)


# (многопоточное) вычисляет среднее количество плохих блоков (веса больше T2) длины N2 в коде длины N с количеством ошибок веса T за count_iter итераций    
def find_avg_bad_blocks_cnt_mp(N, N2, dist, dist_small, count_iter, nCPU):
    
    T = floor((dist - 1)/2)
    T2 = floor((dist_small - 1)/2)
    
    error = []
    for i in range(0, T):
        error.append(1)
        
    for i in range(0, N-T):
        error.append(0)
    
    errors = []
    
    for i in range(nCPU):
        random.shuffle(error)
        errors.append(error)
    
    manager = Manager()
    result = manager.list()
    processes = []
    for error_wrk in errors:

        p = Process(target = calc_bad_blocks_count_worker, args = (error_wrk, N, N2, T2, count_iter, result))
        p.start()
        processes.append(p)
    
    print("All processes started")
    
    for process in processes:
        process.join()
    
    print("All processes ended")
    
    sum=0
    for res in result:
        sum+=res
        
    res=sum*1.0/(count_iter*nCPU)
    
    print(f"avg_bad_blocks_count = {res}")
    return res


# (многопоточное) воркер для calc_avg_bad_blocks_cnt_prob
def calc_avg_bad_blocks_cnt_prob_worker(t, t2, n1, n2, rs, result):
    
    for r in rs:     
        sum=0
        for k in range(r, n1+1):
        
            sum_i = (-1)**(k-r) * binomial(k, r) * binomial(n1, k) * binomial(n2, t2+1)**k * binomial((n1-k)*n2, t-(t2+1)*k)
            
            sum += sum_i   
        
        result[r] = sum


#  (многопоточное) вычисляет среднее количество плохих блоков, используюя формулу включения-исключения
def calc_avg_bad_blocks_cnt_prob(n1, n2, dist, dist_small, nCPU):
    
    
    t = floor((dist - 1)/2)
    t2 = floor((dist_small - 1)/2)
    
    #print(f"t = {t}, t2 = {t2}")
    
    rs_list = []
    
    count_r_for_worker = floor(1.0*(n1+1) / nCPU)

    
    for i in range(nCPU):
        r = []
        for j in range(count_r_for_worker):
        
            cur_r = i * count_r_for_worker + j
            
            if cur_r <= n1:
                r.append(cur_r)
                
        rs_list.append(r)
    
    
    cur_r = count_r_for_worker * nCPU
    for i in range(nCPU):
        if cur_r <= n1:
            rs_list[i].append(cur_r)
        else: 
            break
        cur_r += 1

    
    manager = Manager()
    result = manager.dict()
    processes = []
    for rs in rs_list:

        p = Process(target = calc_avg_bad_blocks_cnt_prob_worker, args = (t, t2, n1, n2, rs, result))
        p.start()
        processes.append(p)
    
    #print("All processes started")
    
    for process in processes:
        process.join()
    
    #print("All processes ended")
    
    
    all_comb = binomial(n1*n2, t)
    
    sum=0
    for r, p in result.items():
        sum += r*p
        #print(f"Prob of {r} = {p*1.0/all_comb}")
        
    res=sum*1.0/all_comb
    
    #print(f"avg_bad_blocks_count_prob = {res}")
    return ceil(res)


# вычисляет вероятность модифицированной атаки для заданного D-кода
def print_d_code_attack_prob(d_code, count_iter = 0):
    
    print(f"D-code: {d_code}")
    
    matr = get_d_code_gen_matr(d_code)
    dim = dim_d_code(d_code)
    dist = dist_d_code(d_code)
    
    if d_code[0][0] >= M_1 / 2:
        reverse_flg = false
    else:
        reverse_flg = true
                    
    print(f"reverse_flg is {reverse_flg}")
    
    if not reverse_flg:
        r_small = d_code[0][1]
        bad_blocks_count = floor((dist - 1) / (2 ** (M_2 - r_small) + 1))
        good_block_count = 2 ** M_1 - bad_blocks_count
        dim_small = dim_RM(r_small, M_2)
        all_blocks_count = 2 ** M_1
        block_length = 2 ** M_2
        dist_small = 2 ** (M_2 - r_small)
    else:
        r_small = d_code[-1][0]
        bad_blocks_count = floor((dist - 1) / (2 ** (M_1 - r_small) + 1))
        good_block_count = 2 ** M_2 - bad_blocks_count
        dim_small = dim_RM(r_small, M_1)
        all_blocks_count = 2 ** M_2
        block_length = 2 ** M_1
        dist_small = 2 ** (M_1 - r_small)
       
    print(f"dim = {dim}, dist = {dist}, dist_small = {dist_small}, r_small = {r_small}, bad_blocks_count = {bad_blocks_count}, good_block_count = {good_block_count}, dim_small = {dim_small}, all_blocks_count = {all_blocks_count}, block_length = {block_length}")
        
    if count_iter > 0:
        act_count_iter = count_iter
        blocks_count = good_block_count
    else:
        act_count_iter = 1000
        blocks_count = find_blocks_count(matr, dim, dim_small, all_blocks_count, bad_blocks_count, block_length, reverse_flg)
        
    
    
    P1_min = get_correct_prob(all_blocks_count, blocks_count, bad_blocks_count)
    P2 = float(count_full_rank_d_code(matr, dim, act_count_iter, all_blocks_count, blocks_count, block_length, false, reverse_flg)/act_count_iter)
    P_attack_min = P1_min * P2
    print(f"bad_blocks_count = {bad_blocks_count}, good_block_count = {good_block_count}, K = {blocks_count}, P1_min = {'{0:.3E}'.format(P1_min)}, P2 = {'{0:.3E}'.format(P2)}, P_attack_min = {P_attack_min}")


#Для заданного кода Рида-Маллера ноходим корректирующую способность декодера Сидельникова-Першакова
def print_sid_persh_decoder(r, m, init_lmda = 0, t_sid = 0, print_mode = True):
    
    dim = dim_RM(r, m)
    C = 0.5
    
    if init_lmda == 0 and t_sid == 0:    
        lmda = 2.0 * sqrt((2 ** r - 1) * ln(2) * dim)
        t_sid = floor((2 ** m - lmda * sqrt(2 ** m))/2)
        error_prob = C * 1.0 /(2 ** dim)
        sec_lvl = get_sec_lvl(2 ** m, dim, 2 * t_sid + 1)
        
    elif t_sid != 0:
        N = 2 ** m
        lmda = (N - 2*t_sid)/sqrt(N)
        eps = 1 / (2 ** r - 1)
        first_add = - (lmda ** 2) * eps / (2 * ln(2))
        #print(float(first_add))
        error_prob = float(C * 2 ** (first_add + dim))
        #error_prob = 1
        
        t_sid = floor((2 ** m - lmda * sqrt(2 ** m))/2)
        sec_lvl = get_sec_lvl(2 ** m, dim, 2 * t_sid + 1)
    else:
        lmda = init_lmda
        eps = 1 / (2 ** r - 1)
        first_add = - (lmda ** 2) * eps / (2 * ln(2))
        #print(float(first_add))
        error_prob = float(C * 2 ** (first_add + dim))
        #error_prob = 1
        
        t_sid = floor((2 ** m - lmda * sqrt(2 ** m))/2)
        sec_lvl = get_sec_lvl(2 ** m, dim, 2 * t_sid + 1)
        

    t_gar = floor((2 ** (m - r) - 1) / 2)
    
    # + 2720 * 2720 + 3488 * 3488
    # + 2 ** m * 2 ** m + dim * dim
    
    if print_mode:
        print(f"RM key size: {2 ** m * dim}")
        print(f"max_lambda = {floor(sqrt(2 ** m))}")
        print(f"RM({r}, {m}) - [{2 ** m}, {dim}, {2 ** (m - r)}]-code: t_gar = {t_gar}, t_sid = {t_sid} (lambda = {ceil(lmda)}), error_prob = {error_prob}, sec_lvl = {sec_lvl}")

    return t_sid, error_prob


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


# Для D-кодов из списка выводит корректирующую способность с использованием декодера Сидельникова-Першакова 
# для разных вероятностей ошибочного декодирования
def print_dcodes_sid_decoder():
    m1_list = [7,8]
    m2_list = [7,8]
    dcodes = [[(0, 3), (1, 2), (2, 1)],
              [(1, 2), (2, 1), (3, 0)],
              [(1, 2), (2, 1)]]
    
    #dcodes = [[(0, 2), (1, 1)], [(1, 1), (2, 0)]]
    
    r = 3
    min_error_prob =  0.001
    
    
    for m1 in m1_list:
        for m2 in m2_list:
            print("\n****************************")
            print(f"m1 = {m1}, m2 = {m2}")
            print("****************************\n")
            
            for dcode in dcodes:
                
                dim_dcode = dim_d_code(dcode, m1, m2)
                dist_dcode = dist_d_code(dcode, m1, m2)
                t_gar = floor((dist_dcode-1)/2)
                i_lambda = floor(sqrt(2 ** (m1 + m2)))
                t_sid, error_prob = print_sid_persh_decoder(3, m1+m2, i_lambda, False)
                key_size = dim_dcode * (2 ** (m1 + m2))
                
                print("\t==============================================================================================")
                print(f"\t\tD-code: {dcode}, dim = {dim_dcode}, t_gar = {t_gar}, key_size = {key_size}")
                print("\t==============================================================================================")
                
                while error_prob <= min_error_prob:
    
                    sec_lvl = get_sec_lvl(2 ** (m1+m2), dim_dcode, 2*t_sid+1)
        
                    if t_sid > t_gar and sec_lvl <= error_prob*(10 ** 7):
                        print(f"\tLambda: {i_lambda}, t_sid = {t_sid}, sec_lvl = {sec_lvl}, error_prob = {error_prob}")
                    
                    i_lambda -= 1
                    t_sid, error_prob = print_sid_persh_decoder(3, m1+m2, i_lambda, False)
                    
                print('\t----------------------------------------------------------------------------------------------')


# Для заданных параметров выводит список D-кодов, отсортированный по стойкости с использованием блочного декодера
def print_d_codes_block_decoder(r1, r2, m1, m2, t2, calc_blocks_flg=False):
    
    print(f"m1 = {m1}, m2 = {m2}, r2 = {r2}, t2 = {t2}")
    if calc_blocks_flg:
        print(f"All blocks count = {2^m1}")
    print('-------------------------------------')
    
    length = 2^(m1+m2)
    t_block = floor(t2*2^m1)
    
    dcodes = get_d_codes(r1=r1, r2=r2, m1=m1, m2=m2, dots_count=0, strength=True, sec_lvl_flg=False, print_mode=False)
    dcodes = sorted(list(dcodes), key = lambda element : get_sec_lvl(length, dim_d_code(element, m1, m2), t_block*2+1), reverse=True)
    
    for d_code in dcodes:
        
        dim = dim_d_code(d_code, m1, m2)
        dist = dist_d_code(d_code, m1, m2)
        t_gar = floor((dist-1)/2)
        pk_size = round(float((length * dim) / 8 / 1024 / 1024),3)
        sec_lvl_gar = get_sec_lvl(length, dim, dist)
        sec_lvl_block = get_sec_lvl(length, dim, t_block*2+1)
        print(f"d_code = {d_code}, [{length},{dim},{dist}], t_gar = {t_gar}, t_block = {t_block}, pk_size_MB = {pk_size}, \n sec_lvl_gar = {sec_lvl_gar}, sec_lvl_prob = {sec_lvl_block}")
        if calc_blocks_flg:
            min_cnt_blocks = find_blocks_count(matr=get_d_code_gen_matr(d_code, m1, m2, False), dim_full=dim_d_code(d_code, m1, m2), 
                            dim_small=dim_RM(r2,m2), all_blocks_count=2^m1, bad_blocks_count=2, 
                            block_length=2^m2, reverse_flg=False, print_mode=False)
            print(f"Min count blocks = {min_cnt_blocks}")
        #print(dcode_degrees(d_code, m1, m2, True))
        print('-------------------------------------')


# Вычисление вероятности P = P1*P2 для вероятностного декодера        
def get_block_decoder_probability(total_blocks, blocks_cnt, max_bad_blocks, p2):
    print(f"K={blocks_cnt} (P2 = {p2})")
    for j in range (1, max_bad_blocks+1):
        p1 = float(binomial(total_blocks-j, blocks_cnt)/binomial(total_blocks, blocks_cnt))
        print(f"({j} бл.) P1 = {p1}, P = {p1*p2}")


# Вычисление DFR вероятностного декодера
def get_block_decoder_dfr(K, blocks_cnt, p2, dfr_block, max_iter = 1000000):
    
    DFR_ISD_ITER = 1
    i=0
    while DFR_ISD_ITER > 10^(-9) and i <= max_iter:
        
        summ=0
        for j in range(1, blocks_cnt+1):
            p1 = binomial(blocks_cnt - j, K) / binomial(blocks_cnt, K)
            p_isd_j = p1*p2
            p_block_j = binomial(blocks_cnt, j) * dfr_block^j * (1 - dfr_block)^(blocks_cnt - j)
            summ += (1 - (1 - p_isd_j)^i)*p_block_j
        
        DFR_ISD_ITER = 1 - ((1 - dfr_block)^blocks_cnt + summ)

        if DFR_ISD_ITER < 10^(-8.95) or i == 0 or i == 3 or i == 7:
            print(f"i={i}\nDFR={DFR_ISD_ITER}")
        i+=1


# Вычисление вероятности появления плохих блоков
def get_prob_of_blocks(p, blocks_cnt, min_blocks_cnt):
    p_final = 0
    for i in range(min_blocks_cnt, blocks_cnt + 1):
        p_final += binomial(blocks_cnt, i) * p^i * (1 - p)^(blocks_cnt - i)
        
    print(f"DFR = {p_final}")        


# Вычисление времени работы вероятностных декодеров
def get_avg_clocks(p_bad, iter_max, p2, K, n1, c_dec, c_isd):
    
    cl=0
    
    for i in range(0, iter_max+1):
        p_iter = 0
        
        for b in range(0, n1):
            p1 = binomial(n1 - b, K) / binomial(n1, K)
            p_dec = binomial(n1, b) * p_bad ^ b * (1 - p_bad)^(n1 - b)
            p_iter += (1 - p1 * p2)^(i-1) * p1 * p2 * p_dec
        
        #print(f"i = {i}, p_iter = {p_iter}")
        cl += i * c_isd * p_iter
        
    cl += c_dec
    
    return ceil(cl)
    

