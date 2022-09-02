

def cut_x(sequence, cut_len): # алгоритм достаточно быстрой предварительной нарезки олигов
    seq_len = len(sequence)

    # подсчет длинны олигов
    # выполняется за 1 или 2 итерации, в основном числе случаев
    count_list = []
    while True:
        oligos_count = seq_len // cut_len
        tail_lenght = seq_len % cut_len
        c = tail_lenght / oligos_count
        if c == 0: # тот случай когда длинна олига cut_len без остатка укладывается на строке sequence
            count_dict = {'oligo_lenght': cut_len, 'oligo_count': oligos_count}
            count_list.append(count_dict)
            break
        elif c < 1 and c > 0: # тот случай когда остаток от целого числа уложенных олигов delta меньше числа олигов
                            # в этом случае увеличиваем на 1 длинну нарезки у delta олигов
            delta = tail_lenght % oligos_count
            count_dict = {'oligo_lenght': cut_len + 1, 'oligo_count': delta}
            count_list.append(count_dict)
            count_dict = {'oligo_lenght': cut_len, 'oligo_count': oligos_count - delta}
            count_list.append(count_dict)
            break
        else: # случай когда остаток от целого числа уложенных олигов больше числа олигов
            # в этом случае увеличиваем длинну нарезки и повторяем рассчет
            cut_len += 1

    #быстрая разбивка всей строки на подстроки
    if len(count_list) == 1:
        out_list = list(map(''.join, zip(*[iter(sequence)] * count_list[0]['oligo_lenght'])))
        return out_list
    else:
        len_seq1 = count_list[0]['oligo_lenght'] * count_list[0]['oligo_count']
        s1, s2 = sequence[:len_seq1], sequence[len_seq1 :]
        out_list = list(map(''.join, zip(*[iter(s1)] * count_list[0]['oligo_lenght'])))
        out_list.extend(list(map(''.join, zip(*[iter(s2)] * count_list[1]['oligo_lenght']))))
        return out_list

def test_cut_x():
    seq = 'ACTGCTGTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCC'
    print(seq)
    print(cut_x(seq, 11))

if __name__ == '__main__':
    test_cut_x()