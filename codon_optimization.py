#Наработки по кодоновой оптимизации генов белков эукариот для синтеза в E. coli

#######################################################################################

import random
from Bio.Seq import Seq

######################################################################################

#Функция, разбивающая код на кодоны (возвращает строку)
def codone_splitter(s):
    return ' '.join([s[i:i + 3] for i in range(0, len(s), 3)])
#Функция, заменяющая редкие для E. coli кодоны на частые. 
#Для нее важно наличие пробелов в строке (строка, выдаваемая codone_splitter), так как иначе собьется рамка считывания
def codone_optimizer(s):
    #leu=['TTA', 'TTG', 'CTT', 'CTG']
    s=s.replace('CTC', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'],  weights=[0.21, 0.16, 0.18, 0.45])))
    s=s.replace('CTA', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'],  weights=[0.21, 0.16, 0.18, 0.45])))
    #ser=['TCT', 'TCA']
    s=s.replace('TCC', ''.join(random.choices(['TCT', 'TCA'],  weights=[0.5, 0.5])))
    s=s.replace('TCG', ''.join(random.choices(['TCT', 'TCA'],  weights=[0.5, 0.5])))
    #pro=['CCT', 'CCA', 'CCG']
    s=s.replace('CCC', ''.join(random.choices(['CCT', 'CCA', 'CCG'],  weights=[0.29, 0.27, 0.44])))
    #his
    s=s.replace('CAC', 'CAT')
    #arg=['CGT', 'CGC']
    s=s.replace('CGA', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
    s=s.replace('CGG', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
    s=s.replace('AGA', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
    s=s.replace('AGG', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
    return s
#Функция для подсчета GC состава
def GC_calculate(s):
    return (s.count('G')+s.count('C'))/len(s)
#Функция для поиска мест вхождений редких кодонов
def find_rare_codone(s):
    rare=['CTC', 'CTA', 'TCC', 'TCG', 'CCC', 'CAC', 'CGA', 'CGG', 'AGA', 'AGG']
    for codone in rare:
        res = [i for i in range(len(s)) if s.startswith(codone, i)]
        print(res)
    return res
#Функция для трансляции и записи в строку и вывода
def translate(s):
    seq = Seq(s)
    aa = ' '+'   '.join(list((seq.translate())))+' ' 
    #А зачем я здесь сначала преобразовываю seq в список, а потом список в строку?
    dna = codone_splitter(s)
    for i in range (0, len(dna), 64):
        print(dna[i:i+64])
        print(aa[i:i+64])
#Функция для выведения двух последовательнойстей друг под другом и выеделения замен цветом
def compare(s1,s2):
    k=1
    for i in range (0, len(s1), 64):
        print(k, ': ', s1[i:i+64])
        print(k, ': ', s2[i:i+64])
        k+=1


############################################################################################

with open('sequence.txt', 'r') as f:
    sequence = ''.join(f.readlines()).replace('\n', '') #cчитываем строку из файла без \n

############################################################################################
print('Здравствуйте! Последовательность mRNA белка эукариот (в кодировке A,T,G,C без пробелов), записанная в файле sequence.txt будет оптимизирована для экспресии в E. coli')
print('Аминокислотная последовательность оригинальной ДНК:')
print(translate(sequence))
print()
sequence_optimized= codone_optimizer(codone_splitter(sequence)).replace(' ','')
print('Аминокислотная последовательность оптимизированной ДНК:')
print(translate(sequence_optimized))
print()
print('Сравнение старой и новой последовательности:')
compare(codone_splitter(sequence), codone_splitter(sequence_optimized))
print()
print('GC состав до:', "{0:.2f}".format(GC_calculate(sequence)))
print('GC состав после:', "{0:.2f}".format(GC_calculate(sequence_optimized)))
