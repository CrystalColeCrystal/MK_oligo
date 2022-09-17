#  В модуле представлена система классов для реализации алгоритмов рассчета ТД параметров вторичных структур ДНК олигов
#  Реализован и валидирован метод рассчета температур плавления олигов методом ближайших соседей
#  В коментариям к классам указаны ссылки на публикации, из которых взяты экспериментальные ТД параметры

#import oligoUtils
import math

class oligoSequenceBase():

    # класс нужен для манипуляции с сиквенсами
    # реализованы такие операции как:
    # реверс, комплемент, проверка на самокомплементарность

    def __init__(self, sequence):
        self.sequence = sequence

    def reverse(self):
        return self.sequence[::-1]

    def complement(self):
        complement = self.sequence
        complement = complement.replace('A', 't')
        complement = complement.replace('T', 'a')
        complement = complement.replace('C', 'g')
        complement = complement.replace('G', 'c')
        return complement.upper()

    def reverse_complement(self):
        return oligoSequenceBase(self.reverse()).complement()

    def is_self_complement(self):
        return self.sequence == self.reverse_complement()


class oligoNearestNeighborTDParams():

    #  родительский класс для семейства классов интерфейсов с эксперементальными данными
    #   реализованы базовые методы, которые наследуются потомками класса
    #

    def __init__(self):
        self._init_params()

    def _init_params(self):
        pass

    def _getKey(self, NNkey):
        key = NNkey
        if key not in self._keys:
            key = f'{NNkey}|{oligoSequenceBase(NNkey).complement()}'
            if key not in self._keys:
                _NNseq = oligoSequenceBase(NNkey).reverse()
                key = f'{oligoSequenceBase(_NNseq).complement()}|{_NNseq}'
                if key not in self._keys:
                    key = 'error'
        return key

    def deltaH(self, NNseq):
        return self._params[self._getKey(NNseq)]['dH']

    def deltaS(self, NNseq):
        return self._params[self._getKey(NNseq)]['dS']

    def deltaG37(self, NNseq):
        return self._params[self._getKey(NNseq)]['dG37']



class oligos_Nearest_neighbor_Watson_Crick_SantaLucia(oligoNearestNeighborTDParams):

    #  сюда загружаются ТД параметры для рассчета температуры плавления олига методом ближ.соседей

    def _init_params(self):
        self._params = {   # kcal / mol
                            # John SantaLucia, Jr. and Donald Hicks // doi: 10.1146/annurev.biophys.32.110601.141800
                            # Nearest-neighbor thermodynamic parameters for DNA
                            # Watson-Crick pairs in 1 M NaCl
                            'AA|TT': {'dH': -7.6, 'dS': -21.3, 'dG37': -1.00},
                            'AT|TA': {'dH': -7.2, 'dS': -20.4, 'dG37': -0.88},
                            'TA|AT': {'dH': -7.2, 'dS': -21.3, 'dG37': -0.58},
                            'CA|GT': {'dH': -8.5, 'dS': -22.7, 'dG37': -1.45},
                            'GT|CA': {'dH': -8.4, 'dS': -22.4, 'dG37': -1.44},
                            'CT|GA': {'dH': -7.8, 'dS': -21.0, 'dG37': -1.28},
                            'GA|CT': {'dH': -8.2, 'dS': -22.2, 'dG37': -1.30},
                            'CG|GC': {'dH': -10.6, 'dS': -27.2, 'dG37': -2.17},
                            'GC|CG': {'dH': -9.8, 'dS': -24.4, 'dG37': -2.24},
                            'GG|CC': {'dH': -8.0, 'dS': -19.9, 'dG37': -1.84},
                            'initiation': {'dH': 0.2, 'dS': -5.7, 'dG37': 1.96},
                            'terminal AT': {'dH': 2.2, 'dS': 6.9, 'dG37': 0.05},
                            'Symmetry correction': {'dH': 0.0, 'dS': -1.4, 'dG37': 0.43},
                            'error': {'dH': 0, 'dS': 0, 'dG37': 0}
                         }
        self._keys = list(self._params.keys())


class oligos_Nearest_neighbor_Watson_Crick_SantaLucia_2(oligoNearestNeighborTDParams):

    #  сюда загружаются ТД параметры для рассчета температуры плавления олига методом ближ.соседей

    def _init_params(self):
        self._params = {   # kcal / mol
                            # John SantaLucia, Jr.,* Hatim T. Allawi, and P. Ananda Seneviratne //  DOI: 10.1021/bi951907q
                            # Biochemistry 1996, 35, 3555-3562 Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability
                            # Nearest-neighbor thermodynamic parameters for DNA
                            # Watson-Crick pairs in 1 M NaCl
                            'AA|TT': {'dH': -8.4, 'dS': -23.6, 'dG37': -1.02},
                            'AT|TA': {'dH': -6.5, 'dS': -18.8, 'dG37': -0.73},
                            'TA|AT': {'dH': -6.3, 'dS': -18.5, 'dG37': -0.60},
                            'CA|GT': {'dH': -7.4, 'dS': -19.3, 'dG37': -1.38},
                            'GT|CA': {'dH': -8.6, 'dS': -23.0, 'dG37': -1.43},
                            'CT|GA': {'dH': -6.1, 'dS': -16.1, 'dG37': -1.16},
                            'GA|CT': {'dH': -7.7, 'dS': -20.3, 'dG37': -1.46},
                            'CG|GC': {'dH': -10.1, 'dS': -25.5, 'dG37': -2.09},
                            'GC|CG': {'dH': -11.1, 'dS': -28.4, 'dG37': -2.28},
                            'GG|CC': {'dH': -6.7, 'dS': -15.6, 'dG37': -1.77},
                            'initiation': {'dH': 0.0, 'dS': -5.9, 'dG37': 1.82},
                            'terminal AT': {'dH': 0.4, 'dS': 0.0, 'dG37': 0.4},
                            'Symmetry correction': {'dH': 0.0, 'dS': -1.4, 'dG37': 0.4},
                            'init GC': {'dH': 0.0, 'dS': -5.9, 'dG37': 1.82},
                            'init AT': {'dH': 0.0, 'dS': -9.0, 'dG37': 2.8},
                            'error': {'dH': 0, 'dS': 0, 'dG37': 0}
                         }
        self._keys = list(self._params.keys())


class dnaSeq2NNKeys():

    # Данный класс требуется для генерации правильных ключей
    # при обращении к табличным параметрам
    # Класс выполняет роль посредника между табличными данными и классом, который вычисляет непосредственно температуру плавления
    # Данный класс будет расширятся или дифференцироваться в дальнейшем при реализации рассчетов других элементов вторичной структуры

    def __init__(self, sequence, sequence2=''):
        self.sequence = sequence
        self.sequence2 = sequence2

    def __getSimplekeys(self):
        keys = {}
        for i in range(1, len(self.sequence)):
            pair = self.sequence[i-1:i+1]
            if pair not in list(keys.keys()):
                keys[pair] = 1
            else:
                keys[pair] += 1

        keys['initiation'] = 0
        keys['terminal AT'] = 0
        keys['Symmetry correction'] = 0
        keys['init GC'] = 0
        keys['init AT'] = 0
        return keys

    def __getitem__(self, item):
        if item in ['simpleKey', 'SimpleKeys', 'simpleKeys', 'simple', 'Simple']:
            return self.__getSimplekeys()
        else:
            return []


class dnaOligoNNTm():
    # NN_SantaLucia, NN_SantaLucia_2

    # Класс вычисляет температуру плавления олига по сиквенсе
    # пока 2 параметра метода рассчета: NN_SantaLucia, NN_SantaLucia_2 соответствуют таблицам параметров из 2-х публикаций
    # эти 2 метода хорошо коррелируют между собой.
    # хотя NN_SantaLucia более точный в сравнении с экспериментом

    def __init__(self, method='NN_SantaLucia_2'):
        self.method = method

    def _getSantaLucia_key_1(self, sequence, keys):
        keys['initiation'] = 1
        if sequence[0] in ['A', 'T']:
            keys['terminal AT'] = 1
        if sequence[-1] in ['A', 'T']:
            keys['terminal AT'] += 1

    def _getSantaLucia_key_2(self, sequence, keys):
        if sequence.count('C') > 0 or sequence.count('G') > 0:
            keys['init GC'] = 1
        else:
            keys['init AT'] = 1
        if sequence[0] in ['A', 'T']:
            keys['terminal AT'] = 1
        elif sequence[-1] in ['A', 'T']:
            keys['terminal AT'] += 1

    def _santaLucia_NN(self, sequence):

        NN_Keys = dnaSeq2NNKeys(sequence) # раскладываем исходную сиквенсу на пары букв и другие ключи
        keys = NN_Keys['simple']

        if self.method == 'NN_SantaLucia':
            self._getSantaLucia_key_1(sequence, keys)
            NN_santaLucia = oligos_Nearest_neighbor_Watson_Crick_SantaLucia()
        elif self.method == 'NN_SantaLucia_2':
            self._getSantaLucia_key_2(sequence, keys)
            NN_santaLucia = oligos_Nearest_neighbor_Watson_Crick_SantaLucia_2()

        if oligoSequenceBase(sequence).is_self_complement():
            Cm = 1e-4 # если олиг самокомплементарен
            keys['Symmetry correction'] = 1
        else:
            Cm = 4e-4 # если олиг не самокомплементарен
            keys['Symmetry correction'] = 0

        R = 1.9872 # универсальная газовая постоянная

        deltaH = sum([NN_santaLucia.deltaH(key) * number for key, number in zip(keys.keys(), keys.values())])
        deltaS = sum([NN_santaLucia.deltaS(key) * number for key, number in zip(keys.keys(), keys.values())])

        Tm = deltaH * 1000 / (deltaS + R * math.log(Cm)) - 273.15

        return Tm

    def Tm(self, sequence):
        return self._santaLucia_NN(sequence)


class dnaOligoTm():

    # Простые модели рассчета температур плавления
    # "плохо" коррелируют с экспериментом

    def __init__(self, method='second'):
        self.method = method

    def __set_count(self, sequence):
        self.wA = sequence.count('A')
        self.xT = sequence.count('T')
        self.yG = sequence.count('G')
        self.zC = sequence.count('C')

    def __Tm_first(self, sequence):
        self.__set_count(sequence)
        if len(sequence) < 14:
            Tm = 2 * (self.wA + self.xT) + 4 * (self.yG + self.zC)
        else:
            Tm = 64.9 + 41 * (self.yG + self.zC - 16.4) / (self.wA + self.xT + self.yG + self.zC)
        return Tm

    def __Tm_second(self, sequence):
        self.__set_count(sequence)
        Z = self.wA + self.xT + self.yG + self.zC
        if len(sequence) < 14:
            Tm = 2 * (self.wA + self.xT) + 4 * (self.yG + self.zC)
        else:
            Tm = 100.5 + (41 * (self.yG + self.zC) / Z) - 820 / Z + 16.6 * math.log10(0.05)
        return Tm

    def Tm(self, sequence):
        if self.method == 'first':
            return self.__Tm_first(sequence)
        elif self.method == 'second':
            return  self.__Tm_second(sequence)
        else:
            return 0.


def test():
    NN_santaLucia = oligos_Nearest_neighbor_Watson_Crick_SantaLucia()

    print(NN_santaLucia.deltaG37('Symmetry correction'))
    print(NN_santaLucia.deltaG37('AG'))

    # CGTTGA

    dG = NN_santaLucia.deltaG37('initiation') + NN_santaLucia.deltaG37('Symmetry correction') + \
         NN_santaLucia.deltaG37('CG') + NN_santaLucia.deltaG37('GT') + NN_santaLucia.deltaG37('TT') + \
         NN_santaLucia.deltaG37('TG') + NN_santaLucia.deltaG37('GA') +  NN_santaLucia.deltaG37('terminal AT')

    print(dG)

def test2():
    #NNseq = dnaSeq2NNKeys('CCGCGG')
    #NNseq = dnaSeq2NNKeys('CGATCG')
    #NNseq = dnaSeq2NNKeys('CCGG')
    #NNseq = dnaSeq2NNKeys('GCGC')
    #NNseq = dnaSeq2NNKeys('CAAGCTTG')
    #NNseq = dnaSeq2NNKeys('GTTGCAAC')
    #NNseq = dnaSeq2NNKeys('CATATGGCCATATG')
    #NNseq = dnaSeq2NNKeys('CTGACAAGTGTC')
    NNseq = dnaSeq2NNKeys('CGTCGACG')
    print(NNseq['simple'])

    NN_santaLucia = oligos_Nearest_neighbor_Watson_Crick_SantaLucia()

    keys = NNseq['simple']
    deltaH = sum([NN_santaLucia.deltaH(key) for key in keys])
    deltaS = sum([NN_santaLucia.deltaS(key) for key in keys])
    deltaG37 = sum([NN_santaLucia.deltaG37(key) for key in keys])
    print(deltaH, deltaS, deltaG37)
    R = 1.9872
    Tm = deltaH * 1000 / (deltaS + R * math.log(0.0001)) - 273.15
    print(Tm)

def test3():
    #seq = 'CGATCG'
    seq = 'GCTAGC'
    #seq = 'CGTTGA'
    seq = 'CGCGTACGCGTACGCG'
    seq = 'CCGG'

    NN_Tm = dnaOligoNNTm(method='NN_SantaLucia_2')
    print('NN Tm', NN_Tm.Tm(seq))

    base_Tm = dnaOligoTm(method='first')
    print('base Tm', base_Tm.Tm(seq))

def test_dataset():
    # data from: Table 2: Experimental and Predicted Thermodynamic Parameters of Duplex Formation
    # John SantaLucia, Jr.,* Hatim T. Allawi, and P. Ananda Seneviratne //  DOI: 10.1021/bi951907q
    data = [
        {'seq': 'CCGG', 'Tm_exp': 16.6},
        {'seq': 'CGCG', 'Tm_exp': 23.7},
        {'seq': 'GCGC', 'Tm_exp': 27.5},
        {'seq': 'CCGCGG', 'Tm_exp': 55.2},
        {'seq': 'CGATCG', 'Tm_exp': 34.3},
        {'seq': 'CGCGCG', 'Tm_exp': 55.7},
        {'seq': 'CGGCCG', 'Tm_exp': 59.6},
        {'seq': 'CGTACG', 'Tm_exp': 35.0},
        {'seq': 'GACGTC', 'Tm_exp': 36.1},
        {'seq': 'GCATGC', 'Tm_exp': 36.5},
        {'seq': 'GCCGGC', 'Tm_exp': 57.7},
        {'seq': 'GCGAGC', 'Tm_exp': 33.2},
        {'seq': 'GCGCGC', 'Tm_exp': 56.1},
        {'seq': 'GCTAGC', 'Tm_exp': 34.3},
        {'seq': 'GGATCC', 'Tm_exp': 30.8},
        {'seq': 'GGCGCC', 'Tm_exp': 53.9},
        {'seq': 'GGGACC', 'Tm_exp': 44.9},
        {'seq': 'GTGAAC', 'Tm_exp': 33.2},
        {'seq': 'CAAAAAG', 'Tm_exp': 31.5},
        {'seq': 'CAAAAAAG', 'Tm_exp': 36.9},
        {'seq': 'CAAGCTTG', 'Tm_exp': 44.6},
        {'seq': 'CATCGATG', 'Tm_exp': 47.4},
        {'seq': 'CGATATCG', 'Tm_exp': 44.1},
        {'seq': 'CGTCGACG', 'Tm_exp': 58.3},
        {'seq': 'GGAGCTCC', 'Tm_exp': 54.4},
        {'seq': 'CAACTTGATATTATTA', 'Tm_exp': 58.8},
        {'seq': 'AAAAAAAA', 'Tm_exp': 31.2},

        {'seq': 'GTATACCGGTATAC', 'Tm_exp': 61.9},
        {'seq': 'CATATTGGCCAATATG', 'Tm_exp': 65.3},
        {'seq': 'GTATAACCGGTTATAC', 'Tm_exp': 65.9},
    ]

    return data

def validate():
    data = test_dataset()
    results = []

    NN_Tm1 = dnaOligoNNTm(method='NN_SantaLucia')
    NN_Tm2 = dnaOligoNNTm(method='NN_SantaLucia_2')
    base_Tm = dnaOligoTm(method='first')

    for d in data:
        dd = {}
        dd['seq'] = d['seq']
        dd['exp_Tm'] = d['Tm_exp']
        dd['SL_1_Tm'] = NN_Tm1.Tm(d['seq'])
        dd['SL_2_Tm'] = NN_Tm2.Tm(d['seq'])
        dd['base_Tm'] = base_Tm.Tm(d['seq'])
        results.append(dd)

    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.metrics import r2_score

    df = pd.DataFrame(results)

    #print(df['SL_2_Tm'].corr(df['exp_Tm']))
    print(r2_score(df['exp_Tm'], df['SL_2_Tm'])) # 0.8477589738328025
    #print(df['SL_1_Tm'].corr(df['exp_Tm']))
    print(r2_score(df['exp_Tm'], df['SL_1_Tm'])) # 0.861240183665382
    #print(df['SL_2_Tm'].corr(df['SL_1_Tm']))
    print(r2_score(df['SL_1_Tm'], df['SL_2_Tm'])) # 0.9745704234917622
    #print(df['base_Tm'].corr(df['exp_Tm']))
    print(r2_score(df['exp_Tm'], df['base_Tm'])) # -1.8855817791494016

    plt.scatter(df['SL_2_Tm'], df['exp_Tm'])
    #plt.scatter(df['base_Tm'], df['exp_Tm'])
    #plt.scatter(df['SL_2_Tm'], df['SL_1_Tm'])
    #plt.scatter(df['SL_2_Tm'], df['SL_1_Tm'])
    plt.show()



if __name__ == '__main__':
    #test3()
    validate()