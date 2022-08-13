# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 18:06:27 2022

@author: user
"""

#Наработки по расчету температуры отжига через Nearest Neigbour
#встает вопрос, какие значения энтальпий и энтропий использовать
import math

def Tm_calculation_NN(comp_seq):
    deltaH=0
    deltaS=0
    conc=0.6*10**(-12)
    NN_deltaH_values={'AA': -9.1, 'TT': -9.1, 'AT': -8.6, 'TA': -6.0, 'CA': -5.8, 'AC': -6.5, 'GT': -6.5, 'TG': -5.8, 'CT': -7.8, 'TC': -5.6, 'GA': -5.6, 'AG': -7.8, 'CG': -11.9, 'GC': -11.1, 'GG': -11.0, 'CC': -11.0}
    NN_deltaS_values={'AA': -24.0, 'TT': -24.0, 'AT': -23.9, 'TA': -16.9, 'CA': -12.9, 'AC': -17.3, 'GT': -17.3, 'TG': -12.9, 'CT': -20.8, 'TC': -13.5, 'GA': -13.5, 'AG': -20.8, 'CG': -27.8, 'GC': -26.7, 'GG': -26.6, 'CC': -26.6}
    for i in range (1, len(comp_seq)):
        deltaH += NN_deltaH_values[comp_seq[i-1:i+1]]
        deltaS += NN_deltaS_values[comp_seq[i-1:i+1]]
    Tm=(-deltaH*1000-3400)/(-deltaS+1.9859*math.log(1/conc))+16.6*math.log10(0.05)-273
    return Tm


