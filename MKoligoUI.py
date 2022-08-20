#-*- coding: utf-8 -*-
"""
Created on Sun Jul 10 10:34:43 2022

@author: user
"""

#Программа для разбивки гена на олигонуклеотиды для синтеза при помощи ПЦР assembly

#######################################################################################


import math, io, os, random
from Bio.Seq import Seq
from PyQt5 import QtWidgets, uic #, QtCore, QtGui
from PyQt5.QtWidgets import QFileDialog, QMessageBox#, QCheckBox

class oligoPool():
    def __init__(self, oligopool):
        self.oligopool = oligopool.replace(" ", "")
        self.oligos = []

    def create(self):
        #oligos=[]
        list1=[]
        #i - длина первого олига, x+20 - длина всех остальных (20- комплементарный участок)
        for i in range (40, 61):
            for x in range (30, 41):
                if (len(self.oligopool)-i)%x == 0 and ((len(self.oligopool)-i)/x)%2!=0:
                    list1.append((x,i))
        #print(list1)
        x,i=list1[len(list1)-1]
        self.oligos.append(self.oligopool[0:i])
        k=i #положение в строке
        while k != len(self.oligopool):
            self.oligos.append(self.oligopool[k:k+x])
            k+=x
        for i in range (2, len(self.oligos), 2):
            self.oligos[i]=self.oligos[i-1][len(self.oligos[i-1])-20:len(self.oligos[i-1])]+self.oligos[i]
        for i in range (1, len(self.oligos), 2):
            self.oligos[i]=str(Seq(self.oligos[i-1][len(self.oligos[i-1])-20:len(self.oligos[i-1])]).complement()+Seq(self.oligos[i]).complement())[::-1]
        return self.oligos
        #return self.oligos

    #Функция для создания олигонуклеотидов
    '''def create(self):
        for i in range (2, len(self.oligos), 2):
            self.oligos[i]=self.oligos[i-1][len(self.oligos[i-1])-20:len(self.oligos[i-1])]+self.oligos[i]
        for i in range (1, len(self.oligos), 2):
            self.oligos[i]=str(Seq(self.oligos[i-1][len(self.oligos[i-1])-20:len(self.oligos[i-1])]).complement()+Seq(self.oligos[i]).complement())[::-1]
        return self.oligos'''

class Tm():
    def __init__(self, oligos):
        self.oligos = oligos.replace(" ", "")

    def Tm_calculation(self):   #ПОЧЕМУ-ТО НЕ РАБОТАЕТ!!! НАДО РАЗАБРАТЬСЯ!!!
        Tms = []
        Tm = 0
        print(self.oligos)
        for i in range(0, len(self.oligos)-1):
            if i%2 == 0:
                wA = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('A')
                xT = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('T')
                yG = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('G')
                zC = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('C')
            else:
                self.oligos[i] = self.oligos[i][::-1]
                wA = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('A')
                xT = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('T')
                yG = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('G')
                zC = self.oligos[i][len(self.oligos[i])-20:len(self.oligos[i])].count('C')
                self.oligos[i] = self.oligos[i][::-1]
            Tm = 100.5+(41*(yG+zC)/(wA+xT+yG+zC))-(820/(wA+xT+yG+zC)) + 16.6*math.log10(0.05)
            Tms.append(Tm)
            print(Tms)
            #print('Tm = ', "{:.1f}".format(Tm))
        #print(Tms)  
        #return(Tms)
        
    def Tm_calculation_NN(self):
        #print('Температуры отжига комплементарных участков: ')
        Tms = []
        deltaH = 0
        deltaS = 0
        conc = 0.6*10**(-12)
        NN_deltaH_values = {'AA': -9.1, 'TT': -9.1, 'AT': -8.6, 'TA': -6.0, 'CA': -5.8,
                            'AC': -6.5, 'GT': -6.5, 'TG': -5.8, 'CT': -7.8, 'TC': -5.6,
                            'GA': -5.6, 'AG': -7.8, 'CG': -11.9, 'GC': -11.1, 'GG': -11.0, 'CC': -11.0}
        NN_deltaS_values = {'AA': -24.0, 'TT': -24.0, 'AT': -23.9, 'TA': -16.9, 'CA': -12.9,
                            'AC': -17.3, 'GT': -17.3, 'TG': -12.9, 'CT': -20.8, 'TC': -13.5,
                            'GA': -13.5, 'AG': -20.8, 'CG': -27.8, 'GC': -26.7, 'GG': -26.6, 'CC': -26.6}
        for i in range (1, len(self.oligos)):
            deltaH += NN_deltaH_values[self.oligos[i-1:i+1]]
            deltaS += NN_deltaS_values[self.oligos[i-1:i+1]]
        Tm = (-deltaH*1000-3400)/(-deltaS+1.9859*math.log(1/conc))+16.6*math.log10(0.05)-273
        #print(Tm)
        return Tm

class codonOpt():
    def __init__(self, seq):
        self.seq = seq

    def optimaze(self):
        if len(self.seq) < 3:
            return("Не задана последовательность олигонуклеотида")
        else:
            self.seq = self.seq.replace(" ", "")
            s = ' '.join([self.seq[i:i + 3] for i in range(0, len(self.seq), 3)])
            s=s.replace('CTC', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'],  weights=[0.21, 0.16, 0.18, 0.45])))
            s=s.replace('CTA', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'],  weights=[0.21, 0.16, 0.18, 0.45])))
            s=s.replace('TCC', ''.join(random.choices(['TCT', 'TCA'],  weights=[0.5, 0.5])))
            s=s.replace('TCG', ''.join(random.choices(['TCT', 'TCA'],  weights=[0.5, 0.5])))
            s=s.replace('CCC', ''.join(random.choices(['CCT', 'CCA', 'CCG'],  weights=[0.29, 0.27, 0.44])))
            s=s.replace('CAC', 'CAT')
            s=s.replace('CGA', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
            s=s.replace('CGG', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
            s=s.replace('AGA', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
            s=s.replace('AGG', ''.join(random.choices(['CGT', 'CGC'],  weights=[0.53, 0.47])))
            return(s)
        

class seqAnalysis():
    def __init__(self, GCcontent):
        self.GCcontent = GCcontent.replace(" ", "")

    def obrabotka(self):
        if len(self.GCcontent)==0:
            return("Не задана последовательность олигонуклеотида")
        else:
            return("GC состав: " + "{0:.2f}".format((self.GCcontent.count('G')+self.GCcontent.count('C')*100)/len(self.GCcontent)) + "%")
        
class Window(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        uic.loadUi("MainWindow.ui", self)
        self.pushButton_6.clicked.connect(self.loadFile)
        self.pushButton.clicked.connect(self.codonOptimize)
        self.analiseseq_pushButton.clicked.connect(self.analiseSeq)
        self.run_pushButton.clicked.connect(self.runScript)

    def popupwin(self, text, title):
        infoBox = QMessageBox()
        infoBox.setIcon(QMessageBox.Information)
        infoBox.setText(text)
        infoBox.setWindowTitle(title)
        infoBox.setStandardButtons(QMessageBox.Ok)
        infoBox.exec_()
        
    def loadFile(self):
        self.way = QFileDialog.getOpenFileName(None, 'Open File', './', "Fasta (*.fa *.fasta);;CSV Files (*.csv)")
        self.textEdit.setText(self.way[0])
        if self.way[0][-2:] == 'fa' or self.way[0][-5:] == 'fasta':
            self.sequence = io.open(self.way[0], mode='r').read().split()
        else:
            return 0
        self.plainTextEdit.clear()
        self.plainTextEdit.insertPlainText(self.sequence[1].upper())
        self.checkOligo()
        

    def checkOligo(self):
        if len(self.plainTextEdit.toPlainText()) < 120:
            self.popupwin("Lenght sequence is to shirt","Info")
        elif len(self.plainTextEdit.toPlainText()) > 800:
            self.popupwin("Lenght sequence is to shirt","Info")
        else:
            return 0
       
    '''def markOligo(self):
        #def cut(s):
        #self.checkOligo()
        self.oligos = []
        list1=[]
        #i - длина первого олига, x+20 - длина всех остальных (20- комплементарный участок)
        for i in range (40, 61):
            for x in range (30, 41):
                if (len(s)-i)%x == 0 and ((len(s)-i)/x)%2!=0:
                    list1.append((x,i))
        #print(list1)
        x,i=list1[len(list1)-1]
        oligos.append(s[0:i])
        k=i #положение в строке
        while k != len(s):
            oligos.append(s[k:k+x])
            k+=x
        print(self.oligos)
        for i in range (2, len(oligos), 2):
            oligos[i]=oligos[i-1][len(oligos[i-1])-20:len(oligos[i-1])]+oligos[i]
        for i in range (1, len(oligos), 2):
            oligos[i]=str(Seq(oligos[i-1][len(oligos[i-1])-20:len(oligos[i-1])]).complement()+Seq(oligos[i]).complement())[::-1]
        return oligos'''
        
    def codonOptimize(self):
        self.popupwin(codonOpt(self.plainTextEdit.toPlainText()).optimaze(),"Info")
        optimazeSeq = codonOpt(self.plainTextEdit.toPlainText()).optimaze()
        self.plainTextEdit.clear()
        self.plainTextEdit.insertPlainText(optimazeSeq)
        pass

    def analiseSeq(self):
        print(Tm(self.plainTextEdit.toPlainText()).Tm_calculation_NN())
        self.popupwin(seqAnalysis(self.plainTextEdit.toPlainText()).obrabotka(),"Info")
        seqAnalysis(self.plainTextEdit.toPlainText()).obrabotka()
        pass

    def runScript(self):
        oligos = oligoPool(self.plainTextEdit.toPlainText()).create()
        print(oligos)
        for oligo in oligos:
            if (oligos.index(oligo)+1)%2 !=0:
                print(oligos.index(oligo)+1, 'f: ', oligo) 
            else:
                print(oligos.index(oligo)+1, 'b: ', oligo)
        #Tms=Tm_calculation(oligos)
        #pretty_output(oligos, Tms)
        print("Олигонуклеотиды с областями комплиментарности и температурами отжига представлены в файле pretty_output.txt")
        
        #oligos=[hairpin(oligo, oligos.index(oligo), sequence, 5) for oligo in oligos]
        #oligos=[hairpin(oligo[::-1], oligos.index(oligo), sequence, 3)[::-1] for oligo in oligos]
        print('')
        print('Олигонуклеотиды, исправленные с учетом шпилек:')
        for oligo in oligos:
            if (oligos.index(oligo)+1)%2 !=0:
                print(oligos.index(oligo)+1, 'f: ', oligo) 
            else:
                print(oligos.index(oligo)+1, 'b: ', oligo)
        pass




if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = Window()
    w.show()
    sys.exit(app.exec_())
