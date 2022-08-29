#-*- coding: utf-8 -*-
"""
Created on Sun Jul 10 10:34:43 2022

@author: user
"""

#Программа для разбивки гена на олигонуклеотиды для синтеза при помощи ПЦР assembly

#######################################################################################


import math, io, os, random
from Bio.Seq import Seq
from PyQt5 import QtWidgets, uic , QtCore, QtGui
from PyQt5.uic import *
from PyQt5.QtWidgets import * 
from PyQt5.QtCore import *
from PyQt5.QtGui import *

class menuAbout(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = uic.loadUi("about.ui") 
        self.ui.exec()

class baseCodes(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = uic.loadUi("base.ui") 
        self.ui.exec()

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


'''def loop_increment(loop_size):
    loop_increments={3: 7.4, 4: 5.9, 5: 4.4, 6: 4.3, 7: 4.1, 8: 4.1, 9: 4.2, 10: 4.3, 11: 4.6, 12: 4.9, 13: 5.25, 14: 5.6, 15: 5.85, 16: 6.1, 17: 6.4, 18: 6.7, 19: 6.9, 20: 7.1}
    if loop_size<21: 
        return loop_increments[loop_size]
    else:
        return 4.2545*math.log(loop_size)-5.6243

#Функция для рассчета температуры плавления шпильки
def hairpin_Tm(stem_sequence, loop_size):
    deltaH=0
    deltaS=0
    NN_deltaH_values={'AA': -9.1, 'TT': -9.1, 'AT': -8.6, 'TA': -6.0, 'CA': -5.8, 'AC': -6.5, 'GT': -6.5, 'TG': -5.8, 'CT': -7.8, 'TC': -5.6, 'GA': -5.6, 'AG': -7.8, 'CG': -11.9, 'GC': -11.1, 'GG': -11.0, 'CC': -11.0}
    NN_deltaS_values={'AA': -24.0, 'TT': -24.0, 'AT': -23.9, 'TA': -16.9, 'CA': -12.9, 'AC': -17.3, 'GT': -17.3, 'TG': -12.9, 'CT': -20.8, 'TC': -13.5, 'GA': -13.5, 'AG': -20.8, 'CG': -27.8, 'GC': -26.7, 'GG': -26.6, 'CC': -26.6}
    for i in range (1, len(stem_sequence)):
        deltaH += NN_deltaH_values[stem_sequence[i-1:i+1]]
        deltaS += NN_deltaS_values[stem_sequence[i-1:i+1]]
    deltaG_loop_incriment=loop_increment(loop_size)
    Tm=1000*((deltaH+deltaG_loop_incriment)/deltaS)-273.15
    return Tm

def hairpin(oligo, index, sequence, flag):
    there_is_a_hairpin=0
    for i in range (10, len(oligo)+1):
        if str(Seq(oligo[:3][::-1]).complement()) == oligo[i-3:i]:
            stem=oligo[:3]
            loop_size=i-6
            for k in range (3, i-8):
                if  str(Seq(oligo[k]).complement()) == oligo[i-k-1]:
                    stem+=oligo[k]
                    loop_size-=2
                else:
                    break
            if hairpin_Tm(stem, loop_size)>19:
                there_is_a_hairpin=1
                print('Обнаружена шпилька в олигонуклеотиде номер ', index+1, ', температура плавления: ' "{:.1f}".format(hairpin_Tm(stem, loop_size)), ',проводится корректировка')
                if (flag==5 and index%2==0):
                    oligo= sequence[sequence.index(oligo)-1]+oligo
                    #print(1, index)
                elif (flag==5 and index%2!=0):
                    oligo=str(Seq(sequence[sequence.index(str(Seq(oligo[::-1]).complement()))+len(oligo)]).complement())+oligo
                    #print(2, index)
                elif (flag==3 and index%2==0):
                    oligo=oligo[::-1]+sequence[sequence.index(oligo[::-1])+len(oligo)]
                    #print(3, index)
                elif (flag==3 and index%2!=0):
                    oligo=oligo[::-1]+str(Seq(sequence[sequence.index(str(Seq(oligo).complement()))-1]).complement())
                    #print(4, index)
    if (flag==5 and there_is_a_hairpin==1):
        return oligo
    elif (flag==3 and there_is_a_hairpin==1):
        return oligo[::-1]
    else:
        return oligo'''


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
    def __init__(self, seq, species):
        self.seq = seq.upper().replace(" ", "")
        self.species = species

    def optimaze(self):
        if len(self.seq) < 3:
            return("Не задана последовательность гена")
        elif self.species == 'E.coli':
            s = ' '.join([self.seq[i:i + 3] for i in range(0, len(self.seq), 3)])
            s=s.replace('CTC', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'], weights=[0.21, 0.16, 0.18, 0.45])))
            s=s.replace('CTA', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'], weights=[0.21, 0.16, 0.18, 0.45])))
            s=s.replace('TCC', ''.join(random.choices(['TCT', 'TCA'], weights=[0.5, 0.5])))
            s=s.replace('TCG', ''.join(random.choices(['TCT', 'TCA'], weights=[0.5, 0.5])))
            s=s.replace('CCC', ''.join(random.choices(['CCT', 'CCA', 'CCG'], weights=[0.29, 0.27, 0.44])))
            s=s.replace('CAC', 'CAT')
            s=s.replace('CGA', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            s=s.replace('CGG', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            s=s.replace('AGA', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            s=s.replace('AGG', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            return(s)
        elif self.species == 'Mouse':
            s = ' '.join([self.seq[i:i + 3] for i in range(0, len(self.seq), 3)])
            s=s.replace('CTC', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('CTA', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('TCC', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('TCG', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('CCC', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('CAC', 'CAT')
            s=s.replace('CGA', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('CGG', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('ATA', ''.join(random.choices(['TTT'], weights=[1.0])))
            s=s.replace('TAT', ''.join(random.choices(['TTT'], weights=[1.0])))
            return(s)
        elif self.species == 'Drosofila':
            s = ' '.join([self.seq[i:i + 3] for i in range(0, len(self.seq), 3)])
            s=s.replace('CTC', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('CTA', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('TCC', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('TCG', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('CCC', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('CAC', 'CAT')
            s=s.replace('CGA', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('CGG', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('ATA', ''.join(random.choices(['AAA'], weights=[1.0])))
            s=s.replace('TAT', ''.join(random.choices(['AAA'], weights=[1.0])))
            return(s)
        else:
            s = ' '.join([self.seq[i:i + 3] for i in range(0, len(self.seq), 3)])
            s=s.replace('CTC', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'], weights=[0.21, 0.16, 0.18, 0.45])))
            s=s.replace('CTA', ''.join(random.choices(['TTA', 'TTG', 'CTT', 'CTG'], weights=[0.21, 0.16, 0.18, 0.45])))
            s=s.replace('TCC', ''.join(random.choices(['TCT', 'TCA'], weights=[0.5, 0.5])))
            s=s.replace('TCG', ''.join(random.choices(['TCT', 'TCA'], weights=[0.5, 0.5])))
            s=s.replace('CCC', ''.join(random.choices(['CCT', 'CCA', 'CCG'], weights=[0.29, 0.27, 0.44])))
            s=s.replace('CAC', 'CAT')
            s=s.replace('CGA', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            s=s.replace('CGG', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            s=s.replace('AGA', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            s=s.replace('AGG', ''.join(random.choices(['CGT', 'CGC'], weights=[0.53, 0.47])))
            return(s)
        

class seqAnalysis():
    def __init__(self, GCcontent):
        self.GCcontent = GCcontent.replace(" ", "")

    def obrabotka(self):
        if len(self.GCcontent)==0:
            return("Не задана последовательность гена")
        else:
            # Molecular Weight = (An x 313.21) + (Tn x 304.2) + (Cn x 289.18) + (Gn x 329.21) - 61.96
            return("Mr: " + "{0:.2f}".format(self.GCcontent.count('G')*329.21 + 
            self.GCcontent.count('C')*289.18 + 
            self.GCcontent.count('A')*313.21 + 
            self.GCcontent.count('T')*304.2 - 61.96) + '\n' + 
            "GC состав: " + "{0:.2f}".format((self.GCcontent.count('G')+self.GCcontent.count('C')*100)/len(self.GCcontent)) + "%")

'''class prettyOutput():
    def __init__(self, oligopool):
        self.oligopool = oligopool.replace(" ", "")
        self.oligos = []

    def outputFile(oligos, Tms):
        with open('pretty_output.txt', 'w') as f:
            str1 = oligos[0]
            str2 = ' '*20
            str3 = ' '*(len(oligos[0])-20)
            for Tm in Tms:
                str3 = str3 + 'Tm = ' + "{:.1f}".format(Tm) + ' '*(len(oligos[1])-29)
            for i in range (1,len(oligos)):
                if i%2!=0:
                    str2 = str2 + ' '*(len(oligos[i-1])-40) + oligos[i][::-1]
                else:
                    str1 = str1 + ' '*(len(oligos[i-1])-40) + oligos[i]
            f.write(str1+'\n')
            f.write(str2+'\n')
            f.write(str3+'\n')
            f.write('\n')
            for oligo in oligos:
                if (oligos.index(oligo)+1)%2 !=0:
                    f.write(str(oligos.index(oligo)+1) + ' f: ' + oligo + '\n') 
                else:
                    f.write(str(oligos.index(oligo)+1) + ' b: ' + oligo + '\n') 
        f.close()'''

class Window(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self)
        uic.loadUi("MainWindow.ui", self)
        self.pushButton_6.clicked.connect(self.loadFile)
        self.pushButton.clicked.connect(self.codonOptimize)
        self.analiseseq_pushButton.clicked.connect(self.analiseSeq)
        self.run_pushButton.clicked.connect(self.runScript)
        self.actionAbout.triggered.connect(self.menuAbout)
        self.actionBase_codes.triggered.connect(self.menuBaseCodes)
        self.actionLoad_File.triggered.connect(self.loadFile)
        self.actionExit.triggered.connect(self.closedApp)

    def closedApp(self):
        reply = QMessageBox.question(self, 'Close App', 'You progress wil be lost! Continue?')
        if reply == QMessageBox.Yes:
            self.close()
        if reply == QMessageBox.No:
            pass
                
    def menuAbout(self):
        menuAbout(self)

    def menuBaseCodes(self):
        baseCodes(self)
        

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
       
            
    def codonOptimize(self):
        self.popupwin(codonOpt(self.plainTextEdit.toPlainText(), self.comboBox.currentText()).optimaze(),"Codon Optimisation")
        optimazeSeq = codonOpt(self.plainTextEdit.toPlainText(), self.comboBox.currentText()).optimaze()
        self.plainTextEdit.clear()
        if len(optimazeSeq) > 3 and optimazeSeq != "Не задана последовательность гена":
            self.plainTextEdit.insertPlainText(optimazeSeq)
        else:
            pass

    def analiseSeq(self):
        #print(Tm(self.plainTextEdit.toPlainText()).Tm_calculation_NN())
        self.popupwin(seqAnalysis(self.plainTextEdit.toPlainText()).obrabotka(),"Analise Sequence")
        seqAnalysis(self.plainTextEdit.toPlainText()).obrabotka()
        pass

    def runScript(self):
        if len(self.plainTextEdit.toPlainText()) > 3:
            oligos = oligoPool(self.plainTextEdit.toPlainText()).create()
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
        else:
            self.popupwin("Не задана последовательность гена","Run App")
            pass
        

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = Window()
    w.show()
    sys.exit(app.exec_())
