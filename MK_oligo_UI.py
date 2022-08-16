# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 10:34:43 2022

@author: user
"""

#Программа для разбивки гена на олигонуклеотиды для синтеза при помощи ПЦР assembly

#######################################################################################


import math, io, os, random
from Bio.Seq import Seq
from PyQt5 import QtWidgets, uic #, QtCore, QtGui
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QCheckBox

class codonOpt():
    def __init__(self, seq):
        self.seq = seq

    def optimaze(self):
        if len(self.seq) < 3:
            return("Не заденна последовательность олигонуклеотида")
        else:
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
        self.GCcontent = GCcontent

    def obrabotka(self):
        if len(self.GCcontent)==0:
            return("Не заденна последовательность олигонуклеотида")
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
       
    def markOligo(self):
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
        return oligos
        
    def codonOptimize(self):
        self.popupwin(codonOpt(self.plainTextEdit.toPlainText()).optimaze(),"Info")
        pass

    def analiseSeq(self):
        self.popupwin(seqAnalysis(self.plainTextEdit.toPlainText()).obrabotka(),"Info")
        #seqAnalysis(self.plainTextEdit.toPlainText()).obrabotka()
        pass

    def runScript(self):
        print("runScript")
        '''oligos=create(cut(sequence))

        for oligo in oligos:
            if (oligos.index(oligo)+1)%2 !=0:
                print(oligos.index(oligo)+1, 'f: ', oligo) 
            else:
                print(oligos.index(oligo)+1, 'b: ', oligo)
        Tms=Tm_calculation(oligos)
        pretty_output(oligos, Tms)

        print("Олигонуклеотиды с областями комплиментарности и температурами отжига представлены в файле pretty_output.txt")
        
        oligos=[hairpin(oligo, oligos.index(oligo), sequence, 5) for oligo in oligos]
        oligos=[hairpin(oligo[::-1], oligos.index(oligo), sequence, 3)[::-1] for oligo in oligos]
        print('')
        print('Олигонуклеотиды, исправленные с учетом шпилек:')

        for oligo in oligos:
            if (oligos.index(oligo)+1)%2 !=0:
                print(oligos.index(oligo)+1, 'f: ', oligo) 
            else:
                print(oligos.index(oligo)+1, 'b: ', oligo)'''
        pass




if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = Window()
    w.show()
    sys.exit(app.exec_())
