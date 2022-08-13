# -*- coding: utf-8 -*- #
"""
Created on Sun Aug  7 13:17:51 2022

@author: user
"""

#MK_oligo.py реализованная в формате телеграм бота

import re
import telebot
import math
from Bio.Seq import Seq

bot = telebot.TeleBot('5510594397:AAFLNzH-F9QRAJYyExc2Af2QgcfOCPedzxU')

@bot.message_handler(commands=['start'])
def start(message):
    bot.send_message(message.chat.id, 'Привет! В ответ на это сообщение пришлите мне ген для обработки в формате FASTA')
    
@bot.message_handler()
def get_user_text(message):
    if re.search(r'[^AaGgTtCc\n]',message.text ):
        bot.send_message(message.chat.id, 'Это не строка в формате FASTA')
    else:
        bot.send_message(message.chat.id, 'Отлично, провожу расчеты...')
        sequence=str(message.text).replace('\n', '').upper()
        oligos=create(cut(sequence))
        Tms=Tm_calculation(oligos)
        with open('pretty_output.txt', 'w') as f:
            str1=oligos[0]
            str2=' '*20
            str3=' '*(len(oligos[0])-20)
            for Tm in Tms:
                str3=str3+'Tm = '+"{:.1f}".format(Tm)+' '*(len(oligos[1])-29)
            for i in range (1,len(oligos)):
                if i%2!=0:
                    str2=str2+' '*(len(oligos[i-1])-40)+oligos[i][::-1]
                else:
                    str1=str1+' '*(len(oligos[i-1])-40)+oligos[i]
            f.write(str1+'\n')
            f.write(str2+'\n')
            f.write(str3+'\n')
            f.write('\n')
            for oligo in oligos:
                if (oligos.index(oligo)+1)%2 !=0:
                    f.write(str(oligos.index(oligo)+1)+' f: '+oligo+'\n') 
                else:
                    f.write(str(oligos.index(oligo)+1)+' b: '+oligo+'\n') 
        f.close()
        file = open('pretty_output.txt', 'r')
        bot.send_message(message.chat.id, 'В этом файле представлены олигонуклеотиды для синтеза вашего гена. Под областями комплементарности указаны температуры отжига:')
        bot.send_document(message.chat.id, file)
        bot.send_message(message.chat.id, "Избавляюсь от тугоплавких шпилек на 5' и 3' концах...")
        oligos=[hairpin(oligo, oligos.index(oligo), sequence, 5) for oligo in oligos]
        oligos=[hairpin(oligo[::-1], oligos.index(oligo), sequence, 3)[::-1] for oligo in oligos]
        with open('no_hairpins.txt', 'w') as f:
            for oligo in oligos:
                if (oligos.index(oligo)+1)%2 !=0:
                    f.write(str(oligos.index(oligo)+1)+' f: '+oligo+'\n') 
                else:
                    f.write(str(oligos.index(oligo)+1)+' b: '+oligo+'\n') 
        f.close()
        file = open('no_hairpins.txt', 'r')
        bot.send_message(message.chat.id, 'Шпильки скорректированы!')
        bot.send_document(message.chat.id, file)
        
        
###############################################################################


#Функция для разбивки гена на куски одной цепи
def cut(s):
    oligos=[]
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
    return oligos

#Функция для создания олигонуклеотидов
def create(oligos):
    for i in range (2, len(oligos), 2):
        oligos[i]=oligos[i-1][len(oligos[i-1])-20:len(oligos[i-1])]+oligos[i]
    for i in range (1, len(oligos), 2):
        oligos[i]=str(Seq(oligos[i-1][len(oligos[i-1])-20:len(oligos[i-1])]).complement()+Seq(oligos[i]).complement())[::-1]
    return oligos
            
#Функция для расчета температуры отжига олигов - !!! желательно переделать на nearest neighbour !!!
def Tm_calculation(oligos):
    Tms=[]
    for i in range (0, len(oligos)-1):
        if i%2==0:
            wA=oligos[i][len(oligos[i])-20:len(oligos[i])].count('A')
            xT=oligos[i][len(oligos[i])-20:len(oligos[i])].count('T')
            yG=oligos[i][len(oligos[i])-20:len(oligos[i])].count('G')
            zC=oligos[i][len(oligos[i])-20:len(oligos[i])].count('C')
        else:
            oligos[i]=oligos[i][::-1]
            wA=oligos[i][len(oligos[i])-20:len(oligos[i])].count('A')
            xT=oligos[i][len(oligos[i])-20:len(oligos[i])].count('T')
            yG=oligos[i][len(oligos[i])-20:len(oligos[i])].count('G')
            zC=oligos[i][len(oligos[i])-20:len(oligos[i])].count('C')
            oligos[i]=oligos[i][::-1]
        Tm=100.5+(41*(yG+zC)/(wA+xT+yG+zC))-(820/(wA+xT+yG+zC)) + 16.6*math.log10(0.05)
        Tms.append(Tm)
        #print('Tm = ', "{:.1f}".format(Tm))
    return Tms

#Функция для поиска шпилек 
#flag =5 - поиск шпилек с 5' конца
#flag =3 - поиск шпилек с 3' конца
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
                #print('Обнаружена шпилька в олигонуклеотиде номер ', index+1, ', температура плавления: ' "{:.1f}".format(hairpin_Tm(stem, loop_size)), ',проводится корректировка')
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
        return oligo
    

# Из GeneRunner: Tm is the thermodynamic melting temperature for the structure displayed in the window. 
#For hairpins it is the temperature where the ΔG value is 0, 
#the temperature above which the hairpin is no longer stable. 
#These values only indicate the relative stability of the secondary structures and are not exact Tms. 
#They should only be used to compare the relative stability of the structures.

#dG is the free energy of the dimer or hairpin. 
#The ΔG for a hairpin is determined by the free energy values of Freier at al. 
#delta G loop incriment - Freier, Susan M. et al. “Improved free-energy parameters for predictions of RNA duplex stability.” Proceedings of the National Academy of Sciences of the United States of America 83 24 (1986): 9373-7 .
#The values are listed in the appendix A. 
#NN delta G и delta S - Breslauer KJ, Frank R, Blöcker H, Marky LA. Predicting DNA duplex stability from the base sequence. Proc Natl Acad Sci U S A. 1986;83(11):3746-3750. doi:10.1073/pnas.83.11.3746  Table 2

#Функция посчета delta G loop incriment (увеличение дельты Жэ в зависимости от длины петли шпильки см. ссылки перед hairpin_Tm(stem_sequence, loop_size))                
def loop_increment(loop_size):
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



#################################################################################
bot.polling(none_stop=True)