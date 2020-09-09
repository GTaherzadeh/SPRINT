###################################
#!/usr/bin/env python
#coding=utf-8
# SPRINT.py
# Ghazaleh Taherzadeh
###################################
"""predict Peptide-binding residues for a new chain"""
import sys
import string
import math
import os
import numpy
import scipy.io

# function
def encode_restype(res):
    """binary encoding residue type."""
    AAs = 'ARNDCQEGHILKMFPSTWYV'
    code = []
    for a in AAs:
        if res == a:
            code.append('1')
        else:
            code.append('0')

              
    return code
    
def compute_entropy(dis_list):
    """compute shannon entropy for a distribution.
    base = len(dis_list) is the base of log function 
    to make entropy between 0 and 1."""
    
    if sum(dis_list) == 0:
        return 0.0
    prob_list = map(lambda x:(x+0.0)/sum(dis_list),dis_list)
    ent = 0.0
    for prob in prob_list:
        if prob != 0:
            ent -= prob*math.log(prob,len(dis_list))
    return ent

def compute_ss_content(ss_seq_win):
    """compute ss content in a window."""
    con_C = con_H = con_E = 0
    for ss in ss_seq_win:
        if ss == 'C':
            con_C += 1
        elif ss == 'H':
            con_H += 1
        elif ss == 'E':
            con_E += 1
        else:
            print('X')
             
    act_len = con_C+con_H+con_E+0.0
    return ['%.3f'%(con_C/act_len),'%.3f'%(con_H/act_len),'%.3f'%(con_E/act_len)]

def compute_pcc(x,y):
    """compute the PCC between vector x and y"""
    mean_x = (sum(x)+0.0)/len(x)
    dev_x = sum([(i-mean_x)**2 for i in x])
    mean_y = (sum(y)+0.0)/len(y)
    dev_y = sum([(i-mean_y)**2 for i in y])
    if dev_x == 0 or dev_y == 0:
        return 0.0
    ret = 0.0
    for i in xrange(len(x)):
        ret += (x[i]-mean_x)*(y[i]-mean_y)
    return ret/math.sqrt(dev_x*dev_y)

def ss_to_num(sin_ss):
    """C->0,H->1,E->2,'$'->-1"""
    if sin_ss == 'C':
        return 0
    elif sin_ss == 'H':
        return 1
    elif sin_ss == 'E':
        return 2
    else:
        return -1

def tri_ss_to_num(tri_ss):
    """CCC->0,EEE->26,$XX->-1,XX$->-1
    return a 27-dimensional vector"""
    tri_num = []
    for ts in tri_ss:
        tri_num.append(ss_to_num(ts))
    num = -1
    if tri_num[0] != -1 and tri_num[2] != -1:
        num = 9*tri_num[0]+3*tri_num[1]+tri_num[2]
    ret = []
    for i in xrange(27):
        if i == num:
            ret.append('1')
        else:
            ret.append('0')
    return ret

def seg_bound(s_win):
    """Two boundaries of a segment"""
    c_ss = s_win[len(s_win)/2]
    l_len = r_len = 0
    i = len(s_win)/2 - 1
    while i >= 0:
        if s_win[i] != c_ss:
            break
        l_len += 1
        i -= 1
    i = len(s_win)/2 + 1
    while i < len(s_win):
        if s_win[i] != c_ss:
            break
        r_len += 1
        i += 1
    return (l_len,r_len)

def com_mea(tp,fp,tn,fn):
    """compute sen,spe,acc,mcc,pre."""
    sen = (tp+0.0)/(tp+fn)
    spe = (tn+0.0)/(tn+fp)
    acc = (tp+tn+0.0)/(tp+fp+tn+fn)
    mcc = 0
    if (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) != 0:
        mcc = (tp*tn-fp*fn+0.0)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    pre = 0
    if tp+fp != 0:
        pre = (tp+0.0)/(tp+fp)
    return (sen,spe,acc,mcc,pre)

def com_acc(pre,act):
    """compute the accuracy.
    pre and act are list of integers(1:positive;-1:negetive)."""
    tp = fp = tn = fn = 0
    for i in xrange(len(pre)):
        if pre[i] == 1:
            if act[i] == 1:
                tp += 1
            else:
                fp += 1
        else:
            if act[i] == 1:
                fn += 1
            else:
                tn += 1
    (sen,spe,acc,mcc,pre) = com_mea(tp,fp,tn,fn)
    return (tp,fp,tn,fn,sen,spe,acc,mcc,pre)


ext_state = {'A':110.2,'D':144.1,'C':140.4,'E':174.7,'F':200.7,\
             'G':78.7,'H':181.9,'I':185.0,'K':205.7,'L':183.1,\
             'M':200.1,'N':146.4,'P':141.9,'Q':178.6,'R':229.0,\
             'S':117.2,'T':138.7,'V':153.7,'W':240.5,'Y':213.7}

#pid = r'1a81A' ./SPRINT.py 1a81A
pid = sys.argv[1] 
G1_win = 1
G2_win = 9
G3_win = 2
G4_win = 2 
base_path = 'path_to_your_sequence_dir'
res_path = base_path+'/Results/'
ss_path = base_path+'/SPIDER_1/'
rsa_path = base_path+'/SPIDER_1/'
pssm_path = base_path+'/pssm_pssm/'
fasta_path = base_path+'/fasta/'
info_path = base_path+'info/'
fea_path = base_path+'fea/'
error_file = base_path+'error.log'
fname = file(base_path+'fn.txt','w')
fname.write('%s\n' %pid)
fname.close()
#/*************************************************/
# build info file
fin = file(ss_path+pid+'.fasta.spd3','r')
ss = fin.readlines()[1:] 
fin.close()
ss_res = map(lambda x:x.split()[1],ss)
ss_pre = map(lambda x:x.split()[2],ss)
ss_prob_c = map(lambda x:x.split()[3],ss)
ss_prob_h = map(lambda x:x.split()[5],ss)
ss_prob_e = map(lambda x:x.split()[4],ss)
fin = file(rsa_path+pid+'.fasta.spd3','r')
sp = fin.readlines()[1:]
fin.close()
rsa_res = ''.join([x.split()[1] for x in sp])
rsa_pre = [string.atof(x.split()[6]) for x in sp]
fin = file(pssm_path+pid+'.pssm','r')
pssm = fin.readlines()
fin.close()
if len(pssm[-6].split()) != 0 or pssm[3].split()[0] != '1': 
    print 'error on reading pssm, line -6 is not a spare line;\
     or line 3 is not the first line'
    sys.exit(1)
pssm = pssm[3:-6]
fin = file(fasta_path+pid+'.fasta','r')
ann = fin.readlines()
fin.close()
if len(ann) != 2:
    print 'check sequence',pid
    sys.exit(1)
fastaseq = ann[1].split()[0]
if not fastaseq == rsa_res == ''.join(ss_res):
    print 'Sequence inconsistent!'
    print 'fasta: ',fastaseq
    print '   ss: ',''.join(ss_res)
    print '  rsa: ',rsa_res
    exit(1)
fout = file(info_path+pid+'.info','w')
fout.write('>%s\n' %pid)
pos = 0
for i in xrange(len(fastaseq)):
    res = fastaseq[i]
    fout.write('%5d%5s%5s'%(i+1,res,res))
    if pssm[pos].split()[1] == res:
        for p_e in pssm[pos].split()[2:22]:
            fout.write(':%2s' %p_e)
        for p_e in pssm[pos].split()[22:42]:
            fout.write(':%3s' %p_e)
        fout.write(':%5s' %pssm[pos].split()[42])
    else:
        print 'Error reading pssm file!'
        flog = file(error_file,'a')
        flog.write(pid+': error on writing pssm, %s:%s\n' \
        %(pssm[pos].split()[1],res))
        flog.close()
        sys.exit(1)
    if ss_res[pos] == res:
        fout.write(':%s' %ss_pre[pos])
        fout.write(':%s' %ss_prob_c[pos])
        fout.write(':%s' %ss_prob_h[pos])
        fout.write(':%s' %ss_prob_e[pos])
    else:
        print 'Error reading ss file!'
        flog = file(error_file,'a')
        flog.write(pid+': error on writing ss, %s:%s\n' %(ss_res[pos],res))
        flog.close()
        sys.exit(1)
    if rsa_res[pos] == res:
        fout.write(':%5.3f' %rsa_pre[pos])
    else:
        print 'Error reading rsa file!'
        flog = file(error_file,'a')
        flog.write(pid+': error on writing rsa, %s:%s\n' %(rsa_res[pos],res))
        flog.close()
        sys.exit(1)
    pos += 1
    fout.write('\n')
fout.close()
#/*************************************************/
# build feature file
fin = file(info_path+'%s.info'%pid,'r')
info = fin.readlines()[1:]
fin.close()
output = file(fea_path+'%s.fea'%pid,'w')
seq_len = len(info)
out_list = []
for i in xrange(len(info)):
    out_list.append([])
rt = map(lambda x:x.split(':')[0][-1],info)

for i in xrange(G1_win):
    rt.insert(0,'X')
    rt.append('X')
for i in xrange(G1_win,seq_len+G1_win):
    for j in xrange(i-G1_win,i+G1_win+1):
        out_list[i-G1_win].append(':'.join(encode_restype(rt[j])))
pssm = map(lambda x:map(lambda y:'%7.5f' %(1/(1+math.pow(math.e,-string.atoi(y)))),x.split(':')[1:21]),info)
pssm_t = []
for i in xrange(20):
    pssm_t.append('%7.5f' %(1/(1+math.e**0)))
for i in xrange(G2_win):
    pssm.insert(0,pssm_t)
    pssm.append(pssm_t)
for i in xrange(G2_win,seq_len+G2_win):
    for j in xrange(i-G2_win,i+G2_win+1):
        out_list[i-G2_win].append(':'.join(pssm[j]))
wop = ['%7.5f' %compute_entropy(z) for z in map(lambda x:map(lambda y:string.atoi(y),x.split(':')[21:41]),info)]
for i in xrange(G2_win):
    wop.insert(0,'%7.5f' %(0))
    wop.append('%7.5f' %(0))
for i in xrange(G2_win,seq_len+G2_win):
    for j in xrange(i-G2_win,i+G2_win+1):
        out_list[i-G2_win].append(wop[j])
pssm = [[string.atoi(y) for y in x.split(':')[1:21]] for x in info]
pssm_t = []
for i in xrange(20):
    pssm_t.append(0)
for i in xrange(G2_win):
    pssm.insert(0,pssm_t)
    pssm.append(pssm_t)
for i in xrange(G2_win,seq_len+G2_win):
    for j in xrange(i-G2_win,i+G2_win+1):
        if j != i:
            out_list[i-G2_win].append('%.4f'%(compute_pcc(pssm[i],pssm[j])))
ss_prob = map(lambda x:x.split(':')[43:46],info)
for i in xrange(G3_win):
    ss_prob.insert(0,['0.000','0.000','0.000'])
    ss_prob.append(['0.000','0.000','0.000'])
for i in xrange(G3_win,seq_len+G3_win):
    for j in xrange(i-G3_win,i+G3_win+1):
        out_list[i-G3_win].append(':'.join(ss_prob[j]))
ss_seq = map(lambda x:x.split(':')[42],info)
for i in xrange(G3_win):
    ss_seq.insert(0,'$')
    ss_seq.append('$')
for i in xrange(G3_win,seq_len+G3_win):
    out_list[i-G3_win].append(':'.join(compute_ss_content(ss_seq[i-G3_win:i+G3_win+1])))
for i in xrange(G3_win,seq_len+G3_win):
    out_list[i-G3_win].append(':'.join(tri_ss_to_num(ss_seq[i-1:i+2])))
for i in xrange(G3_win,seq_len+G3_win):
    [l_b,r_b] = seg_bound(ss_seq[i-G3_win:i+G3_win+1])
    for j in xrange(3):#0-C,1-H,2-E
        if j == ss_to_num(ss_seq[i]):
            out_list[i-G3_win].append('%.3f'%((l_b+r_b+1.0)/(2*G3_win+1)))
        else:
            out_list[i-G3_win].append('0.000')
    for j in xrange(3):
        if j == ss_to_num(ss_seq[i]):
            out_list[i-G3_win].append('%.3f'%((min(l_b,r_b)+0.0)/G3_win))
        else:
            out_list[i-G3_win].append('0.000')
    for j in xrange(3):
        if j == ss_to_num(ss_seq[i]):
            out_list[i-G3_win].append('%.3f'%((max(l_b,r_b)+0.0)/G3_win))
        else:
            out_list[i-G3_win].append('0.000')
rsa = map(lambda x:x.split(':')[46].split()[0],info)
for i in xrange(G4_win):
    rsa.insert(0,'1.000')
    rsa.append('1.000')
for i in xrange(G4_win,seq_len+G4_win):
    for j in xrange(i-G4_win,i+G4_win+1):
        out_list[i-G4_win].append(rsa[j])
rsa = [string.atof(x.split(':')[46]) for x in info]
for i in xrange(G4_win):
    rsa.insert(0,1.0)
    rsa.append(1.0)
for i in xrange(G4_win,seq_len+G4_win):
    for j in xrange(1,G4_win+1):
        out_list[i-G4_win].append('%.4f'%(sum(rsa[i-j:i+j+1])/(2*j+1)))


myList = []
for i in xrange(len(out_list)):
    output.write('-1')
    output.write('\t')
    out_list1 = []
    out_list1.append(':'.join(out_list[i]))
    myList = [t.split(':') for t in out_list1]
    for k in myList:
        for p in enumerate(k):
            D = (p)
	    output.write(':'.join([str(D[0]+1),D[1]]))
	    output.write('\t')
    output.write('\n')
output.close()
#/*************************************************/
#Classifier
outfile = file(res_path+'%s.prob'%pid,'w')
sys.path.append ('path_to_libsvm/libsvm.20/libsvm-3.20/python')
import svm
from svm import *
from svmutil import *
#y, x = svm_read_problem('./train_feature_data.txt')
#m = svm_train(y, x, '-s 0 -t 2 -g 0.05 -c 1 -b 1')
#svm_save_model('Fullmodel.model', m)
m = svm_load_model(base_path+'Fullmodel.model')
a, b = svm_read_problem(fea_path+'%s.fea'%pid)
p_label, p_acc, p_val = svm_predict(a, b, m,'-b 1')
for item in p_val:
  outfile.write("%s\n" % item[0])
outfile.close()
#/*************************************************/
