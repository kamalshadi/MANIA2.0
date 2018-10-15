import json
import numpy as num
from py2neo import Node, Relationship, Graph
from curve import *
from utils import *

with open('config.json') as f:
    config = json.load(f)

def getData(sub,hem):
    # get subject raw data for a hemisphere
    query = '''match (n:ROI)-[r]->(m:ROI)
    where n.hemisphere = "ZZZ" and r.SUBJECT=SSS
    return n.name as n1,m.name as n2,r._length as _length,r._weight as _weight, r.weight as weight, r.length as length'''
    query = query.replace('ZZZ',hem)
    query = query.replace('SSS',str(sub))
    graph = Graph(password='1234')
    A = graph.run(query).data()
    return A

def getRawMatrix(sub,hem):
    # get subject raw weighted matrix for a hemisphere
    query = '''match (n:ROI)-[r]->(m:ROI)
    where n.hemisphere = "ZZZ" and r.SUBJECT=SSS
    return n.name as n1,m.name as n2, r.weight as weight, r.length as length'''
    query = query.replace('ZZZ',hem)
    query = query.replace('SSS',str(sub))
    graph = Graph(password='1234')
    A = graph.run(query)
    W = num.zeros((180,180))
    F = num.zeros((180,180))
    for w in A:
        n1 = w['n1']
        n2 = w['n2']
        if(n1==n2):continue
        i1 = int(n1[1:]) - 1
        i2 = int(n2[1:]) - 1
        W[i1][i2] = w['weight']
        F[i1][i2] = w['length']
    return (W,F)

def getAllRawMatrix():
    Z = {}
    for sub in config['SUBS']:
        for hem in ['L','R']:
            print(sub,hem)
            W, F = getRawMatrix(sub,hem)
            Z[(sub,hem)] = (W,F)
    return Z




def getLocalRegressors(A):
    # get brain regressor
    R = []
    for w in A:
        if(w['n1'] == w['n2'] ):continue
        t1 = sorted(zip(w['_length'],w['_weight']))
        a = [xx[0] for xx in t1]
        b = [num.log(max(xx[1],1.0)/config['N_STREAMLINES']) for xx in t1]
        D = regress_noplot(a,b)
        D['a'] = a
        D['b'] = b
        D['n1'] = w['n1']
        D['n2'] = w['n2']
        R.append(D)
    return R

def correcttion(R,A):
    # get corrected data from regressor
    i1 = 0
    i2 = 0
    i3 = 0
    A_dic = {}
    C = {}
    NC = {}
    for i,w in enumerate(A):
        A_dic[w['n1']+w['n2']] = i
    for i,z in enumerate(A):
        try:
            C[z['n2']+z['n1']]
            continue
        except KeyError:
            j = A_dic[z['n2']+z['n1']]
            w = A[i]
            v = A[j]
            e1 = w['envelope']
            e2 = v['envelope']
            wbi = num.argmax(w['b'])
            vbi = num.argmax(v['b'])
            wb = num.max(w['b'])
            vb = num.max(v['b'])
            if ( len(e1)>0 and len(e2)>0 ):
                i1 = i1 + 1
                # both direction have envelope points
                R1 = R['s-'+w['n1']]
                R2 = R['t-'+w['n2']]
                m = []
                for p in e1:
                    m.append(p[1]-R1[1]*p[0])
                    m.append(p[1]-R2[1]*p[0])
                C[w['n1']+w['n2']] = num.median(m)
                NC[w['n1']+w['n2']] = max(w['b'])
                R1 = R['s-'+w['n2']]
                R2 = R['t-'+w['n1']]
                m = []
                for p in e2:
                    m.append(p[1]-R1[1]*p[0])
                    m.append(p[1]-R2[1]*p[0])
                C[w['n2']+w['n1']] = num.median(m)
                NC[w['n2']+w['n1']] = max(v['b'])
            elif (wb>noise_floor and vb>noise_floor):
                i2 = i2 + 1
                # both direction have points above noise floor
                R1 = R['s-'+w['n1']]
                R2 = R['t-'+w['n2']]
                m = []
                p = [w['a'][wbi],wb]
                m.append(p[1]-R1[1]*p[0])
                m.append(p[1]-R2[1]*p[0])
                C[w['n1']+w['n2']] = num.median(m)
                NC[w['n1']+w['n2']] = max(w['b'])
                R1 = R['s-'+w['n2']]
                R2 = R['t-'+w['n1']]
                m = []
                p = [v['a'][vbi],vb]
                C[w['n2']+w['n1']] = num.median(m)
                NC[w['n2']+w['n1']] = max(v['b'])
            else:
                i3 = i3 + 1
                NC[w['n1']+w['n2']] = max(w['b'])
                NC[w['n2']+w['n1']] = max(v['b'])
                C[w['n1']+w['n2']] = max(w['b'])
                C[w['n2']+w['n1']] = max(v['b'])
    return (C,NC)

def getNetByDensity(M,d):
    for i in range(-1,config['N_STREAMLINES']):
        net = num.zeros((config['N_ROI'],config['N_ROI']))
        net[M>i] = 1
        if density(net)<d:
            return net

def getMatrix(C,hem):
    nroi = config['N_ROI']
    z = num.zeros((nroi,nroi))
    for i in range(nroi):
        for j in range(nroi):
            if(i==j):continue
            r1 = hem+str(i+1)
            r2 = hem+str(j+1)
            try:
                if(not num.isnan(C[r1+r2])):
                    z[i][j]= min(num.exp(C[r1+r2]),1)*config['N_STREAMLINES']
            except KeyError:
                pass
    return z

def getMANIA(M):
    for i in range(config['CLIP'],config['N_STREAMLINES']-config['CLIP']):
        netc = num.zeros((config['N_ROI'],config['N_ROI']))
        netc[M>i] = 1
        q = i - config['CLIP']
        _minc = NAR(netc)
        if(q==0):
            out_net = netc
            minc = _minc
        else:
            if(_minc<=minc):
                out_net = netc
                minc = _minc
    return out_net

def randM(d):
    M = num.zeros((180,180))
    q = 0
    m = num.floor(d*180)
    while True:
        i = num.random.randint(0,180)
        j = num.random.randint(0,180)
        if(i==j):continue
        if(M[i][j]==1):continue
        M[i][j]=1
        q = q + 1
        if(q==m):break
    return M

def augment(net,m):
    q = 0
    out = set([])
    while True:
        i = num.random.randint(0,180)
        j = num.random.randint(0,180)
        if(i==j):continue
        if(net[i][j]==1):continue
        net[i][j]=1
        q = q + 1
        out.add((i,j))
        if(q==m):break
    return out

def randomizeTest(R,L,dR,dL):
    nR = dR*180*180 - density(R)*180*180
    nL = dL*180*180 - density(L)*180*180
    s1 = augment(R,nR)
    s2 = augment(L,nL)
    return (s1,s2)

def getNetwroks():
    Z = {}
    E = []
    for sub in config['SUBS']:
        for hem in ['L','R']:
            print(sub,hem)
            try:
                A = getData(sub,hem)
                R = getLocalRegressors(A)
                RL = roi_regressors(R,hem)
                C,NC = correcttion(RL,R)
                M = getMatrix(NC,hem)
                MN1 = getMANIA(M)
                M = getMatrix(C,hem)
                MN2 = getMANIA(M)
                Z[(sub,hem)] = (MN1,MN2)
            except:
                print('<><><><ERR><><><><>')
                print(sub,hem)
                print('<><><><><><><><><><>')
                E.append((sub,hem))
    return (Z,E)
