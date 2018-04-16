# read matrix_seed_to_all_target
import numpy as num
from pandas import HDFStore,DataFrame, read_hdf, Series
import pickle as pk
from regDiff import TVRegDiff
from scipy import interpolate
from tqdm import tqdm
def AR(A):
	# Computes Assymetry ratio
	B=A-A.transpose()
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	eu=float(sum(sum(abs(B))))/2
	if total_edges==0:
		return 1.0
	return eu/total_edges

def NAR(A):
	# Computes the normalized asymmetry ratio
	l1,l2=num.shape(A)
	if l1!=l2:
		print('Error: Input matrix not square')
		return None
	l=l1
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	AR_rand=1-(float(total_edges)/(l*(l-1)))
	if AR_rand==0:
		return float('inf')
	return AR(A)/AR_rand

def sim(A1,A2,mode='jac'):
	# Jaccard similarity on edges
	l1,l2=num.shape(A1)
	if l1!=l2:
		print('Error: Input matrix not square')
		return None
	l=l1
	l1,l2=num.shape(A1)
	if l1!=l or l2!=l:
		print('Error: Matrices must be of the same size')
		return None
	if mode=='jac':
		P=A1+A2
		P=P-num.diag(num.diag(P))
		P[P>0]=1
		S=A1*A2
		S=S-num.diag(num.diag(S))
		return float(sum(sum(S)))/sum(sum(P))
	else:
		S=A1*A2
		S=S-num.diag(num.diag(S))
		A1=A1-num.diag(num.diag(A1))
		A2=A2-num.diag(num.diag(A2))
		denom=min(sum(sum(A1)),sum(sum(A2)))
		return float(sum(sum(S)))/denom

def density(A):
	# returns the density of edges in A
	l1,l2=num.shape(A)
	if l1!=l2:
		print('Error: Input matrix not square')
		return None
	l=l1
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	return float(total_edges)/(l*(l-1))


def readS2R(fn = "strengthL1"):
    with open(fn) as f:
        for numSeeds,line in enumerate(f):
            if numSeeds == 0:
                w = [int(float(xx)) for xx in line.strip().split()]
                numROIs = len(w)
            else:
                pass
    Z = num.zeros((numSeeds+1,numROIs))
    with open(fn) as f:
        for i,line in enumerate(f):
            w = [int(float(xx)) for xx in line.strip().split()]
            Z[i,:] = w
    return Z

def readS2R_L(fn = "distanceL1"):
    with open(fn) as f:
        for numSeeds,line in enumerate(f):
            if numSeeds == 0:
                w = [float(xx) for xx in line.strip().split()]
                numROIs = len(w)
            else:
                pass
    Z = num.zeros((numSeeds+1,numROIs))
    with open(fn) as f:
        for i,line in enumerate(f):
            w = [float(xx) for xx in line.strip().split()]
            Z[i,:] = w
    return Z

def pdf(z):
    a1, a2 = num.percentile(z, [0,99])
    a = [xx for xx in z if (xx>a1 and xx<a2)]
    a = z
    hist, bin_edges = num.histogram(a, bins=20, density=True)
    return (hist, bin_edges[1:])

def randomDisPDF():
    # Get the PDF for 9 random connection
    D = num.zeros((9,20)) # Distances
    C = [0.0]*9 # connections
    Z = [(0.0,0.0,0.0)]*9 # connections
    X = num.zeros((9,20))
    A = readS2R()
    L = readS2R_L() # Read all distances
    i = 0
    S = set([])
    while True:
        if i==9:
            break
        cur = num.random.randint(1,179)
        if i == 0:
            cur = 3
        if cur in S:
            continue
        else:
            S.add(cur)
        NOS = num.max(A[:,cur])
        can = int(num.argmax(A[:,cur]))
        if NOS < 100:
            continue
        C[i] = (cur+1,NOS)
        L1 = L[:,cur]; # Get connected ROI distances (L1 to L?) - python indecies start from 0
        L2 = L1.flatten()
        L3 = L2[num.nonzero(L2)] # Get distances from all connected seeds
        mn = num.mean(L3)
        md = num.median(L3)
        ma = L[can,cur]
        if len(L3) < 5:
            continue
        y, x = pdf(L3)
        X[i,:] = x
        D[i,:] = y
        Z[i] = (mn,md,ma)
        i = i + 1
    return (D,X,C,Z)

def distanceBias():
    D = readS2R_L() # distances
    P = readS2R() # Probabilities
    ns, nr = D.shape
    s2r = {}
    for r in range(nr):
        if r == 0:
            continue
        for s in range(ns):
            if P[s,r]>2000:
                try:
                    s2r[r].append((P[s,r],D[s,r]))
                except KeyError:
                    s2r[r] = [(P[s,r],D[s,r])]
    return s2r

def sample_csv():
    D1 = readS2R_L("distanceL1") # distances
    D4 = readS2R_L("distanceL4") # distances
    P1 = readS2R("strengthL1") # Probabilities
    P4 = readS2R("strengthL4") # Probabilities
    L1ToL4_d = D1[:,3];
    L4ToL1_d = D4[:,0];
    L1ToL4_p = P1[:,3];
    L4ToL1_p = P4[:,0];
    a = num.zeros((len(L1ToL4_d),2))
    a[:,0] = L1ToL4_d
    a[:,1] = L1ToL4_p
    b = num.zeros((len(L4ToL1_d),2))
    b[:,0] = L4ToL1_d
    b[:,1] = L4ToL1_p
    num.savetxt("L1.csv", a, delimiter=",")
    num.savetxt("L4.csv", b, delimiter=",")

def create_store(sub):
    hdf =HDFStore('all.h5')
    d = DataFrame(columns=['SUB','SEED','SEED ROI','TARGET ROI','HEMISPHERE','DISTANCE','STRENGTH','CAT1','CAT2','CAT3'])
    for i in range(1,181):
        LSfname = '../'+sub+'/out/L'+str(i)+'/matrix_seeds_to_all_targets'
        LDfname = '../'+sub+'/out/L'+str(i)+'/matrix_seeds_to_all_targets_lengths'
        RSfname = '../'+sub+'/out/R'+str(i)+'/matrix_seeds_to_all_targets'
        RDfname = '../'+sub+'/out/R'+str(i)+'/matrix_seeds_to_all_targets_lengths'
        ls = readS2R(LSfname)
        rs = readS2R(RSfname)
        ld = readS2R_L(LDfname)
        rd = readS2R_L(RDfname)
        numSeeds ,numROIs = ls.shape
        for j in tqdm(range(numSeeds),total=numSeeds):
            for q in range(numROIs):
                tmp = Series([sub,j+1,i+1,q+1,'L',ld[j,q],ls[j,q],'','',''])
                d = d.append(tmp,ignore_index=True)
        # numSeeds ,numROIs = rs.shape
        # for j in range(numSeeds):
        #     for q in range(numROIs):
        #         tmp = Series([sub,j+1,i+1,q+1,'R',rd[j,q],rs[j,q],'','',''])
        #         d = d.append(tmp,ignore_index=True)
        if i == 1: break
    hdf.put(sub,d)



    # a1 = readS2R()
    # a2 = readS2R_L()
    # hdf.put('sub/L1/p', a2, format='table', data_columns=True)
    # hdf.put('sub/L1/c', a1, format='table', data_columns=True)

def basic_results(sub):
    with open('../'+sub+'.res','rb') as f:
        D = pk.load(f)
    den = [xx[2] for xx in D]
    T = [xx[3] for xx in D]
    A = [xx[1] for xx in D]
    net = [xx[0] for xx in D]
    return (den[50:-50],A[50:-50],T[50:-50],net[50:-50])

def wideFlat(sub):
    with open('../'+sub+'.res','rb') as f:
        D = pk.load(f)
    NET = [xx[0] for xx in D]



def smoothers(sub):
    d,a,t = basic_results(sub)
    d = list(reversed(d))
    t = list(reversed(t))
    a = num.array(list(reversed(a)))
    deriv_sm = TVRegDiff(100*a, 100, 10e10, dx=1, ep=1e-1, scale='small', plotflag=0,diagflag=False)
    deriv_sm = abs(deriv_sm)
    return (deriv_sm,d,a,t)
