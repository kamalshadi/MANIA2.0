# Curve fitting for distance bias
import numpy as num
from scipy.optimize import curve_fit
from scipy.stats import norm
from sklearn.linear_model import TheilSenRegressor
# N is number of bins
# m is number of selected points per bin
import pickle as pk
from bokeh.plotting import output_notebook,figure, show
from bokeh.layouts import row,column,gridplot
from bokeh.models import Label

estimator = TheilSenRegressor(random_state=42)

# Base Model > Log(NOS) = linear(distance)
def func(x, a, b):
    return a*x + b


def regress(a,b, tit = ''):
    c = sorted(zip(a,b))
    d = [xx for xx in c if (xx[1]>100 and xx[0]>1)]
    a = [xx[0] for xx in d]
    b = [num.log(xx[1]) for xx in d]
    D = maxBin4(a,b,10)
    # D1 = maxBin3(a,b,10)
    outs = D['outs']
    bins = D['bins']
    _intercept = D['popt'][1] - num.log(5000)
    intercept = "{:5.2f}".format(_intercept)
    slope = num.ceil(-100*D['popt'][0])
    params = 'intercept ='+str(intercept)+' , '+'slope = '+str(slope)+'%'
    # normalize function
    fn = lambda x : x - num.log(5000)
    # outs1 = D1['outs']
    # bins1 = D1['bins']
    l = len(outs)
    # l1 = len(outs1)
    title = tit+'(R2%='+str(D['r'])+','+params+')'
    p1 = figure(plot_width=500, plot_height=500, title = title,y_axis_label = "Probability(log scale)",x_axis_label = "Distance")
    p1.circle(a,list(map(fn,b)))
    p1.line(D['x'],list(map(fn,D['z'])),legend="Max",color="red")
#     p1.line(D1['x'],D1['z'],legend="Soft-Max",color="blue")
#     p1.line(D1['x'],D1['z_pred'],legend="thiel-soft-Max",color="black")
#     p1.line(D['x'],D['z_pred'],legend="thiel-Max",color="cyan")
# #     p1.line(D['cxw'],D['czw'],legend="Donahue - weighted regression",color="blue")
    for i in range(l):
            p1.circle_cross(outs[i][0][0],fn(outs[i][0][1]),color='green',size=8)
    # for i in range(l1):
    #     p1.square(outs1[i][0][0],outs1[i][0][1],color='yellow',size=4,alpha=.8)
    return p1

max_nos = 5000    # probtrackx run parameter
noise_floor = num.log(15.0/max_nos)
min_r2 = 75
min_envelope_points = 5
min_points_on_right = 5


# sorted by x
# Local regressor using all points
# Both envelope detection and regressor
def pre_mania2(x,y):
    D = {} # regression model
    D['isSuccess'] = False
    D['reason'] = None
    D['popt'] = [0.0,0.0] # regression parameters
    D['slope'] = 0.0
    D['intercept'] = 0
    D['envelope'] = [] # envelope points
    D['r2'] = 0.0
    D['x'] = []
    D['z'] = []
    D['std'] = 0.0

    outs = [] # output points
    stopping = False
    l = len(x)
    # w = int(num.ceil(l*wp))
    tail = False
    for i,cur in enumerate(x):
        if i >= (l - min_points_on_right):
            break
        if y[i] <= noise_floor:
            continue
        ql = y[(i+1):]
        ax = num.max(ql)
        if y[i] >= ax:
            outs.append((cur,y[i]))

    D['envelope'] = outs

    if (len(outs) < min_envelope_points ):
        D['reason'] = 'envelope'
        return D

    at = [xx[0] for xx in outs]
    D['x'] = at
    bt = [xx[1] for xx in outs]

    popt, pcov = curve_fit(func, at, bt,maxfev = 10000)
    D['popt'] = popt
    D['std'] = num.sqrt(pcov[0,0])
    # a = sorted(list(set(a)))
    z = [func(xx,*popt) for xx in at]
    D['z'] = z
    residuals = [(bt[q] - z[q])**2 for q in range(len(z))]
    ss_res = num.sum(residuals)
    ss_tot = num.sum((num.array(bt)-num.mean(bt))**2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared = r_squared * 100
    D['r2'] = r_squared
    D['intercept'] = popt[1]
    D['slope'] = popt[0]
    # estimator.fit(num.array(at).reshape(-1, 1), num.array(bt).reshape(-1, 1))
    # y0 = estimator.predict(0)[0]
    # y1 = estimator.predict(1)[0]
    # D['intercept'] = y0
    # D['slope'] = y1-y0
    if (r_squared <= min_r2):
        D['reason'] = 'r2'
        return D
    D['isSuccess'] = True
    return D

# without envelope detection
# x need not to be sorted
def lslinear(x,y):
    D = {}
    w = sorted(zip(x,y))
    x = [xx[0] for xx in w]
    y = [xx[1] for xx in w]
    D['x'] = x
    D['y'] = y
    popt, pcov = curve_fit(func, x, y,maxfev = 10000)
    D['popt'] = popt
    D['std'] = num.sqrt(pcov[0,0])
    # a = sorted(list(set(a)))
    z = [func(xx,*popt) for xx in x]
    D['z'] = z
    residuals = [(y[q] - z[q])**2 for q in range(len(z))]
    ss_res = num.sum(residuals)
    ss_tot = num.sum((num.array(y)-num.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared = r_squared * 100
    D['r2'] = r_squared
    D['intercept'] = popt[1]
    D['slope'] = popt[0]
    return D

def regress_noplot(a,b):
    c = sorted(zip(a,b))
    a = [xx[0] for xx in c]
    b = [xx[1] for xx in c]
    D = pre_mania2(a,b)
    D['a'] = a
    D['b'] = b
    # D1 = maxBin3(a,b,10)
    return D



# hard max
def correct(x,y,N,m=1):

    run = True
    # if len(y) < N*m: # minimum number of points
    #     return (x,[num.log(max(xx,1)) for i,xx in enumerate(z)])


    elit = 1
    q = 0

    while True:
        outs = [] # output points
        l = len(x)
        tail = False
        for i,cur in enumerate(x):
            if i >= (l - 3):
                break
            ql = y[i:]

            tmp = 100 - q*elit
            ax = num.percentile(ql,tmp)
            if y[i] >= ax:
                outs.append([(cur,y[i]),(i,l-1)])

        q = q + 1
        if len(outs)>= 5:
            break

    at = [xx[0][0] for xx in outs]
    bt = [num.log(xx[0][1]) for xx in outs]
    popt, pcov = curve_fit(func, at, bt,maxfev = 10000)
    # a = sorted(list(set(a)))
    y = [num.log(xx) for xx in y]
    zf = zip(x,y)
    tmp = [xx[1]-popt[1]*xx[0] for xx in zf]
    return tmp


def correcttion(R,A):
    A_dic = {}
    C = {}
    for i,w in enumerate(A):
        A_dic[w['n1']+w['n2']] = i
    for i,w in enumerate(A):
        try:
            C[w['n2']+w['n1']]
            continue
        except KeyError:
            j = A_dic[w['n2']+w['n1']]
            v = A[j]
        e1 = w['envelope']
        e2 = v['envelope']
        R1 = R['s-'+w['n1']]
        R2 = R['t-'+w['n2']]
        if ( len(e1)>0 and len(e2)>0 ):
            m = []
            for p in e1:
                m.append(p[1]-R1[1]*p[0])
                m.append(p[1]-R2[1]*p[0])
            C[w['n1']+w['n2']] = num.median(m)
    return C




    # correct all data

    # both pairs have envelope points
    #
    # One direction envelope point - one direction max point
    #
    # One pair below noise
    #
    #
    pass


def roi_regressors(A,hem):
    l = len(A)
    D = {}
    for i in range(l):
        D[i] = {}
        D[i]['n1'] = A[i]['n1']
        D[i]['n2'] = A[i]['n2']
        if(not A[i]['isSuccess']):
            D[i]['x'] = []
            D[i]['y'] = []
            continue
        tmp = A[i]['envelope']
        ce = num.sum([xx[1] for xx in A[i]['envelope']])/len(A[i]['envelope'])
        cex = num.sum([xx[0] for xx in A[i]['envelope']])/len(A[i]['envelope'])
        D[i]['x'] = [xx[0]-cex for xx in tmp]
        D[i]['y'] = [xx[1]-ce for xx in tmp]
        D[i]['z'] = [xx-ce for xx in A[i]['z']]
        D[i]['slope'] = A[i]['slope']
        D[i]['r2'] = A[i]['r2']
    p = []
    C = {}
    N = [hem+str(i+1) for i in range(180)]
    R = {}
    for roi in N:
        # source taget
        datax = [xx['x'] for xx in D.values() if xx['n1'] == roi]
        datay = [xx['y'] for xx in D.values() if xx['n1'] == roi]
        fx = [item for sublist in datax for item in sublist]
        fy = [item for sublist in datay for item in sublist]
        F = lslinear(fx,fy)
        # estimator.fit(num.array(fx).reshape(-1, 1), num.array(fy).reshape(-1, 1))
        # y0 = estimator.predict(0)[0]
        # y1 = estimator.predict(1)[0]
        # z_pred = estimator.predict(alo)
        R["s-"+roi]=(F['r2'],F['slope'],F['intercept'])


        # target
        datax = [xx['x'] for xx in D.values() if xx['n2'] == roi]
        datay = [xx['y'] for xx in D.values() if xx['n2'] == roi]
        fx = [item for sublist in datax for item in sublist]
        fy = [item for sublist in datay for item in sublist]
        F = lslinear(fx,fy)
        # estimator.fit(num.array(fx).reshape(-1, 1), num.array(fy).reshape(-1, 1))
        # y0 = estimator.predict(0)[0]
        # y1 = estimator.predict(1)[0]
        # z_pred = estimator.predict(alo)
        R["t-"+roi]=(F['r2'],F['slope'],F['intercept'])
    return R

def mixed_regressor(R,local,r2,std,s,t):
    S = R["s-"+s]
    T = R["t-"+t]
    # max r2
    r = [r2,S[0],T[0]]
    sl = [local,S[1],T[1]]
    stds = [std,S[-1],T[-1]]

    # avg_local
    tmp = num.sum([sl[i]*r[i] for i in range(1,3)])/num.sum(r[1:])
    fx = (sl[0]*r[0]+(100-r[0])*tmp)/100
    a = r[0]/100
    den = num.sum(r[1:])
    b = (r[1]/den)*((100-r[0])/100)
    c = (r[2]/den)*((100-r[0])/100)
    cur = num.sqrt(a**2 * stds[0]**2 + b**2 * stds[1]**2 + c**2 * stds[2]**2)
    model = (fx,cur)
    return model

def all_regressors(R,A):
    l = len(A)
    D = {}
    # (slope, intercept, r2, std, kind=0,1,2,3)
    for i in range(l):
        n1 = A[i]['n1']
        n2 = A[i]['n2']
        if (A[i]['isSuccess']):
            D[n1+n2] = (A[i]['slope'],A[i]['intercept'],A[i]['std'],0)
        else:
            if (A[i]['reason'] == 'r2'):
                v1, v2 = mixed_regressor(R, A[i]['slope'],A[i]['r2'],A[i]['std'],n1,n2)
                D[n1+n2] = (v1,A[i]['intercept'],v2,1)
            else:
                v1, v2 = mixed_regressor(R, 0,0,0,n1,n2)
                D[n1+n2] = (v1,None,v2,2)
    for i in range(1,181):
        for j in range(1,181):
            if (i==j):continue
            n1 = 'L'+str(i)
            n2 = 'L'+str(j)
            try:
                D[n1+n2]
            except:
                v1, v2 = mixed_regressor(R, 0,0,0,n1,n2)
                D[n1+n2] = (v1,None,v2,3)
    return D

# a = [1,5.6,6,7,7,8,9,11,11,12.5,13,13.2,13.3,24,24,31,31.2,34,41,51.1,51.2]
# b = [10,43,60,17,98,80,91,1,2,12.5,13,13.2,13.3,24,24,31,31.2,34,41,51.1,51.2]
# outs,bins = maxBin(a,b,10)
# print(outs)
# print(bins)
# b=[0,2,21,3,-1,2]
# D = binCollapse(a,b,3,2)
# print(D)


    # collapse x,y
    # return x,y collapsed and candidate
