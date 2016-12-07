#!/usr/bin/python

'''
python multiprocessing scan offsets for 5 chambers
'''

import multiprocessing
from multiprocessing import Process, Value, Array
from ROOT import TH1F, TFile, TF1, TH2F, TMath
import os, time


file = TFile("res.root", "recreate")

#########################################################
z = [0, 132.487, 266.47, 397.00, 529]

NN = 100
scan_range = 100.0 #mm
z_step = [ -scan_range/2 + scan_range/NN*float(i) for i in range(NN) ]
x_step = [ -scan_range/2 + scan_range/NN*float(i) for i in range(NN) ]

PI = 3.14159265358793
theta_range = 50.0*PI/180.0 #rad
theta_step = [ -theta_range/2 + theta_range/float(NN)*float(i) for i in range(NN)]

h_scan = TH2F("h_scan", "h_scan", NN, -scan_range/2, scan_range/2, NN, -scan_range/2, scan_range/2)


#########################################################
# offset results
# z_offset = [0, 0, -15.5, 0, -30.3]
z_offset = [0, 0, 0, 0, 0]
x_offset = [0, 0, 0, 0, 0]


#########################################################
# load file to memory
x0 = []
x2 = []
x4 = []
y0 = []
y2 = []
y4 = []

def load_file(name):
    try:
        file = open(name, 'r')
	print name, 'openned.'
    except Exception, e:
        print e

    for line in file:
        ll = line.split()
	x0.append( float(ll[0]) )
	x2.append( float(ll[1]) )
	x4.append( float(ll[2]) )

	y0.append( float(ll[3]) )
	y2.append( float(ll[4]) )
	y4.append( float(ll[5]) )

    print 'length:', len(x0), len(x2), len(x4), len(y0), len(y2), len(y4)


#########################################################
# for printing, no useful
def apply_offset(x0, x4):
    return float(x0)+x_offset[0] + ( float(x4)+x_offset[4]-float(x0)+x_offset[0] )*(z[2]+z_offset[2])/(z[4]+z_offset[4])

def get_resolution_after_correction(name):
    load_file(name)

    hname = 'h_resolution'
    h = TH1F(hname,hname,1000, -10, 10);

    entry= len(x0)
    index = 0
    for index in range(entry):
        x_project = apply_offset(x0[index], x4[index])
        xm = x2[index] - x_project
	h.Fill(xm)

    _ax = h.GetMaximumBin()
    ax = h.GetXaxis().GetBinCenter(_ax)
    h.Fit("gaus","Q0", "", ax-3.0, ax+3.0)
    h.Write()


'''
---------------------------------------------------------
functioning core start
---------------------------------------------------------
'''
#########################################################
# core fuc (single loop)
def stat_rms(z2offset, z4offset, x2offset, x4offset):
    z2 = z[2] + z2offset
    z4 = z[4] + z4offset

    hname = 'h'
    h = TH1F(hname,hname,10000, -100, 100);

    entry= len(x0)
    index = 0
    for index in range(entry):
        x_project =x0[index] + (x4[index]+x4offset-x0[index])*z2/z4
        xm = x2[index]+x2offset - x_project
	h.Fill(xm)

    _ax = h.GetMaximumBin()
    ax = h.GetXaxis().GetBinCenter(_ax)
    h.Fit("gaus","Q0", "", ax-3.0, ax+3.0)
    res = h.GetFunction("gaus").GetParameter(2)
    return res

def stat_theta(theta2offset, theta4offset):
    theta2 = theta2offset
    theta4 = theta4offset

    hname = 'h_'
    h = TH1F(hname,hname,10000, -100, 100);

    entry= len(x0)
    index = 0
    for index in range(entry):
        tx4 = (x4[index])*TMath.Cos(theta4)
	tz4 = z[4] + (x4[index])*TMath.Sin(theta4)
	tx2 = (x2[index])*TMath.Cos(theta2)
	tz2 = z[2] + (x2[index])*TMath.Sin(theta2)
	tx0 = x0[index]

        x_project =tx0 + (tx4 - tx0)*tz2/tz4
        xm = tx2 - x_project
	h.Fill(xm)

    _ax = h.GetMaximumBin()
    ax = h.GetXaxis().GetBinCenter(_ax)
    h.Fit("gaus","Q0", "", ax-3.0, ax+3.0)
    res = h.GetFunction("gaus").GetParameter(2)
    return res

def get_rms_single_process(start, end, arr):
    t1 = 0.0
    t2 = 0.0
    for i in range(start, end):
        for j in range(NN):
	    # scan z
	    #res = stat_rms(z_step[i], z_step[j], t1, t2)

	    # scan x
	    #res = stat_rms(t1, t2, x_step[i], x_step[j])

	    # scan kx
	    res = stat_theta(theta_step[i], theta_step[j])

            index = (i-start)*NN + j
	    arr[index] = res

def main_process(name):
    load_file(name)

    quantity = int(NN/4) * int(NN)
    d = [Array('d', quantity) for i in range(4)]

    job = [0, int(NN/4), int(NN/4)*2, int(NN/4)*3, NN]

    start_time = time.time()
    p = [ Process(target = get_rms_single_process, args=(job[i], job[i+1], d[i])) for i in range(4) ]
    for sp in p: sp.start()
    for sp in p: sp.join()

    # fill results to histogram
    for i in range(4):
        start = i * int(NN/4)
	index = 0
	for res in d[i]:
	    y = index % NN
	    x = int(index/NN) + start
	    h_scan.SetBinContent(x, y, res)
	    index += 1

    h_scan.Write()
    elapsed_time = time.time() - start_time
    print 'total elapsed time:', elapsed_time
'''
---------------------------------------------------------
functioning core end
---------------------------------------------------------
'''

#########################################################
# multiprocessing scan dz or dx or dy
def get_rms(ii, jj, aa, bb, d):
    '''
    if scan z, put aa = NN/2, bb = NN/2 (NN=100 or 1000)
    if scan x, put ii = NN/2, jj = NN/2 (NN= 100 or 100)
    '''
    i = int(ii)
    j = int(jj)
    z2 = z[2] + z_step[i]
    z4 = z[4] + z_step[j]

    a = int(aa)
    b = int(bb)
    x2offset = x_step[a]
    x4offset = x_step[b]

    hname = 'h_'+str(i)+'_'+str(j)+'_'+str(aa)+'_'+str(bb)
    h = TH1F(hname,hname,10000, -100, 100);

    entry= len(x0)
    index = 0
    for index in range(entry):
        x_project =x0[index] + (x4[index]+x4offset-x0[index])*z2/z4
        xm = x2[index]+x2offset - x_project
	h.Fill(xm)

    _ax = h.GetMaximumBin()
    ax = h.GetXaxis().GetBinCenter(_ax)
    h.Fit("gaus","Q0", "", ax-3.0, ax+3.0)
    h.SetDirectory(file)
    res = h.GetFunction("gaus").GetParameter(2)
    d.value = res

#######################################################
# multiprocessing scan theta
def rms_scan_theta(ii, jj, d):
    i = int(ii)
    j = int(jj)

    theta2 = theta_step[i]
    theta4 = theta_step[j]

    hname = 'h_'+str(i)+'_'+str(j)
    h = TH1F(hname,hname,10000, -100, 100);

    entry= len(x0)
    index = 0
    for index in range(entry):
        tx4 = (x4[index])*TMath.Cos(theta4)
	tz4 = z[4] + (x4[index])*TMath.Sin(theta4)
	tx2 = (x2[index])*TMath.Cos(theta2)
	tz2 = z[2] + (x2[index])*TMath.Sin(theta2)
	tx0 = x0[index]

        x_project =tx0 + (tx4 - tx0)*tz2/tz4
        xm = tx2 - x_project
	h.Fill(xm)

    _ax = h.GetMaximumBin()
    ax = h.GetXaxis().GetBinCenter(_ax)
    h.Fit("gaus","Q0", "", ax-3.0, ax+3.0)
    #h.SetDirectory(file)
    #h.Write()
    res = h.GetFunction("gaus").GetParameter(2)
    d.value = res


########################################################
# single core, no multiprocessing
def get_rms_single(ii,jj):
    i = int(ii)
    j = int(jj)
    z2 = z[2] + z_step[i]
    z4 = z[4] + z_step[j]

    hname = 'h_'+str(i)+'_'+str(j)
    h = TH1F(hname,hname,1000, -10, 10);

    entry= len(x0)
    index = 0
    for index in range(entry):
        x_project =x0[index] + (x4[index]-x0[index])*z2/z4
        xm = x2[index] - x_project
	h.Fill(xm)

    _ax = h.GetMaximumBin()
    ax = h.GetXaxis().GetBinCenter(_ax)
    h.Fit("gaus","Q0", "", ax-3.0, ax+3.0)
    res = h.GetFunction("gaus").GetParameter(2)
    return res


#########################################################
def main(name):
    load_file(name)

    d = [ Value('d', -99.0) for i in range(6) ]

    for i in range(NN):
        start_time = time.time()
        for j in range(0, NN, 4):
	    # scan z
            p = [ Process(target = get_rms, args=(i, j+k, NN/2, NN/2, d[k])) for k in range(4) ]

	    # scan x
            #p = [ Process(target = get_rms, args=(NN/2, NN/2, i, j+k,  d[k])) for k in range(4) ]

	    # scan theta
	    #p = [ Process(target = rms_scan_theta, args=(i, j+k, d[k])) for k in range(4) ]

	    for sp in p: sp.start()
	    for sp in p: sp.join()

	    for k in range(4):
	        h_scan.SetBinContent(i, j+k, d[k].value)
        print i,'/',NN,'elapsed time:', time.time() - start_time
    h_scan.Write()


def single(name):
    load_file(name)
    for i in range(10):
        for j in range(10):
	    res = get_rms_single(i, j)
	    print res
	    h_scan.SetBinContent(i, j, res)

    h_scan.Write()

if __name__=='__main__':
    #main('MinInput_288_2.dat')
    #single('MinInput_2.dat')
    #get_resolution_after_correction('MinInput_2.dat')
    main_process('MinInput_288_2.dat')
