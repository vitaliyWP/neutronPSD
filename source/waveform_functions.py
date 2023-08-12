import numpy as np
from numpy import max, min, argmax, sqrt, exp, pi, argmin

c1 = '#dc2f18'
c2 = '#2d82b7'
c3 = '#07004d'

nice_fonts = {
    "text.usetex": True,
    "font.size" : 18,
    "axes.formatter.limits" :( -2, 4),
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    }


cosmic6runs = ['11639_']
cosmic3runs = ['11667_']

neutruns6  = ['11624_']
neutruns3  = ['11655']
neutruns0  = ['11771']

def extract_data(h1):
    a = np.array(h1.values(), dtype='int64')
    b = np.array(h1.axis().edges())
    b = (b[:-1] + (b[1] - b[0])/2) * 0.00665
    a = (a - np.bincount(a).argmax())
    a = a[:-1]
    b = b[:-1]

    grad1 = np.gradient(a)
    x1 = np.where(grad1 > 5)
    if (len(x1[0]) != 0 ):
        x = x1[0][0]
    else:
        return [0], [0]

    a = a[x-1:] 
    b = b[x-1:] - (x-1)*0.00665

    return a, b

def get_data(h1):
    a = np.array(h1.values(), dtype='int64')
    b = np.array(h1.axis().edges())
    b = (b[:-1] + (b[1] - b[0])/2) #* 0.00667
    a = (a - np.argmax(np.bincount(a)))
    a = a[:-1]
    b = b[:-1]

    return a, b

def get_hists(f, s, arapu):
    hists = list(f.keys())
    i=0
    tmph = []
    while i <len(hists):
        ch  = hists[i].strip().split('_')[6]
        #print(ch)
        if s in hists[i] and int(ch) in arapu:
            tmph.append(hists[i])
        i += 1
    #print(tmph)
    return tmph

'''def get_hists_nobeam(f, s, arapu):
    hists = list(f.keys())
    i=0
    tmph = []
    events = []
    while i <len(hists):
        ch  = hists[i].strip().split('_')[3]
        events.append(hists[i].strip().split('_')[1])
        #print(ch)
        if s in hists[i] and int(ch) in arapu:
            tmph.append(hists[i])
        i += 1
    events = set(events)
    #print(tmph)
    return tmph'''

def get_hists_nobeam2(f, s, arapu):
    hists = list(f.keys())
    i=0
    tmph = []
    events = []
    while i <len(hists):
        ch  = hists[i].strip().split('_')[3]
        events.append(hists[i].strip().split('_')[1])
        #print(ch)
        if s in hists[i] and int(ch) in arapu:
            tmph.append(hists[i])
        i += 1
    events = set(events)
    #print(tmph)
    return tmph, events

def get_hists4event(f, tag_event, channels):
    hists = list(f.keys())
    tmph = []
    i=0
    while i <len(hists):
        #print(hists[i].strip().split('_'))
        ch      = hists[i].strip().split('_')[3]
        event   = hists[i].strip().split('_')[1]
        
        if (int(ch) in channels) and (int(event) == int(tag_event)):
            tmph.append(hists[i])
            #print(ch, event)
        i += 1
    return tmph

'''def get_event_waveforms(f, s, arapu, beam):
    if beam:
        hists = get_hists(f, s, arapu)
        ch_position = 4
    else:
        hists, events = get_hists_nobeam(f, s, arapu)
        ch_position = 1
    muhis = []
    i = 0
    event0  = hists[i].strip().split('_')[ch_position]
    while i < len(hists):
        tmplist = []
        if i+1 == len(hists):
            break
        event1  = hists[i+1].strip().split('_')[ch_position]
        while event0 == event1:
            a, b = get_data(f[hists[i]])
            if np.max(a) > 100:
                tmplist.append(a)
            i += 1
            event0 = event1
            if i < len(hists):
                event1 = hists[i].strip().split('_')[ch_position]
            else:
                break
        event0 = event1
        i += 1
        if len(tmplist) > 1 :
            muhis.append(np.sum(tmplist, axis=0))
    return muhis'''

def get_average(hists1, loc, maxcut):
    hists = hists1[:]
    locmax = np.argmax(hists, axis=0)
    i = 0
    l1 = len(hists)
    while i < len(hists) or i < len(locmax):
        if locmax[i] > loc[1] or locmax[i]<loc[0] or np.max(hists[i]) < maxcut:
            hists.pop(i)
            locmax = np.delete(locmax, i)
        i += 1
    l     = len(hists)
    if l == 0:
        return
    thist = np.sum(hists, axis=0)
    print(f"Histos difference -- {l1 - l}")
    return thist/l

def align_histos(hists, msize, grad):
    aligned      = []
    gradProblems = 0

    for i in range(len(hists)):
        grad1 = np.gradient(hists[i])
        x1    = np.where(grad1 > grad)

        if len(x1[0]) == 0:
            gradProblems += 1
            continue

        y = list(hists[i][x1[0][0]-2:])
        if len(y)>msize:
            aligned.append(np.array(y[:msize+1]))
    print(len(aligned))
    print(gradProblems)

    return aligned

def get_risetime(aligned):
    risetimes = []

    for h in aligned:
        risetimes.append( (h[np.argmax(h)]-1 )*0.00665)

    return risetimes

def get_fwhm(aligned):
    fwhm = []

    for h in aligned:
        maxval = max(h)
        maxloc = argmax(h)

        y1 = np.array(h[:maxloc])/maxval
        y2 = np.array(h[maxloc:])/maxval
        if len(y1) < 1 or len(y2) < 1:
            print(y1, y2)
            continue
        x1 = argmin(abs(y1 - 0.5))
        x2 = argmin(abs(y2 - 0.5))

        
        fwhm.append( (maxloc + x2 - x1)*0.00665)

    return fwhm

#def get_fps_aligned(aligned):


def get_fps_bins(y):
    grad1 = np.gradient(y)
    x1 = np.where(grad1 > 10)
    if len(x1[0]) == 0:
        print(x1)
        return 0, 0, 0
    y = list(y[x1[0][0]-2:])

    maxy    = np.max(y)
    indmax  = np.argmax(y)
    #print(indmax)
    y1 = np.array(y[y.index(maxy):])/maxy
    y2 = np.array(y[:y.index(maxy)])/maxy

    idx = (np.abs(y1 - 0.85)).argmin()
    #idx2 = (np.abs(y2 - 0.1)).argmin()
    #idx3 = (np.abs(y2 - 0.9)).argmin()
    #idx4 = (np.abs(y1 - 0.5)).argmin()
    #idx5 = (np.abs(y2 - 0.5)).argmin()
    fp = np.sum(y[:indmax+idx])/np.sum(y)
    #fp2 = np.sum(y[:y.index(maxy)])/np.sum(y)
    fp3 = np.sum(y[indmax-50:indmax+100])/np.sum(y)
    return fp, fp3, np.sum(y)


def get_events(f):
    hists = list(f.keys())
    i=0
    events = []
    for i in range(len(hists)):
        events.append(hists[i].strip().split('_')[1])
    events = list(set(events))

    return events

def get_event_sums(f, events, cut):
    hists = list(f.keys())
    hists_cut = []
    for ev in events:
        hists2 = [hists[i] for i in range(len(hists)) if ev in hists[i]]
        #print(hists2)
        evsum = 0
        for i in range(len(hists2)):
            a, b = get_data(f[hists[i]])
            if np.max(a) > cut:
                a = a[np.where(a>0)]
                evsum += np.sum(a)
        hists_cut.append(evsum)
    return hists_cut
        
def cut_only_neutrons(f, events, lcut):
    hists = list(f.keys())
    revents = []
    for ev in events:
        hists2 = [hists[i] for i in range(len(hists)) if ev in hists[i]]
        evsum = 0
        argm = []
        for i in range(len(hists2)):
            a, b = get_data(f[hists2[i]])
            if np.max(a) > 100:
                argm.append(np.argmax(a))
                evsum += np.sum(a)
        #print(np.mean(argm))
        if (np.mean(argm) < 800 or np.mean(argm) > 1000) and evsum > lcut:
            #print(np.mean(argm))
            revents.append(ev)    
    return revents

def cut_only_muons(f, events, lcut):
    hists = list(f.keys())

    revents = []
    
    for ev in events:

        hists2  = [hists[i] for i in range(len(hists)) if ev in hists[i]]
        evsum   = 0
        argm    = []

        for i in range(len(hists2)):
            a, b = get_data(f[hists2[i]])
            if np.max(a) > 100:
                argm.append(np.argmax(a))
                evsum += np.sum(a)
        #print(np.mean(argm))
        if np.mean(argm) > 800 and np.mean(argm) < 1000 and evsum > lcut:
            #print(np.mean(argm))
            revents.append(ev)
    return revents


def sort_munu(f, events, lcut):
    hists = list(f.keys())

    mevents, nevents = [], []
    msum, nsum       = [], []
    maxlocs          = []
    for ev in events:

        hists2  = [hists[i] for i in range(len(hists)) if ev in hists[i]]
        evsum   = 0
        argm    = []

        for i in range(len(hists2)):
            a, b = get_data(f[hists2[i]])
            if np.max(a) > 100:
                argm.append(np.argmax(a))
                evsum += np.sum(a[np.where(a>0)])
                maxlocs.append(np.argmax(a))

        if evsum > lcut:
            if np.mean(argm) > 800 and np.mean(argm) < 1000 and evsum > lcut:
                mevents.append(ev)
                msum.append(evsum)
            else:
                nevents.append(ev)
                nsum.append(evsum)

    return mevents, nevents, msum, nsum, maxlocs



def get_hists_nobeam(f, s, arapu):
    hists = list(f.keys())

    i       = 0
    tmph    = []
    events  = []

    while i <len(hists):
        ch  = hists[i].strip().split('_')[3]
        events.append(hists[i].strip().split('_')[1])
        #print(ch)
        if s in hists[i] and int(ch) in arapu:
            tmph.append(hists[i])
        i += 1
    return tmph

def get_event_waveforms(f, s, arapu, beam):
    if beam:
        hists = get_hists(f, s, arapu)
        ch_position = 4
    else:
        hists= get_hists_nobeam(f, s, arapu)
        ch_position = 1
    muhis = []
    i = 0
    event0  = hists[i].strip().split('_')[ch_position]
    while i < len(hists):
        tmplist = []
        if i+1 == len(hists):
            break
        event1  = hists[i+1].strip().split('_')[ch_position]
        while event0 == event1:
            a, b = get_data(f[hists[i]])
            if np.max(a) > 200:
                tmplist.append(a)
            i += 1
            event0 = event1
            if i < len(hists):
                event1 = hists[i].strip().split('_')[ch_position]
            else:
                break
        event0 = event1
        i += 1
        if len(tmplist) > 1 :
            muhis.append(np.sum(tmplist, axis=0))
    return muhis


def get_amp_array(file, thisev, chns):

    rarray = np.zeros(30).reshape(10, 3)

    for a in range(3):

        channels = chns[a]
        for i in range(len(channels)):

            hists = get_hists4event(file, thisev, channels[i])

            if hists == []:
                continue
        
            for h in range(len(hists)):
                dat, xx = get_data(file[hists[h]])
                if np.max(dat) < 100:
                    continue
                dat = dat[np.where(dat>0)]
                rarray[i][a] += np.sum(dat)
    return np.array(rarray)

def get_hist_fps(y1):
    
    y = y1[np.where(y1>0)]
    
    grad1 = np.gradient(y)
    x1 = np.where(grad1 > 12)
    
    if len(x1[0]) == 0:
        print(x1)
        return 0

    y = list(y[x1[0][0]-2:])

    indmax  = np.argmax(y)

    fp1 = np.sum(y[:indmax])/np.sum(y)
    fp3 = np.sum(y[indmax-80:indmax+100])/np.sum(y)
    
    return fp3


def get_fps_bars(file, thisev, chns):

    rarray = np.zeros(30).reshape(10, 3)

    for a in range(3):

        channels = chns[a]
        for i in range(len(channels)):

            hists = get_hists4event(file, thisev, channels[i])

            if hists == []:
                continue
            
            hsum = []
            for h in range(len(hists)):
                dat, xx = get_data(file[hists[h]])
                if np.max(dat) < 100:
                    continue
                hsum.append(dat)
            
            if len(hsum) < 1:
                continue
            hsum = np.sum(hsum, axis=0)
            
            rarray[i][a] = get_hist_fps(hsum)
                
    return rarray
