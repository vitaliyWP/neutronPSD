#! python

import numpy as np
from channel_map import arapu, arachns
from uproot import open


class waveformsFile:
    """
    Class for the custom *.root file containing histograms with waveforms;
    Histogram name example : 'waveform_raw_mu_evt_391_ch_264_track_1;1'
    """
    EVENT_POSITION      = 4 # in histogram name 
    CHANNEL_POSITION    = 6 # in histogram name

    def __init__(self, filename):
        self.file   = open(filename)
        self.events = self.get_event_list()
        #self.nevent = self.get_neut_events()


    def get_hists(self, event=None, channels=[]):
        all = list(self.file.keys())

        if event == None and channels == []:
            return all
        
        elif event != None:
            if channels != []:
                #print(hist.split('_')[self.CHANNEL_POSITION])
                selected = [hist for hist in all if (event == hist.split('_')[self.EVENT_POSITION] and 
                            int(hist.split('_')[self.CHANNEL_POSITION]) in channels)]
                #print(selected)
            else:
                selected = [hist for hist in all if event == hist.split('_')[self.EVENT_POSITION]]
        else:
            selected = [hist for hist in all if hist.split('_')[self.CHANNEL_POSITION] in channels]
        
        return selected
    
    def get_event_list(self):
        hists = self.get_hists()

        events = []
        for i in range(len(hists)):
            events.append(hists[i].strip().split('_')[self.EVENT_POSITION])
    
        return set(events)
    
    def get_data(self, hname):

        hist = self.file[hname]
        waveform  = np.array(hist.values(), dtype='int32')
        waveform = (waveform - np.argmax(np.bincount(waveform)))[:-1]

        return waveform
    
    def get_bar_waveforms(self, event, chs):
        
        hists = self.get_hists(event=event, channels=chs)
        waves = []
        for hist in hists:
            waves.append(self.get_data(hist))

        return waves

    def get_bar_sum(self, event, chs):
        return np.sum(self.get_bar_waveforms(event, chs), axis=0)
    
    # 
    # 
    def get_neut_events(self, chs):
        nevents = []
        for ev in self.events:
            sumw = self.get_bar_sum(ev, chs)
            #print(sumw)
            locmax = np.argmax(sumw)
            if locmax > 250 and locmax<1700 and np.max(sumw)>400:
                nevents.append(ev)
        
        return nevents
        

class Waveform:
    """
    Main class for SiPM waveforms analysis

    """
    def __init__(self, data):
        self.data = data
    
    def max(self):
        return np.max(self.data)
    
    def locmax(self):
        return np.argmax(self.data)
    
    def get_integral(self):
        return np.sum(self.data[ np.where(self.data > 0) ])
    
    def get_fps(self):
        ...
        #gradself.data


if __name__ == '__main__':

    filename = '/Users/vitaliy/Desktop/MicroBoone/neutronData/output_beam_waveforms_filtered11657.root'
    r57 = waveformsFile(filename=filename)
    eventList = r57.get_event_list()
    #print(eventList)
    print(type(None))
    #print(r57.get_hists(eventList[0], arapu))
