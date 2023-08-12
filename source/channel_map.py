#! python
"""
File contains definitions for various groups of channels;

Ex.:
By sensor type : sensel or hamamatsu
By location:     Anaode plane assembleys (apa1, apa2, ...)


"""

import numpy as np


chan5 = [np.arange(216, 220), np.arange(220, 224),
         np.arange(224, 228), np.arange(228, 232),
         np.arange(192, 196), np.arange(232, 236),
         np.arange(196, 200), np.arange(200, 204),
         np.arange(236, 240), np.arange(204, 208)]

chan6 = [np.arange(240, 244), np.arange(244, 248),
         np.arange(248, 252), np.arange(252, 256),
         np.arange(256, 260), np.arange(264, 276),
         np.arange(260, 264), np.arange(276, 280),
         np.arange(280, 284), np.arange(284, 288)]

chan4 = [np.arange(144, 148), np.arange(148, 152),
         np.arange(152, 156), np.arange(156, 160),
         np.arange(160, 164), np.arange(164, 168),
         np.arange(168, 172), np.arange(172, 176),
         np.arange(176, 180), np.arange(180, 184)]

chan1 = [np.arange(0, 4), np.arange(4, 8),
         np.arange(8, 12), np.arange(12, 16),
         np.arange(16, 20), np.arange(20, 24),
         np.arange(24, 28), np.arange(28, 32),
         np.arange(32, 36), np.arange(36, 40)]

chan2 = [np.arange(48, 52), np.arange(52, 56),
         np.arange(56, 60), np.arange(60, 64),
         np.arange(64, 68), np.arange(68, 72),
         np.arange(72, 76), np.arange(76, 80),
         np.arange(80, 84), np.arange(84, 88)]

chan3 = [np.arange(96, 100), np.arange(100, 104),
         np.arange(104, 108), np.arange(132, 144),
         np.arange(108, 112), np.arange(112, 116),
         np.arange(116, 120), np.arange(120, 124),
         np.arange(124, 128), np.arange(128, 132)]

chns_daqside  = [chan5, chan6, chan4]
chns_rackside = [chan3, chan2, chan1]

apa1 = np.concatenate(chan1)
apa2 = np.concatenate(chan2)
apa3 = np.concatenate(chan3)
apa4 = np.concatenate(chan4)
apa5 = np.concatenate(chan5)
apa6 = np.concatenate(chan6)

arachns = [np.arange(132, 144), np.arange(264, 276)]

allCh = np.concatenate((apa1, apa2, apa3, apa4, apa5, apa6), axis=0)
apas = [apa1, apa2, apa3, apa4, apa5, apa6]

sensL = np.concatenate((np.arange(40), np.arange(48, 88),
                        np.arange(96, 108), np.arange(108, 132), np.arange(144, 184),
                        np.arange(192, 204)), axis=0)
arapu = np.concatenate((np.arange(132, 144), np.arange(264, 276)), axis=0)








if __name__ == '__main__':

    print("\nLight bars in the first Anode Assembley Plane")
    print(chan1)
    print("\nChannels in the first Anode Assembley Plane")
    print(apa1)