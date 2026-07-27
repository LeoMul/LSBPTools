'''
Reads AUTOSTRUCTURE outputs and generates R-matrix input file 
To be used for streamlining of calculations 
'''

from rmatrix import *


import argparse

parser = argparse.ArgumentParser()


parser.add_argument('--minst', default=MINST_DEFAULT,  help='Min value of 2S+1',type=int)
parser.add_argument('--maxst', default=MAXST_DEFAULT,  help='Max value of 2S+1',type=int)
#parser.add_argument('--minlt', default=MINLT_DEFAULT,  help='Min value of L',type=int)
parser.add_argument('--j2min', default=J2MIN_DEFAULT,  help='Min value of 2J',type=int)
parser.add_argument('--j2max', default=J2MAX_DEFAULT,  help='Max value of 2J',type=int)
parser.add_argument('--maxc',  default=MAXCDEFAULT,  help='Number of continuum orbitals.',type=int)
parser.add_argument('--maxe',  default=MAXEDEFAULT,  help='Target max scattering energy.',type=int)

args = parser.parse_args()



def main():
    
    calc = rmatrix(maxc=args.maxc,
                   maxe=args.maxe,
                   J2MIN=args.j2min,
                   J2MAX=args.j2max,
                   MINST=args.minst,
                   MAXST=args.maxst)
    
    calc.getConfigs()
    calc.generateCorrelationConfigs()
    calc.setScatteringMomenta()
    calc.getLSPi()
    calc.writeDSTG1()
    calc.writeDSTG2()
    calc.getJPi()
    calc.writeDSTGJK()

    return 0 

main()