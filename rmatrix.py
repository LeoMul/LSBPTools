STANDARD_ORDER     = [ (1,0), (2,0), (2,1), (3,0), (3,1), (3,2), (4,0), (4,1), (4,2), (4,3)]
STANDARD_ORDER_STR = ['  1 0','  2 0','  2 1','  3 0','  3 1','  3 2','  4 0','  4 1','  4 2','  4 3']
STANDARD_ORDER_MAX = [    2,     2,     6,     2,     6,     10,    2,     6,     10,    14]

MAXCDEFAULT = 10 
MAXEDEFAULT =  3.0

J2MIN_DEFAULT =  1 
J2MAX_DEFAULT = 33
 
MINST_DEFAULT =  2 
MAXST_DEFAULT =  6 

MINLT_DEFAULT =  0


import numpy as np 

class rmatrix:
    def __init__(self,maxc = MAXCDEFAULT, maxe = MAXEDEFAULT, J2MIN=-1,J2MAX = -1, MINST=-1, MAXST=-1) -> None:
        
        self.termsS   = []
        self.termsL   = []
        self.termsP   = []
        self.levelsJ   = []
        self.levelsP   = []        
        
        self.orbitals = []
        self.termstrings = []
        self.levelstrings = []

        self.maxc     = maxc
        self.maxe     = maxe 



        self.MINST    = MINST
        self.MAXST    = MAXST 
        self.J2MIN    = J2MIN 
        self.J2MAX    = J2MAX            
        pass

    def getConfigs(self):
        
        file = open('CONFIG.DAT','r')
        
        firstLine = file.readline().split()
        
        kcor1 = 1 
        kcor2 = int(firstLine[0])
        totalNumOrbs = int(firstLine[1])
        orbitalsLine = file.readline()
        orbitalsLineSplit = orbitalsLine.split()
        
        self.orbitalString =  orbitalsLine
        
        self.kcor1 = kcor1
        self.kcor2 = kcor2
        self.orbitalValenceNum = totalNumOrbs - kcor2
        self.orbitalTotalNum   = totalNumOrbs
        self.valenceorbitals   = [(int(a), int(b)) for a, b in zip(orbitalsLineSplit[0::2], orbitalsLineSplit[1::2])]
        
        self.hardMax = []
        for (_,b) in self.valenceorbitals:
            self.hardMax.append(4*b+1)
        
        
        self.numConfigs = int(file.readline().split()[0])
        
        self.minoccline = file.readline()
        self.maxoccline = file.readline()

        self.configsMin = [int(x) for x in self.minoccline.split()]
        self.configsMax = [int(x) for x in self.maxoccline.split()]

        self.configstrings = []
        self.configsNumpy = np.zeros( [ self.numConfigs,self.orbitalValenceNum+1],dtype=int )
        
        
        for ii in range(0, self.numConfigs):
            line = file.readline()
            self.configstrings.append(line)
            self.configsNumpy[ii,:] = np.array([int(x) for x in line.split()])
        
        #print(self.configsNumpy[:, 0:self.orbitalValenceNum ])
        self.nelecsum = np.sum ( self.configsNumpy[:, 0:self.orbitalValenceNum ],axis=1  )
        
        if np.any(self.nelecsum != self.nelecsum[0]):
            print("Error in CONFIG.DAT - inconsistent number of electrons.")
            import sys 
            sys.exit()
        
        self.NELEC = self.nelecsum[0] + np.sum( STANDARD_ORDER_MAX[0:self.kcor2])
        self.corestrin = "".join(STANDARD_ORDER_STR[0:self.kcor2])
        self.coreoccstring = "".join( [' '+ str(x) for x in STANDARD_ORDER_MAX[0:self.kcor2]])
        return None 
    
    def generateCorrelationConfigs(self):
        
        #self.configsCorrelationNumpy = np.zeros_like(self.configsNumpy)
        
        self.correlationConfigs = generate_correlation_configurations(self.configsNumpy[:, 0:self.orbitalValenceNum ], self.hardMax)
        
        self.correlationConfigsNumpy = np.array(self.correlationConfigs)
        print(self.correlationConfigsNumpy)
        
        self.correlationConfigsMaxs = np.max(self.correlationConfigsNumpy,axis=0)
        self.correlationConfigsMins = np.min(self.correlationConfigsNumpy,axis=0)

        print(self.correlationConfigsMins)
        print(self.correlationConfigsMaxs)
        
        return None 
    

    def setScatteringMomenta(self):
        
        if self.nelecsum[0] % 2 == 1:
            #Odd number of electrons
            #N+1 is then even. 
            #scattered 2S+1 is odd 
            #And scattered 2J is even.
            
            if self.J2MIN == J2MIN_DEFAULT:
                self.J2MIN = J2MIN_DEFAULT-1
            if self.J2MAX == J2MAX_DEFAULT:
                self.J2MAX = J2MAX_DEFAULT+1
            if self.MAXST == MAXST_DEFAULT:
                self.MAXST = MAXST_DEFAULT-1
            if self.MINST == MINST_DEFAULT:
                self.MINST = MINST_DEFAULT-1
                
        else:
            #Even number of electrons 
            #N+1 is then odd.  
            #scattered 2S+1 is even 
            #And scattered 2J is odd.
            self.J2MAX = J2MAX_DEFAULT

            
        
        '''
        triangle inequality, for a requested J2MAX with max multiplicity
        MAXST -  we need to go up to at least this L. 
        '''
        
        sarray = [self.MINST]
        
        current = self.MINST
        while (current!= self.MAXST):
            current +=2 
            sarray.append(current)
        
        sarray = np.array(sarray)
        sarray = (sarray-1)
        print(sarray)
        
        #minltcheck = int( min( abs(  self.J2MIN - sarray) / 2 ) ) 
        #print(minltcheck)
        self.MINLT = int( min( abs(  self.J2MIN - sarray) / 2 ) ) 
        self.MAXLT = int ( ( self.J2MAX + self.MAXST - 1 ) / 2 ) 

        self.numcpus_stg2 = (self.MAXLT - self.MINLT + 1) * (self.MAXST - self.MINST + 2)

        print('Input:')
        print(' NELEC      = ',self.NELEC)
        print(' J2MIN      = ',self.J2MIN)
        print(' J2MAX      = ',self.J2MAX)
        print(' MINST      = ',self.MINST)
        print(' MAXST      = ',self.MAXST)

        print('Setting:')
        
        print(' MINLT      = ',self.MINLT)
        print(' MAXLT      = ',self.MAXLT)
        print(' NUMCPUSTG2 = ',self.numcpus_stg2)

        
        if ( (self.J2MIN % 2 != ( (self.NELEC+ 1) % 2 )) or (self.J2MAX % 2 != ( (self.NELEC+ 1) % 2 )) or (self.MAXST % 2 != self.NELEC % 2 ) or (self.MINST % 2 != self.NELEC % 2 ) ):
            print("Invalid parity in 2S+1 or 2J.")
            #print(self.J2MIN % 2 != ((self.NELEC+ 1) % 2 ))
            import sys 
            sys.exit()
            
        return None 
    
    
    def getJPi(self):
        file = open('LEVELS','r')
        
        file.readline()
        unFinished = True 
        while (unFinished):
            line = file.readline()
            
            
            
            J = int(line[0:2])
            P = int(line[2:4])

            if int(line[7]) == 0: 
                unFinished = False 
            else: 
                self.levelsJ.append(J)
                self.levelsP.append(P)
                self.levelstrings.append(line[0:4])

        return None
    
    def getLSPi(self):
        
        file = open('TERMS','r')
        
        file.readline()
        unFinished = True 
        while (unFinished):
            line = file.readline()
            
            
            S = int(line[0:2])
            L = int(line[2:4])
            P = int(line[4:6])

            if int(S) == 0: 
                self.lamax = max(self.termsL)
                unFinished = False 
            else: 
                self.termsS.append(S)
                self.termsL.append(L)
                self.termsP.append(P)
                self.termstrings.append(line[0:6])

        return None 
    
    def writeDSTG1(self):
        
        dstg1 = open('dstg1','w')
        dstg1.write('S.S.\n')
        dstg1.write(f"&STG1A RELOP='YES' &END \n")
        dstg1.write(f"&STG1B MAXLA={self.lamax} MAXLT={self.MAXLT} MAXC={self.maxc} MAXE={int(round(self.maxe,0))} &END \n")

        dstg1.close()
        
        return None
    
    def writeDSTG2(self):

        
        dstg2 = open('dstg2','w')
        dstg2.write('S.S.\n')
        dstg2.write(f"&STG2A RELOP='YES' isort=1 &END \n")
        dstg2.write(f"&STG2B MAXORB={self.orbitalTotalNum} NELC={self.NELEC} NAST={len(self.termsL)} INAST=0 MINLT={self.MINLT} MAXLT={self.MAXLT} MINST={self.MINST} MAXST={self.MAXST}  &END \n")
        dstg2.write(self.corestrin+self.orbitalString)
        dstg2.write(f'{self.numConfigs}\n')
        dstg2.write(self.coreoccstring+self.minoccline)
        dstg2.write(self.coreoccstring+self.maxoccline)
        for config in self.configstrings:
            dstg2.write(self.coreoccstring+config)

        for term in self.termstrings:
            dstg2.write(term+'\n')
        
        #print(len(self.correlationConfigs),len(self.correlationConfigs[0]))
        
        numcorrelationConfigs = len(self.correlationConfigs)
        
        dstg2.write('{:3}\n'.format(numcorrelationConfigs))
        
        string = ''
        for jj in self.correlationConfigsMins:
            string += '{:>3}'.format(jj)
        string = self.coreoccstring + string
        string += '\n'
        dstg2.write(string)
        string = ''
        for jj in self.correlationConfigsMaxs:
            string += '{:>3}'.format(jj)
        string = self.coreoccstring + string

        string += '\n'
        dstg2.write(string)
        
        for ii in range(0,numcorrelationConfigs):
            string = ''
            for jj in range(0,len(self.correlationConfigs[ii])):
                string += '{:>3}'.format(self.correlationConfigs[ii][jj])
            string += '  0\n'
            string = self.coreoccstring + string
            dstg2.write(string)
            
        dstg2.close()
        
        return None
    
    def writeDSTGJK(self):
        
        dstgjk = open('dstgjk','w')
        dstgjk.write('S.S.\n')
        dstgjk.write(f"&STGJA RELOP='YES' &END \n")
        dstgjk.write(f"&STGJB JNAST={len(self.levelstrings)} IJNAST=0 J2MIN={self.J2MIN} J2MAX={self.J2MAX} &END \n")

        for level in self.levelstrings:
            dstgjk.write(level+'\n')
        
        dstgjk.close()
        
        return None 
    

def generate_correlation_configurations(input_configs, max_capacities=None):
    #gemini generated
    """
    Generates unique (N+1)-electron correlation configurations by adding 
    a single electron to each N-electron target configuration.
    
    Parameters:
    -----------
    input_configs : list of lists/tuples
        The N-electron base configurations (e.g., [[2, 2, 0]])
    max_capacities : list, optional
        The maximum electron occupancy for each slot. If None, it defaults 
        to the standard atomic shell sequence: s=2, p=6, d=10, f=14, g=18...
        
    Returns:
    --------
    list of lists
        A sorted list of all unique (N+1)-electron configurations.
    """
    #if not input_configs:
    #    return []
        
    # Determine the number of subshells from the first configuration
    num_subshells = len(input_configs[0])
    
    # Default to standard capacities: 2 for s (idx 0), 6 for p (idx 1), 10 for d (idx 2)...
    if max_capacities is None:
        import sys
        print('No max occs for correlation configs') 
        sys.exit()
        
    unique_correlation_set = set()
    
    for config in input_configs:
        for i in range(num_subshells):
            # Check if the subshell has room for another electron
            if config[i] < max_capacities[i]:
                # Create a copy, add the electron, and store as a hashable tuple
                new_config = list(config)
                new_config[i] += 1
                unique_correlation_set.add(tuple(new_config))
                
    # Convert back to a sorted list of lists for clean output
    return sorted([list(cfg) for cfg in unique_correlation_set])

