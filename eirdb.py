# Script to create en EIRRUN database of eirene cases
# Based on database.py by holm10
# Created May 29th 2020


def natsort(l): 
    from re import split
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)



class SETUP():
    ''' Class containing a list of CASE objects '''
    def __init__(self,caselist):
        ''' Stores the list sorted by midplane separatrix electron density '''
        self.cases=caselist
        self.ev=1.602e-19



    '''==========================================
    Handle the case list
    =========================================='''


    def sort_ind(self,var,ind, **kwargs):
        ''' Sorts list by value of var at index (ix,iy) 
        sort_mp(var,ix,iy,**keys)

        Variables:
        var:        String of variable to be used in sorting, e.g. 'bbb.ne'
        ix:         Polidal index for sorting
        iy:         Radial index for sorting

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0

        Returns:
        Void
        '''
        # Check whether the requested parameter has a species index
        if isinstance(ind,int): # EIRENE index
            ind = ind # Use as is
        elif len(ind) == 3: # UEDGE index with triangle specified
            ind = self.cases[0].UEarr[ind[0], ind[1], ind[2]]
        else:
            print('Unspecified index!')
            return
        self.cases.sort(key=lambda case: case.get(var, **kwargs)[ind])


    def sort_rowmax(self,var,row,s=0,supress=False):
        ''' Sorts list by value of var at index (ix,iy) 
        sort_mp(var,row,**keys)

        Variables:
        var:        String of variable to be used in sorting, e.g. 'bbb.ne'
        Row:        Poloidal row along which the max value is found, and used to sort the cases

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0

        Returns:
        Void
        '''
        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            if supress is False:
                print('WARNING! Requested row is guard cell row')

        # Check whether the requested parameter has a species index
        self.cases.sort(key=lambda case: max(case.get(var)[row,:]))

    def get_closest_index(self,var,val,ix,iy,s=0):
        ''' Get the index of the case with var closest to val at location specified by index'''
        return abs(self.index(var,ix,iy,s=s)-val).argmin()


    '''==========================================
    Get parameter
    =========================================='''

    def get(self, var, **kwargs):
        ''' Returns a list of arrays containing var
        get(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        if 'processing' in kwargs:
            if kwargs.get('processing') in []:
                kwargs['processing'] = 'raw' # Set the processing type to be passed

        return asarray([case.get(var, **kwargs) for case in self.cases])
        # Implement checks for further operations

    def get_P_throughflux(self, istra, species):
        from numpy import asarray
        return asarray([case.P_throughflux(istra, species) for case in self.cases])

    def get_E_throughflux(self, istra, species):
        from numpy import asarray
        return asarray([case.E_throughflux(istra, species) for case in self.cases])

    def getZ(self,var):
        from numpy import asarray
        return asarray([case.getZ(var) for case in self.cases])
        

    def getR(self,var):
        from numpy import asarray
        return asarray([case.getR(var) for case in self.cases])
        

    def getN(self,var,s=0):
        from numpy import asarray
        return asarray([case.getN(var,s=s) for case in self.cases])
        
    def getRC(self):
        ''' Gets recombination from eirene out '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),2))
        for i in range(len(self.cases)):
            c=self.cases[i]
            with open(c.path+'/eirene.out','rt') as outfile:
                contents=outfile.read()
            contents=contents[contents.find("SUM OVER THE STRATA"):]
            reca=contents.find("PPATI  = \n")
            try:
                ret[i,0]=float(contents[reca+13:reca+33].strip())
            except: pass
            recm=contents.find("PPMTI  = \n")
            try:
                ret[i,1]=float(contents[recm+13:recm+33].strip())
            except: pass
        return ret


    def get_STRATASUM(self,path,text='SUM OVER THE STRATA',linelist=False):

        with open(path+'/eirene.out','rt') as outfile:
            if linelist is True:
                contents=[l.strip() for l in outfile] 
                indices = [i for i, elem in enumerate(contents) if text in elem]
                contents=contents[indices[0]:]
            else:
                contents=outfile.read()
                contents=contents[contents.find(text):]
        
        return contents
        
    def get_runtime(self):
        from numpy import zeros

        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            ret[i]=float(self.get_STRATASUM(self.cases[i].path,'TOTAL CPU_TIME OF THIS RUN',True)[0].split(' ')[-1].strip())
        return ret

    def get_text(self,text,x,l):
        ''' Gets recombination from eirene out '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=self.get_STRATASUM(self.cases[i].path)
            nhist=contents.find(text)
            ret[i]=float(contents[nhist+x:nhist+x+l].strip())
        
        return ret

    def get_incidentflux(self,NSPEZ=2):
        return self.get_pflux('FLUX INCIDENT ON SURFACE:',NSPEZ=2)


    def get_total_reemittedflux(self,NSPEZ=2,REC=False):
        from numpy import zeros
        ret=self.get_reemittedflux('ATOMS')
        ret+=self.get_reemittedflux('MOLECULES')
        ret+=self.get_reemittedflux('BULK IONS')*False
    
        return ret

    def get_reemittedflux(self,species,NSPEZ=2):
        from numpy import roll
        if species.upper()=='BULK IONS': d=1
        else: d=0
        ret= self.get_pflux('FLUX RE-EMITTED FROM INCIDENT {}:'.format(species.upper()),NSPEZ=2,displace=d)
        if species.upper()=='MOLECULES':
            ret=roll(ret,1,axis=2)
            ret[:,:,0]=0
        return ret
            

    def get_pflux(self,keyword,NSPEZ=2,displace=0,maxlen=100):
        from numpy import zeros,array
        ret=zeros((len(self.cases),maxlen,NSPEZ))
        for i in range(len(self.cases)):
            contents=self.get_STRATASUM(self.cases[i].path,linelist=True)
            retl = [l for l, elem in enumerate(contents) if 'SURFACE AREA' in elem]
            indices = [l for l, elem in enumerate(contents) if keyword in elem]
            for j in range(len(indices)):
                try: juse=next(idx for idx, value in enumerate(retl) if value > indices[j])-1
                except: juse=len(indices)-1

                for k in range(NSPEZ):
                    ind=indices[j]+3*(k+1)+2*k+displace
                    if len(contents[ind])>0:
                        if contents[ind][0] in ['D','H','C']:
                            try: ret[i,juse,k]=float(contents[ind].split(' ')[-1].strip())
                            except: pass
        
        return ret[:,:len(retl),:]


    def get_incidentEflux(self,NSPEZ=2):
        return self.get_Eflux('FLUX INCIDENT ON SURFACE:',NSPEZ)


    def get_total_reemittedEflux(self,NSPEZ=2,REC=False):
        from numpy import zeros
        ret=self.get_reemittedEflux('ATOMS')
        ret+=self.get_reemittedEflux('MOLECULES')
        ret+=self.get_reemittedEflux('BULK IONS')*False
    
        return ret

    def get_reemittedEflux(self,species,NSPEZ=2):
        from numpy import roll
        if species.upper()=='BULK IONS': d=1
        else: d=0
        ret= self.get_Eflux('FLUX RE-EMITTED FROM INCIDENT {}:'.format(species.upper()),NSPEZ=2,displace=d)
        if species.upper()=='MOLECULES':
            ret=roll(ret,1,axis=2)
            ret[:,:,0]=0
        return ret

    def get_Eflux(self,keyword,NSPEZ=2,displace=0,maxlen=100):
        from numpy import zeros,array
        ret=zeros((len(self.cases),maxlen,NSPEZ))
        for i in range(len(self.cases)):
            contents=self.get_STRATASUM(self.cases[i].path,linelist=True)
            retl = [l for l, elem in enumerate(contents) if 'SURFACE AREA' in elem]
            indices = [l for l, elem in enumerate(contents) if keyword in elem]
            for j in range(len(indices)):
                try: juse=next(idx for idx, value in enumerate(retl) if value > indices[j])-1
                except: juse=len(indices)-1

                for k in range(NSPEZ):
                    ind=2+indices[j]+3*(k+1)+2*k+displace
                    if len(contents[ind])>0:
                        if contents[ind][0] in ['D','H','C']:
                            try: ret[i,juse,k]=float(contents[ind].split(' ')[-1].strip())
                            except: pass
        
        return ret[:,:len(retl),:]





    def getNHIST(self):
        ''' Gets recombination from eirene out '''
        return self.get_text('NHIST= ',7,20)

    
    def getInfluxa(self):
        ''' Gets atomic eirene influx from surfaces '''
        #return self.get_text('NHIST= ',20,10)
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("PRFAAI = ")
            try: ret[i]+=float(contents[nhist+20:nhist+33].strip())
            except: pass
            nhist=contents.find("PRFAMI = ")
            try: ret[i]+=float(contents[nhist+20:nhist+33].strip())
            except: pass
        return ret

 
    def getISTRAfluxa(self,ISTRA=1):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("PRFAAI = ")
            try: ret[i]+=float(contents[nhist+20:nhist+33].strip())
            except: pass
            nhist=contents.find("PRFAMI = ")
            try: ret[i]+=float(contents[nhist+20:nhist+33].strip())
            except: pass
        return ret


    
    def getISTRAfluxm(self,ISTRA=1):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("PRFMMI = ")
            try: ret[i]=float(contents[nhist+20:nhist+33].strip())
            except: pass
        
        return ret
    
    def getInfluxm(self):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("PRFMMI = ")
            try: ret[i]=float(contents[nhist+20:nhist+33].strip())
            except: pass
        
        return ret
 
    def geteffluxm(self):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("POTMLI = ")
            try: ret[i]=float(contents[nhist+20:nhist+33].strip())
            except: pass
        
        return ret


    def getISTRAeffluxm(self,ISTRA=1):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("POTMLI = ")
            try: ret[i]=float(contents[nhist+20:nhist+33].strip())
            except: pass
        
        return ret

 
    def geteffluxa(self):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            contents=get_STRATASUM(self.cases[i].path)
            nhist=contents.find("POTATI = ")
            try: ret[i]=float(contents[nhist+20:nhist+33].strip())
            except: pass
        
        return ret



 
    def getISTRAeffluxa(self,ISTRA=1):
        ''' Gets atomic eirene influx from surfaces '''
        from numpy import zeros
        
        ret=zeros((len(self.cases),))
        for i in range(len(self.cases)):
            c=self.cases[i]
            with open(c.path+'/eirene.out','rt') as outfile:
                contents=outfile.read()
            contents=contents[contents.find("THIS IS STRATUM NUMBER ISTRA=     {}".format(ISTRA)):]
            nhist=contents.find("POTATI = ")
            try: ret[i]=float(contents[nhist+20:nhist+33].strip())
            except: pass
        
        return ret





    '''==========================================
    Get locations
    =========================================='''

    def row(self,var,row,s=0,supress=False):
        ''' Returns a list of 1D arrays containing var for row
        row(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        row:        Row along which to return variable

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

      
        return asarray([case.row(var,row,s=s) for case in self.cases])



    def ft(self,var,ft,s=0,supress=False):
        ''' Returns a list of 1D arrays containing var along flux-tube ft
        ft(var,ft,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ft:         Flux tube along which to return variable

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray
        return asarray([case.ft(var,ft,s=s) for case in self.cases])

    def index(self,var,idx,s=None,supress=False):
        ''' Returns a list values of var at location (ix,iy)
        index(var,ix,iy,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ix:         Polidal index to use
        iy:         Radial index to use

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        return asarray([case.get_ind(var,idx) for case in self.cases])


    '''==========================================
    Get min/max
    =========================================='''

    def row_min(self,var,row,s=None,supress=False):
        ''' Returns a list values of the minimum value of var the specified row
        row_min(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        row:        Row along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            if supress is False:
                print('WARNING! Requested row is guard cell row')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([min(case.get_row(var,row)) for case in self.cases])
        else: 
           return asarray([min(case.get_row(var,row,s=suse)) for case in self.cases])



    def row_max(self,var,row,s=None,supress=False):
        ''' Returns a list values of the maximum value of var the specified row
        row_max(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        row:        Row along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            if supress is False:
                print('WARNING! Requested row is guard cell row')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([max(case.get_row(var,row)) for case in self.cases])
        else: 
           return asarray([max(case.get_row(var,row,s=suse)) for case in self.cases])


        
    def ft_min(self,var,ft,s=None,supress=False):
        ''' Returns a list values of the minimum value of var the specified flux tube
        ft_min(var,ft,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ft:         Flux tube along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if ft in [0,-1,self.cases[0].nx()+1]:   
            if supress is False:
                print('WARNING! Requested flux tube is guard cell flux tube')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([min(case.get_ft(var,ft)) for case in self.cases])
        else: 
           return asarray([min(case.get_ft(var,ft,s=suse)) for case in self.cases])


    def ft_max(self,var,ft,s=None,supress=False):
        ''' Returns a list values of the maximum value of var the specified flux tube
        ft_max(var,ft,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ft:         Flux tube along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if ft in [0,-1,self.cases[0].nx()+1]:
            if supress is False:
                print('WARNING! Requested flux tube is guard cell flux tube')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([max(case.get_ft(var,ft)) for case in self.cases])
        else: 
           return asarray([max(case.get_ft(var,ft,s=suse)) for case in self.cases])



    '''==========================================
    Get locations
    =========================================='''

    def get_maxlocation(self, var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns a list of (R,Z) coordinates of max of var in xind/yind index interval 
        get_maxlocation(var,**keys)

        Variables:
        var:            String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        xind[=(1,-1)]   Poloidal index space to search
        yind[=(1,-1)]   Radial index space to search
        ''' 
        from numpy import asarray
        return asarray([case.get_maxlocation(var,s,xind,yind) for case in self.cases])

    def get_minlocation(self, var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns a list of (R,Z) coordinates of min of var in xind/yind index interval 
        get_maxlocation(var,**keys)

        Variables:
        var:            String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        xind[=(1,-1)]   Poloidal index space to search
        yind[=(1,-1)]   Radial index space to search
        ''' 
        from numpy import asarray
        return asarray([case.get_minlocation(var,s,xind,yind) for case in self.cases])
    
    '''==========================================
    Get case indices
    =========================================='''
    def get_closest(self,val,var,ix,iy,s=None):
        ''' Get index of  case with closest value to var at (ix,iy) 
        get_closest(vel,var,ix,iy,**keys)

        Variables:
        val:            Value to match
        var:            String of variable to be returned, e.g. 'bbb.ne'
        ix:             Poloidal location index
        iy:             Radial location index

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        '''
        ind=abs(self.index(var,ix,iy,s)-val).argmin()
        print('Closest value is: '+str(self.index(var,ix,iy,s)[ind]))
        return ind
       





def create_database(savename=None,sortind=(1,1,0),sortvar='ne',outpath='.',path='.',ret=True):
    ''' Creates a database
        Parameters:
            savename        If set, saves dict as pickle named 'savename'
            sortlocation    Location for sorting by te: 'core' or 'mp'
            outpath         Path to location where pickle is saved
            path            Path to parent directory
            subpath         Path to input.py within child directories of path: 
                            path/*/supath/input.py
            commands        List of commands to be executed before restoring solution
            variables       List of all variable names, including package, to be stored
            ret             Boolean whether to return dict or not
    '''
    from os import getcwd,chdir,remove,walk
    from os.path import abspath  
    from pickle import dump
    from eirplot import EIRRUN

    outpath=abspath(outpath)    # Get absolute path of out directory
    chdir(path)                 # Go to the parent directory
    parent=getcwd()               # Get path to parent
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    # Omit supporting file directories

    try:
        dirs.remove('ignore')
    except:
        pass
    try:
        dirs.remove('grid')
    except:
        pass
    try:
        dirs.remove('eirin')
    except:
        pass
    
    if len(dirs)==0:
        return 'No directories found! Aborting...'


    # Empty list to return
    retl=[]
    

    for child in dirs:      # Loop through all directories
        print('******************************')
        print('*** Directory: '+child+' ***')
        print('******************************')
        retl.append(EIRRUN(path+'/'+child))        

        chdir(parent)
    
    
    lst=SETUP(retl) 
    # Get the sep and xpt locations
    lst.sort_ind(sortvar,sortind)
    chdir(outpath)
    # Check if save requested
    if savename is not None:
        with open(savename,'wb') as f:
            dump(lst,f)
    if ret:
        return lst



def restore_database(name):    
    ''' Restores pickled case '''
    from pickle import load
        
    with open(name,'rb') as f:
        ret=load(f)

    try:
        with open(name,'rb') as f:
            ret=load(f)
    except:
        print('ERROR! Could not find file "'+name+'"! Aborting...')
        ret=[]
    return ret
