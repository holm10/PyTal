# TODO: store data to memory rather than read from file
# Script created to read EIRENE grids and parameters from fort.3[3-5] files
# Created from scratch by holma2 on 191219
# Changelog
# 191219 -  Created
from matplotlib.pyplot import ion
ion()


class TALLY:
    def __init__(self, fname, path='.', altname='', ISTRA=0, nx=0, ny=0, ueproj=None):
        ''' Reads the data from tally, returns list with values '''
        from numpy import zeros, transpose,array

        self.fname = fname # Tally name
        self.altname = altname # Variable shorthand name

        # Read the file line-by line
        intal = fname.split("_")[0]=='intal' # Flag to mark if we are in the correct strata
        strata = True*intal # Flag to mark if we are in the correct strata
        celldata = False # Flag to mark whether cell data is being read
        speciesline = False # Flag whether next line is species handles
        unitline = False # Flag whether next line is the variable unit
        totline = False
        meanline = False
        valbuff = []
        with open('{}/{}'.format(path,fname),'r') as f:
            # Loop through the file looking for the requested strata:
            for l in f:
                # Strata detected
                if l[:5] == 'ISTRA':
                    if int(l[8:]) == ISTRA:
                        strata = True
                    else: 
                        strata = False                   
                elif l[0] == '=':
                    if strata and (not celldata):
                        celldata = True
                    else:
                        celldata = False
                elif strata and celldata:
                    valbuff.append([int(l.split()[0]),[float(x) for x in l.split()[1:]]])
                elif strata and (not celldata):
                    if l[:6] == 'NCELLS':
                        self.N = int(l.split()[-1]) # Number of cells
                    elif l[:8] == 'NSPECIES':
                        self.Ns = int(l.split()[-1]) # Number of species
                        if self.Ns == 0:
                            return
                    elif l[:7] == 'SPECIES':
                        speciesline = True
                    elif speciesline:
                        self.species = l.split() # List with species handles
                        speciesline = False
                    elif l[:5] == 'UNITS':
                        unitline = True
                    elif unitline:
                        self.unit = l.split()[0] # String with EIRENE unit
                        unitline = False
                    elif l[:5] == 'TOTAL':
                        totline = True
                    elif totline:
                        self.volsum = array([float(x) for x in l.split()]) # Sum of product of cell value and volume
                        totline = False
                    elif l[:4] == 'MEAN':
                        meanline = True
                    elif meanline:
                        self.mean = array([float(x) for x in l.split()]) # Mean
                        meanline = False

        # TODO: Figure out if index 1 is the GC or first real cell (i.e. Fotran or Python indexing.. Assume Fortran)
        # Create numpy array with data
        self.data = {}
        # Create EIRENE data array
        eirbuff = zeros((self.N, self.Ns))
        for i in valbuff[:len(valbuff)-2*intal]:
            eirbuff[i[0]] = i[1]

        self.data['eir'] = eirbuff.transpose()


        # Create UEDGE data array
        uebuff = zeros((nx+2, ny+2, self.Ns, 2))
        
        for s in range(self.Ns): # Loop through each species
            for i in range(1,self.N-1):
                uebuff[ueproj[i,0], ueproj[i,1], :, ueproj[i,2]] = self.data['eir'][:,i]

        self.data['ue'] = transpose(uebuff, (2,0,1,3))


        
        
            
                    
        
class NODES:
    def __init__(self,path='.',npco=True):
        # Read the node data
        tria=False
        nodex,nodey=[],[]
        self.nodes=[]
        if npco is True: 
            nodefile='triang_new.npco_char'
        else:
            nodefile='fort.33'

        with open('{}/{}'.format(path,nodefile),'r') as f:
            Nnodes=int(f.readline().strip())
            for i in range(Nnodes):
                # Read line-wise
                temp=[float(x) for x in f.readline().split()]
                # If there are less than Nnodes lines, abort
                # This is the case for tria-style files
                if temp==[]: 
                    for i in range(Nnodes):
                        self.nodes.append([nodex[i],nodey[i]])
                    return
            
                # Check data format on first line and choose read algorithm
                if len(temp)==4 and i==0: 
                    tria=True

                # Read node-by-node
                if tria is False:
                    self.nodes.append(temp[1:])
                # Read x-and-y's
                else:
                    for x in temp:
                        if len(nodex)<=Nnodes-1:
                            nodex.append(x)
                        else:
                            nodey.append(x) 
        
        if tria is True: 
            for i in range(Nnodes):
                self.nodes.append([nodex[i],nodey[i]])
        



    def get(self,N):
        ''' Get the coordinates of the Nth fortran index node '''
        return self.nodes[N-1]

class TRIMAP():
    def __init__(self,path='.'):
        # Read the trimap-file
        self.tria=[]
        self.data=[]
        try:
            with open('{}/eirene.trimap'.format(path),'r') as f:
                self.N=int(f.readline().strip())
                for i in range(self.N):
                    line=[float(x) for x in f.readline().split()[1:]]
                    self.data.append(line[:3])
                    self.tria.append([line[3:5],line[5:7],line[7:9]])
        finally:
            f.close

    def get(self,N):
        return self.tria[N-1]
    
    def getdata(self,N):
        return self.data[N-1]

    def getN(self):
        return self.N

    def plottria(self,linewidth=0.5,color='k',aspect=True):
        ''' Plots the trimap grid '''
        from matplotlib.pyplot import figure
        from numpy import array

        fig=figure()
        ax=fig.add_subplot(111)
        for p in self.tria:
            temp=array(p+[p[0]])
            ax.plot(temp[:,0],temp[:,1],linewidth=linewidth,color=color)

        if aspect is True:
            ax.set_aspect('equal')
        return fig
            
        


class TRIANGLES():
    def __init__(self,path='.'):
        # Read the triangle nodes
        self.tria=[]
        try:
            with open('{}/fort.34'.format(path),'r') as f:
                self.N=int(f.readline().strip())
                for i in range(self.N):
                    line=[int(x) for x in f.readline().split()[1:]]
                    if len(line)==11:
                        self.tria.append(line[:3])
                    else:
                        self.tria.append(line)
        finally:
            f.close

    def getN(self):
        return self.N

    def get(self,N):
        return self.tria[N-1]

class NEIGHBORS():
    def __init__(self,path='.'):
        from numpy import array, amax

        self.neigh=[]
        self.ind=[]
        with open('{}/fort.35'.format(path),'r') as f:
            N=int(f.readline().strip())
            for i in range(N):
                line=[int(x) for x in f.readline().split()]
                temp1=[]
                for j in range(3):
                    temp2=[]
                    for k in range(3):
                        temp2.append(line[1+j*3+k])
                    temp1.append(temp2)
                self.neigh.append(temp1) 
                self.ind.append(line[-2:])
        self.ind = array(self.ind) 

    def get(self,N):
        ''' Get the neighbors of the Nth fortran index node '''
        return self.neigh[N-1]
    
    def get_ind(self,N):
        ''' Get the indices of the Nth fortran index node '''
        return self.ind[N-1]

    def get_N(self,N):
        ''' Returns the array of all triangle indices and vertices that neighbor surface with index N '''
        ret=[]
        for t in range(len(self.neigh)):
            for v in range(3):
                if self.neigh[t][v][-1]==N:
                    ret.append([t+1,v])
        return ret

    
    def get_dimensions(self):
        from numpy import amax
        return amax(self.ind, axis=0)    
        


class FORT30:
    ''' Class plotting fort.30 data '''
    def __init__(self,path='.'):
        self.cells=[]
        try:
            with open('{}/fort.30'.format(path),'r') as f:
                isx=False
                start=False
                x,y=[],[]
                for line in f:
                     if start is False:
                            if len(line.split())==4: 
                                x.append([float(i) for i in line.split()])
                                start=True
                     else:
                        if isx is True:
                            x.append([float(i) for i in line.split()])
                            isx=False
                        elif isx is False:
                            y.append([float(i) for i in line.split()])
                            isx=True
            for i in range(len(x)):
                self.cells.append([x[i],y[i]])
        finally:
            f.close

    def plot_cell(self, N, ax, color='k',linewidth=0.5):
        ''' Plots the Nth cell on ax '''
        ax.plot(self.cells[N][0]+[self.cells[N][0][0]],self.cells[N][1]+[self.cells[N][1][0]],color=color,linewidth=linewidth)
    
    def plot_grid(self):
        ''' Plots the plasma grid '''
        from matplotlib.pyplot import figure
        fig=figure()
        ax=fig.add_subplot(111)
        for i in range(len(self.cells)):
            self.plot_cell(i,ax)

class EIRRUN:
    def __init__(self,path='.',ISTRA=0):
        from numpy import linspace,reshape,transpose,pad,zeros,sum
        from os.path import abspath
        self.path=abspath(path)


        self.nodes=NODES(path=path,npco=False)
        self.triangles=TRIANGLES(path=path) 
        self.neigh=NEIGHBORS(path=path) 
        nx, ny = self.neigh.get_dimensions()


        self.nx=nx
        self.ny=ny
        self.N=self.triangles.getN()
        # Create a set of mapping arrays between UEDGE and EIRENE
        self.UEarr = pad( transpose( reshape( 
                    linspace(1,self.N,self.N).astype(int), (ny,nx,2)),
                    axes=(1,0,2)), ((1,1),(1,1),(0,0)), constant_values=0)
        # The UEDGE array is a 0:nx+1 by 0:ny+1 by 2 array.
        # The 3rd dimension contains the index of each of the two triangles
        # that the UEDGE grid is split into. 
        self.EIRarr=zeros((self.N+1,3))
        for i in range(nx+1):
            for j in range(ny+1):
                for k in range(2):
                    self.EIRarr[self.UEarr[i,j,k]]=[i,j,k]
        self.EIRarr=self.EIRarr.astype(int)
        # The EIRENE array is an ordered array according to the EIRENE 
        # triangle indices, where each index contains the coordinates
        # to the cell in the UEDGE grid (index 0 is a placeholder to
        # recreate the fortran indexing. The third index indicates whether
        # it is the first or second traingle in the sell


        self.data={}
        for var, tally in self.tallies().items():
            try:
                self.data[tally] = TALLY(tally, path=path, altname=var, ISTRA=ISTRA, nx=nx, ny=ny, ueproj=self.EIRarr)
            except:
                self.data[tally] = None

    def get_bounds(self,mult=1e-2):
        ''' Returns the min and max bounds in each direction

            Returns: [(xmin, xmax), (ymin, ymax)]
        '''
        from numpy import array
        n=array(self.nodes.nodes)*mult
        # TODO: Subtract min x-bound on flag: for simple geos
        return [(n[:,0].min(),n[:,0].max()),(n[:,1].min(),n[:,1].max())]

    def get_tria(self,N):
        ''' Gets the nodes of the triangle with frotran index N '''
        from numpy import array
        nodes=[]
        for i in self.triangles.get(N):
            nodes.append(self.nodes.get(i))
        return array(nodes)

    def get_vertex(self,N,V):
        ''' Gets vertex V of triangle with Fortran index N '''
        t=self.get_tria(N)
        return t[V:V+2]

    def plot_tria(self, N, ax,label=True,V=True,color='b',scale=1,zerox=0):
        ''' Plots the triangle with Nth Fortran index on axes '''
        from matplotlib.pyplot import text
        nodes=self.get_tria(N)
        x,y=[],[]
        for i in [0,1,2,0]:
            x.append(nodes[i][0]*scale-zerox)
            y.append(nodes[i][1]*scale)
        xc=(1/3)*sum(x[1:])        
        yc=(1/3)*sum(y[1:])        
        if V is True:
            ax.plot(x,y,'k-',linewidth=0.5,color=color,alpha=0.3)
        else:
            ax.plot(x[V:V+2],y[V:V+2],color=color,linewidth=5,alpha=0.3)
        if label is True and V is True:
            text(xc,yc,str(N),fontsize=4)

    def plot_bound(self,N,ax,color='b',linewidth=2):
        ''' Plots the boundary with index N on ax '''
        x,y=[],[]
        Vlist=self.neigh.get_N(N)
        for V in Vlist:
            self.plot_tria(V[0],ax,V=V[1],color=color)
            '''
            v=self.get_vertex(V[0],V[1])
            for i in v:
                x.append(i[0])
                y.append(i[1])

            ax.plot(x,y,color=color,linewidth=linewidth)      
            '''  

    def plot_npco(self):
        from matplotlib.pyplot import figure
        
        f=figure()
        ax=f.add_subplot(111)

        with open('{}/triang_new.npco_char'.format(self.path),'r') as f:
            for l in f:
                try:
                    [_,x,y]=l.strip().split('  ')
                    ax.plot(float(x),float(y),'k.')
                except:
                    pass

    def getR(self,ft,zerox=True):  
        from numpy import zeros,mean
        [xlim,_]=self.get_bounds()
        ret=zeros((self.nx+1,))
        for i in range(self.nx+1):
            x=self.UEarr[i,ft]
            for j in x:
                ret[i]+=mean(self.get_tria(j)[:,0])-zerox*xlim[0]*100
            ret[i]/=2
        return ret[1:]*1e-2

    def getZarr(self):
        from numpy import array
        ret=[]
        for i in range(self.nx):
            ret.append(self.getZ(i))
        return array(ret)

    def getZ(self,row):
        from numpy import zeros,mean

        ret=zeros((self.ny+1,))
        for i in range(self.ny+1):
            x=self.UEarr[row,i]
            for j in x:
                ret[i]+=mean(self.get_tria(j)[:,1])
            ret[i]/=2
        return ret[1:]*1e-2


    def get(self, var, code='eir', processing='raw', s=None):
        from numpy import zeros,sum,mean, transpose
        
        try:
            ret = self.data[self.tallies()[var]] # Try to access the data by the aternative handle
        except:
            ret = self.data[var] # If not, assume the variable to the tally name

        if processing == 'volsum': # EIRENE tally volumetric sum
            return ret.volsum
        elif processing == 'meansum': # EIRENE tally mean
            return ret.mean
        elif (ret.Ns !=1) and (s is not None): # Else, pick the right species
            ret = ret.data[code][s]
        elif (ret.Ns != 1) and (s is None):
            if code == 'ue':
                ret = transpose(ret.data[code], (1,2,3,0))
            elif code == 'eir':
                ret = transpose(ret.data[code], (1,0))
        else:
            ret = ret.data[code][0]
            
        # Check if UEDGE mean is requested
        if (processing == 'mean') and (code == 'ue'):
            ret = sum(ret, axis=2)
        else:
            return ret

        # TODO: get meansum - ret.mean
        # TODO: get volsum - ret.colsum
        # TODO: get mean - ret cell mean
        # TODO: get weighted - ret cellval*cellvol/totvol



    def getN(self,var,s=0):
        from numpy import array

        # Conversion from energy density to temperature
        if var[0]=='t':
            Ze=array(self.vars['e'+var[1:]])
            Zn=array(self.vars['n'+var[1:]])*1e-6
            Zn[Zn<=0]=1e-8
            Ze[Ze<=0]=1e-8
            Z=Ze/Zn
        elif var[0]=='p':
            Z=array(self.vars['e'+var[1:]])*1e6*1.602e-19
        else:
            Z=self.vars[var]

        return Z

        

    def ft(self,var,ft,typ='mean',ISTRA=0):
        from numpy import zeros,sum,mean
        ret=zeros((self.nx+1,2)) 

        # Conversion from energy density to temperature
        if var[0]=='t':
            Ze=array(self.vars['e'+var[1:]])
            Zn=array(self.vars['n'+var[1:]])*1e-6
            Zn[Zn<=0]=1e-8
            Ze[Ze<=0]=1e-8
            Z=Ze/Zn
        elif var[0]=='p':
            Z=array(self.vars['e'+var[1:]])*1e6*1.602e-19
        else:
            Z=self.vars[var]
        val=Z
        
        if typ=='weighted':
            vol=self.vars['vol']

        for i in range(self.nx+1):
            tvol=0
            for k in range(2):
                ret[i,k]=val[self.UEarr[i,ft,k]]
                try: 
                    ret[i,k]=val[self.UEarr[i,ft,k]]*vol[self.UEarr[i,ft,k]]
                    tvol+=vol[self.UEarr[i,ft,k]] 
                except:     
                    pass
            if typ=='weighted':
                ret[i,0]=sum(ret[i])/tvol
        
        if typ=='weighted':
            return ret[1:,0]
        elif typ=='mean':
            return mean(ret[1:],axis=1)
        elif typ=='sum':
            return sum(ret[1:],axis=1)
        else:
            print('Type {} not recognized! Aborting...'.format(typ))
            return


    def row(self,var,row,typ='mean',s=0):
        from numpy import zeros,sum,mean
        ret=zeros((self.nx+1,2)) 
        
        # Conversion from energy density to temperature
        if var[0]=='t':
            Ze=array(self.vars['e'+var[1:]])
            Zn=array(self.vars['n'+var[1:]])*1e-6
            Zn[Zn<=0]=1e-8
            Ze[Ze<=0]=1e-8
            Z=Ze/Zn
        elif var[0]=='p':
            Z=array(self.vars['e'+var[1:]])*1e6*1.602e-19
        else:
            Z=self.vars[var]
        val=Z
        

        if typ=='weighted':
            vol=self.vars['vol']

        for i in range(self.ny+1):
            tvol=0
            for k in range(2):
                ret[i,k]=val[self.UEarr[row,i,k]]
                try: 
                    ret[i,k]=val[self.UEarr[row,i,k]]*vol[self.UEarr[row,i,k]]
                    tvol+=vol[self.UEarr[row,i,k]] 
                except:     
                    pass
            if typ=='weighted':
                ret[i,0]=sum(ret[i])/tvol
        
        if typ=='weighted':
            return ret[1:,0]
        elif typ=='mean':
            return mean(ret[1:],axis=1)
        elif typ=='sum':
            return sum(ret[1:],axis=1)
        else:
            print('Type {} not recognized! Aborting...'.format(typ))
            return



    def plot_ft(self,var,ft,ax=None,pl='plot',zerox=0):
        from matplotlib.pyplot import subplots
        from numpy import sum,zeros
        if ax is None:
            fig, ax = subplots()
        sx = []
        for x in range(self.ny*ft + 1, self.ny*(ft+2) + 1, 2):
            sx.append(sum(self.get_tria(x), axis=0)[0]/3 - zerox)
        ax.plot(sx, self.ft(var,ft),'.-')
        try:
            fig.show()
            return fig
        except:
            pass

    def plot_grid(self,label=True,bounds=[],linewidth=2,showgrid=True,aspect=False,zerox=True):
        ''' Plots the EIRENE grid '''
        from matplotlib.pyplot import figure

        [xlim,ylim]=self.get_bounds()
        # TODO: create functionality for plotting other ISTRA
        colors=['b','darkgreen','red','teal','gold','c','m']
        fig=figure()
        ax=fig.add_subplot(111)
        if showgrid is True:
            for i in range(1,self.triangles.getN()+1):
                self.plot_tria(i,ax,label=label,zerox=xlim[0]*zerox*100)
        for i in range(len(bounds)):
            self.plot_bound(bounds[i],ax,color=colors[i],linewidth=linewidth)
        if aspect is True:
            ax.set_aspect('equal')
        return fig



    def UEcomp( self,var,zUE,zrange=False,zaxis='lin',cmap='bwr',grid=False,ISTRA=0,
                title=None,units='',zerox=True,NOM='U-E',DENOM='U',scal=1,UEscal=1,Escal=1,
                xlabel='Poloidal position [m]',ylabel='Radial position [m]',ax=None):

        from numpy import array,zeros
        

        # Conversion from energy density to temperature
        if var[0]=='t':
            ZE=array(self.vars['e'+var[1:]])
            Zn=array(self.vars['n'+var[1:]])*1e-6
            Zn[Zn<=0]=1e-8
            ZE/=Zn
        elif var[0]=='p':
            ZE=array(self.vars['e'+var[1:]])*1e6*1.602e-19
        else:
            ZE=array(self.vars[var])

        ZU=zeros(ZE.shape)

        for i in range(len(ZE)):
            ZU[i]=zUE[tuple(self.EIRarr[i])]


        ZE*=Escal
        ZU*=UEscal
            
        if NOM=='U-E':
            Z=ZU-ZE
        elif NOM=='E-U':
            Z=ZE-ZU
        elif NOM=='E':
            Z=ZE
        elif NOM=='U':
            Z=ZU

        if DENOM=='U':
            Z/=ZU
        elif DENOM=='E':
            Z/=ZE
        else:
            Z/=DENOM

        Z*=scal
            
        self.plot_var(Z,zrange=zrange,zaxis=zaxis,cmap=cmap,grid=grid,units=units,title=title,zerox=zerox,xlabel=xlabel,ylabel=ylabel,ax=ax)

        

        


    def heatmap(self,var,zrange=False,zaxis='lin',cmap='magma',grid=False,ISTRA=0,
                title=None,units='',zerox=True,ax=None,
                xlabel='Poloidal position [m]',ylabel='Radial position [m]'):
        from numpy import array

        # Conversion from energy density to temperature
        if var[0]=='t':
            Ze=array(self.vars['e'+var[1:]])
            Zn=array(self.vars['n'+var[1:]])*1e-6
            Zn[Zn<=0]=1e-8
            Ze[Ze<=0]=1e-8
            Z=Ze/Zn
        elif var[0]=='p':
            Z=array(self.vars['e'+var[1:]])*1e6*1.602e-19
        else:
            Z=self.vars[var]

        self.plot_var(Z,zrange=zrange,zaxis=zaxis,cmap=cmap,grid=grid,units=units,title=title,zerox=zerox,ax=ax,xlabel=xlabel,ylabel=ylabel)
        


        

    def plot_var(   self,Z,zrange=False,zaxis='lin',cmap='magma',grid=False,
                    title=None,units='',zerox=True,ax=None,
                    xlabel='Poloidal position [m]',ylabel='Radial position [m]'):

        ''' Plots the values of var on the EIRENE grid '''
        from numpy import array,transpose,log10,floor,ceil
        from matplotlib.pyplot import get_cmap,colorbar,subplots
        from matplotlib.patches import Polygon
        from matplotlib.colors import Normalize,LogNorm
        from matplotlib.cm import ScalarMappable


        Z=Z[1:]

        if str(type(ax))!="<class 'matplotlib.axes._subplots.AxesSubplot'>": f,ax=subplots()

        [xlim,ylim]=self.get_bounds()

        
        if zaxis=='log':
            Z[Z<=0]=1e-8

        # Set heatmap limits if requested
        if zrange is False:
            Zmax=Z.max()
            Zmin=Z.min()
        else:
            if isinstance(zrange[0],str): # Choosing to limit only one boundary
                if zrange[0]=="min":
                    Zmin,Zmax=zrange[1],Z.max()
                elif zrange[0]=="max":
                    Zmin,Zmax=Z.min(),zrange[1]
                else:
                    print("zrange can only be defined as 'min'/'max'! Terminating...")
                    exit(0)
            else:
                Zmin=zrange[0]
                Zmax=zrange[1]


        if zaxis=="lin":
            Z=(Z-Zmin)/(Zmax-Zmin)
        elif zaxis=="log":
            Z=((log10(Z)-floor(log10(Zmin)))/(floor(log10(Zmax))-floor(log10(Zmin))))
        else:
            print("Only valid zaxis options are 'lin' and 'log'!")
            return


        ax.set_xlim((xlim)-xlim[0]*zerox)
        ax.set_ylim((ylim))
        cmap=get_cmap(cmap)
        for i in range(len(Z)):
            try:
                xy=array(self.get_tria(i+1))*1e-2
                xy[:,0]-=xlim[0]*zerox
                ax.add_patch(Polygon(xy,closed=True,facecolor=cmap(Z[i]),edgecolor=cmap(Z[i])))
            except:
                pass


        if zaxis=="lin":
            norm = Normalize(vmin=Zmin,vmax=Zmax)
        elif zaxis=="log":
            norm = Normalize(vmin=floor(log10(Zmin)),vmax=floor(log10(Zmax)))
            norm = LogNorm(vmin=Zmin,vmax=Zmax)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    
        # See if colorbar needs to be extended
        if zrange is False:
            extend="neither"
        elif Zmin>Z.min() and Zmax<Z.max():
            extend="both"
        elif Zmin>Z.min():
            extend="min"
        elif Zmax<Z.max():
            extend="max"
        else:
            extend="neither"
        cbar=colorbar(sm,ax=ax,extend=extend)
        
#        if units is not None:

        cbar.set_label(units)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    
        if grid is True:
            for i in range(1,self.triangles.getN()+1):
                self.plot_tria(i,ax,label=False,color='grey',scale=1e-2,zerox=zerox*xlim[0])


        return Z
                    
        
    def tallies(self):
        return {
                    'te':'intal_1',
                    'ti':'intal_2',
                    'ne':'intal_3',
                    'ni':'intal_4',
                    'uup':'intal_5',
                    'vy':'intal_6',
                    'vzd':'intal_7',
                    'bx':'intal_8',
                    'by':'intal_9',
                    'bz':'intal_10',
                    'b':'intal_11',
                    'I12':'intal_12',
                    'eikind':'intal_13',
                    'vol':'intal_14',
                    'I15':'intal_15',
                    'I16':'intal_16',
                    'I17':'intal_17',
                    'I18':'intal_18',
                    'I19':'intal_19',
                    'I20':'intal_20',
                    'I21':'intal_21',
                    'I22':'intal_22', 
                    'na':'outtal_1',
                    'nm':'outtal_2',
                    'nmi':'outtal_3',
                    'ng':'outtal_4',
                    'ea':'outtal_5',
                    'em':'outtal_6',
                    'emi':'outtal_7',
                    'eg':'outtal_8',
                    'apsore':'outtal_9',
                    'apsora':'outtal_10',
                    'apsorm':'outtal_11',
                    'apsormi':'outtal_12',
                    'apsorg':'outtal_13',
                    'apsori':'outtal_14',
                    'mpsore':'outtal_15',
                    'mpsora':'outtal_16',
                    'mpsorm':'outtal_17',
                    'mpsormi':'outtal_18',
                    'mpsorg':'outtal_19',
                    'mpsori':'outtal_20',
                    'ipsore':'outtal_21',
                    'ipsora':'outtal_22',
                    'ipsorm':'outtal_23',
                    'ipsormi':'outtal_24',
                    'ipsormg':'outtal_25',
                    'ipsori':'outtal_26',
                    'gpsore':'outtal_27',
                    'gpsora':'outtal_28',
                    'gpsorm':'outtal_29',
                    'gpsormi':'outtal_30',
                    'gpsormg':'outtal_31',
                    'gpsori':'outtal_32',
                    'aesore':'outtal_33',
                    'aesora':'outtal_34',
                    'aesorm':'outtal_35',
                    'aesormi':'outtal_36',
                    'aesorg':'outtal_37',
                    'aesori':'outtal_38',
                    'mesore':'outtal_39',
                    'mesora':'outtal_40',
                    'mesorm':'outtal_41',
                    'mesormi':'outtal_42',
                    'mesorg':'outtal_43',
                    'mesori':'outtal_44',
                    'iesore':'outtal_45',
                    'iesora':'outtal_46',
                    'iesorm':'outtal_47',
                    'iesormi':'outtal_48',
                    'iesormg':'outtal_49',
                    'iesori':'outtal_50',
                    'gesore':'outtal_51',
                    'gesora':'outtal_52',
                    'gesorm':'outtal_53',
                    'gesormi':'outtal_54',
                    'gesormg':'outtal_55',
                    'gesori':'outtal_56',
                    'O57':'outtal_57',
                    'O58':'outtal_58',
                    'O59':'outtal_59',
                    'O60':'outtal_60',
                    'O61':'outtal_61',
                    'O62':'outtal_62',
                    'O63':'outtal_63',
                    'O64':'outtal_64',
                    'O65':'outtal_65',
                    'O66':'outtal_66',
                    'O67':'outtal_67',
                    'O68':'outtal_68',
                    'O69':'outtal_69',
                    'O70':'outtal_70',
                    'O71':'outtal_71',
                    'O72':'outtal_72',
                    'O73':'outtal_73',
                    'O74':'outtal_74',
                    'O75':'outtal_75',
                    'O76':'outtal_76',
                    'O77':'outtal_77',
                    'O78':'outtal_78',
                    'O79':'outtal_79',
                    'O80':'outtal_80',
                    'O81':'outtal_81',
                    'O82':'outtal_82',
                    'O83':'outtal_83',
                    'O84':'outtal_84',
                    'O85':'outtal_85',
                    'O86':'outtal_86',
                    'O87':'outtal_87',
                    'O88':'outtal_88',
                    'O89':'outtal_89',
                    'O90':'outtal_90',
                    'O91':'outtal_91',
                    'O92':'outtal_92',
                    'O93':'outtal_93',
                    'O94':'outtal_94',
                    'O95':'outtal_95',
                    'O96':'outtal_96',
                    'O97':'outtal_97',
                    'O98':'outtal_98',
                    'O99':'outtal_99',
                    'O100':'outtal_100',
                }



