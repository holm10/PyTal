class SONNET():
    def __init__(self,path):
        from re import findall
        from numpy import zeros

        elements=[]
        with open(path) as sfile:
            for line in sfile:

                if line[0] in ['-','=']:
                    element={'idx':[],'nodes':[],'field':[],'ratio':[]}
                    for i in range(3):
                        try:
                            line=next(sfile)
                        except:
                            break
                        buff=findall('\(([^)]+)', line)
                        if i==0:
                            element['idx']=[int(x) for x in buff[0].split(',')]
                            for j in range(1,3):
                                element['nodes'].append([float(x) for x in buff[j].split(',')])
                        elif i==1:
                            element['ratio']=float(findall('\=([^(]+)', line)[0])
                            element['field']=[float(x) for x in buff[0].split(',')]
                        else:
                            for j in range(2):
                                element['nodes'].append([float(x) for x in buff[j].split(',')])
                    elements.append(element)


        self.dim=(elements[-2]['idx'][0]+1,elements[-2]['idx'][1]+1)

        self.x=zeros(self.dim+(4,))
        self.y=zeros(self.dim+(4,))
        self.bx=zeros(self.dim)
        self.by=zeros(self.dim)
        self.ratio=zeros(self.dim)

        for e in elements[:-1]:
            [ix,iy]=e['idx']
            for i in range(4):
                self.x[ix,iy,i]=e['nodes'][i][0]
                self.y[ix,iy,i]=e['nodes'][i][1]
            self.bx[ix,iy]=e['field'][0]
            self.by[ix,iy]=e['field'][1]
            self.ratio[ix,iy]=e['ratio']


    def plot_grid(self):
        from matplotlib.pyplot import figure

        f=figure()
        ax=f.add_subplot(111)

        for ix in range(self.dim[0]):
            for iy in range(self.dim[1]):
                lx,ly=[],[]
                for i in [0,1,3,2,0]:
                    lx.append(self.x[ix,iy,i])
                    ly.append(self.y[ix,iy,i])
                ax.plot(lx,ly,'k-',linewidth=0.5)

        f.show()



    def heatmap(self,Z,zrange=None,zaxis='lin'):
        from matplotlib.patches import Polygon
        from matplotlib.colors import Normalize,LogNorm
        from matplotlib.pyplot import get_cmap,colorbar,figure
        from numpy import shape,zeros,log10,floor,ceil
        from matplotlib.cm import ScalarMappable

        f=figure()
        ax=f.add_subplot(111)
        rm=self.x
        zm=self.y

        xlim,ylim= [rm.min(),rm.max()],[zm.min(),zm.max()]



        # Set heatmap limits if requested
        if zrange is None:
            Zmax=Z.max()
            Zmin=Z.min()
        else:
            if isinstance(zrange[0],str): # Choosing to limit only one boundary
                if zrange[0]=="min":
                    Zmin,Zmax=zrange[1],Z.max()
                elif zrange[0]=="max":
                    Zmin,Zmax=Z.min(),zrange[1]
            else:
                Zmin=zrange[0]
                Zmax=zrange[1]

        if zaxis=="lin":
            Zcol=(Z-Zmin)/(Zmax-Zmin)
        elif zaxis=="log":
            Zcol=((log10(Z)-floor(log10(Zmin)))/(floor(log10(Zmax))-floor(log10(Zmin))))
        else:
            print("Only valid zaxis options are 'lin' and 'log'!")
            return
            
        # Set colormap
        cmap=get_cmap('magma')
        

        # Plot heatmap using polygons
        for ix in range(self.dim[0]):    
            for iy in range(self.dim[1]):    
                xy=zeros((4,2))
                # Create polygon for each grid cell
                lx,ly=[],[]
                for i in [0,1,3,2]:
                    lx.append(self.x[ix,iy,i])
                    ly.append(self.y[ix,iy,i])
                xy[:,0]=lx
                xy[:,1]=ly
                # Set color based on Z-value
    #            col=cmap((Z[i,j]-Zmin)/Zmax)
                col=cmap(Zcol[ix,iy])
                # Plot the patch
                ax.add_patch(Polygon(xy,closed=True,facecolor=col,edgecolor=col))
                

        for ix in range(self.dim[0]):    
            for iy in range(self.dim[1]): 
                lx,ly=[],[]
                for i in [0,1,3,2,0]:
                    lx.append(self.x[ix,iy,i])
                    ly.append(self.y[ix,iy,i])
                ax.plot(lx,ly,'-',linewidth=0.5,color='grey')
        

        # Set colorbar if requested
        if zaxis=="lin":
            norm = Normalize(vmin=Zmin,vmax=Zmax)
        elif zaxis=="log":
            norm = Normalize(vmin=floor(log10(Zmin)),vmax=floor(log10(Zmax)))
            norm = LogNorm(vmin=Zmin,vmax=Zmax)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        # See if colorbar needs to be extended
        if zrange is None:
            extend="neither"
        elif zrange[0]>Z.min() and zrange[1]<Z.max():
            extend="both"
        elif zrange[0]>Z.min():
            extend="min"
        elif zrange[1]<Z.max():
            extend="max"
        else:
            extend="neither"
        cbar=colorbar(sm,ax=ax,extend=extend)
        

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    

    def plot_btor(self):
        self.heatmap(self.bx)

    def plot_bpol(self):
        self.heatmap(self.by)

    def plot_ratio(self):
        self.heatmap(self.ratio)


