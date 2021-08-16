import numpy as np
import params as par
from itertools import product
import matplotlib.pyplot as plt
from matplotlib  import cm
from params import *
import json
from scipy.optimize import minimize
from scipy.fftpack import fft,fft2,fftshift
from scipy.ndimage import zoom
import time


def generate_Hk(kpoint, Ggrid, hdim, gpot, mpot, emass_field_inv, gdim, include_emass_correction):
     
    for i in range(len(gpot[:,0])):
        if np.linalg.norm(gpot[i,:]) <= par.eps:
           emass_inv_0 = np.real(emass_field_inv[i]) 
           break

    Hk = np.zeros((hdim, hdim), dtype = 'cfloat')
    for i,gi in enumerate(Ggrid):
        for j,gj in enumerate(Ggrid):

            if i == j:
               Hk[i,i] = (np.linalg.norm(kpoint + gi)**2) * par.hb2m0 * emass_inv_0        
            else:
               gd = gi - gj 

               for k,gp in enumerate(gpot):

                   if np.linalg.norm(gd - gp) <= par.eps:
                      Hk[i,j] += mpot[k] 
                      if include_emass_correction:
                         Hk[i,j] += (np.linalg.norm(gd)**2) * par.hb2m0 * emass_field_inv[k]
                      break

    return Hk        
 

def generate_G_grid(gdim):

    Ggrid = []
    gnorm = par.gdim*np.linalg.norm(par.G1m)

    for i in range(3*par.gdim):
        for j in range(3*par.gdim):
         
          ishift = i - int(3*par.gdim/2) 
          jshift = j - int(3*par.gdim/2) 
          GGl = par.G1m * ishift + par.G2m * jshift

          if np.linalg.norm(GGl) <= gnorm + par.eps:
             Ggrid.append(GGl)
     
    Ggrid = np.asarray(Ggrid)

    return Ggrid, len(Ggrid)        
             
def plot_moire(mat, name, cell, dtype = 'field', colormap= cm.jet, repeat = 2):
    m1 = len(mat)
    m2 = repeat*m1
    X=np.zeros([m2+1,m2+1])
    Y=np.zeros([m2+1,m2+1])
    Z=np.zeros([m2+1,m2+1])
    for i,j in product(range(m2+1),range(m2+1)):
        if(cell=='unit'):
            r=par.R1*i/m1+par.R2*j/m1
        elif(cell=='moire'):
            r=par.R1m*i/m1+par.R2m*j/m1
        X[i][j]=r[0]                    # x coordinate of displacement
        Y[i][j]=r[1]                    # y coordinate of displacement
        Z[i][j]=np.real(mat[np.mod(i,m1),np.mod(j,m1)])
        
    if dtype == 'emass':
       Z = 1. / Z
    elif dtype == 'field':
       Z = Z -  np.min(Z) 
             
    Z = zoom(Z, 1)
    Y = zoom(Y, 1)
    X = zoom(X, 1)
    plt.figure(figsize=[20,14])
    plt.contourf(X,Y,Z,30,cmap=colormap)   

    plt.plot([X[m1][m1]],[Y[m1][m1]],c='white')
    plt.axis("equal")
    cbar = plt.colorbar()
    if dtype == 'ldos':
       cbar.set_ticks([])
    plt.axis('off')
    plt.savefig(name)
    plt.close()
       
  
def plotterk(PATH, field, ftype, frelax):

    X, Y, Z = [], [], []
    par.ngrid2 = len(field)
    
    for ishift in range(-par.ngrid2, par.ngrid2):
        for jshift in range(-par.ngrid2, par.ngrid2):
          
          GG =  par.G1m * ishift  + par.G2m * jshift

          if np.linalg.norm(GG) <= (par.Gcut + par.eps) * np.linalg.norm(par.G1m):

             X.append(par.G1m[0] * ishift + par.G2m[0] * jshift)
             Y.append(par.G1m[1] * ishift + par.G2m[1] * jshift) 
             if ishift  == 0 and jshift == 0 and ftype == 'GSFE':
                Z.append(0.)              
             else:
                Z.append(field[np.mod(ishift, par.ngrid2),np.mod(jshift,par.ngrid2)])            

    fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
    for i, ax in enumerate(axes):

        ax.set_xlabel("kx",fontsize=12)
        ax.set_ylabel("ky",fontsize=12)
        ax.grid(True,linestyle='-',color='0.5', linewidth=0.3 )

        # scatter with colormap mapping to z value
        if i == 0:
           Zplot = np.abs(Z)
           plt.title('Modulus')
           cmap = cm.plasma
        else:
           Zplot = np.degrees(np.angle(Z))  
           cmap = cm.jet
           plt.title('Phase') 
       
        ax.set(adjustable='box', aspect='equal')
        im = ax.scatter(X,Y,s=200,c=Zplot, marker = 'o', cmap = cmap );
      
        plt.colorbar(im, ax=ax);
          
    plt.savefig(PATH + '_' + frelax + '_' + 'FT.pdf');
    plt.close();
    
    np.savetxt(PATH + '_FT.dat', np.c_[X,Y,Z])   

def symmetrize_c3(mat):
    mat = np.array(mat)
    out=np.zeros(np.shape(mat));
    m0=np.shape(out)[0]
    if(np.dot(par.R1,par.R2)>0):
        for i,j in product(range(m0),range(m0)):
            out[i,j]=(mat[i, j]+mat[np.mod(-i-j,m0), i]+mat[j, np.mod(-i-j,m0)])/3;
    else:
        for i,j in product(range(m0),range(m0)):
            out[i,j]=(mat[i, j]+mat[-j, i-j]+mat[j-i, -i])/3;
    return out;

def symmetrize_invert(mat):
    mat = np.array(mat)
    out=np.zeros(np.shape(mat));
    m0=np.shape(out)[0]
    for m1 in range(m0):
        m1_tar = np.mod(m0-m1,m0)
        for m2 in range(m0):
            m2_tar = np.mod(m0-m2,m0)
            out[m1,m2] = mat[m1_tar,m2_tar]
    return out;

def generate_kpath(hspoints_frac, N_points):

    hspoints_frac = np.asarray(hspoints_frac)
    
    hspoints = np.zeros(hspoints_frac.shape)
    for i in range(len(hspoints_frac)):
        hspoints[i,:] = hspoints_frac[i,0] * par.G1m + hspoints_frac[i,1] * par.G2m
    
    kpath = []
    for i in range(len(hspoints) -1):

        for j in range(N_points):
              kpoint = hspoints[i,:] + j*(hspoints[i+1,:] - hspoints[i,:])/N_points
              kpath.append(kpoint)

    kpath.append(hspoints[-1,:])

    return kpath

def checkfolders():

    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

    if not os.path.exists(ROOT_DIR + '/Bilayer'):
           os.makedirs(ROOT_DIR + '/Bilayer')

    if not os.path.exists(ROOT_DIR + '/CM_Data'):
           os.makedirs(ROOT_DIR + '/CM_Data')
           
    if not os.path.exists(ROOT_DIR + '/Monolayer/' + '/Layer_1/' + 'G'):
           os.makedirs(ROOT_DIR + '/Monolayer/' + '/Layer_1/' +  'G')
       
    if not os.path.exists(ROOT_DIR + '/Monolayer/' + '/Layer_1/' +  'K'):
           os.makedirs(ROOT_DIR + '/Monolayer/' + '/Layer_1/' + 'K')
           
    if not os.path.exists(ROOT_DIR + '/Monolayer/' + '/Layer_2/' +  'G'):
           os.makedirs(ROOT_DIR + '/Monolayer/' + '/Layer_2/' +  'G')
       
    if not os.path.exists(ROOT_DIR + '/Monolayer/' + '/Layer_2/' +  'K'):
           os.makedirs(ROOT_DIR + '/Monolayer/' + '/Layer_2/' +  'K')
    
    return ROOT_DIR
    
    
def test_func(x, A0, A1):
    return A0 + A1 * (x**2)    
    
    
def compute_bands_bilayer():
    directory = 'bands'
    PATH_Bilayer = os.path.join(ROOT_DIR, 'Bilayer', directory)
    PATH_SCF = os.path.join(ROOT_DIR, 'input_files', 'par.qe_scf.in')
    PATH_BANDS = os.path.join(ROOT_DIR, 'input_files', 'par.qe_bands.in')
    PATH_BANDSX = os.path.join(ROOT_DIR, 'input_files', 'par.qe_bandsx.in')
    PATH_QQEB = os.path.join(ROOT_DIR, 'input_files', 'par.qpar.qeb.job')
    shutil.copyfile(PATH_SCF, PATH_Bilayer)
    shutil.copyfile(PATH_BANDS, PATH_Bilayer)
    shutil.copyfile(PATH_BANDSX, PATH_Bilayer)
    shutil.copyfile(PATH_QQEB, PATH_Bilayer)
    os.system('sbatch ' + 'par.qpar.qeb.job' + '\n')
    

def strainanalysis(route, verbose = 'True'):
    os.chdir(route)
    Klist={}
    Glist={}
    Elist=[]
    for i in range(11):
        name=str((i-int(par.nstrain/2))*3)
        expansion=(i-int(par.nstrain/2)) * par.dstrain
        E1=lattice_energy(route, 'K/'+name)
        E2=lattice_energy(route, 'G/'+name)
        if((E1!='fail')&(E2!='fail')):
            Elist.append([expansion,E1,E2])
        Klist[name]=E1
        Glist[name]=E2
    
    A=[[u[0]**2,u[0],1] for u in Elist]
    bK=[u[1]  for u in Elist]
    bG=[u[2]  for u in Elist]
    explist = [u[0] for u in Elist]

    Ks1=np.linalg.lstsq(A,bK,rcond=None)[0]
    Gs1=np.linalg.lstsq(A,bG,rcond=None)[0]
    a1=par.alat*(1-Ks1[1]/2/Ks1[0])
    K=Ks1[0]/2*1000
    G=Gs1[0]/2*1000
    strainElist={'isotropic':Klist,'anisotropic':Glist,'K':K,'G':G,'a':a1}
    
    if verbose:
       print('the computed lattice constant is:', a1)
       print('K = ', K, ' G = ', G)
       
    return strainElist, K, G

def lattice_energy(route, directory):   #extract the DFT energy in folder 
    PATH = os.path.join(route, directory)
    os.chdir(PATH)
    
    out = 'fail'
    with open('relax.out', 'r') as infile:
        lines = infile.readlines()
    for line in lines:
        string = line.split()        
        if 'Final' in string and 'energy' in string:
            out = float(string[3]) * par.Ry_to_eV    
    if out == 'fail':
       print('ERROR: No final energy found in \n', PATH)                    
    return out


def readGSFE(PATH, symm_c3 = True, plot_GSFEk = True):

    GSFE = [[lattice_energy(PATH, str(i) + '_' + str(j)) for j in range(par.ngrid)] for i in range(par.ngrid)]
    energy={}
    for i,j in product(range(par.ngrid+1),range(par.ngrid+1)):
        shift= list((i/par.ngrid) * par.lattice_matrix[0] + (j/par.ngrid) * par.lattice_matrix[1])[:2]
        name = str(i)+ '_' + str(j)
        if(i == par.ngrid):
            i=0
        if(j == par.ngrid):
            j=0
        energy[name]=[ shift + [GSFE[i][j]]]  
        
    GSFEdic = {'GSFElist':energy}   
    GSFE = np.zeros((par.ngrid,par.ngrid))    
    for i,j in product(range(par.ngrid),range(par.ngrid)):
        u = np.mod(i,par.ngrid)
        v = np.mod(j,par.ngrid)
        name = str(u) + '_' + str(v)
        GSFE[i][j] = GSFEdic["GSFElist"][name][0][2]
        
    GSFE = GSFE*1000.
    GSFE -= np.min(GSFE)*np.ones([par.ngrid,par.ngrid])
    
    if symm_c3:            
       GSFE=np.array(symmetrize_c3(GSFE))
    
    GSFE_k=fft2(GSFE)/par.ngrid**2   
    
    if plot_GSFEk:
       ROOT_DIR = checkfolders()
       PATH = os.path.join(ROOT_DIR, 'CM_Data', 'GSFE_k')
       plotterk(PATH, GSFE_k, 'GSFE', 'relax')
       file=open(ROOT_DIR + '/GSFE.json','w')
       json.dump(GSFEdic, file )
    
    np.save(os.path.join(ROOT_DIR,'GSFE_k'), GSFE_k)
    return GSFE, GSFE_k, GSFEdic

def GSFE_u(u, GSFE_k): 
    klist=[np.array([k1,k2]) for k1,k2 in product(range(-par.ngrid,par.ngrid),range(-par.ngrid,par.ngrid))
          if np.linalg.norm(k1 * par.gv1 + k2 * par.gv2) < par.kmax * np.linalg.norm(par.gv1)]
    Energy = 0
    for k in klist:
        core = 0,   
        for i,j in product(range(par.q),range(par.q)):
            pos = 2*(i*par.q+j)
            u1 = np.array([u[pos],u[pos+1]])
            core += np.exp(2*np.pi*1j*np.dot(k,np.array([i,j])/par.q+2*u1))
        Energy += GSFE_k[k[0],k[1]]*core
    return np.real(Energy)

def GSFE_der(u, GSFE_k):
    klist=[np.array([i,j]) for i,j in product(range(-par.ngrid,par.ngrid),range(-par.ngrid,par.ngrid))
            if np.linalg.norm(i*par.gv1+j*par.gv2)<par.kmax*np.linalg.norm(par.gv1)];
    Du = np.zeros(2*par.q**2);
    for i,j in product(range(par.q),range(par.q)):
        core = 0
        pos = 2*(i*par.q+j);
        u1 = np.array([u[pos],u[pos+1]]);
        for k in klist:
            core += GSFE_k[k[0],k[1]]*2*np.pi*2j*k*np.exp(2*np.pi*1j*np.dot(k,np.array([i,j])/par.q+2*u1))
        Du[pos] = np.real(core[0]);
        Du[pos+1]= np.real(core[1]);
    return Du;
    
def ltf(mat):
    Uint = np.linalg.inv(par.MSL)*mat*par.lattice_matrix[:2,:2]
    return (Uint+Uint.T)/2
    
def strain(u,K,G):     
    E = 0
    for i,j in product(range(par.q),range(par.q)):
        ip,jp = np.mod(i+1,par.q),np.mod(j+1,par.q)
        u12 = (u[2*(i*par.q+jp)]-u[2*(i*par.q+j)])*par.q
        u11 = (u[2*(ip*par.q+j)]-u[2*(i*par.q+j)])*par.q
        u22 = (u[2*(i*par.q+jp)+1]-u[2*(i*par.q+j)+1])*par.q
        u21 = (u[2*(ip*par.q+j)+1]-u[2*(i*par.q+j)+1])*par.q
        U2 = ltf(np.matrix([[u11,u21],[u12,u22]]))
        E += (K+G)/2*np.trace(U2)**2-2*G*np.linalg.det(U2)
    return np.sum(E)
       
def strain_der(u,K,G):
    dE = np.zeros(2*par.q**2)
    for i,j in product(range(par.q),range(par.q)):
        ip,jp = np.mod(i+1,par.q),np.mod(j+1,par.q)        
        u12 = (u[2*(i*par.q+jp)]-u[2*(i*par.q+j)])*par.q
        u11 = (u[2*(ip*par.q+j)]-u[2*(i*par.q+j)])*par.q
        u22 = (u[2*(i*par.q+jp)+1]-u[2*(i*par.q+j)+1])*par.q
        u21 = (u[2*(ip*par.q+j)+1]-u[2*(i*par.q+j)+1])*par.q           
        U0 = np.matrix([[u11,u21],[u12,u22]])
        dU1 = ltf(np.matrix([[-1,0],[-1,0]])*par.q)
        dU2 = ltf(np.matrix([[0,-1],[0,-1]])*par.q)
        dUi1 = ltf(np.matrix([[1,0],[0,0]])*par.q)
        dUi2 = ltf(np.matrix([[0,1],[0,0]])*par.q)
        dUj1 = ltf(np.matrix([[0,0],[1,0]])*par.q)
        dUj2 = ltf(np.matrix([[0,0],[0,1]])*par.q)
        U2 = ltf(U0)
        C1 = 2*(K+G)/2*np.trace(U2)
        dE[2*(i*par.q+j)] += C1*np.trace(dU1)-2*G*(U2[0,0]*dU1[1,1]+U2[1,1]*dU1[0,0]-2*U2[0,1]*dU1[0,1])
        dE[2*(i*par.q+j)+1] += C1*np.trace(dU2)-2*G*(U2[0,0]*dU2[1,1]+U2[1,1]*dU2[0,0]-2*U2[0,1]*dU2[0,1])
        dE[2*(i*par.q+jp)] += C1*np.trace(dUj1)-2*G*(U2[0,0]*dUj1[1,1]+U2[1,1]*dUj1[0,0]-2*U2[0,1]*dUj1[0,1])
        dE[2*(i*par.q+jp)+1] += C1*np.trace(dUj2)-2*G*(U2[0,0]*dUj2[1,1]+U2[1,1]*dUj2[0,0]-2*U2[0,1]*dUj2[0,1])
        dE[2*(ip*par.q+j)] += C1*np.trace(dUi1)-2*G*(U2[0,0]*dUi1[1,1]+U2[1,1]*dUi1[0,0]-2*U2[0,1]*dUi1[0,1])
        dE[2*(ip*par.q+j)+1] += C1*np.trace(dUi2)-2*G*(U2[0,0]*dUi2[1,1]+U2[1,1]*dUi2[0,0]-2*U2[0,1]*dUi2[0,1])
    return dE

def relax(GSFE_k, Kt, Gt, Kb, Gb): 
    ROOT_DIR = checkfolders()   
    u_rel = np.zeros([par.q,par.q,2])
    time_start = time.time()
    history = []
    
    def f(u):
        ut = np.array([u[i] for i in range(len(u)) if np.mod(i,4)<2])
        ub = np.array([u[i] for i in range(len(u)) if np.mod(i,4)>=2])
        Energy = float(GSFE_u(ut-ub, GSFE_k) + strain(ut,Kt,Gt) + strain(ub,Kb,Gb))
        print('GSFE + strain energy:', Energy, 'meV')
        history.append(Energy)
        return Energy 

    def df(u):
        ut = np.array([u[i] for i in range(len(u)) if np.mod(i,4)<2])
        ub = np.array([u[i] for i in range(len(u)) if np.mod(i,4)>=2])
        Dt = GSFE_der(ut-ub, GSFE_k) + strain_der(ut,Kt,Gt)
        Db = -GSFE_der(ut-ub, GSFE_k) + strain_der(ub,Kb,Gb)
        res = np.array([Dt[int(i/2)] for i in range(2*len(ut))])
        for i in range(par.q**2):
            res[4*i+2]=Db[2*i]
            res[4*i+3]=Db[2*i+1]
        return res

    init = np.zeros(4*par.q**2)
    res = minimize(f,init ,jac = df, method='BFGS', options = {'gtol':5e-4})
    time_end = time.time()
    
    print('Lattice relaxation took: ', time_end - time_start, ' s')
    np.savetxt(ROOT_DIR + '/Relaxation_history.txt', history)
    
    for i,j in product(range(par.q),range(par.q)):
        u_rel[i,j,0] = res.x[4*(i*par.q+j)]-res.x[4*(i*par.q+j)+2]
        u_rel[i,j,1] = res.x[4*(i*par.q+j)+1]-res.x[4*(i*par.q+j)+3]

    file = ROOT_DIR + '/Bilayer_relaxed'
    np.save(file, u_rel)
    return u_rel


def include_relax(field, u_raw):

   def d0(r, theta):
      return theta*np.matrix([[-r[1]],[r[0]]]);

   def u1(u,v):      
       u_re=par.uv1[:2]*u_raw[u,v,0]+par.uv2[:2]*u_raw[u,v,1];
       out = np.matrix([[u_re[0]],[u_re[1]]]);
       return out
       
   def M_realspace(mat):   

      def f(r1):
          r=[r1[0,0],r1[1,0]]
          u1=-np.cross(par.uv2[:2],r)/par.Scross_norm
          u2=np.cross(par.uv1[:2],r)/par.Scross_norm     
          klist=[[i,j] for i,j in product(range(-par.ngrid,par.ngrid),range(-par.ngrid,par.ngrid)) if np.linalg.norm(i*par.gv1+j*par.gv2)<par.kmax*np.linalg.norm(par.gv1)];
          real=np.sum([mat[k[0],k[1]]*np.exp(2*np.pi*1j*(u1*k[0]+u2*k[1])) 
                         for k in klist])
          return real
            
      out = np.zeros([par.q,par.q])*1j;
      for u,v in product(range(par.q),range(par.q)):
          r = d0(u*par.uv1m/par.q+v*par.uv2m/par.q, np.radians(par.theta)) + u1(u,v)
          out[u][v] = f(r)  

      return np.real(out)
   
   field_real_relax = M_realspace(field)  
   field_real_relax = symmetrize_c3(field_real_relax)
   field_k_relax = fft2(field_real_relax)/par.q**2 
   
   return field_real_relax, field_k_relax



