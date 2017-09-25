from dolfin import *
#from cbcpost import *
#from cbcpost.utils import create_submesh, create_slice
import numpy as np


TR    = 5.12*10**-3
TE    = 2.29*10**-3
theta = 8.0
r1,r2 = 3.6 , 6.3
R1 = 0.92
k1 = TR*r1
k2  = TE*r2 
aux0 = np.cos(np.deg2rad(theta)) 
aux1= 1.-aux0
aux2 = 1.+aux0  
G =  np.exp(-TR*R1) 

def BPII (c ):
    return  ( (1.-aux0*G)*(1.-G*np.exp(-k1*c*10**-3))*np.exp(-k2*c*10**-3)/((1.-aux0*G*np.exp(-k1*c*10**-3))*(1.-G)) -1 )*100
   

def boundary(x, on_boundary): 
  return on_boundary


if __name__=='__main__':
	import argparse
        import sys
	parser = argparse.ArgumentParser(prog='mri-contrast.py')
	parser.add_argument('--dt', default=0.1, type=float, help='')
	parser.add_argument('--D',  default=1.0, type=float , help='')
 	parser.add_argument('--C', default=1.0,  type=float ,help='')
	parser.add_argument('--mesh',  default="bad_mesh.xml" , help='')   
        parser.add_argument('--out',  default="mri-contrast",type=str , help='')   
        parser.add_argument('--init',  default="" , type=str, help='')  
        parser.add_argument('--Tend',  default=12.0 , type=float, help='')  
	Z = parser.parse_args()
        print Z
        mesh = Mesh()
        if Z.mesh.endswith("xml"):
           mesh = Mesh(Z.mesh)
        elif Z.mesh.endswith("xdmf"):
		xdmf =XDMFFile(mesh.mpi_comm(),Z.mesh)
		xdmf.read(mesh)
        else :
           print "undetermined filetype"
           sys.exit()



####################################################
	V = FunctionSpace(mesh, "CG",1)

	u = TrialFunction(V)
	v = TestFunction(V)

	U = Function(V)   # current U
	U_ = Function(V)  # previous U 

	U0 = Function(V)  # initial U
       

        bpii = Function(V)
        
##################################################   

        dt_ = float(Z.dt)
        print dt_
	
        C = Constant(float(Z.C))
	dt = Constant(dt_)
        Tend = Z.Tend
    
      

        if Z.init.endswith("h5"):
          in0= HDF5File(mesh.mpi_comm(),Z.init,"r") # basename
	  in0.read(U0,"/mric")
          in0.close()
          C = Constant ( U0.vector().array().max() ) #5000 # parallel
          bc  = DirichletBC(V, C, boundary ) 
         
        else :
          C = Constant(float(Z.C))
          bc  = DirichletBC(V, C, boundary ) 
          U0.vector()[:]=0

        U_.assign(U0)
#############################################################      
	m =u*v*dx 
	k = dt*Constant(Z.D)*dot(grad(v),grad(u))*dx

        L =U_*v*dx  

       
        print dt
###########################################################
        dia =  assemble (action(m, Constant(1)))
	M= assemble(m)
        M.zero()
        M.set_diagonal(dia)
          
        K = assemble(k)
        A = M+K
      
        a = u*v*dx + dt*Constant(Z.D)*dot(grad(v),grad(u))*dx
	A =assemble(a)
   
	bc.apply(A)


	t = 0.0
	time_step =0
        
        save_resolution=0.01      

        outfile = File("results%sD%sdt%s/mri-contrast.pvd"%(Z.out,Z.D,Z.dt)) 
        
        save_step = save_resolution*(Tend/dt_)

     
        print "solving"
#################################################################
	while t < Tend+dt_/2:


	  t += dt_

          b = assemble(L)
          bc.apply(b)
          

	  solve(A, U.vector(), b)

          U_.assign(U)

################################################################
         
	  if time_step%save_step==0: 
               bpii.vector()[:] = BPII(U.vector().array())
               outfile << U
               print time_step
               
          time_step+=1
	 
        

