from dolfin import *
#from cbcpost import *
#from cbcpost.utils import create_submesh, create_slice


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
       
        
##################################################   

        dt_ = float(Z.dt)

	D = Constant(float(Z.D))
        C = Constant(float(Z.C))
	dt = Constant(dt_)
        Tend = Z.Tend
    
        print Z.D, Tend   


        if Z.init.endswith("h5"):
          in0= HDF5File(mesh.mpi_comm(),Z.init,"r") # basename
	  in0.read(U0,"/mric")
          in0.close()
          C = U0.vector().array().max() #5000 # parallel
          bc  = DirichletBC(V, C, boundary ) 
          print C
        else :
          bc  = DirichletBC(V, C, boundary ) 
          U0.vector()[:]=0

#############################################################      
	m =u*v*dx 
	k = dt*D*inner(grad(v),grad(u))*dx

        L =U_*v*dx  



###########################################################
        action_form = action(m, Constant(1))
	#A = assemble(a)
        M= assemble(m)
        M.zero()
        M.set_diagonal(assemble(action_form))
          
        K = assemble(k)
        A = M+K
        U_.assign(U0)
  
	
   
	

   

          
	t = 0.0
	time_step =0
        
        save_resolution=0.1      
	#writer = XDMFFile(mesh.mpi_comm(), Z.out)
        outfile = File("resultsD%sdt%s/mri-contrast.pvd"%(Z.D,Z.dt)) 
        
        save_step = save_resolution*(Tend/dt_)
        print save_step 
     
        print "solving"
#################################################################
	while t < Tend-dt_/2:


	  t += dt_

          b = assemble(L)
          bc.apply(b)
          bc.apply(A)

	  solve(A, U.vector(), b, "gmres" , "amg" )

          print assemble((U-U_)**2*dx)
          U_.assign(U)
          U.vector()[:]
################################################################
       
	  if time_step%save_step==0:
               outfile << U 
          time_step+=1
	 

