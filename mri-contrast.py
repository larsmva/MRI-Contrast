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
        parser.add_argument('--save_step',  default=20, type=int , help='') 
        parser.add_argument('--out',  default="mri_contrast.xdmf",type=str , help='')   
        parser.add_argument('--init',  default="" , type=str, help='')  
        parser.add_argument('--Tend',  default=8.0 , type=float, help='')  
	Z = parser.parse_args()
    
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
        Ub = Function(V)  # initial U

##################################################   

        dt_ = float(Z.dt)

	D = Constant(float(Z.D))
        C = Constant(float(Z.C))
	dt = Constant(dt_)
        Tend = Z.Tend

#############################################################      
	a =u*v*dx+ D*dt*inner(grad(u), grad(v))*dx 
	L =U_*v*dx  
############################################################

        if Z.init=="" :
          bc  = DirichletBC(V, C, boundary ) 
          U0.vector()[:]=0
        elif Z.mesh.endswith("h5"):
          in0= HDF5File(mesh.mpi_comm(),Z.init,"r") # basename
	  in0.read(U0,"/mric")
          in0.close()
        else :
           print "undetermined filetype" 
           sys.exit()

###########################################################

	A = assemble(a)
        bc.apply(A)
	b = assemble(L)
        bc.apply(b)
	U_.assign(U0)

   
	t = 0.0
	time_step =0
        
        
	writer = XDMFFile(mesh.mpi_comm(), Z.out)
        """
          Diffusion on surface :: 

          use prevoius solution boundary as dirichlet  ? 

          bc = DirichletBC(V,U, boundary)
          check Abel
        """
#################################################################
	while t < Tend-dt_/2:


	  t += dt_
	  solve(A, U.vector(), b, "gmres", "amg")
	 

##################################################################
	  a = U.vector()[:]
          b=  U.vector()[:]
          b[a>b]  = a[a>b]
          Ub.vector()[:] = b

#################################################################

          U_.assign(U)
	  if (time_step%Z.save_step==0):
	     writer.write(U, t) 
	  
	  time_step+=1
          print time_step	  

