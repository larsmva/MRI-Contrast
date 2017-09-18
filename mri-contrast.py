from dolfin import *
#from cbcpost import *
#from cbcpost.utils import create_submesh, create_slice

#default values
No_refinements=0 
dt_val=0.1 
D_val=1.0
C0=1 
meshfile = "bad_mesh.xml"
#hack to read input arguments 
import sys

def boundary(x, on_boundary): 
  return on_boundary


if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser(prog='mri-contrast.py')
	parser.add_argument('--dt', default=0.1, help='')
	parser.add_argument('--D',  default=1.0 , help='')
 	parser.add_argument('--C', default=1.0, help='')
	parser.add_argument('--mesh',  default="bad_mesh.xml" , help='')   
        parser.add_argument('--save_step',  default=20 , help='') 
        parser.add_argument('--out',  default="mri_contrast.xdmf" , help='')   
	Z = parser.parse_args()
    
        mesh = Mesh()
        if Z.mesh.endswith("xml"):
           mesh = Mesh(Z.mesh)
        elif Z.mesh.endswith("xdmf"):
		xdmf =XDMFFile(mesh.mpi_comm(),Z.mesh)
		xdmf.read(mesh)
        else :
           print "undetermined filetype"

	V = FunctionSpace(mesh, "CG",1)

	u = TrialFunction(V)
	v = TestFunction(V)

	U = Function(V)   # current U
	U_ = Function(V)  # previous U 

	U0 = Function(V)  # initial U
	U1 = Function(V) 

	
	#in0= HDF5File(mesh.mpi_comm(),"MRI0.h5","r") # basename
	#in0.read(U0,"/mri0")
	#in0.close()

	#in1 = HDF5File(mesh.mpi_comm(),"MRI1.h5","r")
	#in1.read(U1,"/mri1")
	#in1.close()





	Tend = 18.0
	t = 0.0
	D = Constant(Z.D)
        C = Constant(Z.C)
	dt = Constant(Z.dt)

	a = D*dt*inner(grad(u), grad(v))*dx 
	L =U_*v*dx  

	A = assemble(a)


	U_.assign(U0)


	time_step =0


	writer = XDMFFile(mesh.mpi_comm(), Z.out)

	while t < Tend-dt_val/2:
	 

	  bc  = DirichletBC(V, C, boundary)  
	  t += dt_val 
	  bc.apply(A)
	  b = assemble(L)
	  bc.apply(b)
	  solve(A, U.vector(), b, "gmres", "amg")

	  U_.assign(U)
	  
	  if (time_step%Z.save_step==0):
	     writer.write(U,time_step)
	  
	  time_step+=1
          print time_step	  

