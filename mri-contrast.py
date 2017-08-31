from dolfin import *
from cbcpost import *
from cbcpost.utils import create_submesh, create_slice

#default values
No_refinements=0 
dt_val=0.1 
D_val=1.0
C0=1 
#hack to read input arguments 
import sys
for s in sys.argv[1:]:
  exec(s)
  print s
print "dt_val ", dt_val, " No_ref ", No_refinements , "D" , D_val

def boundary(x, on_bounday): 
  return on_bounday

mesh = Mesh()
xdmf =XDMFFile(mesh.mpi_comm(),"pial_mesh.xdmf")
xdmf.read(mesh)


V = FunctionSpace(mesh, "Lagrange", 1)

u = TrialFunction(V)
v = TestFunction(V)

U = Function(V)   # current U
U_ = Function(V)  # previous U 

U0 = Function(V)  # initial U
U1 = Function(V) 


    
	
in0= HDF5File(mesh.mpi_comm(),"MRI0.h5","r") # basename
in0.read(U0,"MRI0")
in0.close()

in1 = HDF5File(mesh.mpi_comm(),"MRI1.h5","r")
in1.read(U1,"MRI1")
in1.close()




Tend = 18.0
t = 0.0
D = Constant(D_val)
dt = Constant(dt_val)

a = u*v*dx + D*dt*inner(grad(u), grad(v))*dx 
L = U_*v*dx 

A = assemble(a)


#slicemesh_x = create_slice(mesh,(-9.80,4.98,27.24) , (1,0,0))
#slicemesh_y = create_slice(mesh,(-7.49,4.98,27.24)  ,( 0,1,0))
#pp =  PostProcessor(dict(casedir="nscmD%sDt%s"%(D_val,dt_val)))
#pp.add_fields( [SolutionField("Error",dict(save=False)) ,Norm("Error",dict(save=True,plot=True))])
#pp.add_field(SolutionField("C_X"))
#pp.add_fields([SubFunction("C_X",slicemesh_x,dict(save=True))])
#pp.add_field(SolutionField("C_Y"))
#pp.add_fields([SubFunction("C_Y",slicemesh_y,dict(save=True))])

time_step =0
_v = Function(V)
U_.assign(U1)

save_step=10
writer = XDMFFile(mesh.mpi_comm(), "mri-contrast.xdmf")
I = 200

while t < Tend-dt_val/2:
 

  bc  = DirichletBC( V, I, boundary)  
  t += dt_val 
  bc.apply(A)
  b = assemble(L)
  bc.apply(b)
  solve(A, U.vector(), b, "gmres", "amg")
  
  E = assemble((U-U3)**2*dx)
  U_.assign(U)

  if (t%save_step==0):
     writer.write(U)

  time_step+=1
  
  #pp.update_all({"C_X": lambda: U} , t ,time_step ) 
  #pp.update_all({"C_Y": lambda: U} , t ,time_step )
  #pp.update_all({"Error":lambda:  E}, t ,time_step )
  
#pp.finalize_all()

