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
U2 = Function(V)
U3 =Function(V)

in0= HDF5File(mesh.mpi_comm(),"C00.h5","r") # basename
in0.read(U0,"C00")
in0.close()

in1 = HDF5File(mesh.mpi_comm(),"C01.h5","r")
in1.read(U1,"C01")
in1.close()

in2 =HDF5File(mesh.mpi_comm(),"C02.h5","r")
in2.read(U2,"C02")
in2.close()

in3 = HDF5File(mesh.mpi_comm(),"C03.h5","r")
in3.read(U3,"C03")
in3.close()


#w1  = XDMFFile(mesh.mpi_comm(),"nscm_C01.xdmf")
#w1.write(U1)
#w2 = XDMFFile(mesh.mpi_comm(),"nscm_C02.xdmf")
#w2.write(U2)
#w3 = XDMFFile(mesh.mpi_comm(),"nscm_C03.xdmf")
#w3.write(U3)


Tend = 5.0
t = 0.0
D = Constant(D_val)
dt = Constant(dt_val)

a = u*v*dx + D*dt*inner(grad(u), grad(v))*dx 
L = U_*v*dx 

A = assemble(a)


slicemesh_x = create_slice(mesh,(-9.80,4.98,27.24) , (1,0,0))
slicemesh_y = create_slice(mesh,(-7.49,4.98,27.24)  ,( 0,1,0))

#pp =  PostProcessor(dict(casedir="nscmD%sDt%s"%(D_val,dt_val)))

#pp.add_fields( [SolutionField("Error",dict(save=False)) ,Norm("Error",dict(save=True,plot=True))])



#pp.add_field(SolutionField("C_X"))
#pp.add_fields([SubFunction("C_X",slicemesh_x,dict(save=True))])

#pp.add_field(SolutionField("C_Y"))
#pp.add_fields([SubFunction("C_Y",slicemesh_y,dict(save=True))])

time_step =0
_v = Function(V)
U_.assign(U1)

save_step=0
writer6 = XDMFFile(mesh.mpi_comm(), "nscm0_fullD86.xdmf")
writer4 = XDMFFile(mesh.mpi_comm(), "nscm0_fullD84.xdmf")
while t < Tend-dt_val/2:
  
  if t < 3.0 : 
    _v=  U1*(3.-t)/3.0 + U2*t/3.0
  elif t>=3.0 and t<5.0:
    _v= U2*(5.-t)/2.  +U3*(t-3.0)/2.
  else :
    _v= U3
  

  bc  = DirichletBC( V, _v,boundary)  
  t += dt_val 
  bc.apply(A)
  b = assemble(L)
  bc.apply(b)
  solve(A, U.vector(), b, "gmres", "amg")
  
  #E = assemble((U-U3)**2*dx)
  U_.assign(U)

  print "t ", t
  #pp.update_all({"C_X": lambda: U} , t ,time_step ) 
  #pp.update_all({"C_Y": lambda: U} , t ,time_step )
  #pp.update_all({"Error":lambda:  E}, t ,time_step )
  if t>3.-dt_val/2 and t<3. + dt_val/2:  
    print "save"
    writer4.write(U)
  elif t>5.-dt_val/2 and t < 5. +dt_val/2:
    print "save"
    writer6.write(U)
  
  time_step+=1

  
#pp.finalize_all()

