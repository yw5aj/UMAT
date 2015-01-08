os.chdir('x:/workfolder/abaqusfolder/umat')
from abqimport import *
import subprocess, copy
# Taken from abaqus_v6.env
compile_fortran=['ifort',
                 '/c','/DABQ_WIN86_64', '/extend-source',
                 '/iface:cref', '/recursive', '/Qauto-scalar',
                 '/QxSSE3', '/QaxAVX', 
                 '/heap-arrays:1', 
                 # '/Od', '/Ob0'   # <-- Optimization 
                 # '/Zi',          # <-- Debugging
                 '/include:%I']

codename = 'NeoHookeanNumericF90.f90'
# codename = 'NeoHookeanNumericF77.for'
objname = 'NeoHookeanNumeric.obj'
compile_args = copy.deepcopy(compile_fortran)
compile_args.extend([codename, '/object:%s'%objname])
subprocess.call(compile_args)


mdb.jobs['SingleElem'].setValues(
    userSubroutine='X:\\WorkFolder\\AbaqusFolder\\umat\\%s'%objname)
mdb.jobs['SingleElem'].submit(consistencyChecking=OFF)

