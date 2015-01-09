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


def compile_umat(codename_list):
    assert compile_fortran[0] == 'ifort', 'This code only works with IVF.'
    compile_args = copy.deepcopy(compile_fortran)
    compile_args.extend(codename_list)
    assert subprocess.call(compile_args) == 0, \
        'Compilation error. See the command line window for more details.'
    return

def run_umat(codename_list, jobname, objname, wait=False):
    compile_umat(codename_list)
    mdb.jobs[jobname].setValues(userSubroutine=objname)
    mdb.jobs[jobname].submit(consistencyChecking=OFF)
    if wait:
        mdb.jobs[jobname].waitForCompletion()
    return


if __name__ == '__main__':
    # codename_list = ['umat.f90', 'umatutils.f90', 'numerichyper.f90', 'modpsi_nh.f90']
    codename_list = ['NeoHookeanNumeric.f90']
    objname = 'umat.obj'
    jobname = 'SingleElem'
    run_umat(codename_list, jobname, objname)


