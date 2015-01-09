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


def compile_umat(code_list, objname=None):
    if objname is None:
        objname = 'UMAT.obj'
    compile_args = copy.deepcopy(compile_fortran)
    compile_args.extend([code_list, '/o %s'%objname])
    assert subprocess.call(compile_args) == 0, \
        'Compilation error. See the command line window for more details.'
    return objname

def run_umat(code_list, jobname, wait=False):
    objname = compile_umat(code_list)
    mdb.jobs[jobname].setValues(userSubroutine=objname)
    mdb.jobs[jobname].submit(consistencyChecking=OFF)
    if wait:
        mdb.jobs[jobname].waitForCompletion()
    return


if __name__ == '__main__':
    code_list = ['umat_nh.f90', 'umatutils.f90', 'numerichyper.f90', 'modpsi_nh.f90']
    jobname = 'SingleElem'
    run_umat(codelist, jobname)

