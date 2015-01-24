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


def combine_modules(codename_list, objname):
    umatfname = objname.replace('.obj', '.f90')
    with open(umatfname, 'w') as fwrite:
        for codename in codename_list:
            with open(codename, 'r') as fread:
                fwrite.write(fread.read())
                fwrite.write('\n\n')
    return umatfname

    
def compile_umat(umatfname):
    assert compile_fortran[0] == 'ifort', 'This code only works with IVF.'
    compile_args = copy.deepcopy(compile_fortran)
    compile_args.append(umatfname)
    assert subprocess.call(compile_args) == 0, \
        'Compilation error. See the command line window for more details.'
    return


def run_umat(codename_list, jobname, objname, wait=False):
    umatfname = combine_modules(codename_list, objname)
    compile_umat(umatfname)
    mdb.jobs[jobname].setValues(userSubroutine=objname)
    mdb.jobs[jobname].submit(consistencyChecking=OFF)
    if wait:
        mdb.jobs[jobname].waitForCompletion()
    return


if __name__ == '__main__':
    # codename_list = ['umatutils.f90', 'modpsi_neo.f90',  'numerichyper.f90', 'nhinterface.f90']
    # codename_list = ['umatutils.f90', 'modpsi_hgo.f90',  'numerichyper.f90', 'nhinterface.f90']
    # codename_list = ['umatutils.f90', 'modpsi_hg.f90',  'numerichyper.f90', 'nhinterface.f90']
    # codename_list = ['umatutils.f90', 'modpsi_hg.f90',  'numerichyper.f90', 'nhcylinterface.f90']
    codename_list = ['umatutils.f90', 'modpsi_hgo.f90',  'numerichyper.f90', 'nhcylinterface.f90']
    objname = 'umat.obj'
    # jobname = 'SingleElem'
    # jobname = 'SingleElemShearNumeric'
    # jobname = 'HolzapfelNumeric'
    # jobname = 'numeric_circ_k0_displ'
    jobname = 'ArteryInflSymmNumeric'
    run_umat(codename_list, jobname, objname)


