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


def compile_umat(codename, objname=None):
    for fname in os.listdir('.'):
        if fname.endswith('.mod'):
            os.remove(fname)
    if objname is None:
        objname = codename[:codename.index('.')] + '.obj'
    compile_args = copy.deepcopy(compile_fortran)
    compile_args.extend([codename, '/object:%s'%objname])
    assert subprocess.call(compile_args) == 0, \
        'Compilation error. See the command line window for more details.'
    return objname

def run_umat(codename, jobname, wait=False):
    objname = compile_umat(codename)
    mdb.jobs[jobname].setValues(userSubroutine=objname)
    mdb.jobs[jobname].submit(consistencyChecking=OFF)
    if wait:
        mdb.jobs[jobname].waitForCompletion()
    return


if __name__ == '__main__':
    codename = 'NeoHookeanNumeric.f90'
    jobname = 'SingleElem'
    run_umat(codename, jobname)

