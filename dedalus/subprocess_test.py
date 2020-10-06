# import subprocess
# from time import sleep
# import numpy as np 

# nProc = 4
# program = ['mpiexec','-np',str(nProc),'python3','active_matter.py']  
# returnCode = None
# with open('out.txt','w+') as fout:
#     # while not returnCode==0:
#     p = subprocess.call(program,stdout=fout)#(program, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     print(fout.read())
#         # returnCode = None
#         # arr = np.ones(10)
#         # p.stdin.write(arr.tobytes())
#         # while returnCode is None:
#             # message =  p.stdout.readline().rstrip("\n").split()
#             # if message != None: print("\t\t"+" ".join(message)+"\t\t")
#             # sleep(0.1)
#             # out = p.stdout.read()
#             # returnCode = p.poll()
#         # if returnCode != 0: print("\t\t\t\t Error n: ",p.poll()," resetting simulation...")

from mpi4py import MPI    
import sys
client_script = 'active_matter.py'
comm = MPI.COMM_SELF.Spawn(sys.executable, args=[client_script], maxprocs=5)