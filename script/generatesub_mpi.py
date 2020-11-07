import numpy as np

numnode = 7

strongscriptdir = 'strong/'
weakscriptdir = 'weak/'
blockscriptdir = 'block/'
mpscriptdir = 'mp/'
mpiscriptdir = 'mpi/'
mpmpiscriptdir = 'mpmpi/'
strongresultdir = 'results/baseline/strong/'
weakresultdir = 'results/baseline/weak/'
blockresultdir = 'results/block/baseline/'
blocktuneresultdir = 'results/block/tune/'
mpresultdir = 'results/mp/'
mpiresultdir = 'results/mpi/'
mpmpiresultdir = 'results/mpmpi/'
mydir = '$HOME/projects/CS5220/shallow-water-project'

# # generate strong scaling .sub files
# for i in range(0, numnode):
#     filename = strongscriptdir + "shallow"+str(i+1)+".sub"
#     f = open(filename, "w+")
#     f.write('#!/bin/bash\n')
#     f.write('#SBATCH -J shallow\n')
#     f.write('#SBATCH -o '+ strongresultdir + 'shallow_'+str(i+1)+'n.out\n')
#     f.write('#SBATCH -e '+ strongresultdir + 'shallow_'+str(i+1)+'n.err\n')
#     f.write('#SBATCH --nodes=' + str(i+1) + ' \n')
#     f.write('#SBATCH --ntasks=' + str(i+1) + ' \n')
#     f.write('#SBATCH --tasks-per-node=1 \n')
#     f.write('#SBATCH --cpus-per-task=1 \n')
#     f.write('#SBATCH --get-user-env \n')
#     f.write('#SBATCH -t 00:10:00 \n')
#     f.write('#SBATCH --mem-per-cpu=1000 \n')
#     f.write('#SBATCH --partition=cs5220 \n')
#     f.write('\n')
#     f.write('source /etc/profile.d/modules.sh \n')
#     f.write('#module load openmpi-4.0.0 \n')
#     f.write('cd ' + mydir + ' \n')
#     f.write('src/lshallow tests.lua dam 1000 \n')
#     f.write('\n')

# # generate weak scaling .sub files
# for i in range(0, numnode):
#     filename = weakscriptdir + "shallow"+str(i+1)+"_weak.sub"
#     f = open(filename, "w+")
#     f.write('#!/bin/bash\n')
#     f.write('#SBATCH -J shallow\n')
#     f.write('#SBATCH -o '+ weakresultdir + 'shallow_'+str(i+1)+'n.out\n')
#     f.write('#SBATCH -e '+ weakresultdir + 'shallow_'+str(i+1)+'n.err\n')
#     f.write('#SBATCH --nodes=' + str(i+1) + ' \n')
#     f.write('#SBATCH --ntasks=' + str(i+1) + ' \n')
#     f.write('#SBATCH --tasks-per-node=1 \n')
#     f.write('#SBATCH --cpus-per-task=1 \n')
#     f.write('#SBATCH --get-user-env \n')
#     f.write('#SBATCH -t 00:10:00 \n')
#     f.write('#SBATCH --mem-per-cpu=1000 \n')
#     f.write('#SBATCH --partition=cs5220 \n')
#     f.write('\n')
#     f.write('source /etc/profile.d/modules.sh \n')
#     f.write('#module load openmpi-4.0.0 \n')
#     f.write('cd ' + mydir + ' \n')
#     f.write('src/lshallow tests.lua dam ' + str((i+1)*200) + ' \n')
#     f.write('\n')

# generate block param tuning scripts
# (batch time step, 2, 4, 6, 8; block size 32, 64, 128, 256; problem size, 400x400, 800x800, 1600x1600)
# # generate baseline scripts
# grid = [400, 800, 1600]
# for i in range(0, len(grid)):
#     filename = blockscriptdir + "shallow_block_baseline_"+str(i+1)+".sub"
#     f = open(filename, "w+")
#     f.write('#!/bin/bash\n')
#     f.write('#SBATCH -J shallow\n')
#     f.write('#SBATCH -o '+ blockresultdir + 'shallow_'+str(i+1)+'n.out\n')
#     f.write('#SBATCH -e '+ blockresultdir + 'shallow_'+str(i+1)+'n.err\n')
#     f.write('#SBATCH --nodes=1 \n')
#     f.write('#SBATCH --ntasks=1 \n')
#     f.write('#SBATCH --tasks-per-node=1 \n')
#     f.write('#SBATCH --cpus-per-task=1 \n')
#     f.write('#SBATCH --get-user-env \n')
#     f.write('#SBATCH -t 00:10:00 \n')
#     f.write('#SBATCH --mem-per-cpu=1000 \n')
#     f.write('#SBATCH --partition=cs5220 \n')
#     f.write('\n')
#     f.write('source /etc/profile.d/modules.sh \n')
#     f.write('#module load openmpi-4.0.0 \n')
#     f.write('cd ' + mydir + ' \n')
#     f.write('src/lshallow tests.lua dam ' + str(grid[i]) + '\n')
#     f.write('\n')

# generate tuning scripts
grid = [400, 800, 1600]
batch = [1, 2, 4, 6, 8]
block = [32, 64, 128, 256, 512, 1024]
for i in range(0, len(grid)):
    for j in range(0, len(block)):
        if grid[i] > block[j]:
            numblock = int(np.ceil(grid[i]/block[j]))
            mpin = numblock * numblock
            numnode = numblock
            numtasks = mpin
            taskpernode = numblock
            for k in range(0, len(batch)):
                filename = mpiscriptdir + "grid"+str(grid[i])+"_block"+str(block[j])+"_batch"+str(batch[k])+".sub"
                f = open(filename, "w+")
                f.write('#!/bin/bash\n')
                f.write('#SBATCH -J shallow\n')
                f.write('#SBATCH -o '+ mpmpiresultdir + 'grid'+str(grid[i])+'_block'+str(block[j])+'_batch'+str(batch[k])+'.out\n')
                f.write('#SBATCH -e '+ mpiresultdir + 'grid'+str(grid[i])+'_block'+str(block[j])+'_batch'+str(batch[k])+'.err\n')
                f.write('#SBATCH --nodes='+str(numnode)+' \n')
                f.write('#SBATCH --ntasks='+str(numtasks)+' \n')
                f.write('#SBATCH --tasks-per-node='+str(taskpernode)+' \n')
                f.write('#SBATCH --cpus-per-task=6 \n')
                f.write('#SBATCH --get-user-env \n')
                f.write('#SBATCH -t 00:10:00 \n')
                f.write('#SBATCH --mem-per-cpu=1000 \n')
                f.write('#SBATCH --partition=cs5220 \n')
                f.write('\n')
                f.write('source /etc/profile.d/modules.sh \n')
                f.write('#module load openmpi-4.0.0 \n')
                f.write('cd ' + mydir + ' \n')
                f.write('mpirun -n ' + str(mpin) + ' src/lshallow tests.lua dam ' + str(grid[i]) + ' ' + str(batch[k]) + ' ' + str(block[j]) + ' ' + str(block[j]) + '\n')
                f.write('\n')
