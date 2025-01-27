numnode = 7

strongscriptdir = 'strong/'
weakscriptdir = 'weak/'
strongresultdir = 'results/baseline/strong/'
weakresultdir = 'results/baseline/weak/'
mydir = '$HOME/projects/CS5220/shallow-water-project'

# generate strong scaling .sub files
for i in range(0, numnode):
    filename = strongscriptdir + "shallow"+str(i+1)+".sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J shallow\n')
    f.write('#SBATCH -o '+ strongresultdir + 'shallow_'+str(i+1)+'n.out\n')        
    f.write('#SBATCH -e '+ strongresultdir + 'shallow_'+str(i+1)+'n.err\n')
    f.write('#SBATCH --nodes=' + str(i+1) + ' \n')
    f.write('#SBATCH --ntasks=' + str(i+1) + ' \n')
    f.write('#SBATCH --tasks-per-node=1 \n')
    f.write('#SBATCH --cpus-per-task=1 \n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    f.write('\n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('#module load openmpi-4.0.0 \n')
    f.write('cd ' + mydir + ' \n')
    f.write('src/lshallow tests.lua dam 1000 \n')
    f.write('\n')

# generate weak scaling .sub files
for i in range(0, numnode):
    filename = weakscriptdir + "shallow"+str(i+1)+"_weak.sub"
    f = open(filename, "w+")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -J shallow\n')
    f.write('#SBATCH -o '+ weakresultdir + 'shallow_'+str(i+1)+'n.out\n')        
    f.write('#SBATCH -e '+ weakresultdir + 'shallow_'+str(i+1)+'n.err\n')
    f.write('#SBATCH --nodes=' + str(i+1) + ' \n')
    f.write('#SBATCH --ntasks=' + str(i+1) + ' \n')
    f.write('#SBATCH --tasks-per-node=1 \n')
    f.write('#SBATCH --cpus-per-task=1 \n')
    f.write('#SBATCH --get-user-env \n')
    f.write('#SBATCH -t 00:10:00 \n')
    f.write('#SBATCH --mem-per-cpu=1000 \n')
    f.write('#SBATCH --partition=cs5220 \n')
    f.write('\n')
    f.write('source /etc/profile.d/modules.sh \n')
    f.write('#module load openmpi-4.0.0 \n')
    f.write('cd ' + mydir + ' \n')
    f.write('src/lshallow tests.lua dam ' + str((i+1)*200) + ' \n')
    f.write('\n')
