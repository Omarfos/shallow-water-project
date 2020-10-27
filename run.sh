#!/bin/bash

strongdir=script/strong
weakdir=script/weak

# strong scaling analysis (--exclusive flag avoid jobs colliding on nodes)
sbatch --exclusive $strongdir/shallow1.sub
sbatch --exclusive $strongdir/shallow2.sub
sbatch --exclusive $strongdir/shallow3.sub
sbatch --exclusive $strongdir/shallow4.sub
sbatch --exclusive $strongdir/shallow5.sub
sbatch --exclusive $strongdir/shallow6.sub
sbatch --exclusive $strongdir/shallow7.sub

# batch run weak scaling analysis
# sbatch --exclusive $weakdir/shallow1_weak.sub
# sbatch --exclusive $weakdir/shallow2_weak.sub
# sbatch --exclusive $weakdir/shallow3_weak.sub
# sbatch --exclusive $weakdir/shallow4_weak.sub
# sbatch --exclusive $weakdir/shallow5_weak.sub
# sbatch --exclusive $weakdir/shallow6_weak.sub
# sbatch --exclusive $weakdir/shallow7_weak.sub
