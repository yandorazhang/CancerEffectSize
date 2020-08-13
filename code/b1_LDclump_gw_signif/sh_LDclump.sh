for trait in {1..47}
do
qsub -N clump -cwd -l mem_free=20G,h_vmem=20G,h_fsize=200g,chatterjee LDclump.sh $trait
done

