for iter in {1..47}
do
for pic in 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001
do
export iter; export pic;
sbatch fit_GENESIS.sh;
done
done


