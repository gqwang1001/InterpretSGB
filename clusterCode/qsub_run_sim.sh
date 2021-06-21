for i in `seq 1 1`;do
  for j in `seq 1 7`;do 
      qsub -cwd run_gw_sim.sh $i $j
  done
done