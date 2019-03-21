#script="membit.py -f tmp_10000.pdb -simplethickness -n index.ndx"
#script="membit.py -f tmp_1000.pdb -thickness 10 0.1 -n index.ndx"
#script="membit.py -f tmp_1000.pdb -insertion closest -n index.ndx"
#script="membit.py -f tmp_1000.pdb -insertion average -n index.ndx"
script="membit.py -f tmp_1000.pdb -insertion 6 0 0 6 noNaN -n index.ndx"

time python $script 

#time python -m cProfile -o profile.profl -s cumulative $script 
#gprof2dot -f pstats profile.profl | dot -Tsvg -o profile.svg
#google-chrome profile.svg &
