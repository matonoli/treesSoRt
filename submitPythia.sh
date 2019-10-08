#!/bin/bash
echo "How many processes?"
read npr
echo "How many events per process?"
read nev
for ((i = 0 ; i < $npr ; i++)); do
	#make makeTreeSoRt; ./makeTreeSoRt.exe $nev "output_${i}.root" 0 > "output_${i}.log" &
	make makeTreeSoRt; taskset --cpu-list ${i} ./makeTreeSoRt.exe $nev "output_${i}.root" 0 > "output_${i}.log" &
	echo "finished " ${i}
    pids[${i}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done
echo "Finished with pythia jobs, merging now"
touch output.log
for ((i = 0 ; i < $npr ; i++)); do
	echo "output_${i}.log" >> output.log 
	cat output_${i}.log >> output.log
	rm output_${i}.log
done

#hadd output.root output_*.root
#rm output_*.root
