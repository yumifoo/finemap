trait=$1

for i in {1..23}
do

    grun -n "$trait"_ld_"$i" -M 50 bash runLDstore.sh $trait $i

done
