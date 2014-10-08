#/bin/sh
#$ -cwd

rm *aei
for ((i=1 ; i<=2 ; i++))
do 
cp spinp$i.out spinp$i.dat 
cp horb$i.out horb$i.dat 
cp spins.out spins.dat 
done
./element
