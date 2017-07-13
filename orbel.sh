#!/usr/bin/env bash

####Ler o arquivo
# Constants
gm=39.476926421373015 # Factor solar mass
opt="el" # Given orbital elements
ialpha=-1 # Elliptical orbit

while IFS='\n' read -r line || [[ -n "$line" ]]
do
    a=$(echo $line | cut -d" " -f1)
    e=$(echo $line | cut -d" " -f2)
    inc=$(echo $line | cut -d" " -f3)
    capom=$(echo $line | cut -d" " -f4)
    omega=$(echo $line | cut -d" " -f5)
    M=$(echo $line | cut -d" " -f6)

    printf "$opt \n$gm $ialpha $a $e $inc $capom $omega $capm $M\n" | ./orbel.x >> temp.txt
done < $1

mv temp.txt "$(echo $1 | cut -d"_" -f1-2)_xv.txt"


####Calcular posição e velocidade
####Escrever os elementos acrescidos da posição e velocidade