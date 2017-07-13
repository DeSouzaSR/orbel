#!/usr/bin/env bash
clear
echo "ORBEL PROGRAM"
echo "============="
echo

echo "Choose an option"
echo "1 - Orbital Elements to Position and Velocity"
echo "2 - Position and Velocity to Orbital Elements"
echo "3 - Exit"
read opt

if [ $opt == 1 ]
then
    echo
    echo "Enter the orbital elements in this order"
    echo "1) Mass factor - GM"
    echo "2) Type orbit:"
    echo "  - Ellipse   = -1"
    echo "  - Parable   =  0"
    echo "  - Hyperbole = -1"
    echo "3) Eccentricity"
    echo "4) Inclination [deg]"
    echo "5) Long. Nodo - Omega [deg]"
    echo "6) Arg. Periaster - omega [deg]"
    echo "7) Mean anomaly - M [deg]"
    echo 
    echo "Presse Ctrl + D to finish"
    cat > "el.dat" 

    # Adding "el" in the beginner 
    echo "el" | cat - el.dat > temp && mv temp el.dat

    echo
    echo "Orbital elements"
    echo "----------------"
    cat el.dat

    echo
    echo "Position and velocity"
    echo "---------------------"
    ./orbel.x < el.dat

    rm el.dat 
elif [ $opt == 2 ]
then
    echo
    echo "Enter the position and velocity in this order"
    echo "1) Mass factor - GM"
    echo "2) x"
    echo "3) y"
    echo "4) z"
    echo "5) vx"
    echo "6) vy"
    echo "7) vz"
    echo 
    echo "Presse Ctrl + D to finish"
    cat > "xv.dat"

    # Adding "el" in the beginner 
    echo "xv" | cat - xv.dat > temp && mv temp xv.dat

    echo
    echo "Position and velocity"
    echo "---------------------"
    cat xv.dat

    echo
    echo "Orbital elements"
    echo "---------------------"
    ./orbel.x < xv.dat

    rm xv.dat
else 
    exit     
fi 