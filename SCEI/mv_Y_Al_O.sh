#!/bin/bash


#Delet old Al O in 01YY file
rm 01Y/MO_COEFF_Al.txt
rm 01Y/MO_ENERGY_Al.txt
rm 01Y/MO_COEFF_O.txt
rm 01Y/MO_ENERGY_O.txt

#Delet old Y O in 02Al file
rm 02Al/MO_COEFF_Y.txt
rm 02Al/MO_ENERGY_Y.txt
rm 02Al/MO_COEFF_O.txt
rm 02Al/MO_ENERGY_O.txt

#Delet old Y Al in 03O file
rm 03O/MO_COEFF_Y.txt
rm 03O/MO_ENERGY_Y.txt
rm 03O/MO_COEFF_Al.txt
rm 03O/MO_ENERGY_Al.txt

#cp Y file to 02Al 03O
cp 01Y/MO_COEFF_Y.txt 02Al/
cp 01Y/MO_ENERGY_Y.txt 02Al/
cp 01Y/MO_COEFF_Y.txt 03O/
cp 01Y/MO_ENERGY_Y.txt 03O/

#cp Al file to 01Y 03O
cp 02Al/MO_COEFF_Al.txt 01Y/
cp 02Al/MO_ENERGY_Al.txt 01Y/
cp 02Al/MO_COEFF_Al.txt 03O/
cp 02Al/MO_ENERGY_Al.txt 03O/

#cp O file to 01Y 02Al
cp 03O/MO_COEFF_O.txt 01Y/
cp 03O/MO_ENERGY_O.txt 01Y/
cp 03O/MO_COEFF_O.txt 02Al/
cp 03O/MO_ENERGY_O.txt 02Al/