#!/bin/bash


#Delet old Li Y in 01F file
rm 01F/MO_COEFF_Li.txt
rm 01F/MO_ENERGY_Li.txt
rm 01F/MO_COEFF_Y.txt
rm 01F/MO_ENERGY_Y.txt

#Delet old F Y in 02Li file
rm 02Li/MO_COEFF_F.txt
rm 02Li/MO_ENERGY_F.txt
rm 02Li/MO_COEFF_Y.txt
rm 02Li/MO_ENERGY_Y.txt

#Delet old F Li in 03Y file
rm 03Y/MO_COEFF_F.txt
rm 03Y/MO_ENERGY_F.txt
rm 03Y/MO_COEFF_Li.txt
rm 03Y/MO_ENERGY_Li.txt

#cp F file to 02Li 03Y
cp 01F/MO_COEFF_F.txt 02Li/
cp 01F/MO_ENERGY_F.txt 02Li/
cp 01F/MO_COEFF_F.txt 03Y/
cp 01F/MO_ENERGY_F.txt 03Y/

#cp Li file to 01F 03Y
cp 02Li/MO_COEFF_Li.txt 01F/
cp 02Li/MO_ENERGY_Li.txt 01F/
cp 02Li/MO_COEFF_Li.txt 03Y/
cp 02Li/MO_ENERGY_Li.txt 03Y/

#cp Y file to 01F 02Li
cp 03Y/MO_COEFF_Y.txt 01F/
cp 03Y/MO_ENERGY_Y.txt 01F/
cp 03Y/MO_COEFF_Y.txt 02Li/
cp 03Y/MO_ENERGY_Y.txt 02Li/