#!/bin/bash

#region=fuentes
region=dubai

for year in {1979..2019}
do
   ./test_dni.x $region $year
done
