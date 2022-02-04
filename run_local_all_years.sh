#!/bin/bash

for year in {1979..2019}
do
   ./test_dni.x fuentes $year
done
