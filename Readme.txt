This package includes two python files and several data used in the python code

Python file:
1. flux_for_different_runoff_temp_slope_v3_sta.py

This is the file for the first part described in Experiments that can test several different climate and slope datasets to the model
--In this file, we use two different observation dataset to evaluate the model and "sta" variable is the one control the choice
--For the erosion rate revised case, we use "B" and "Be" variable to control and "B" is the erosion rate from TSS and "Be" is the erosion rate from both TSS and Be isotope.

2.new_test_parameter_update_Ec_LAI_global_Etotal.py

This is the file that we redo the parameter selection process and we use 72 cores in this python code
--It will output the .csv result of all the tests which include the basinal river flux and total silicate weathering flux.
--In this code, we don't exclude the possibility that the Ca2++Mg2+ concentration of sediment is higher than metamorphic but when we plot will exclude it.