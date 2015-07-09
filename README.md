# DispersalModel
Contains two versions of the most recent larval dispersal model. One version is for static PLD's, the other for variable (plastic) PLD's. 

# Sample command
```
python DispersalModelPlastic.py -l HIreefsNew.csv -c files_HA8km.csv -d ./hawaii8km_/ -u _100m_u.txt -v _100m_v.txt -o ./out/Dispersal_python_8kmHA_test_site_ -t 2 -m 5 -s 122 -n 130 -r 5 -z 687 -a 438 -b 251 -f 175 -g 210 -j 15 -k 35 -q 0.08 -i 1

python DispersalModelPlastic_MT.py -l HIreefsNew.csv -c files_HA8km.csv -d ./hawaii8km_/ -u _100m_u.txt -v _100m_v.txt -o ./out/Dispersal_python_8kmHA_test_site_ -t 2 -m 5 -s 122 -n 130 -r 5 -z 687 -a 438 -b 251 -f 175 -g 210 -j 15 -k 35 -q 0.08 -i 1 -p 2
```
##Changes to parameters for interpolation code
----------------------------------------------

Many of the parameters are the same in the previous version of the code, and the new version of code.  Only a few parameters have changed or been added, which we will go over.

###-i, --total_steps
The -i parameter is a new parameter.  This is what tells the application how many ways to split the 24 hour period.  For example, if we wanted to have our time steps be every 6 hours, we would use the following: 
```
-i 4 
```
This is determined by 24/4 = 6

If for example, you wanted to have each step be 1 second, you would use -i 86400

---------------------
###-s & -n
These two paramters referred to startarray and endarray.  Due to how the input files are provided in day steps, these two parameters must still reflect julian days.  This limits the simulation to beginning on day boundaries.  Other other calculations for the application to run, (days vs. steps based on -i) are computed internally.

---------------------
###-t & -m
These two parameters, must now reflect the number of steps you want to take, NOT DAYS.  This means, if I had previously used -t 2 and -m 5 for the old code that worked in day steps, but now was moving to try and use 6 hr steps, I would need to multiply the -t and -m values by the value used in -i.  For example, these two would be equivalent, where one uses time steps of a day, while the other users a time step fo 6 hrs:
```      	  	      
-i 1 -t 2 -m 5 (24hrs per step)
-i 4 -t 8 -m 20 (6hrs per step)
```

---------------------
###-p (_MT version only)
In the multi-threaded version of the code, the user may want to specify the number of threads to use.  To provide the user with this option, we supply a command line parameter, -p .  

