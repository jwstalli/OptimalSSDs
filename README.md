# OptimalSSDs
Optimal supersaturated designs with continuous factors using heuristic-initiated lasso sieve.

SIC_Paper_Function_Library.R - which is the library for all of the R functions needed to run the new CS measure code. 

new_cs_measure_demo.R - this is the code for the demonstrations to motivate the move to the new measure. It loads in the function library and some other competing designs. I showed 3 scenarios in there and added a lot of comments so you all can tweak it if needed. 

new_cs_measure_EPCEA_eval.R - this is the script that does the EPCEA algorithm that comes up with designs on the pareto front for n=9, p=10 all positive signs. However, I commented there that these scenarios can be changed. I also added some code to compare the designs on the pareto front as well as some other added in designs. 

In the last two R scripts, you will need to update the filepaths of the designs read in, and the function library. 

I have also added in the designs that are read in, and compared. There are two "d1.txt" designs, but they have different dimensions to cover different scenarios. 
