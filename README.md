# Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft
Collision Probability Calculation and Optimal Avoidance Control for Spacecraft.


## Introduction:
This code package implements the algorithms from the course report "Collision Probability Calculation and Optimal Avoidance Control for Spacecraft".


Before running this program please install gpops tool box first:
Run gpopsSetup.m to install gpops tool box to your matlab. (In gpops folder)

## Run the program. 
(In ECH267_Course_Project folder)

### Step 1:
Run Step1_prob.m.
### Result 1:
A plot of the collision probability before taking the strategy.
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/Prob1.svg">

### Step 2:
Run Step2_OptimalPosition.m
### Result 2:
No plots output.

### Step 3:
Run Step3_spacecraftMain.m
### Result 3:
Five plots. (1)State X values vs time. (2)State Y values vs time. (3)State Z values vs time. (4)State Rs values vs time.
(5)Control values vs time
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/StateX.svg">
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/StateY.svg">
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/StateZ.svg">
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/StateRs.svg">
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/Control.svg">

### Step 4:
Run Step4_ProbAfter.m
### Result 4:
A plot of the collision probability after taking the avoidance strategy.
<img src="https://github.com/ZiquanDeng/Collision-Probability-Calculation-and-Optimal-Avoidance-Control-for-Spacecraft/blob/main/ECH%20Plots/ProbAfter.svg">

