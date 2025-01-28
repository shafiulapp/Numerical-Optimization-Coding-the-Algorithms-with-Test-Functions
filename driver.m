clc; 
clear;
mthd = algorithms();
x0 = [-117; 89];

mthd.run(x0,"rosenbrock", "SDLS", true);
 
mthd.plot_sum.plot_metrics();

% "SDLS" "CG" "NewtonLS" "TRS" "QN" "IQN"
%"rastrigin" "goldsteinprice" "ackley" "rosenbrock" "beale" 
% "booth" "matyas" "himmelblau" "mccormick" "schaffer"