close all
clear
clc

rng(0)
p1 = rand(1000, 300);
p2 = rand(1000, 300);

mu_g1 = mean(p1);
mu_g2 = mean(p2);
s_g1 = cov(p1);
s_g2 = cov(p2);

