clc;
clear all;

%% initialize variables;
global Ybus nbus sys_base

sys_base=100;%MVA
run data_ne;
nbus=39;
% run d9bus;
% nbus=9;
%% solve ybus
Ybus=ybus;

%% solve loadflow

[Vg,thg,Pg,Qg,Pl,Ql,PV,PQ,Vbus,theta]=loadflow;
% theta - bus angle in rad
% Vbus - bus voltage magnitude in per unit

%% Get A, B, C matrix
[A, B] = linearization(Pg,Qg,Vg,thg,Vbus,theta,PV,PQ,Pl,Ql);

%% Starting invariant zero analysis

