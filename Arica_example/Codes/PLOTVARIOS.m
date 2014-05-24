%PLOT EXPERIMENTAL VAROGRAMS

%% IMPORT DATA

clear all

v1=importdata('vario1-10.txt');
v2=importdata('vario110.txt');
v3=importdata('vario100.txt');
v4=importdata('vario010.txt');

%% PLOTS

figure
hold on

plot(v1(:,2),v1(:,3),'r','LineWidth',1.1)
plot(v2(:,2),v2(:,3),'b','LineWidth',1)
plot(v3(:,2),v3(:,3),'g','LineWidth',1)
plot(v4(:,2),v4(:,3),'k','LineWidth',1)