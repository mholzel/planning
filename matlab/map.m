clear;clc;close all;

a = csvread('../data/original.csv');
b = csvread('../data/refined.csv');

figure 
plot(a(:,1),a(:,2), 'x-')
axis off
saveas( gcf, 'original_waypoint_map.jpg' )
        
figure
plot(a(:,1),a(:,2), '-bx', 5 + b(:,1), b(:,2), 'ro-')    
axis off
xlim([150, 250])
ylim([1900, 2000])
saveas( gcf, 'original_vs_smooth_waypoint_map.jpg' )
