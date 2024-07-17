%% clear
clear
close all
clc

%% function
func=@(x) 2*x(1).^2-1.05*x(1).^4+(x(1).^6)/6+x(1)*x(2)+x(2)^2;
X_0_best=[3;3];
X_lb=[-5;-5];
X_ub=[5;5];

N_particles=20;
max_iter=2000;
func_tol=.00001
x_tol=.00001
weight_mat=.1*ones(size(X_lb));
phi_g=2.5
phi_p=2.25
[x_best,f_best,counter,f_conv_history,previous_iters,previous_f]=particleSwarmWithRestarts(func,X_0_best,X_lb,X_ub,N_particles,max_iter,func_tol,x_tol,weight_mat,phi_g,phi_p);

%% plotting
x=linspace(-5,5,100);
y=linspace(-5,5,100);
[X,Y]=meshgrid(x,y);

Z=2*X.^2-1.05*X.^4+(X.^6)/6+X.*Y+Y.^2;
meshc(X,Y,log10(Z));
hold on
colormap jet
colorbar