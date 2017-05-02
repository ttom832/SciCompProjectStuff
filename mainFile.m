clear
clc

Nodes=[10,20,50];
Lambda=1.4;

for p=1:length(Nodes)
[h,...
    X,...
    Y,...
    cycles,...
    percent_error,...
    u_gauss_seidel,...
    u_gauss_seidel_relaxed,...
    u_gauss_seidel_no_force,...
    u_gauss_seidel_no_force_relaxed]=project(Nodes(p),Lambda);


Cycles(p,:)=cycles;
Error(p,:)=percent_error;
H(p,1)=h;
end

disp('Step Size')
disp(H)
disp('Gauss Seidel , Gauss Seidel with Relaxation , Gauss Seidel no Forcing Function , Gauss Seidel with Relaxation no Forcing Function')
Cycles
Error
% figure(1)
% mesh(X,Y,u_gauss_seidel),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel')
% 
% figure(2)
% mesh(X,Y,u_gauss_seidel_relaxed),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel With Relaxation')
% 
% figure(3)
% mesh(X,Y,u_gauss_seidel_no_force),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel No Forcing Function')
% 
% figure(4)
% mesh(X,Y,u_gauss_seidel_no_force_relaxed),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel No Forcing Function With Relaxation')
% 
% figure(5)
% mesh(X,Y,u_gauss_seidel_relaxed-u_gauss_seidel),xlabel('x'),ylabel('y'),zlabel('u'),title('Difference in Gauss Seidel and Gauss Seidel Relaxed')
% 












