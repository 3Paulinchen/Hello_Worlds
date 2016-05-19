function[C,n,h,x,y]=IMEX_Central_Reaktion()
%%  Isotrop mit Reakion (IMEX)

% Number of steps
Steps=10;
% Time step
tau=0.01;
% Length
L=0.002;
% Velocity solid (constant)


% Time step
%delta_t= 1.2500e-06;

% Diffusion coefficient
diffusion_0=1*10^-10;

% F�r Geschwindigkeit
%Viskosit�t
nu=10;
% Permeabilit�t
%k_0=1*10^-15;
k_0=7.5*10^-16;
% Fl�che in Meter
area=10^-12;
% Konstante M
konstant2=4.3;
% Void ratio = Porosit�t
phi_0=0.8;
%Dehnung �ber Zeit tr(e)
e=linspace(0,0.1,Steps);
% Druck gradient
delta_p=linspace(0,2*10^6,Steps);
delta_p_Abstieg=linspace(2*10^6,0,Steps);
% Anstieg
e_Anstieg=linspace(0.1,0,Steps);

a=10^-10;
%a=10^-4;
%Setting up the grid

h=L/50;
x = 0:h:L;
y = 0:h:L;
n=length(x);
ndgrid(n,n); 


%Inital concentration c_0
c_0=1;
%concentration matrix
c=zeros(n^2,1);
c=reshape(c,[n n]);

% (initial concentration left 
c(1,1:n)=c_0;
% and right
%c(n,1:n)=c_0;


% Reshape into vector
c=reshape(c,[n^2 1]);



Ident=speye(n^2,n^2); 
%% Belastung
r=1;
Zyklen=200;
while(r<Zyklen)  
for i=1:Steps
   
    % porosity at step
    phi=1-(1-phi_0)*(1+e(i));
    % velocity darcy law
    k=k_0*exp(konstant2*(phi-phi_0)/(1-phi));    
    k=k/nu;
    Q=k*area*delta_p(i)/L;
    q=Q/(phi*area);
   
    
    
    % diffusion
    diffusion=diffusion_0*exp(10*(phi-phi_0)/(phi*phi_0));
    
    Ident_n=speye(n,n);
      %% Building matrix for Laplace Operator
    T= (1/h^2)*(gallery('tridiag',n,1,-2,1));   
    A=kron(T,Ident_n)+kron(Ident_n,T);
   
   %% building gradient matrix (backwards scheme) upwind
      A1 = (1/(2*h))*(gallery('tridiag',n^2,-1,0,1));
      
    %% Reaction term matrix
%    f=exp(-10^2*x);
 %   f=a*f;
    R= speye(n,n);
 %   for i=1:n
 %       R(i,i)=f(i);
 %   end
    R(1,1)=0;
    R(n,n)=0;
    
    R=kron(R,Ident_n);
    

   % M= (I + tau*a*R) %% steif
     M=  Ident + a*tau/2*R  ;
  %  M= Ident+ (tau/2)*q*A1;
    N = Ident + tau*diffusion*A - tau*q*A1 - a*(tau/2)*R;
   %  N= Ident+ tau*diffusion*A- a*tau*R - (tau/2)*q*A1;
  %% Boundary conditions
  rhs=N*c;
  [c]=Boundary(M,rhs,n,c_0,h,q,a,Q);
     
end
%Static
for i=1:5*Steps
     
    % porosity at step
  %  phi=1-(1-phi_0)*(1+e(i));
    % velocity darcy law
  %  k=k_0*exp(konstant2*(phi-phi_0)/(1-phi));    
  %  k=k/nu;
  %  Q=k*area*delta_p(i)/L;
    q=Q/(phi*area);
  
    
    
    % diffusion
    diffusion=diffusion_0*exp(10*(phi-phi_0)/(phi*phi_0));
    
    Ident_n=speye(n,n);
      %% Building matrix for Laplace Operator
    T= (1/h^2)*(gallery('tridiag',n,1,-2,1));   
    A=kron(T,Ident_n)+kron(Ident_n,T);
   
   %% building gradient matrix (backwards scheme) upwind
      A1 = (1/(2*h))*(gallery('tridiag',n^2,-1,0,1));
      
    %% Reaction term matrix
%    f=exp(-10^2*x);
 %   f=a*f;
    R= speye(n,n);
 %   for i=1:n
 %       R(i,i)=f(i);
 %   end
    R(1,1)=0;
    R(n,n)=0;
    
    R=kron(R,Ident_n);
    

 % M= (I + tau*a*R) %% steif
     M=  Ident + a*tau/2*R  ;
  %  M= Ident+ (tau/2)*q*A1;
    N = Ident + tau*diffusion*A - tau*q*A1 - a*(tau/2)*R;
   %  N= Ident+ tau*diffusion*A- a*tau*R - (tau/2)*q*A1;
  %% Boundary conditions
  rhs=N*c;
  [c]=Boundary(M,rhs,n,c_0,h,q,a,Q);
    
end
%% Entlastung
for i=1:Steps
    % porosity at step
    phi=1-((1-phi_0)*(1+e_Anstieg(i)));
    % velocity darcy law
    k=k_0*exp(konstant2* (phi-phi_0)/(1-phi));    
    k= k/nu;
    Q=k*area*delta_p_Abstieg(i)/L;
    q=Q/(phi*area);
 
    % diffusion
    diffusion=diffusion_0*exp(10*(phi-phi_0)/(phi*phi_0));
    
     %% Building matrix for Laplace Operator
    T= (1/h^2)*(gallery('tridiag',n,1,-2,1));
    Ident_n=speye(n,n);
    A=kron(T,Ident_n)+kron(Ident_n,T);

   %% building gradient matrix (backwards scheme) centered
    A1 = (1/(2*h))*(gallery('tridiag',n^2,-1,0,1));
   %% downwind
   %  T1 = (gallery('tridiag',n,0,1-1));
    
%% Reaction term matrix
%    f=exp(-10^2*x);
 %   f=a*f;
    R= speye(n,n);
 %   for i=1:n
 %       R(i,i)=f(i);
 %   end
    R(1,1)=0;
    R(n,n)=0;
    
    R=kron(R,Ident_n);
    

 % M= (I + tau*a*R) %% steif
     M=  Ident + a*tau/2*R  ;
  %  M= Ident+ (tau/2)*q*A1;
    N = Ident + tau*diffusion*A - tau*q*A1 - a*(tau/2)*R;
   %  N= Ident+ tau*diffusion*A- a*tau*R - (tau/2)*q*A1;
  %% Boundary conditions
  rhs=N*c;
  [c]=Boundary(M,rhs,n,c_0,h,q,a,Q);
 
end
   
    r=r+1;
    r
end
   C=reshape(c,[n,n]);
   [X,Y] = meshgrid(x,y); 

%    Plotting
   figure();
   contourf(X,Y,C');
   xlabel('x','FontSize',14);ylabel('y','FontSize',14);zlabel('Con','FontSize',14); 
   title('IMEX, Isotrop (mit Reaktion)','FontSize',14);

   colorbar;
   % c=reshape(C,[n^2 1]);
end



