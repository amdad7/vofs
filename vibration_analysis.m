clear all
clc
%% VOFS final project work 2020 %%
%% hinged-hinged-beam           18NA10002   /   AMDAD MUHAMMED PV 
%%
% length = room number 
% breadth = last digit of roll number 
% depth = half of breadth

%% DECLARING LENGTH, BREADTH, DEPTH AND DRAFT OF THE BEAM

L = 109;   % length of the beam      FREE AND SLIDING
B = 2;    % breadth of the beam
D = B/2;    % depth of the beam
d = B/4;    % draft of the beam


%% DECLARING THE DENSITIES AND MODULUS OF ELASTICITY

density_steel = 7850; % Kg/m3
Elasticity = 209*10^9; % Pa
density_water = 1023.6; % Kg/m3    
g = 9.81; % m/s2

%% FINDING THE THICKNESS OF THE BEAM

% declaring variable 't' to calculate the thickness of the beam.
syms('t')

% Quadratic equation to find the thickness.
eqn = 4*t^2 - t*(2*B + 2*D) + ((B*d*density_water)/(density_steel));

% Solving the equation to find the thickness.
t = solve(eqn);

% here we are using only suitable values of thickness
t = vpa(t(t<(D/2)))  % thickness
thick=t
%% FINDING MASS AND MOMENT OF INERTIA OF THE BEAM

% mass of beam in kg is given by mass = density of steel*volume of beam
Mass_of_beam = density_steel*L*((B*D) - (B - 2*t)*(D - 2*t))

% moment of inertia in m^4 is given by the equation M = (B*D^3 - B'*D'^3)/12
% where B' is the inside breadth
% where D' is the inside depth
Moment_of_inertia = ((B*D^3) - (B - 2*t)*((D - 2*t)^3))/12




%% FINDING FREQUENCY_PARAMETERS AND WAVE NUMBER
syms x
% a was found using the conditions given for hinged-hinged beam 
a = [1,0,1,0;-1,0,1,0;cos(x),sin(x),cosh(x),sinh(x);-cos(x),-sin(x),cosh(x),sinh(x)]
d=det(a)
%solving===>> -4*sin(x)*sinh(x)=0

% that is sin(x)=0 where x=(beta)*L

% first 10 frequency parameters ((beta)*L) are [n*pi] where n E (1:10)
   
% first 10 wave_numbers is [n*pi/L] where n E (1:10)
%% FINDING CO-EFICIANTS AND PLOTING different MOdeshapes

%declaring the matrix
lis=[pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi,10*pi] %FREQ_PARAMETERS
%co-efficiants b=1,a=0,c=0,d-0 found solving 4 equations where we fix b=1
%and solve the rest
b=1
for j=lis
    
    a= 0;
    c= 0;
    d= 0;
    x=0:0.2:L;
    omega=a*cos(j*x/L)+b*sin(j*x/L)+c*cosh(j*x/L)+d*sinh(j*x/L);
    figure(1)  
    plot(x,omega),grid
    title('plot of first 10 modeshapes')
    xlabel('position(m)')
    ylabel('deflection(m)')
    hold on
end
legend({'first','second','third','fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth', 'tenth'},'Location','southeast');
hold off

%% FINDING DRY NATURAL FREQUENCY

syms x
L=109        %length
         %(beta)*L
E=209*(10^9) %elasticity
m=Mass_of_beam   %mass
I=Moment_of_inertia 
dlis=[]
 %   syms x
 %  num =E*I*(j^4/L^4)*(sin(j*x/L))^2 ;  % PHI(X)=sin(ix/L)
 %  num = int(exp,x,0,209) ;
 %  num=vpa(num);
 % ex=m*I*((sin(j*x/L))^2) ;
 % den= int(ex,x,0,209) ;
 % den=vpa(den);
 % k=num/den;
 % k=sqrt(k)
 % wlis=[wlis k];
  % here we have our integration simplified to (E*I*(i^4))/(m*L^4) where i
  % is the various frequency parameters
 for i=lis
    val=(E*I*(i^4))/(m*L^4);
    val=sqrt(val);
    dlis=[dlis eval(val)];
end

dlis  %dry natural frequency list

%plot([1:10],dlis),grid


%% finding wet natural frequency
%here sin(x)^2 term will again be cancelled in this case
% so we will have to add (B(x)) in numerater 
% and shoud add 80% of total mass in den


wlis=[]
for i=lis
    val=((E*I*(i^4)/L^4)+(density_water*g*(B/4)*B))/(1.8*(m));
    val=sqrt(val);
    wlis=[wlis vpa(val)];
end
%%
wlis %wet natural frequency
figure(2)
plot([1:10],dlis),grid
hold on
plot([1:10],wlis),grid
xlabel('modeshape no')
ylabel('frequency')
legend({'drynatural frequency','wet natural frequency'},'Location','southeast');
hold off
%% orthoganility between two different modes %% And finding [GM] generalized inertia matrix
b=1
phi_lis=[];
lis=[pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi,10*pi] %FREQ_PARAMETERS // beta*L;
for j=lis
    syms x
    a= 0;
    c= 0;
    d= 0;
    omega=a*cos(j*x/L)+b*sin(j*x/L)+c*cosh(j*x/L)+d*sinh(j*x/L);
    phi_lis=[phi_lis omega];

end
phi_lis
inertia_matrix=[]
for i=phi_lis
    res=[];
    for j=phi_lis
        syms x 

        ans= Mass_of_beam*int(i*j,0,L);
        ans=eval(ans);
        res=[res ans]; 
    end
    inertia_matrix=[inertia_matrix;res];
end        
inertia_matrix %here re1 is the generalized inertia matrix (10*10)   
%% generalized stiffness matrix  [GK]

dlis=[]    %creating a list 4th order derivatives of all 10 phi function
for i=phi_lis
    syms x
    fn=diff(i,4);
    dlis=[dlis fn];
end

%taking combination of 10 phi's to each other and integrating 

stiffness_matrix=[]   
for i=phi_lis
    res=[];
    for j=dlis
        syms x
        fn= Moment_of_inertia*Elasticity*int(i*j,0,L);
        fn=eval(fn);
        res=[res fn];
    end
    stiffness_matrix=[stiffness_matrix;res];
end
stiffness_matrix

%% meshing GM
tlis=[]
for i=1:10
   for j=1:10
       inertia_matrix(i,j);
       tlis=[tlis inertia_matrix(i,j)];
   end
end 
tlis

lis2=[]
lis3=[]
for i=1:10
    for j=1:10
        lis2=[lis2 i];
    end
end    
for i=1:10
    for j=1:10
        lis3=[lis3 j];
    end
end
lis2
lis3
figure(3)
plot3(lis2,lis3,tlis)
zlabel('inertia')
title('GENERAL INERTIA MATRIX')

%% meshing stiffness matrix

tlis1=[]
for i=1:10
   for j=1:10
       stiffness_matrix(i,j);
       tlis1=[tlis1 stiffness_matrix(i,j)];
   end
end 
tlis

lis2=[]
lis3=[]
for i=1:10
    for j=1:10
        lis2=[lis2 i];
    end
end    
for i=1:10
    for j=1:10
        lis3=[lis3 j];
    end
end
lis2
lis3
figure(4)
plot3(lis2,lis3,tlis1)
zlabel('inertia')
title('GENERAL STIFNESS MATRIX')



%% finding temporal part #time

%phi0=0   %since initial velocity is 0
%A0=0    &since initial velocit=0 and initial displacement is 1
syms t
fvp=0
A0=0
tp=[]
for i=1:10
    syms x
    syms t
    tw = sqrt(stiffness_matrix(i,i)/inertia_matrix(i,i));
    tw=vpa(tw);
    expp=phi_lis(1)*phi_lis(i)
    A0=(2/L)*int(expp,0,L)
    sp=A0*cos(tw*t) ; 
    tp=[tp sp];
    fvp=fvp+(sp*phi_lis(i))    ;
end
fvp
%% PLOTING TEMPORAL PART OF UNDAMPED
for i=1:10
    syms t
    t=0:0.2:50
    figure(60)
    plot(t,eval(tp(i)))
    xlabel('time')
    ylabel('amp')
    title('temporal part of undamped')
    hold on
end
hold off
 %% MESHING DEFLECTION AS POSITION IN x AXIS AND AND TIME IN Y AXIS 
% % FREE VIBRATION DEFLECTION
% fvp
% X=0:2:L;
% Y=X;
% [x,t]=meshgrid(X)
% f=cos(0.031423048338185370631503445792987.*t)*sin((pi.*x)/109) + cos(0.50276877341096593010405513268779.*t)*sin((2*pi.*x)/109) + cos(8.0443003745754548816648821230046.*t)*sin((4*pi.*x)/109) + cos(128.70880599320727810663811396807.*t)*sin((8*pi.*x)/109) + cos(75.446739059983073616422188933939.*t)*sin((7*pi.*x)/109) + cos(19.639405211365858150429630768485.*t)*sin((5*pi.*x)/109) + cos(314.23048338185373040687409229577.*t)*sin((10*pi.*x)/109) + cos(2.5452669153930149725795217818813*t)*sin((3*pi*x)/109) + cos(40.724270646288239561272348510101.*t)*sin((6*pi.*x)/109) + sin((9*pi.*x)/109)*cos(206.16662014683419101856998167932.*t)
% surf(x,t,f)
% xlabel('displacement')
% ylabel('time')
% zlabel('deflection')

 %% MESHING DEFLECTION AS POSITION IN x AXIS AND AND TIME IN Y AXIS  
exttt=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(fvp,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    exttt=[exttt;extt];
end
exttt
%%
x=0:1:L;
y=x;
[X,Y]=meshgrid(x)
figure(7)
surf(X,Y,exttt)
xlabel("time")
ylabel("position")
zlabel("deflection")
title("undamped free vibration(H-H)")
hold off

%% Meshing INERTIA MATRIX 
x1=1:1:10
y1=x1
[X,Y]=meshgrid(x1)
figure(5)
surf(X,Y,inertia_matrix)
zlabel("Generalised Mass")
xlabel("row")
ylabel("column")
title("Mesh of generalised Mass Matrix")
hold off
%% MESHING STIFFNESS MATRIX
x1=1:1:10
y1=x1
[X,Y]=meshgrid(x1)
figure(6)
surf(X,Y,stiffness_matrix)
zlabel("Generalised Stiffness")
xlabel("row")
ylabel("column")
title("Mesh of generalised Stiffness Matrix")
hold off

%% temporal part of damped vibration

%mx"+ck'+kx=0
%c=0.05*sqrt(m*k)
im=inertia_matrix
sm=stiffness_matrix
xo=1
vo=0
tpdlis=[]
for i=1:10
    syms t
    mt=im(i,i);
    kt=sm(i,i);
    c=0.05*2*sqrt(mt*kt);
    eta=0.05/2;   %c/sqrt(k*M) c=0.05*sqrt(k*m)
    wn=sqrt(kt/mt);
    wd=(wn)*sqrt(1-eta^2);
    expp=phi_lis(1)*phi_lis(i)
    xo=(2/L)*int(expp,0,L)
    fn=exp(-c*t/(2*mt))*((xo*cos(wd*t))+((vo+2*eta*wn*xo)*(sin(wd*t))/wd));
    tpdlis=[tpdlis fn];
end    
tpdlis

fuun=0
for i=1:10
    syms t
    fuun=fuun+tpdlis(i)*phi_lis(i);
end
fuun
% %% ploting temporal part of damped
% syms t
% eval(t)
% t=1:0.5:50
% plot(t,subs(fuun)

%% mesh
% X=0:2:20;
% Y=X;
% [x,t]=meshgrid(X)
% surf(x,t,f)
% %set(gca,'ztick',0:0.1:5)
% xlabel('displacement')
% ylabel('time')
% zlabel('deflection')
%% meshing deflection of damped vibration
extt4=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(fuun,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    extt4=[extt4;extt];
end
extt4
%%
x=0:1:L;
y=0:1:L;
[X,Y]=meshgrid(x,y)
figure(8)
surf(X,Y,extt4)
xlabel("time")
ylabel("position")
zlabel("deflection")
title("damped free vibration(H-H)")
hold off
%% FORCED VIBRATION    % Froude-Krylov forcing, obeying the deep water dispersion relation.
 
%generalized forcing equation and froude-krylov pressure and force per
 %unitlength

 syms t x
k4 = 2*pi/L;
g=9.8
force = ((2*density_water*g)/k4)*.14*L*exp(-k4*D)*sin((k4*B)/2)*cos(sqrt(g*k4)*t) %per_unit_lenfth
gen_force = [];
for i = 1:10
    syms x t
    gen_force = [gen_force int(phi_lis(i)*force, x, 0, L)];
end
gen_force
 %% forcing amplitude
Foo=((2*density_water*g)/k4)*.14*L*exp(-k4*D)*sin((k4*B)/2)
lamda=L

%% generalized forcing plot 

for i=1:10
    syms t
    t=0:0.2:60;
    figure(14)
    plot(t,eval(gen_force(i)));
    xlabel('t')
    ylabel('force')
    hold on
end
legend({'first','second','third','fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth', 'tenth'},'Location','southeast');
hold off
%% FORCED UNDAMPED vibration
%finding complimentory part

cp=[]
g=9.8
cpt=[]
% foo=2*ro*(9.8/k4)*amp*(exp(-k4*B/2))*(sin(k4*B/2))*L
for i=1:10
    syms x t
    wf=sqrt((stiffness_matrix(i,i))/inertia_matrix(i,i));
    m=inertia_matrix(i,i);
    Afo=gen_force(i)/(m*((wf^2)-(g*k4)));
    cpt=[cpt Afo];
    cp=[cp (Afo*phi_lis(i))];
end
cp

%% forced undamped
syms t x
ps=fvp  %purticular solution
tot_for_def=ps+sum(cp)

%% meshing deflection of forced undamped vibration as afunction of space and time
extt2=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(tot_for_def,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    extt2=[extt2;extt];
end
extt2
%%
x=0:1:L;
y=0:1:L;
[X,Y]=meshgrid(x,y)
figure(10)
surf(X,Y,extt2)
xlabel("time")
ylabel("position")
zlabel("deflection")
title("undamped forced vibration(H-H)")
hold off
%% ploting principal_cordinates as a funtion of time including(complementory part)
for i=1:10
    syms t
    t=0:0.2:60
    figure(12)
    plot(t,eval(cpt(i)+tp(i)))
    xlabel('time')
    ylabel('force')
    hold on
end 
legend({'first','second','third','fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth', 'tenth'},'Location','southeast');
hold off
%% BENDING_MOMENT shear_force
syms x t
bending_moment=Elasticity*Moment_of_inertia*diff(tot_for_def,x,2)
shear_force=Elasticity*Moment_of_inertia*diff(tot_for_def,x,3)
%% Meshing bending_moment 
extt3=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(bending_moment,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    extt3=[extt3;extt];
end
extt3
%%
x=0:1:L;
y=0:1:L;
[X,Y]=meshgrid(x,y)
figure(30)
surf(X,Y,extt3)
xlabel("time")
ylabel("position")
zlabel("bending moment")
title("bending moment(H-H)")
hold off
%% MESHING SHEAR FORCE

extt4=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(shear_force,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    extt4=[extt4;extt];
end
extt4
%%
x=0:1:L;
y=0:1:L;
[X,Y]=meshgrid(x,y)
figure(31)
surf(X,Y,extt4)
xlabel("time")
ylabel("position")
zlabel("shear force")
title("shear force(H-H)")
hold off
%% bending_stress 

bending_stress=(bending_moment*(B/4))/Moment_of_inertia

%% meshing bending stress

extt5=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(bending_stress,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    extt5=[extt5;extt];
end
extt5
%% 
x=0:1:L;
y=0:1:L;
[X,Y]=meshgrid(x,y)
figure(31)
surf(X,Y,extt5)
xlabel("time")
ylabel("position")
zlabel("deflection")
title("bending-stress(H-H)")
hold off

%% SHEAR STRESS

shear_stress=(3/2)*(shear_force/((B*B/2)-((B-eval(thick(1,1)*2))*(B/2-eval(thick(1,1)*2)))))



%% meshing bending stress

extt6=[]
for i=0:1:L
    extt=[];
    for j=0:1:L
        syms x t
        exx=subs(shear_stress,{x,t},{i,j});
        exx=eval(exx);
        extt=[extt exx];
    end
    extt6=[extt6;extt];
end
extt6
%% 
x=0:1:L;
y=0:1:L;
[X,Y]=meshgrid(x,y)
figure(32)
surf(X,Y,extt5)
xlabel("time")
ylabel("position")
zlabel("deflection")
title("shear-stress(H-H)")
hold off
%%  finding maxx
%found all using excel

%maximum shear force=>>2138132.0185734
%location=>109m time=>0s

%maximum bending moment=>>[47639885.8727914]
%location=>55m time=>93s

%maximum bending moment=>979885685.1
%location=>55m time 93s

