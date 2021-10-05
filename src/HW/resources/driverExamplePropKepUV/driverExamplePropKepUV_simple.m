clear all
close all


%INITIAL CONDITIONS FOR S/C
ep=1e-5;        %just to make it 3D a bit
vmult=0.9;      %make 1 for near parab, >1 for hyperbola, <1 for hyperb
x0=[1 ep ep ep vmult*sqrt(2) ep]';

%PROPAGATION DETAILS
ni=100;          %number of segments to propagate
dUV=0.5         %for each segment, how much to advance the independent variable (universal variable like Eccentric anomaly)
dvMag=1e-3      %deltaV magnitude at each segment boundary (direction in same direction as current velocity for this demo)

orderCase=2;    %0 for state only propagation, 1 for state plus first order STM, 2 for first and sencond order STMs

bigx=zeros(7,ni);
bigx(1,1:6)=x0;
bigx(1,7)=0;

xnow=x0;

%call it first to intilize Fortran library only
[Yf,dYf,d2Yf,capZ,errFlag]=propKeplerUV( xnow,dUV,orderCase,true,false);

for i=1:ni-1
    
    %propagate plus partials using matlab .m file (no initialization needed)
    [Yf,dYf,d2Yf,capZ,errFlagm]=propKeplerUVmatlab( xnow,dUV,orderCase);
    
    %store outputs for later comparison to fortran call
    Yfm=Yf; dYfm=dYf; d2Yfm=d2Yf; capZm=capZ;
    
    %propagate plus partials using fortran library call
    [Yf,dYf,d2Yf,capZ,errFlag]=propKeplerUV( xnow,dUV,orderCase,false,false);
    
    %----------------------------------------------------------------------
    %store a metric to represent the difference between the fortran and
    %matlab calls for all the computed states/derivs
    qnorm(1)=norm(Yfm-Yf);
    if(orderCase>0)
        qnorm(2)=norm(dYfm-dYf);
    end
    
    if(orderCase>1)
        for qq=1:7
            qnorm(2+qq)=norm(d2Yfm(:,:,qq)-d2Yf(:,:,qq));
        end
    end
    normVec(i)=norm(qnorm);
    %----------------------------------------------------------------------
    
    %update the state for repropagation
    xnow(1:6)=Yf(1:6);
    xnow(4:6)=xnow(4:6)+dvMag*(xnow(4:6)/norm(xnow(4:6)));
    
    %store the values in order to plot the trajectory afterwards
    bigx(i,1:7)=Yf;
    
end

%unload fortran library
[Yf,dYf,d2Yf,capZ,errFlag]=propKeplerUV( xnow,dUV,orderCase,false,true);


figure(2)
plot(normVec)
title(['norm of MATLAB vs Fortran, state and partials, orderCase=',num2str(orderCase)])
xlabel('i^{th} segment')

figure(1)
%plot the trajectory
plot3(bigx(:,1),bigx(:,2),bigx(:,3),'.','linestyle','-')
hold all

%plot the controls
dvdir=bigx(:,4:6)
sz=size(Yf,1)
for j=1:i
    nrm=norm(dvdir(j,1:3));
    dvdir(j,1:3)=dvdir(j,1:3)*(dvMag/nrm);
end
quiver3(bigx(:,1),bigx(:,2),bigx(:,3),dvdir(:,1),dvdir(:,2),dvdir(:,3),'red')

%plot the center
plot3(0,0,0,'o')

xlabel('x')
ylabel('y')
zlabel('z')
title(['Trajectory, velocity pointing control at each node, |\Deltav|=',num2str(dvMag)])


