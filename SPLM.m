function [z,dt,t, upliftSet,U_store,a]=SPLM(varargin)

% Stream Power Law Model: SPLM
%
% Syntax
%
%       [z,dt,t, upliftSet,U_store,a]=SPLM(varargin)
%
% Description
%
%       SPLM solves the 1D stream power law using different numerical
%       methods. Numerical methods are based on TORO et al. Springer (3th
%       edition) Slopes below 1e-4 are set to 0 in order to avoid very
%       small timesteps for n<1
%
% Input arguments
%
%      Input arg 1:
%            Numerical method
%            Case 1 : explicit upwind scheme (see Roberts 2010)
%            Case 2 : implicit upwind scheme (see Braun and Willet 2013)
%            Case 3 : TVD scheme Campforts and Govers 2014
%            Case 4 : comparison of schemes (default)
%
%      Input arg 2:   
%           Parameters=[K m n kappa]
%               default values: [2e-06 .4 1 0].
%      Input arg 3:  
%           Spatial {dx x DA} in m and m²
%               default values: 
%                   dx = 500m x = dx:dx:50E3 DA=x.^2
%      Input arg 4:   
%           timing [t_end dt]: time in years
%               default values:
%                   t_end=5e06y
%                   dt=calcualted automatically using the CFL criterion
%      Input arg 5:
%           UpliftData={UScen maxElevation upliftSet interpolation}
%               UScen:
%                   Case 1 : two instanteneous uplift pulses
%                   Case 2 : 3 Uplift pulses {Default}
%                   Case 3 : No uplift: analytical solution is added.
%                   Case 4 : Use inserted Uplift rate
%               maxElevation: total amount of uplift (m)
%               upliftSet: variation of uplift in time
%               interpolation: if 0, inserted upliftSet is used as absolute uplift
%                   rates. if 1, inserted upliftSet is rescaled according to the given
%                   maximum elevation. {default value: 1}
%      Input arg 6: 
%           oriBed={iniSurf baseLevelDescent}
%           iniSurf: surface before model run
%           baseLevelDescent: lowering of base level during model run
%           default values: iniSurf=x.*0 baseLevelDescent=0
%      Input arg 7:   
%           visibleFlag=display intermediate results or not
%           default value: 0
%      Input arg 8:  
%           plotOut=display output
%           default value: 1 
% 
% Output arguments
%
%      z
%           structure with calculated elevations, depending on selected
%           numerical method for analytical solutions an adiditional field
%           z.an is added
%
%
% Example 1: Linear incision, resolution 100m
%
%   K=5e-6; m=.42; n=1;kappa=0;
%   parameters=[K m n kappa]; 
%   dx=100; x=1:dx:15E3;hackFactor=2;DA=x.^hackFactor;
%   spatial={dx x DA} ; 
%   t_end=1E6;
%   timing=[t_end nan]; 
%   UScen=1; maxElevation=1000;
%   upliftData={UScen maxElevation}; 
%   iniSurf = x.*0;baseLevelDescent=0;
%   oriBed={iniSurf baseLevelDescent}; 
%   visibleFlag=0;
%   plotOut=1; 
%   [z, dt]=SPLM(numM,parameters,spatial,timing,upliftData,oriBed,visibleFlag,plotOut);
%
%
% Example 2 Linear incision, resolution 100m, wit analytical characteristics
%
%   Numerical method
%   numM=4;
%   %Parameters
%   K=5e-6; m=.42; n=1;kappa=0;
%   parameters=[K m n kappa];
%   %Spatial
%   dx=100; x=1:dx:15E3;
%   hackFactor=2;
%   DA=x.^hackFactor;
%   spatial={dx x DA} ;
%   %Timing
%   t_end=1E6;
%   timing=[t_end nan];
%   %Uplift Scen
%   uScen=3; maxElevation=1000;
%   upliftData={uScen maxElevation};
%   %Original bed and baselevel evolution
%   nb =length(x);
%   iniSurf = x.*0;
%   baseLevelDescent=0;
%   oriBed={iniSurf baseLevelDescent};
%   visibleFlag=0;
%   plotOut=1;
%   %Run the model
%   [z, dt]=SPLM(numM,parameters,spatial,timing,upliftData,oriBed,visibleFlag,plotOut);
%
% Example3
%   clearvars; clc;close all force;
%   numM=4;
%   K=5e-6; m=.42; n=.9;kappa=0;
%   parameters=[K m n kappa];
%   dx=100; x=1:dx:15E3;hackFactor=2;DA=x.^hackFactor;
%   spatial={dx x DA} ;
%   t_end=1E6;
%   timing=[t_end nan];
%   UScen=3; maxElevation=1000;
%   upliftData={UScen maxElevation};
%   iniSurf = x.*0;baseLevelDescent=0;
%   oriBed={iniSurf baseLevelDescent};
%   visibleFlag=0;
%   plotOut=1;
%   [z, dt]=SPLM(numM,parameters,spatial,timing,upliftData,oriBed,visibleFlag,plotOut);
%
% Example4
%   clearvars; clc;close all force;
%   numM=4;
%   K=5e-6; m=.42; n=1.1;kappa=0;
%   parameters=[K m n kappa];
%   dx=100; x=1:dx:15E3;hackFactor=2;DA=x.^hackFactor;
%   spatial={dx x DA} ;
%   t_end=1E6;
%   timing=[t_end nan];
%   UScen=3; maxElevation=1000;
%   upliftData={UScen maxElevation};
%   iniSurf = x.*0;baseLevelDescent=0;
%   oriBed={iniSurf baseLevelDescent};
%   visibleFlag=0;
%   plotOut=1;
%   [z, dt]=SPLM(numM,parameters,spatial,timing,upliftData,oriBed,visibleFlag,plotOut);
%
% References
%
% * TVD-FVM: Campforts, B. and Govers, G.:
%           Keeping the edge: A numerical method that avoids knickpoint
%           smearing when solving the stream power law, J. Geophys. Res.
%           Earth Surf., 120(7), 1189–1205, doi:10.1002/2014JF003376, 2015.
%
% * TTLEM:  Campforts B., Schwanghart W., Govers G.: TTLEM 1.0 : a
%           numerical package for accurate simulation of transient landscape
%           evolution in MATLAB. Discussion paper in GMD.
%
% Date: 15. June, 2016


%% Read Input


optargin = size(varargin,2);

%% Paramter initialisation
% Numerical method
if optargin<1
    numM=4;
else
    numM=varargin{1};
    if numM<1 || numM>4
        error('NumM must have a value between 1 and 4')
    end
end

% Model parameters
if optargin<2
    K=2e-06;
    m=.4;
    n=1;
    kappa=0;
else
    K=varargin{2}(1);m=varargin{2}(2);n=varargin{2}(3);kappa=varargin{2}(4);
end
if n<.7
    error('Currently, n cannot be lower than 0.7')
end
    
% Model spatial properties
if optargin<3
    dx=500;
    x=dx:dx:100E3;
    DA=x.^2;
else
    dx=varargin{3}{1}; x=varargin{3}{2}; DA=varargin{3}{3};
end
z.DA=DA;
z.x=x;
x_ori=x;
a = -K.*DA.^m;
cfl = .7;

% Model temporal properties
if optargin<4
    t_end=5e6;
    dt = cfl*dx/max(abs(a));
else
    t_end=varargin{4}(1);
    if isnan(varargin{4}(2))
        if numM~=2
            dt = cfl*dx/max(abs(a));
        else
            dt = cfl*dx/max(abs(a))*1;
        end
    else
        dt = varargin{4}(2);
    end
end

% Model uplift scenario
if optargin<5
    uScen=3;
    maxElevation=1000;
else
    uScen=varargin{5}{1};
    if uScen<0||uScen>4
        error('uScen should be between 1 and 4')
    end
    if numel(varargin{5})==1
        maxElevation=1000;
    elseif  numel(varargin{5})==2
        maxElevation=varargin{5}{2};
    elseif  numel(varargin{5})==3
        maxElevation=varargin{5}{2};
        u_rate_my=varargin{5}{3};
    elseif  numel(varargin{5})==4
        maxElevation=varargin{5}{2};
        u_rate_my=varargin{5}{3};
    end
end

% Model initial surface and baselevel properties
if optargin<6
    iniSurf=x.*0;
    baseLevelDescent=0;
    constantDA=0;
else
    iniSurf=varargin{6}{1};
    baseLevelDescent=varargin{6}{2};
    if size(varargin{6},2)==3
        constantDA=1;
    else
        constantDA=0;
    end
end
if uScen==3 % in case of an analytical solution the initial condition will be set to @shape3
    iniSurf=feval(@shape3,x,x_ori);
    z.iniSurf=iniSurf;
end
if optargin<7
    visible=0;
else
    visible=varargin{7};
end
if visible
    %Preparation movie file
    riverPlot=figure('units','normalized','position',[0.1 0.1 .8 .8],'color',[1 1 1]) ;
    movegui(riverPlot, 'onscreen');
    
    movieFlag=1;
    if movieFlag
        vidObj = VideoWriter('riverModels.avi');
        open(vidObj);
        movegui(riverPlot, 'onscreen');
        rect = get(riverPlot,'Position');
        rect(1:2) = [0 0];
    else
        F(round(t_end/dt)+1) = struct('cdata',[],'colormap',[]);
    end
end
if optargin<8
    plotOut=2;
else
    plotOut=varargin{8};
end

% Uplift
nbScenVar=ceil(t_end/dt);%2E4;
upliftSet=zeros(1,nbScenVar);
upliftSetTime_V=linspace(1,t_end,nbScenVar);
interpolationFlag=1;
if uScen==1
    place1=round(nbScenVar/8);
    upliftSet(1:place1)=1;
    upliftSet(place1*3+1:place1*4)=1;
elseif uScen==2
    place1=round(nbScenVar/10);
    upliftSet(1+place1*.4:place1+place1*.4)=1;
    upliftSet(place1*4+place1*.4+1:place1*5+place1*.4)=1;
    upliftSet(place1*9+1:end)=1;
elseif uScen==3
    upliftSet=0;
elseif uScen==4
    upliftSet=varargin{5}{3};
    interpolationFlag=varargin{5}{4};
end


%Convert UpliftSet to m/y
if uScen<=2 
        upliftSetTime=t_end/numel(upliftSet);
        upliftSet=upliftSet*maxElevation/sum(upliftSet);
        upliftSet=upliftSet/upliftSetTime;% Uin m/y
end

if uScen==4
     if interpolationFlag==1
        disp('Carefull, version dependent!')
        upliftSetTime=t_end/numel(upliftSet);
        upliftSet=upliftSet*maxElevation/sum(upliftSet);
        upliftSet=upliftSet/upliftSetTime;% Uin m/y
    else
         upliftSet=upliftSet;
    end
end
z.upliftSet=upliftSet;

% Contruct Uplift vector
if uScen<=2 || uScen==4
    if interpolationFlag==1
        U=upliftSet*dt;
        U_store=U;
    elseif interpolationFlag==2
        placeInSet=ceil(t./upliftSetTime);
        placeInSet(1)=1;
        U=upliftSet(placeInSet)*dt;
        U_store=U*maxElevation/sum(U);
        U_store(isnan(U_store))=0;
    end
end

% Initialisation of vectors for main loop
t = dt:dt:t_end;
z.ex= iniSurf;
z.im= iniSurf;
z.TVD = iniSurf;
z.iniSurf=iniSurf;

% Define waitbar
if plotOut
    h=waitbar(0,'Time for a coffee');
end

% Extend rivers
x=[x(1) x];
if constantDA
    DA=[DA(1) DA];
else
%   DA=[DA(1) DA(2)+DA];
DA=[DA(1) DA];
end
% Extend inisurf
z.ex= [z.ex(1) z.ex];
z.im=[z.im(1) z.im];
z.TVD = [z.TVD(1) z.TVD];
% 
z.DA=DA;
z.x=x;
nb = length(x);
a = -K.*DA.^m;
a_ori=a;
a_m = min(0,a);
a_p = max(0,a);

% Matrix for Crank-Nicolson if diffusion to river profile is applied
if kappa~=0
    alpha= (kappa)*dt/(dx^2);
    A=zeros(nb); %Creation of matrix
%   r= zeros(1, nbJ);
    %Boundary conditions
    b0=1; c0=0;
    A(1,1)= b0; A(1,2)= c0;
    aJ=0; bJ=1;
    A(nb,nb-1)=aJ; A(nb,nb)=bJ;
    %Create matrix A
    for j=2:nb-1
        k=j-1;
        A(j,k)=-0.5*alpha;
        A(j,k+1)=1+alpha;
        A(j,k+2)=-0.5*alpha;
    end
    At=A^-1;
end

%% Main loop
ui=0;
for k = t 
    
    ui=ui+1;
    if plotOut
        waitbar(ui/numel(t),h);
    end
    
    %% Uplift
    if uScen~=3
        z.im(1:end-1) = z.im(1:end-1)+U_store(ui);
        z.ex(1:end-1)= z.ex(1:end-1)+U_store(ui);
        z.TVD(1:end-1)= z.TVD(1:end-1)+U_store(ui);
    end
    
    zprior.ex=z.ex;
    zprior.TVD=z.TVD;    
    zprior.im=z.im;  
    z.ex(end)=z.ex(end)-baseLevelDescent/numel(t);
    z.im(end)=z.im(end)-baseLevelDescent/numel(t);
    z.TVD(end)=z.TVD(end)-baseLevelDescent/numel(t);
    
    %% Explicit upwind
    if numM==1 || numM==4
        
        dte=dt;
        int_time=dt;
        while int_time>0
            
            z.ex(z.ex<1e-6)=0;
            ex_r=[z.ex(2:end) z.ex(end)];
        %     ex_r(ex_r>z.ex)=z.ex(ex_r>z.ex);
            s1=(z.ex-ex_r)/dx;
            s1(s1<1e-4)=0;
            if n~=1
                exp_f=s1.^(n-1);
                exp_f(isinf(exp_f))=1;
                a=a_ori.*exp_f;
            end
            
        %   Check timestep
            dt_calc = cfl*dx/max(abs(a));
            if dt_calc<dte
                disp(['For stability, dte is set to: ' num2str(dt_calc)]);
                dte=dt_calc;
            end
            int_time=int_time-dte;
            if int_time<0
                dte=dte+int_time;
                int_time=0;
            end
            
            z.ex= a.*(s1.^1)*dte+z.ex;
            if any(~isreal(z.ex))
                z.ex=real(z.ex);
                error('Not real z.ex values')
            end
            
            %Crank-Nicolson
            if kappa~=0
                %create r
                r=z.ex';
                CN=(At*r)';
                z.ex(1:end-1)=CN(1:end-1);
            end
            erosion.ex(ui,:)=zprior.ex-z.ex;
        end
    end
    
    
    %% Implicit upwind
    if numM==2 || numM==4 || numM==5
        zPrev=z.im;
        for j=nb-1:-1:1
            if n==1
                  z.im(j)= (z.im(j)+z.im(j+1)*K*(DA(j)^m)*dt/dx)/(1+(K*DA(j)^m)*dt/dx);
%                 z.im(j)= (z.im(j)-(K*(DA(j)^m)*dt/(2*dx))*(z.im(j)- zPrev(j+1)-z.im(j+1)))...
%                     /(1+(K*DA(j)^m)*dt/(2*dx));
            else
                tempz=z.im(j);
                finalz=tempz-(tempz-z.im(j)+K*DA(j)^m*dt*((tempz-z.im(j+1))/dx)^n)/(1+(n*K*DA(j)^m*dt/dx)*((tempz-z.im(j+1))/dx)^(n-1));
                iter=0;
                while round(tempz*1E2)~=round(finalz*1E2)&&iter<100
                    iter=iter+1;
                    tempz=finalz;
                    finalz=tempz-(tempz-z.im(j)+K*DA(j)^m*dt*((tempz-z.im(j+1))/dx)^n)/(1+(n*K*DA(j)^m*dt/dx)*((tempz-z.im(j+1))/dx)^(n-1));
                end
                z.im(j)=finalz;
            end
            if z.im(j+1)>0 && z.im(j)<z.im(j+1)
            %             disp('Carefull')
                z.im(j)=z.im(j+1);
            end
        end
        %Crank-Nicolson
        if kappa~=0
            r=z.im';
            CN=(At*r)';
            z.im(1:end-2)=CN(1:end-2);
        end
        if any(~isreal(z.im))
            z.im=real(z.im);
%       disp('Not real z.ex values')
        end
    end
   %% TVD, assume a always +, 
   if numM==3 || numM==4 || numM==5
       
       dte=dt;
       int_time=dt;
       while int_time>0
           
           %% TVD
           z.TVD(z.TVD<1e-6)=0;
           if n~=1
               TVD_r_t=[z.TVD(2:end) z.TVD(end)];
               s1=(z.TVD-TVD_r_t)/dx;
               s1(s1<1e-4)=0;
               exp_f=s1.^(n-1);
               exp_f(isinf(exp_f))=1;
               a=a_ori.*exp_f;
               a_m = min(0,a);
               a_p = max(0,a);
           end
               
       %   Check timestep
           dt_calc = cfl*dx/max(abs(a));
           if dt_calc<dte
               disp(['For stability, dte is set to: ' num2str(dt_calc)]);
               dte=dt_calc;
           end
           int_time=int_time-dte;
           if int_time<0
               dte=dte+int_time;
               int_time=0;
           end
           
           TVD_r=[z.TVD(2:end) z.TVD(end)];
           TVD_r=[z.TVD(2:end) z.TVD(end)];
           TVD_r2=[z.TVD(3:end) z.TVD(end) z.TVD(end)];
           TVD_l=[z.TVD(1) z.TVD(1:end-1)];
           r_TVD=(TVD_r2-TVD_r)./(TVD_r-z.TVD);
           r_TVD(diff(z.TVD)==0)=1;
           r_TVD(1) = 1;
%      r_TVD(nb-1) = r_TVD(nb-2);
           r_TVD(nb) = 1;
           
       %   Define Flux Limiter function
           %VANLEER
           phi = (r_TVD + abs(r_TVD))./(1 + abs(r_TVD));
           
       %   Compute fluxes for TVD
           TVD_r=[z.TVD(2:end) z.TVD(end)];
           TVD_l=[z.TVD(1) z.TVD(1:end-1)];
           
           F_rl = a_p.*z.TVD + a_m.*TVD_r;
           F_rh = (1/2)*a.*(z.TVD+TVD_r) - (1/2)*(a.^2).*(dte/dx).*(TVD_r-z.TVD);
           F_ll = a_p.*TVD_l + a_m.*z.TVD;
           F_lh= (1/2)*a.*(TVD_l+z.TVD) - (1/2)*(a.^2).*(dte/dx).*(z.TVD-TVD_l);
       %   Compute nz.ext time step
           phi_l=[phi(1) phi(1:end-1)];
           F_right = F_rl + phi.*(F_rh-F_rl);
           F_left = F_ll+ phi_l.*( F_lh- F_ll);      
           TVD_next= z.TVD-dte*(F_right-F_left)/dx;           
           
       %   UPDATE info
           TVD_next(1) = TVD_next(2);
           TVD_next(nb) = z.TVD(end);
           if any(~isreal(TVD_next))
               error('imaginary number');
           %         real(TVD_next);
           end
           z.TVD = TVD_next;
           erosion.TVD(ui,:)=zprior.TVD-z.TVD;
           
           if kappa~=0
               %diffusive term
               %create r
               r=z.TVD';
               %Crank-Nicolson
               CN=(At*r)';
               z.TVD(1:end-2)=CN(1:end-2);
           end
       end
   end
    %% Plot output
    if visible
        
        figure(riverPlot)
        if uScen==3
            time=t(ui);
            an = feval(@exact3);
            plot(x,an,'color',[0 100 0],'linewidth',2);
            hold on
        end
        
        plot(x,z.ex,'-*r','linewidth',1);
        hold on
        plot(x,z.im,'-*b');
        plot(x,z.TVD,'-*k');
%   plot(x,iniSurf,'--k')
        
        hold off
    %     ylim([0 ,2700]);
        xlabel('Distance, m','fontweight','bold')
        ylabel('Elevation, m','fontweight','bold')
        if uScen==3
            legend('Analytical (n=1)','Explicit FDM','Implicit FDM','TVD\_FVM','Initial','location','northeast')
        else
            legend('Explicit FDM','Implicit FDM','TVD\_FVM','Initial','location','southwest')
        end
        
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        if movieFlag
            movegui(riverPlot, 'onscreen');
            hold all;
        %     datetick;
        %     drawnow;
            writeVideo(vidObj,getframe(gcf));
            hold off
        else
            F(ui)=getframe;
        end
    end
    
end

% Normalise length
x=x(2:end);
DA=DA(2:end);
z.x=x;
z.DA=DA;
z.ex=z.ex(2:end);
z.im=z.im(2:end);
z.TVD=z.TVD(2:end);

% Analytical solution with methods of characteristics
if uScen==3
    time=t(ui);
    z.an = feval(@exact3,x,time, K,m,x_ori);
end

% Close movie if applies
if plotOut
    close(h)
end

if visible
    if movieFlag
        close(vidObj);
    else
        if uScen==3
            legend('Analytical (n=1)','Explicit FDM','Implicit FDM','TVD\_FVM','Initial','location','northeast')
        else
            legend('Explicit FDM','Implicit FDM','TVD\_FVM','Initial','location','northeast')
        end
    end
end


%% Plot output
if plotOut==2
    
    g=figure('units','normalized','outerposition',[0.05 0.05 .55 0.45], 'color', [1 1 1]);
    hold on 
    x_km=x*1E-3;
    if uScen==3
        plot(x_km,iniSurf,'linestyle','-', 'color',[0.4 0.4 0.4],'linewidth',2.5);
        hold on
        plot(x_km,z.an,'color',[0 0.6 0],'linewidth',2);
    end
    if numM==1
        plot(x_km,z.ex,'--r','linewidth',2);
    else
        plot(x_km,z.ex,'--r','linewidth',2);
        plot(x_km,z.im,'--b','linewidth',2);
        plot(x_km,z.TVD,'--k','linewidth',2);
        box on
        hold off
    end
    if ~interpolationFlag
        ylim([-5 ,maxElevation+0.05*maxElevation]);
    end
    xlabel('Distance, km','fontweight','bold')
    ylabel('Elevation, km','fontweight','bold')
    if uScen==3
        legend('Initial',sprintf('Analytical\ncharacteristics'),'Explicit FDM','Implicit FDM','TVD\_FVM','location','northeast')
        legend('boxoff')
        ylim([-5,max(z.an)+max(z.an)*0.01]);
    else
        if numM==1
            legend('Explicit ','location','northeast')
        else
            legend('Explicit FDM','Implicit FDM','TVD\_FVM','location','northeast')
        end
        legend('boxoff')
    end    
    set(findall(gcf,'-property','FontSize'),'FontSize',12)


    
    %% Plot uplift
%     g=figure('units','normalized','outerposition',[0.05 0.05 .55 0.45], 'color', [1 1 1]);
%     subplot(1,2,1)
%     plot(x*1e-3,DA*1e-6,'k','linewidth',2);
%     xlabel('Distance, km','fontweight','bold')
%     ylabel('DA, km^2','fontweight','bold')
%     %
%     if uScen<3 || uScen==4
%         subplot(1,2,2)
%         uplift_mma=upliftSet*1e3;
%         upliftSetTime_V_Ma=upliftSetTime_V*1e-6;
%         ha=area(max(upliftSetTime_V_Ma)-(upliftSetTime_V_Ma),(uplift_mma));
%         set(gca, 'xdir','reverse')
%         set(ha(1),'FaceColor',[.5 0.5 0.5])
%         set(ha(1),'LineStyle','none')
%         xlabel('Age, Ma','fontweight','bold')
%         ylabel('Uplift rate, mm/a','fontweight','bold')
%     %     ylim([0 ,max(uplift_mma)]);
%     %     xlim([0 ,max(t_end)]);
%         box on
%     end    
%     set(findall(gcf,'-property','FontSize'),'FontSize',12)
%     set(findall(gcf,'-property','FontSize'),'Fontname','Arial')

end
