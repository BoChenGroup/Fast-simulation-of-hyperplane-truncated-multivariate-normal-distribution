% Simulation of Multinomial distributed dataset
% Cong Yulai
% 2016 11 10

%% 
clear,clc,close all;

%% 
LineSettings    =  {'--*','-.o',':s','-d','-->','-.h',':x','-p'}    ;
xlabelsize  =   28  ;   ylabelsize  =   28  ;   gcafontsize     =   18  ;
linewidthsize   =   2   ;   Markersize = 10 ;
FontName    =   'Times New Roman'   ;
ColorSet    =   [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1 ; 1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1]     ;
ColorSet    =   [ColorSet ; ColorSet(2:end,:) * 0.5]    ;
LineStyleCellSet    =   {'-' , '--' , ':' , '-.' }  ;
MarkerCellSet       =   {'o' , 'x' , 's' , 'd' , '^' , 'v' , '<' , '>' , 'p' , 'h' , '.' , '+' , '*'}   ;


%%
% V   =   2000    ;
% N   =   1e4     ;
% SparseRate  =   0.02    ;
% PoisRate    =   50      ;
% eta     =   0.01    ;
% 
% phi     =   rand(V,1)   ;
% phi(randperm(V,fix(V*SparseRate)))  =   100     ;
% phi     =   abs(phi) / sum(abs(phi))    ;
% figure(1),hold on,plot(phi,'r-'),title('phi')
% 
% tic
% edges   =   zeros(V+1,1)    ;   edges(1)    =    0  ;
% for i = 2:(V+1)
%     edges(i)    =   edges(i-1) + phi(i-1)     ;
% end
% 
% DataMat     =   sparse(V,N)     ;
% for i = 1:N
%     DataMat(:,i)    =   sparse( histcounts(rand(poissrnd(PoisRate),1),edges) )    ;
% end
% TimeGenData     =   toc     ;
% 
% save('MultDataTmp.mat')     ;

load('MultDataTmp.mat')     ;

% figure(1),hold on,plot(phi,'r-'),title('phi')

tmp     =   sum( DataMat , 2 ) + eta  ;
tmp     =   tmp / sum(tmp)  ;
% figure(1),plot(tmp,'g-'),legend('Truth','MaxLogPost')

ResErr0     =   norm( tmp - phi )   ;

M0   =   sum(DataMat(:))     ;

%% Settings
MBsize  =   10  ;
rho     =   N / MBsize     ;  

tao0        =   0   ;   kappa0  =   0.99     ;     epsi0    =   1     ;
epsit       =   (tao0 + (1:30000)) .^ (-kappa0)    ;     
epsit       =   epsi0 * epsit / epsit(1)    ;

IterMax     =   2000    ;

%% Our solution
if 1
    figure(2),subplot(221),plot(phi,'r-'),title('phi')
    figure(2),subplot(222),plot(epsit,'r-'),title('epsit')

    PhiOur  =   1/V * ones(V,1)     ;
    for epochi  =   1:200000
        idxall     =   randperm(N)  ;  
        for MBt  =   1:floor(rho)
            X   =   DataMat(:,idxall( (MBt-1)*MBsize + (1:MBsize) ))       ;

            MBObserved      =   (epochi-1)*floor(rho) + MBt     ;
            rhon    =   rho * sum(X,2)    ;

            if (MBObserved == 1)
                M     =   sum( rhon , 1 )  ;
            else
                M     =   (1 - epsit(MBObserved)) * M + epsit(MBObserved) * sum( rhon , 1 )   ;
            end

            stepsize    =   epsit(MBObserved) / M   ;

            tic
            tmp     =   PhiOur + stepsize * ((rhon + eta)- sum(rhon + eta) * PhiOur) + sqrt(2 * stepsize * PhiOur) .* randn(size(PhiOur))  ;
            PhiOur  =   ProjSimplexSpecial(tmp , PhiOur , 0)  ; 
            TimeOurCun(MBObserved)  =   toc     ;
            
            ResErr  =   norm( PhiOur - phi )    ;
            ResErrOurCun(MBObserved)  =   ResErr     ;

            figure(2),subplot(223),plot(PhiOur,'g-'),title('PhiOur')

            figure(2),subplot(224),hold on,plot(MBObserved,ResErr,'r*'),
            plot(MBObserved,ResErr0,'gs'),

            pause(0.05)
            
            if MBObserved > IterMax
                break   ;    
            end
        end
        if MBObserved > IterMax
            break   ;    
        end
    end
    
    save('ResultsOur.mat','ResErr0','M0','MBsize','tao0','kappa0',...
        'epsi0','epsit','rho','IterMax','TimeOurCun','ResErrOurCun')
    
end

%% PreGibbs
if 1
    figure(3),subplot(221),plot(phi,'r-'),title('phi')
    figure(3),subplot(222),plot(epsit,'r-'),title('epsit')

    PhiGibbs  =   1/V * ones(V,1)     ;
    PhiBarGibbs     =   PhiGibbs(1:end-1)   ;
    B   =   [-eye(V-1);ones(1,V-1)]  ;
    b   =   [zeros(V-1,1);1]     ;
    SampleEpoch     =   10   ;
    
    for epochi  =   1:200000
        idxall     =   randperm(N)  ;  
        for MBt  =   1:floor(rho)
            X   =   DataMat(:,idxall( (MBt-1)*MBsize + (1:MBsize) ))       ;

            MBObserved      =   (epochi-1)*floor(rho) + MBt     
            rhon    =   rho * sum(X,2)    ;

            if (MBObserved == 1)
                M     =   sum( rhon , 1 )  ;
            else
                M     =   (1 - epsit(MBObserved)) * M + epsit(MBObserved) * sum( rhon , 1 )   ;
            end

            stepsize    =   epsit(MBObserved) / M   ;

            tic
            tmpMu   =   PhiGibbs + stepsize * ((rhon + eta)- sum(rhon + eta) * PhiGibbs)    ;
            tmpSigma    =   2 * stepsize * ( diag(PhiBarGibbs) - PhiBarGibbs * PhiBarGibbs.' )  ;
%             if sum( eig(tmpSigma) < 0)
%                 aaa     =    1  ;
%             end
            PhiBarGibbs     =   TMVN_Gibbs_Cong(tmpMu(1:end-1), tmpSigma, B, b,1,SampleEpoch)  ;
            PhiGibbs    =   max(eps, abs( [PhiBarGibbs ; 1 - sum(PhiBarGibbs)] ) )    ;
            PhiGibbs    =   PhiGibbs / sum(PhiGibbs)   ;
            PhiBarGibbs     =   PhiGibbs(1:end-1)   ;
            TimePhiGibbsCun(MBObserved)  =   toc     ;
            
            ResErr  =   norm( PhiGibbs - phi )    ;
            ResErrPhiGibbsCun(MBObserved)  =   ResErr     ;
            
%             figure(3),subplot(223),plot(PhiGibbs,'g-'),title('PhiGibbs')
% 
%             figure(3),subplot(224),hold on,plot(MBObserved,ResErr,'r*'),
%             plot(MBObserved,ResErr0,'gs'),
%             
%             pause(0.05)

            if MBObserved > IterMax
                break   ;    
            end
        end
        if MBObserved > IterMax
            break   ;    
        end
    end
    
    filename    =   ['ResultsPreGibbs',num2str(SampleEpoch),'.mat']     ;
    save(filename,'SampleEpoch','TimePhiGibbsCun','ResErrPhiGibbsCun')
    
end

%% Figure
if 1
    load('ResultsOur.mat')  ;
    load('ResultsPreGibbs1.mat')  ;
    TimePhiGibbsCun1    =   TimePhiGibbsCun     ;
    ResErrPhiGibbsCun1  =   ResErrPhiGibbsCun   ;
    load('ResultsPreGibbs5.mat')  ;
    TimePhiGibbsCun5    =   TimePhiGibbsCun     ;
    ResErrPhiGibbsCun5  =   ResErrPhiGibbsCun   ;
    load('ResultsPreGibbs10.mat')  ;
    TimePhiGibbsCun10    =   TimePhiGibbsCun     ;
    ResErrPhiGibbsCun10  =   ResErrPhiGibbsCun   ;

    Iteration   =   1:length(TimeOurCun)    ;

    figure(5),hold on;box on;set(gcf, 'Color', 'w');
    plot(Iteration,ResErr0*ones(IterMax+1,1),'Color',ColorSet(1,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrOurCun,'Color',ColorSet(2,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun1,'Color',ColorSet(3,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun5,'Color',ColorSet(4,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun10,'Color',ColorSet(5,:),'LineWidth',linewidthsize)   ;
    legend('Batch posterior mean','SG-MCMC-fast','SG-MCMC-Gibbs 1','SG-MCMC-Gibbs 5','SG-MCMC-Gibbs 10')
    xlim([1 Iteration(end)]),xlabel('Iteration'),ylabel('Residual error'),
    set(gca,'fontsize',gcafontsize-2)
    
    axes('Position',[0.58,0.26,0.28,0.13]),hold on
    plot(Iteration,ResErr0*ones(IterMax+1,1),'Color',ColorSet(1,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrOurCun,'Color',ColorSet(2,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun1,'Color',ColorSet(3,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun5,'Color',ColorSet(4,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun10,'Color',ColorSet(5,:),'LineWidth',linewidthsize)   ;
    xmin    =   1200    ;   xmax    =   2000    ;   ymin    =   1e-3    ;   ymax    =   3e-3    ;
    xlim([xmin xmax]),ylim([ymin ymax]),box on;
    
    
    axes('Position',[0.23,0.55,0.1,0.3]),hold on
    plot(Iteration,ResErr0*ones(IterMax+1,1),'Color',ColorSet(1,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrOurCun,'Color',ColorSet(2,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun1,'Color',ColorSet(3,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun5,'Color',ColorSet(4,:),'LineWidth',linewidthsize)   ;
    plot(Iteration,ResErrPhiGibbsCun10,'Color',ColorSet(5,:),'LineWidth',linewidthsize)   ;
    xlim([1 10]),ylim([0.01 0.06]),box on;
    
    filename    =   ['ResErrIteration',num2str(MBsize)]     ;
    export_fig(filename,'-pdf','-r600','-q101')

    TimeOur     =   TimeOurCun  ;
    TimePhiGibbs1   =   TimePhiGibbsCun1     ;
    TimePhiGibbs5   =   TimePhiGibbsCun5     ;
    TimePhiGibbs10  =   TimePhiGibbsCun10    ;
    for i = Iteration(2:end)
        TimeOur(i)  =   TimeOur(i-1) + TimeOur(i)   ;
        TimePhiGibbs1(i)  =   TimePhiGibbs1(i-1) + TimePhiGibbs1(i)   ;
        TimePhiGibbs5(i)  =   TimePhiGibbs5(i-1) + TimePhiGibbs5(i)   ;
        TimePhiGibbs10(i) =   TimePhiGibbs10(i-1) + TimePhiGibbs10(i)   ;
    end

    figure(6),grid on;box on;set(gcf, 'Color', 'w');
    tmp     =   0.01:max(TimePhiGibbs10)/100:max(TimePhiGibbs10)    ;
    semilogx( tmp , ResErr0*ones(size(tmp)),'Color',ColorSet(1,:),'LineWidth',linewidthsize)   ;hold on;
    semilogx(TimeOur,ResErrOurCun,'Color',ColorSet(2,:),'LineWidth',linewidthsize)   ;
    semilogx(TimePhiGibbs1,ResErrPhiGibbsCun1,'Color',ColorSet(3,:),'LineWidth',linewidthsize)   ;
    semilogx(TimePhiGibbs5,ResErrPhiGibbsCun5,'Color',ColorSet(4,:),'LineWidth',linewidthsize)   ;
    semilogx(TimePhiGibbs10,ResErrPhiGibbsCun10,'Color',ColorSet(5,:),'LineWidth',linewidthsize)   ;
%     legend('Batch posterior mean','SG-MCMC-fast','SG-MCMC-Gibbs 1','SG-MCMC-Gibbs 5','SG-MCMC-Gibbs 10')
    xlim([0.01 max(TimePhiGibbs10)]),xlabel('Time (s)'),ylabel('Residual error'),
    set(gca,'fontsize',gcafontsize-2,'XTick',[0.01 1 100 10000]);
    
    filename    =   ['ResErrTime',num2str(MBsize)]     ;
    export_fig(filename,'-pdf','-r600','-q101')

end
