% Cong Yulai
% 2017 08 06
%% SimuGauLinearEqConst 
clear,clc,close all
LineSettings    =  {'--*','-.o',':s','-d','-->','-.h',':x','-p'}    ;
xlabelsize  =   28  ;   ylabelsize  =   28  ;   gcafontsize     =   18  ;
linewidthsize   =   2   ;   Markersize = 10 ;
FontName    =   'Times New Roman'   ;
ColorSet    =   [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1 ; 1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1]     ;
ColorSet    =   [ColorSet ; ColorSet(2:end,:) * 0.5]    ;
LineStyleCellSet    =   {'-' , '--' , ':' , '-.' }  ;
MarkerCellSet       =   {'o' , 'x' , 's' , 'd' , '^' , 'v' , '<' , '>' , 'p' , 'h' , '.' , '+' , '*'}   ;


%% Figure 1: 2-dimensional x => z
if 0
    mu  =   [1,1.2]   ;
    SIGMA   =   [1,0.3;0.3,1]     ;
    G   =   [1,1]   ;   r   =   1   ;
    thx     =   3   ;

    fourpoint   =   [thx,0 ; -thx,0 ; 0,thx ; 0,-thx]   ;
    tmp     =   ( ones(4,1) * mu + fourpoint * chol(SIGMA,'upper') ) .'     ;
    x1min   =   min(tmp(1,:))  ;   x1max   =   max(tmp(1,:))  ;
    x2min   =   min(tmp(2,:))  ;   x2max   =   max(tmp(2,:))  ;

    x1  =   linspace(x1min,x1max,100)   ;
    x2  =   (linspace(x2min,x2max,100)).'   ;
    x1mat   =    ones(length(x2),1) * x1    ;
    x2mat   =    x2 * ones(1,length(x1))    ;
    X   =   [x1mat(:) , x2mat(:)]       ;
    pX  =   mvnpdf(X,mu,SIGMA)    ;
    pXmat   =   reshape(pX,length(x2),length(x1))   ;

    y   =   (r - G(1) * x1) / G(2)  ;

    H1  =   null(G)     ;
    M   =   eye(1)  ;
    H2  =   SIGMA * G.' * M     ;
    H   =   [H1, H2]    ;
    Hinv    =   inv(H)  ;

    % figure(5),subplot(121),hold on,box on
    figure(5),hold on,box on 
    %text(0,0,'\it O','Color','red','FontSize',28)
    plot([0,100], [0 0],'k:','LineWidth',linewidthsize*1.2),plot([0 0] , [0,100],'k--','LineWidth',linewidthsize*1.2)
    contour(x1,x2,pXmat,'LineWidth',linewidthsize*1.2),plot(x1,y,'k','LineWidth',linewidthsize*1.2)
    axis([x1min,x1max,x2min,x2max]),set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('\it{x}_{\rm{1}}','FontSize',xlabelsize),
    ylabel('\it{x}_{\rm{2}}','FontSize',ylabelsize),

    plot([0 H1(1)*10], [0 H1(2)*10],'r:','LineWidth',linewidthsize*1.2),plot([0 H2(1)*10], [0 H2(2)*10],'r--','LineWidth',linewidthsize*1.2)
    filename    =   'pxOri'    ;
    % saveas(gcf,filename,'pdf')     ;



    muz =   (Hinv * mu.') .'     ;
    SIGMAz  =   Hinv * SIGMA * Hinv.'   ;
    z2constr    =   inv(G * H2) * r     ;

    tmp     =   ( ones(4,1) * muz + fourpoint * chol(SIGMAz,'upper') ) .'     ;
    z1min   =   min(tmp(1,:))  ;   z1max   =   max(tmp(1,:))  ;
    z2min   =   min(tmp(2,:))  ;   z2max   =   max(tmp(2,:))  ;

    z1  =   linspace(z1min,z1max,100)   ;
    z2  =   (linspace(z2min,z2max,100)).'   ;
    z1mat   =    ones(length(z2),1) * z1    ;
    z2mat   =    z2 * ones(1,length(z1))    ;
    Z   =   [z1mat(:) , z2mat(:)]       ;

    pZ  =   mvnpdf(Z,muz,SIGMAz)    ;
    pZmat   =   reshape(pZ,length(z2),length(z1))   ;

    % figure(5),subplot(122),hold on,box on
    figure(6),hold on,box on 
    plot([0,100], [0 0],'r:','LineWidth',linewidthsize*1.2),plot([0 0] , [0,100],'r--','LineWidth',linewidthsize*1.2)
    contour(z1,z2,pZmat,'LineWidth',linewidthsize*1.2),
    plot([z1min z1max], [z2constr z2constr],'k','LineWidth',linewidthsize*1.2)
    axis([z1min,z1max,z2min,z2max]),set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('\it{z}_{\rm{1}}','FontSize',xlabelsize),
    ylabel('\it{z}_{\rm{2}}','FontSize',ylabelsize),

    plot([0 Hinv(1,1)*10], [0 Hinv(2,1)*10],'k:','LineWidth',linewidthsize*1.2),plot([0 Hinv(1,2)*10], [0 Hinv(2,2)*10],'k--','LineWidth',linewidthsize*1.2)
    filename    =   'pzIndep'    ;
    % saveas(gcf,filename,'pdf')     ;
end


%% Figure 2: 2-dimensional x   algorothm 2
if 0
    mu  =   [1,1.2]   ;
    SIGMA   =   [1,0.3;0.3,1]     ;
    G   =   [1,1]   ;   r   =   1   ;
    thx     =   3   ;

    fourpoint   =   [thx,0 ; -thx,0 ; 0,thx ; 0,-thx]   ;
    tmp     =   ( ones(4,1) * mu + fourpoint * chol(SIGMA,'upper') ) .'     ;
    x1min   =   min(tmp(1,:))  ;   x1max   =   max(tmp(1,:))  ;
    x2min   =   min(tmp(2,:))  ;   x2max   =   max(tmp(2,:))  ;

    x1  =   linspace(x1min,x1max,100)   ;
    x2  =   (linspace(x2min,x2max,100)).'   ;
    x1mat   =    ones(length(x2),1) * x1    ;
    x2mat   =    x2 * ones(1,length(x1))    ;
    X   =   [x1mat(:) , x2mat(:)]       ;
    pX  =   mvnpdf(X,mu,SIGMA)    ;
    pXmat   =   reshape(pX,length(x2),length(x1))   ;

    y   =   (r - G(1) * x1) / G(2)  ;


    yp  =   [1;2]   ;
    ypvec   =   SIGMA * G.' * (G * SIGMA * G.')^(-1) * (r - G * yp)     ;
    xp  =   yp + ypvec  ;

    SIGMA * G.' * (G * SIGMA * G.')^(-1) * r
    eye(length(mu)) - SIGMA * G.' * (G * SIGMA * G.')^(-1) * G


    figure(7),hold on,box on 
    plot([0,100], [0 0],'k:','LineWidth',linewidthsize*1.2),plot([0 0] , [0,100],'k--','LineWidth',linewidthsize*1.2)
    contour(x1,x2,pXmat,'LineWidth',linewidthsize*1.2),plot(x1,y,'k','LineWidth',linewidthsize*1.2)
    axis([x1min,x1max,x2min,x2max]),set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('\it{x}_{\rm{1}}','FontSize',xlabelsize),
    ylabel('\it{x}_{\rm{2}}','FontSize',ylabelsize),

    plot(yp(1),yp(2),['r',MarkerCellSet{1}],'MarkerSize',Markersize,'LineWidth',linewidthsize*1.2),
    text(yp(1)+0.2,yp(2)+0.1,'y','FontSize',30,'FontName',FontName,'FontWeight','bold','FontAngle','italic')
    plot(xp(1),xp(2),['r',MarkerCellSet{5}],'MarkerSize',Markersize,'LineWidth',linewidthsize*1.2),
    text(xp(1)-0.42,xp(2),'x','FontSize',30,'FontName',FontName,'FontWeight','bold','FontAngle','italic')
    plot([yp(1) xp(1)] , [yp(2) xp(2)] , 'r:','LineWidth',linewidthsize*1.2),

    filename    =   'pxOriProjYtoX'    ;
    % saveas(gcf,filename,'pdf')     ;
end


%% Figure 3 and 4: [Alg 1 vs Alo 2]  Time vs Dimensions
if 0
    folder  =   'Figures/'  ;
    SampleNum   =   10000    ;  TrialNum    =   5  ;
    kset    =   [50 200 500 1000 3000 5000]    ;

    DataOut     =   []  ;     DataOut.kset  =   kset   ;  DataOut.TrialNum  =   TrialNum    ;

    TimeAlg1Set     =   zeros(1,TrialNum)   ;   
    TimeAlg2Set     =   zeros(1,TrialNum)   ;

%     k2settingsSet  =   [20,0.1:0.2:0.9]    ;
%     for k2i     =   1:length(k2settingsSet)    ;
%         k2settings  =   k2settingsSet(k2i)  ;
%         for DiagSigma   = [0 1]
%             if k2settings >= 1
%                 k2set   =   k2settings * ones(1,length(kset))  ;
%             else
%                 k2set   =   fix(k2settings * kset)  ;
%             end
%     
%             AverTimeAlg1    =   []  ;   AverTimeAlg2    =   []  ;
%     
%             for i   =   1:length(kset)
%                     k   =   kset(i)     ;
%                     k2      =   k2set(i)    ;  
%                     k1  =   k - k2      ;
%     
%                     for trial = 1:TrialNum
%                         %%%%%%%%%%%%%%%%%%%%     Generate Mu,Sigma,G,r   %%%%%%%%%%%%%%%%%%
%                         Mu  =   randn(k,1)     ;
%     
%                         if DiagSigma
%                             Sigma   =   0.05 + rand(k,1)  ;
%                         else
%                             U   =   orth(rand(k,k))   ;
%                             Sigma   =   U.' * diag(0.05 + rand(k,1)) * U    ;
%                         end
%     
%                         G   =   randn(k2,k)     ;
%                         [U,S,V]     =   svd(G)  ;
%                         S((S<=0.05) & (S>0))    =   0.05    ;
%                         G   =  U * S * V.'  ; 
%     
%                         r   =   randn(k2,1)     ;
%     
%                         %%%%%%%%%%%%%%%%%%%%     Algorithm 1   %%%%%%%%%%%%%%%%%%
%                         tic
%     
%                         [U,S,V]     =   svd(G)  ;   
%                         H1  =   V(:,(k2+1):k)   ;   H2  =   V(:,1:k2)   ;   H   =   [H1,H2]     ;
%     
%                         z2  =   inv(G * H2) * r     ;
%                         if DiagSigma
%                             tmp    =   bsxfun(@times, H1.', 1./Sigma.' )  ;
%                         else
%                             tmp    =   H1.' * inv(Sigma)  ;
%                         end
%                         Lambda11    =   tmp * H1    ;
%                         Lambda12    =   tmp * H2    ;
%     
%                         Lambda11inv     =   inv(Lambda11)   ;
%                         tmp     =   inv(H) * Mu    ;
%                         Muz1    =   tmp(1:k1) - Lambda11inv * Lambda12 * (z2 - tmp((k1+1):k,:) )  ;
%     
%                         R   =   chol(Lambda11inv,'lower')   ;
%                         z1setAlg1   =   bsxfun(@plus, Muz1, R * randn(k1,SampleNum)  )    ;
%                         xsetAlg1    =   bsxfun(@plus, H1 * z1setAlg1, H2 * z2 )    ;
%     
%                         TimeAlg1Set(trial)    =    toc    ;
%     
%                         %%%%%%%%%%%%%%%%%%%%     Algorithm 2   %%%%%%%%%%%%%%%%%%
%                         tic
%     
%                         if DiagSigma
%                             yset    =   bsxfun(@plus, Mu, bsxfun(@times, sqrt(Sigma), randn(k,SampleNum) )  )    ;
%                             SigmaGT     =   bsxfun(@times, Sigma , G.' )    ;
%                             tmp     =   SigmaGT * inv(G * SigmaGT)  ;
%                             xsetAlg2    =   bsxfun(@plus, tmp * r , (eye(k) - tmp * G) * yset )   ;
%                         else
%                             R   =   chol(Sigma,'lower')   ;
%                             yset    =   bsxfun(@plus, Mu, R * randn(k,SampleNum)  )    ;
%                             SigmaGT     =   Sigma * G.'     ;
%                     %         alphaset    =   inv( G * SigmaGT ) * bsxfun(@minus, r, G * yset)     ;
%                     %         xsetAlg2    =   yset + SigmaGT * alphaset   ;
%                             tmp     =   SigmaGT * inv(G * SigmaGT)  ;
%                             xsetAlg2    =   bsxfun(@plus, tmp * r , (eye(k) - tmp * G) * yset )   ;
%                         end
%     
%                         TimeAlg2Set(trial)    =   toc    ;
%     
%                         %%%%%%%%%%%%%%%%%%%%     Patches for correction    %%%%%%%%%%%%%%%%%%
%                         if 0 & (k2settings == k2settingsSet(1)) & (k == 5000) & (DiagSigma==1)
%                             folderPatch     =   'Figures/Patches/'  ;
%                             for ii  =   1:50
%                                 dim     =   randperm(k,2)   ;   nbins   =   [10 10]     ;
%                                 xmin    =	mean(xsetAlg1(dim(1),:)) - 4 * std(xsetAlg1(dim(1),:))    ;
%                                 xmax    =   mean(xsetAlg1(dim(1),:)) + 4 * std(xsetAlg1(dim(1),:))    ;
%                                 ymin    =   mean(xsetAlg1(dim(2),:)) - 4 * std(xsetAlg1(dim(2),:))    ;
%                                 ymax    =   mean(xsetAlg1(dim(2),:)) + 4 * std(xsetAlg1(dim(2),:))    ;
%                                 if xmax - xmin > ymax - ymin
%                                     tmp     =   (ymax+ymin)/2 - (xmax - xmin)/2     ;
%                                     tmp1    =   (ymax+ymin)/2 + (xmax - xmin)/2     ;
%                                     ymin    =    tmp    ;   ymax    =   tmp1    ;
%                                 else
%                                     tmp     =   (xmax+xmin)/2 - (ymax - ymin)/2     ;
%                                     tmp1    =   (xmax+xmin)/2 + (ymax - ymin)/2     ;
%                                     xmin    =    tmp    ;   xmax    =   tmp1    ;
%                                 end
% %                                 figure(1),subplot(121),plot(xsetAlg1(dim(1),:),xsetAlg1(dim(2),:),'r.'),axis([xmin xmax ymin ymax])
% %                                 figure(1),subplot(122),plot(xsetAlg2(dim(1),:),xsetAlg2(dim(2),:),'r.'),axis([xmin xmax ymin ymax]),pause(0.3)
%                                 [Z,Xedges,Yedges]   =   histcounts2(xsetAlg1(dim(1),:),xsetAlg1(dim(2),:),nbins)    ;
%                                 Xedges  =   Xedges(1:(end-1)) + 0.5 * (Xedges(2) - Xedges(1))   ;
%                                 Yedges  =   Yedges(1:(end-1)) + 0.5 * (Yedges(2) - Yedges(1))   ;
%                                 [X,Y]   =   meshgrid(Yedges,Xedges)     ;       
%                                 figure(9),axes('Position',[0,0,1,1]),contour(X,Y,Z,'LineWidth',linewidthsize),%axis([xmin xmax ymin ymax])
%                                 axis([-3 3 -3 3]),
%                                 filename    =   ['ZA1A2Trial',num2str(trial),'D',num2str(dim(1)),'D',num2str(dim(2)),'A1']     ;
%                                 saveas(gcf,[folderPatch,filename],'pdf')     ; saveas(gcf,[folderPatch,filename],'jpeg')     ;    
%                                 close ;
%                                 [Z,Xedges,Yedges]   =   histcounts2(xsetAlg2(dim(1),:),xsetAlg2(dim(2),:),nbins)    ;
%                                 Xedges  =   Xedges(1:(end-1)) + 0.5 * (Xedges(2) - Xedges(1))   ;
%                                 Yedges  =   Yedges(1:(end-1)) + 0.5 * (Yedges(2) - Yedges(1))   ;
%                                 [X,Y]   =   meshgrid(Yedges,Xedges)     ;      
%                                 figure(9),axes('Position',[0,0,1,1]),contour(X,Y,Z,'LineWidth',linewidthsize),%axis([xmin xmax ymin ymax])
%                                 axis([-3 3 -3 3]),
%                                 filename    =   ['ZA1A2Trial',num2str(trial),'D',num2str(dim(1)),'D',num2str(dim(2)),'A2']     ;
%                                 saveas(gcf,[folderPatch,filename],'pdf')     ; saveas(gcf,[folderPatch,filename],'jpeg')     ;    
%                                 close ;
%                             end
%                         end
% %                         a=1;
%                     end
%     
%                     AverTimeAlg1    =   [AverTimeAlg1 , mean(TimeAlg1Set)]  ;
%                     AverTimeAlg2    =   [AverTimeAlg2 , mean(TimeAlg2Set)]  ;
%     
%             %     end
%             end
%     
%             figure(8),hold on;box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
%             plot(kset,AverTimeAlg1,'r*','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
%             plot(kset,AverTimeAlg2,'bo','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
%             xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
%             hh = legend('Algorithm 1','Algorithm 2');    set(hh,'Location','northwest','FontSize',23)
%     
%             pNaive  =   polyfit(kset,AverTimeAlg1,4)    ;  pMargCut    =   polyfit(kset,AverTimeAlg2,4)    ;
%             x   =   kset(1):kset(end)   ;
%             yNaive  =   polyval(pNaive,x)   ;   yMargCut  =   polyval(pMargCut,x)   ;  
%             plot(x,yNaive,'r-','LineWidth',linewidthsize*1.2),
%             plot(x,yMargCut,'b-','LineWidth',linewidthsize*1.2),
%     
%             if DiagSigma
%                 tmp1     =   num2str(k2settings)     ;   tmp1(tmp1 == '.')     =   []  ;
%                 filename    =   ['A1A2K2',tmp1,'DiagCov']     ;
%                 tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA1 = AverTimeAlg1 ;']   ;
%                 eval(tmp)   ;
%                 tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA2 = AverTimeAlg2 ;']   ;
%                 eval(tmp)   ;            
%             else
%                 tmp1     =   num2str(k2settings)     ;   tmp1(tmp1 == '.')     =   []  ;
%                 filename    =   ['A1A2K2',tmp1,'GenCov']     ;
%                 tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA1 = AverTimeAlg1 ;']   ;
%                 eval(tmp)   ;
%                 tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA2 = AverTimeAlg2 ;']   ;
%                 eval(tmp)   ;  
%             end
%             saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ;  
%             close ;
%         end
%     end
%     
%     filename    =   ['DimTimeAlg1Alg2.mat']     ;
%     save([folder,filename],'DataOut')     ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;
    load('E:\0Research\[201606]SimuGauLinearEqConst\Matlab Code\Figures\DimTimeAlg1Alg2.mat')
    k2settingsSet  =   [20,0.1:0.2:0.9]    ;    folder  =   'Figures/'  ;

    kset    =   DataOut.kset    ;

    AverTimeAlg1    =   DataOut.TimeDimK220GenCovA1     ;
    AverTimeAlg2    =   DataOut.TimeDimK220GenCovA2     ;
    figure(8),hold on;box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    plot(kset,AverTimeAlg1,'r*','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    plot(kset,AverTimeAlg2,'bo','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    hh = legend('Algorithm 1','Algorithm 2');    set(hh,'Location','northwest','FontSize',23)
    ylim([0 20])
    pNaive  =   polyfit(kset,AverTimeAlg1,4)    ;  pMargCut    =   polyfit(kset,AverTimeAlg2,4)    ;
    x   =   kset(1):kset(end)   ;
    yNaive  =   polyval(pNaive,x)   ;   yMargCut  =   polyval(pMargCut,x)   ;  
    plot(x,yNaive,'r-','LineWidth',linewidthsize*1.2),
    plot(x,yMargCut,'b-','LineWidth',linewidthsize*1.2),
    tmp1     =   num2str(k2settings)     ;   tmp1(tmp1 == '.')     =   []  ;
    filename    =   ['A1A2K220GenCov']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ;  
    close ;

    AverTimeAlg1    =   DataOut.TimeDimK220DiagCovA1     ;
    AverTimeAlg2    =   DataOut.TimeDimK220DiagCovA2     ;
    figure(8),hold on;box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    plot(kset,AverTimeAlg1,'r*','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    plot(kset,AverTimeAlg2,'bo','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    hh = legend('Algorithm 1','Algorithm 2');    set(hh,'Location','northwest','FontSize',23)
    ylim([0 20])
    pNaive  =   polyfit(kset,AverTimeAlg1,4)    ;  pMargCut    =   polyfit(kset,AverTimeAlg2,4)    ;
    x   =   kset(1):kset(end)   ;
    yNaive  =   polyval(pNaive,x)   ;   yMargCut  =   polyval(pMargCut,x)   ;  
    plot(x,yNaive,'r-','LineWidth',linewidthsize*1.2),
    plot(x,yMargCut,'b-','LineWidth',linewidthsize*1.2),
    tmp1     =   num2str(k2settings)     ;   tmp1(tmp1 == '.')     =   []  ;
    filename    =   ['A1A2K220DiagCov']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ;  
    close ;

    figure(10),%subplot(221),title('GenA1'),
    hold on,box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA1']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        plot(kset,GainVSkGenCov,MarkerCellSet{i},'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    end,ylim([0 50])
    hh = legend('\it{k_2=0.1k}','\it{k_2=0.3k}','\it{k_2=0.5k}','\it{k_2=0.7k}','\it{k_2=0.9k}'); 
    set(hh,'Location','northwest','FontSize',23)
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA1']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        pNaive  =   polyfit(kset,GainVSkGenCov,4)    ;  
        x   =   kset(1):kset(end)   ;   yNaive  =   polyval(pNaive,x)   ;
        plot(x,yNaive,'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),    
    end
    filename    =   ['Alg1GenCovVaryK2']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ; 
    close ;

    figure(10),%subplot(222),title('DiagA1'),
    hold on,box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA1']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        plot(kset,GainVSkGenCov,MarkerCellSet{i},'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    end,ylim([0 50])
    hh = legend('\it{k_2=0.1k}','\it{k_2=0.3k}','\it{k_2=0.5k}','\it{k_2=0.7k}','\it{k_2=0.9k}'); 
    set(hh,'Location','northwest','FontSize',23)
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA1']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        pNaive  =   polyfit(kset,GainVSkGenCov,4)    ;  
        x   =   kset(1):kset(end)   ;   yNaive  =   polyval(pNaive,x)   ;
        plot(x,yNaive,'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),    
    end
    filename    =   ['Alg1DiagCovVaryK2']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ; 
    close ;

    figure(10),%subplot(223),title('GenA2'),
    hold on,box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA2']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        plot(kset,GainVSkGenCov,MarkerCellSet{i},'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    end,ylim([0 50])
    hh = legend('\it{k_2=0.1k}','\it{k_2=0.3k}','\it{k_2=0.5k}','\it{k_2=0.7k}','\it{k_2=0.9k}'); 
    set(hh,'Location','northwest','FontSize',23)
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA2']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        pNaive  =   polyfit(kset,GainVSkGenCov,4)    ;  
        x   =   kset(1):kset(end)   ;   yNaive  =   polyval(pNaive,x)   ;
        plot(x,yNaive,'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),    
    end
    %small figure
    xmin    =   3000    ;   xmax    =   5000    ;   ymin    =   3   ;   ymax    =   16  ;
    plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],'k:','LineWidth',linewidthsize*1.2)
    axes('Position',[0.5,0.5,0.3,0.2]),hold on,box on
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'GenCovA2']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        plot(kset,GainVSkGenCov,MarkerCellSet{i},'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),
        pNaive  =   polyfit(kset,GainVSkGenCov,4)    ;  
        x   =   kset(1):kset(end)   ;   yNaive  =   polyval(pNaive,x)   ;
        plot(x,yNaive,'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),    
    end,axis([xmin,xmax,ymin,ymax])
    filename    =   ['Alg2GenCovVaryK2']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ; 
    close ;

    figure(10),%subplot(224),title('DiagA2'),
    hold on,box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA2']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        plot(kset,GainVSkGenCov,MarkerCellSet{i},'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    end,ylim([0 50])
    hh = legend('\it{k_2=0.1k}','\it{k_2=0.3k}','\it{k_2=0.5k}','\it{k_2=0.7k}','\it{k_2=0.9k}'); 
    set(hh,'Location','northwest','FontSize',23)
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA2']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        pNaive  =   polyfit(kset,GainVSkGenCov,4)    ;  
        x   =   kset(1):kset(end)   ;   yNaive  =   polyval(pNaive,x)   ;
        plot(x,yNaive,'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),    
    end
    %small figure
    xmin    =   3000    ;   xmax    =   5000    ;   ymin    =   3   ;   ymax    =   16  ;
    plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],'k:','LineWidth',linewidthsize*1.2)
    axes('Position',[0.5,0.5,0.3,0.2]),hold on,box on
    for i   =   1:(length(k2settingsSet)-1)
        tmp1     =   num2str(k2settingsSet(i+1) )     ;   tmp1(tmp1 == '.')     =   []  ;
        tmp     =   ['DataOut.TimeDimK2',tmp1,'DiagCovA2']   ;
        GainVSkGenCov   =   eval(tmp)   ;
        plot(kset,GainVSkGenCov,MarkerCellSet{i},'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),
        pNaive  =   polyfit(kset,GainVSkGenCov,4)    ;  
        x   =   kset(1):kset(end)   ;   yNaive  =   polyval(pNaive,x)   ;
        plot(x,yNaive,'Color',ColorSet(i,:),'LineWidth',linewidthsize*1.2),    
    end,axis([xmin,xmax,ymin,ymax])
    filename    =   ['Alg2DiagCovVaryK2']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ; 
    close
end


%% Figure 6: Time vs Dimension 
if 0
    VSet   =   [100 1000 5000 10000]    ;
    TimeNaive   =   zeros(size(VSet))   ;   TimeMargCut     =   zeros(size(VSet))   ;
    for i   =   1:length(VSet)
        V   =   VSet(i)     ;
        folder  =   ['Figures/V',num2str(V),'/']  ;
        filename    =   ['DimTimeV',num2str(V),'.mat']     ;
        load([folder,filename])     ; 

        tmp     =   mean(TimeSet,2)     ;
        TimeNaive(i)    =   tmp(1)      ;
        TimeMargCut(i)  =   tmp(2)      ;
    end
    pNaive  =   polyfit(VSet,TimeNaive,2)    ;  pMargCut    =   polyfit(VSet,TimeMargCut,2)    ;
    x   =   VSet(1):VSet(end)   ;
    yNaive  =   polyval(pNaive,x)   ;   yMargCut  =   polyval(pMargCut,x)   ;  

    figure(4),hold on;box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    plot(VSet,TimeNaive,'r*','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    plot(VSet,TimeMargCut,'bo','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    xlabel('Dimension \it{k}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    hh = legend('Cholesky','Truncated-MVN');    set(hh,'Location','northwest','FontSize',23)
    plot(x,yNaive,'r-','LineWidth',linewidthsize*1.2),
    plot(x,yMargCut,'b-','LineWidth',linewidthsize*1.2),

    filename    =   'TimeDimension'    ;
    % saveas(gcf,filename,'pdf')     ;
end


%% Figure 7: V = 10000 dimensions
if 0
    SampleNum   =   10000    ;  fignum  =   10  ;
    DimSet  =   []  ;   TimeSet     =   []  ;

    folder  =   'Figures/V10000/'  ;

    for i   =   1:100
        %%%%%%%%%%%%%%%%%%%%     Naive simulation  %%%%%%%%%%%%%%%%%%
        N   =   9999   ;   V   =   N + 1   ;
        Mu  =   ones(V,1) / V   ;   
        BarMu   =   Mu(1:N,1)   ;
        Phi     =   randg( ones(V,1) )  ;   Phi     =   Phi / sum(Phi,1)    ;
        BarPhi  =   Phi(1:N,1)   ;
        b   =   0.5     ;

        tic
        BarSigmaMat     =   b * ( diag(BarPhi) - BarPhi * BarPhi.' )    ;
        R   =   chol(BarSigmaMat,'lower')   ;
        SampleNaive     =   bsxfun(@plus, BarMu, R * randn(N,SampleNum)  )    ;
        time(1)     =    toc    ;

        dim     =   randperm(N,2)   ;   DimSet  =   [DimSet,dim']   ;
        nbins   =   [10 10]     ;
        [Z,Xedges,Yedges]   =   histcounts2(SampleNaive(dim(1),:),SampleNaive(dim(2),:),nbins)    ;
        Xedges  =   Xedges(1:(end-1)) + 0.5 * (Xedges(2) - Xedges(1))   ;
        Yedges  =   Yedges(1:(end-1)) + 0.5 * (Yedges(2) - Yedges(1))   ;
        [X,Y]   =   meshgrid(Yedges,Xedges)     ;       
        figure(2),subplot(121),contour(X,Y,Z),axis([-0.03 0.03 -0.03 0.03])

        filename    =   ['NaiveContour',num2str(i)]     ;
        figure(3),axes('Position',[0,0,1,1]),contour(X,Y,Z,'LineWidth',linewidthsize),axis([-0.03 0.03 -0.03 0.03])
    %     saveas(gcf,[folder,filename],'pdf')     ; saveas(gcf,[folder,filename],'jpeg')     ;  
        close ;

        %%%%%%%%%%%%%%%%%%%%   Marg-Cut simulation  %%%%%%%%%%%%%%%%%%
        V   ;
        Mu  ;   
        Phi ;
        b 	;

        SigmaMat        =   b * diag(Phi)   ;

        tic
        SampleY      =   bsxfun(@plus, Mu,   bsxfun(@times, sqrt(b * Phi), randn(V,SampleNum)  )   )  ;
        SampleX      =   SampleY   +   Phi * (1 - sum(SampleY,1))    ;
        SampleMargCut   =   SampleX(1:N,:)  ;
        time(2)     =   toc     ;   TimeSet  =   [TimeSet,time']   ;


        [Z,Xedges,Yedges]   =   histcounts2(SampleMargCut(dim(1),:),SampleMargCut(dim(2),:),nbins)    ;
        Xedges  =   Xedges(1:(end-1)) + 0.5 * (Xedges(2) - Xedges(1))   ;
        Yedges  =   Yedges(1:(end-1)) + 0.5 * (Yedges(2) - Yedges(1))   ;
        [X,Y]   =   meshgrid(Yedges,Xedges)     ;      
        figure(2),subplot(122),contour(X,Y,Z),axis([-0.03 0.03 -0.03 0.03]);

        filename    =   ['MargCutContour',num2str(i)]     ;
        figure(3),axes('Position',[0,0,1,1]),contour(X,Y,Z,'LineWidth',linewidthsize),axis([-0.03 0.03 -0.03 0.03])
    %     saveas(gcf,[folder,filename],'pdf')     ;  saveas(gcf,[folder,filename],'jpeg')     ;  
        close ;

        fprintf('V: %d, Trial: %d, Dim: %d-%d, Time: %1.2d------%1.2d. \n',V,i,dim(1),dim(2),time(1),time(2))  ; 

        pause(0.01)
    end

    filename    =   ['DimTimeV',num2str(V),'.mat']     ;
    % save([folder,filename],'DimSet','TimeSet')     ; 
end


%% Figure 8: [Alg 3]  Time vs Dimensions
if 0
    folder  =   'Figures/'  ;
    SampleNum   =   1    ;  TrialNum    =   50  ;

    k1      =   4000    ;
    k2set   =   [1,500,1000,2000,3000,4000,5000,6000,7000,8000]     ;%

%     TimeAlg3NaiveSet    =   zeros(TrialNum,length(k2set))   ;   
%     TimeAlg3Set         =   zeros(TrialNum,length(k2set))   ;
%     
%     for k2i     =   1:length(k2set)    ;
%         k2      =   k2set(k2i)  ;
%         k   =   k1 + k2;
%         
%         [U,S,V]     =   svd(randn(k,k))     ;
% 
%         for trial = 1:TrialNum
% 
%             fprintf('Alg 3: K2: %d, Trial: %d\n',k2,trial)  ;
%             
%             %%%%%%%%%%%%%%%%%%%%     Generate Mu1,Sigma11,Sigma12,Sigma22,Sigma21   %%%%%%%%%%%%%%%%%%
%             Mu1     =   randn(k1,1)     ;
% %             U   =   orth(rand(k,k))   ;
% %             U1  =   U(:,1:k1)   ;   U2  =   U(:,(k1+1):end)     ;
% %             S   =   diag( max(0.1, rand(k,1) ) )     ;
% %             Sigma11     =   U1.' * S * U1  ;    
% %             Sigma22     =   U2.' * S * U2  ;
% %             [u,s,v]     =   svd(Sigma11)    ;   
% %             U1  =   U1 * u  ;
% %             [u,s,v]     =   svd(Sigma22)    ;   
% %             U2  =   U2 * u  ;
% %             Sigma11     =   U1.' * S * U1  ;
% %             Sigma12     =   U1.' * S * U2  ;
% %             Sigma21     =   U2.' * S * U1  ;
% %             Sigma22     =   U2.' * S * U2  ;
% %                 A       =   max(0.05, rand(k1,1) )    ;
% %                 C       =   max(0.05, rand(k2,1) )    ;
% %                 U       =   randn(k1,k2)  ;
% %                 Sigma11     =   1./A ;
% %                 Sigma22     =   bsxfun(@plus, 1./C , U.' * diag(1./A) * U )  ;
% %                 Sigma12     =   bsxfun(@times, 1./A , U) ;
% %                 Sigma21     =   Sigma12.'  ;
%             iii     =   1   ;
%             while 1
% %                 [U,S,V]     =   svd(randn(k,k))     ;
%                 U1   =   U(:,randperm(k,k1))   ;   V1   =   V(:,randperm(k,k2))   ;
%                 Sigma11     =   U1.' * U1     ;   
%                 Sigma11     =   diag(Sigma11) + rand(k1,1) * 0.1  ;
%                 Sigma22     =   V1.' * V1     ;   
%                 Sigma22     =   diag(Sigma22) + rand(k2,1) * 0.1  ;
%                 Sigma12     =   U1.' * V1     ;
%                 Sigma21     =   V1.' * U1     ;                
% 
%                 [L,p1]   =   chol( diag(Sigma11) - Sigma12 * inv(diag(Sigma22)) *  Sigma21 ,'lower')   ;
%                 [L,p2]   =   chol( diag(Sigma22) - Sigma21 * inv(diag(Sigma11)) *  Sigma12  , 'lower')  ;
%                 if (p1 == 0) & (p2 == 0)
%                     break ;
%                 end
%                 iii     =   iii + 1
%             end
%             %%%%%%%%%%%%%%%%%%%%     Algorithm 3 Naive   %%%%%%%%%%%%%%%%%%
%             tic
%             tmp     =   diag(Sigma11) - bsxfun(@times, Sigma12, 1./Sigma22.') *  Sigma21   ;
%             R   =   chol(tmp,'lower')   ;
%             x1setAlg3Naive  =   bsxfun(@plus, Mu1, R * randn(k1,SampleNum)  )    ;
%             TimeAlg3NaiveSet(trial,k2i)    =    toc    ;
% 
%             %%%%%%%%%%%%%%%%%%%%     Algorithm 3   %%%%%%%%%%%%%%%%%%
%             tic
%             y1set   =   bsxfun(@times, sqrt(Sigma11), randn(k1,SampleNum) )     ;
%             tmp1    =   bsxfun(@times, Sigma21, 1./Sigma11.')   ;
%             y2set   =   chol( diag(Sigma22) - tmp1 *  Sigma12  , 'lower')  * randn(k2,SampleNum)      ;
%             x1setAlg3   =   bsxfun(@plus, Mu1, y1set - bsxfun(@times, Sigma12, 1./Sigma22.') * (tmp1*y1set + y2set ) )   ;
%             TimeAlg3Set(trial,k2i)    =   toc    ;
% 
%             %%%%%%%%%%%%%%%%%%%%     Patches for correction    %%%%%%%%%%%%%%%%%%
%             if 0 & (SampleNum > 1000)
%                 for ii  =   1:50
%                     dim     =   randperm(k1,2)   ;   nbins   =   [10 10]     ;
%                     xmin    =	mean(x1setAlg3Naive(dim(1),:)) - 4 * std(x1setAlg3Naive(dim(1),:))    ;
%                     xmax    =   mean(x1setAlg3Naive(dim(1),:)) + 4 * std(x1setAlg3Naive(dim(1),:))    ;
%                     ymin    =   mean(x1setAlg3Naive(dim(2),:)) - 4 * std(x1setAlg3Naive(dim(2),:))    ;
%                     ymax    =   mean(x1setAlg3Naive(dim(2),:)) + 4 * std(x1setAlg3Naive(dim(2),:))    ;
%                     if xmax - xmin > ymax - ymin
%                         tmp     =   (ymax+ymin)/2 - (xmax - xmin)/2     ;
%                         tmp1    =   (ymax+ymin)/2 + (xmax - xmin)/2     ;
%                         ymin    =    tmp    ;   ymax    =   tmp1    ;
%                     else
%                         tmp     =   (xmax+xmin)/2 - (ymax - ymin)/2     ;
%                         tmp1    =   (xmax+xmin)/2 + (ymax - ymin)/2     ;
%                         xmin    =    tmp    ;   xmax    =   tmp1    ;
%                     end
%                     figure(1),subplot(121),plot(x1setAlg3Naive(dim(1),:),x1setAlg3Naive(dim(2),:),'r.'),axis([xmin xmax ymin ymax])
%                     figure(1),subplot(122),plot(x1setAlg3(dim(1),:),x1setAlg3(dim(2),:),'r.'),axis([xmin xmax ymin ymax]),pause(0.4)
%                 end
%             end
%             a   =   1   ;
%         end
%     end
%     AverTimeAlg3Naive   =   mean(TimeAlg3NaiveSet,1)  ;
%     AverTimeAlg3        =   mean(TimeAlg3Set,1)  ;
%     filename    =   ['DimTimeAlg3NaiveOur.mat']     ;
%     save([folder,filename],'SampleNum','TrialNum','k1','k2set','TimeAlg3NaiveSet','TimeAlg3Set','AverTimeAlg3Naive','AverTimeAlg3')     ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;
    load('E:\0Research\[201606]SimuGauLinearEqConst\Matlab Code\Figures\DimTimeAlg3NaiveOur.mat')
    
    figure(11),hold on;box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    plot(k2set,AverTimeAlg3Naive,'r*','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    plot(k2set,AverTimeAlg3,'bo','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    xlabel('Dimension \it{k_2}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    hh = legend('Naive','Algorithm 3');    set(hh,'Location','northwest','FontSize',23)
    ylim([0 6])
    
    pNaive  =   polyfit(k2set,AverTimeAlg3Naive,4)    ;  pMargCut    =   polyfit(k2set,AverTimeAlg3,4)    ;
    x   =   k2set(1):k2set(end)   ;
    yNaive  =   polyval(pNaive,x)   ;   yMargCut  =   polyval(pMargCut,x)   ;  
    plot(x,yNaive,'r-','LineWidth',linewidthsize*1.2),
    plot(x,yMargCut,'b-','LineWidth',linewidthsize*1.2),

    filename    =   ['A3FixK1VaryK2DiagCov']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ;  
    close ;
end


%% Figure 9: [Alg 4]  Time vs Dimensions
if 0
    folder  =   'Figures/'  ;
    SampleNum   =   1    ;  TrialNum    =   50  ;

    n       =   4000    ;
    pset    =   [1,500,1000,2000,3000,4000,5000,6000,7000,8000]     ;
    
    TimeAlg4NaiveSet    =   zeros(TrialNum,length(pset))   ;   
    TimeAlg4Set         =   zeros(TrialNum,length(pset))   ;
    
    for pi	=   1:length(pset)    ;
        p	=   pset(pi)  ;
        
        for trial = 1:TrialNum
            
            fprintf('Alg 4: P: %d, Trial: %d\n',p,trial)  ;
   
            %%%%%%%%%%%%%%%%%%%%     Generate Mubeta,A,Omega,Phi   %%%%%%%%%%%%%%%%%%
            Mubeta     =   randn(p,1)     ;
            while 1
                A       =   max(0.05, rand(p,1) )    ;
                Omega   =   max(0.05, rand(n,1) )    ;
                Phi     =   randn(n,p)  ;
                
                [L,p1]   =   chol( inv( diag(A) + Phi.' * diag(Omega) * Phi) ,'lower')   ;
                if p1 == 0
                    break ;
                end
            end
            %%%%%%%%%%%%%%%%%%%%     Algorithm 4 Naive   %%%%%%%%%%%%%%%%%%
            tic
            tmp     =   inv( diag(A) + Phi.' * diag(Omega) * Phi)   ;
            R   =   chol(tmp,'lower')   ;
            betasetAlg4Naive  =   bsxfun(@plus, Mubeta, R * randn(p,SampleNum)  )    ;
            TimeAlg4NaiveSet(trial,pi)    =    toc    ;

            %%%%%%%%%%%%%%%%%%%%     Algorithm 4   %%%%%%%%%%%%%%%%%%
            tic
            Ainv    =   1./A    ;   Omegainv    =   1./Omega    ;
            y1set   =   bsxfun(@times, sqrt(Ainv), randn(p,SampleNum) )     ;
            y2set   =   bsxfun(@times, sqrt(Omegainv), randn(n,SampleNum) )     ;
            tmp     =   bsxfun(@times, Ainv, Phi.' )    ;
            betasetAlg4    =   bsxfun(@plus, Mubeta, y1set - tmp*inv(diag(Omegainv)+Phi*tmp)*(Phi*y1set+y2set) )  ;
            TimeAlg4Set(trial,pi)    =   toc    ;

            %%%%%%%%%%%%%%%%%%%%     Patches for correction    %%%%%%%%%%%%%%%%%%
            if 0 & (SampleNum > 1000)
                for ii  =   1:50
                    dim     =   randperm(n,2)   ;   nbins   =   [10 10]     ;
                    xmin    =	mean(betasetAlg4Naive(dim(1),:)) - 4 * std(betasetAlg4Naive(dim(1),:))    ;
                    xmax    =   mean(betasetAlg4Naive(dim(1),:)) + 4 * std(betasetAlg4Naive(dim(1),:))    ;
                    ymin    =   mean(betasetAlg4Naive(dim(2),:)) - 4 * std(betasetAlg4Naive(dim(2),:))    ;
                    ymax    =   mean(betasetAlg4Naive(dim(2),:)) + 4 * std(betasetAlg4Naive(dim(2),:))    ;
                    if xmax - xmin > ymax - ymin
                        tmp     =   (ymax+ymin)/2 - (xmax - xmin)/2     ;
                        tmp1    =   (ymax+ymin)/2 + (xmax - xmin)/2     ;
                        ymin    =    tmp    ;   ymax    =   tmp1    ;
                    else
                        tmp     =   (xmax+xmin)/2 - (ymax - ymin)/2     ;
                        tmp1    =   (xmax+xmin)/2 + (ymax - ymin)/2     ;
                        xmin    =    tmp    ;   xmax    =   tmp1    ;
                    end
                    figure(1),subplot(121),plot(betasetAlg4Naive(dim(1),:),betasetAlg4Naive(dim(2),:),'r.'),axis([xmin xmax ymin ymax])
                    figure(1),subplot(122),plot(betasetAlg4(dim(1),:),betasetAlg4(dim(2),:),'r.'),axis([xmin xmax ymin ymax]),pause(0.4)
                end
            end
            a   =   1   ;
        end
    end
    AverTimeAlg4Naive   =   mean(TimeAlg4NaiveSet,1)  ;
    AverTimeAlg4        =   mean(TimeAlg4Set,1)  ;
    filename    =   ['DimTimeAlg4NaiveOur.mat']     ;
    save([folder,filename],'SampleNum','TrialNum','n','pset','TimeAlg4NaiveSet','TimeAlg4Set','AverTimeAlg4Naive','AverTimeAlg4')     ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;
    load('E:\0Research\[201606]SimuGauLinearEqConst\Matlab Code\Figures\DimTimeAlg4NaiveOur.mat')
    
    figure(11),hold on;box on;set(gca,'fontsize',gcafontsize,'FontName',FontName)
    plot(pset,AverTimeAlg4Naive,'r*','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    plot(pset,AverTimeAlg4,'bo','LineWidth',linewidthsize*1.2,'MarkerSize',Markersize),
    xlabel('Dimension \it{p}','FontSize',23),    ylabel('Average Time (s)','FontSize',23),
    hh = legend('Naive','Algorithm 4');    set(hh,'Location','northwest','FontSize',23)

    pNaive  =   polyfit(pset,AverTimeAlg4Naive,4)    ;  pMargCut    =   polyfit(pset,AverTimeAlg4,4)    ;
    x   =   pset(1):pset(end)   ;
    yNaive  =   polyval(pNaive,x)   ;   yMargCut  =   polyval(pMargCut,x)   ;  
    plot(x,yNaive,'r-','LineWidth',linewidthsize*1.2),
    plot(x,yMargCut,'b-','LineWidth',linewidthsize*1.2),

    filename    =   ['A4FixNVaryPDiagCov']     ;
    saveas(gcf,[folder,filename],'pdf')     ;  % saveas(gcf,[folder,filename],'jpeg')     ;  
    close ;
end


