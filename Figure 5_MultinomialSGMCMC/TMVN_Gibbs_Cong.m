function [Xset]     =   TMVN_Gibbs_Cong(Mu, Sigma, B, b,Nsample,epoch)
% x ~ N_{T}(Mu,Sigma),T={x|Bx \le b}
% Effcient Gibbs Sampling of Truncated Multivariate Normal with Application
% to Constrained Linear Regression 
% By Gabriel Rodriguez-Yam, Richard A. Davis, and Louis L. Scharf March, 2004
% Cong Yulai
% 2016 11 07
if ~exist('Nsample','var')
    Nsample     =   1   ;
end
if ~exist('epoch','var')
    epoch   =   1   ;
end

P   =   size(Sigma,2)  ;

invA    =   chol(Sigma,'lower')     ;
A       =   inv(invA)   ;

Muz     =   A * Mu  ;
Dz      =   B * invA    ;

z   =   Dz \ (b-0.001) * ones(1,Nsample)  ;
for ii  =   1: epoch
    for i   =   1:P
        di  =   Dz(:,i)     ;
        tmp     =   z   ;   tmp(i,:)  =   0   ;
        bi  =   bsxfun(@minus, b , Dz * tmp )   ;

        tmp     =   bsxfun(@times, bi , 1./di )   ;
        tmp1    =   tmp     ;   tmp1(di <= 0,:)   =   inf     ;
        zmax    =   min( tmp1 , [] , 1)     ;
        tmp1    =   tmp     ;   tmp1(di >= 0,:)   =   -inf     ;
        zmin    =   max( tmp1 , [] , 1)     ;

        if sum(zmin > zmax) | sum(bi(di==0) < 0)
            zmax    =    zmin   ;
            warning('Null truncated region of z!')    ;
%             error('Null truncated region of z!')    ;
        end

        z(i,:)    =   Muz(i) + trandn(zmin - Muz(i) , zmax - Muz(i))  ;

    end
end
Xset   =   invA * z    ;






