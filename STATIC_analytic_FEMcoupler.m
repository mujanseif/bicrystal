function [uhat,fend,Ubar, f_dln] = STATIC_analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
    gamma_disp, gammat, gamma_mixed, fixedDofs,freeDofs,dx,t,...
    gamma_dln, x3x6, n_nodes, n_nodes_t, n_se, idxi, ...
    f_dln_node, f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, eps)

%Coupling of FEM and DDD
% u = uhat + utilda
% f = f_hat + ftilda

% segments = constructsegmentlist(rn,links);
% Udot = 1*1E3*dx*(1E-4/160E9); %for tungsten...
% Udot = 100*1E3*dx*(1E-4/160E9); %for tungsten...
%Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle
Udot = 0;
Ubar = Udot*t;
%Ubar = 0.1*1E4; for debucontourfggin
u=zeros(3*(mno),1);

% u(3*gamma_mixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes

uhat=zeros(3*mno,1);
% utilda=zeros(3*mno,1);

% gn = gamma_disp(:,1); % global node number
% x0 = xnodes(gn,1:3); % field point
% point_array_length = size(x0,1);
% segments_array_length = size(segments,1);

%Matlab wrapper
% tic;
% displacements = displacement_fivel(x0,segments,NU); %utilda array ordered by x0
% toc;

%Full C version (including combinatorials, checks, corrections etc.)
% tic;
% [Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
%                        segments(:,3), segments(:,4), segments(:,5),... %burgers vector
%                        segments(:,6), segments(:,7), segments(:,8),... %start node segs
%                        segments(:,9), segments(:,10), segments(:,11),... %end node segs
%                        segments(:,12), segments(:,13), segments(:,14),... %slip plane
%                        NU,point_array_length,segments_array_length);
% displacements = horzcat(Ux,Uy,Uz);
% toc;
% disp(displacementsMEX-displacements);
%  [Ux, Uy, Uz]  = displacement_fivel(x0,segments,NU);

% utilda(3*gn -2) = Ux;
% utilda(3*gn -1) = Uy;
% utilda(3*gn   ) = Uz;
% 
% if any(isnan(utilda))
%     disp('some entries of utilda are NaN')
%     pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
% end

% uhat(:) = u(:);% - utilda(fixedDofs);
uhat(fixedDofs) = u(fixedDofs);% - utilda(fixedDofs);

[x1x2, b, n_dln] = extract_dislocation_nodes(rn, links);
f_dln(:,1) = 0;
f_hat(:,1) = 0;
[f_dln,~] = analytic_traction(x3x6 , x1x2, b, n_nodes, n_nodes_t,...
                        n_se, n_dln, 3*gamma_dln, idxi,...
                        f_dln_node, f_dln_se, f_dln,...
                        MU, NU, a, use_gpu, n_threads, para_scheme, eps);

% Apply all of these to a single surface
% f_hat(gamma_disp) = -f_dln(gamma_disp)';% no applied forces
f_hat(freeDofs) = -f_dln(freeDofs)';% no applied forces

% f    = f_hat-kg(:,gamma_disp)*uhat(gamma_disp);
f    = f_hat-kg(:,fixedDofs)*uhat(fixedDofs);

bcwt = mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

% f(gamma_disp) = bcwt*uhat(gamma_disp);
% f(fixedDofs) = bcwt*uhat(fixedDofs);
uhat = U\(L\f); %using LU decomposition

rhat=kg*uhat;

% fend = rhat(gamma_disp)+f_dln(gamma_disp);
fend = rhat;
% fend = sum(fend);
end
