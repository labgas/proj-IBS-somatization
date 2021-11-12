function C = LCN_clustering_coef_wu(W,method)

% LCN_clustering_coef_wu.m     
%
% this function calculates the clustering coefficient for weighted
% undirected networks using different generalizations which are reported in
% the literature.
%
%   format: C = LCN_clustering_coef_wu_Yu(W,method);
%
%   The weighted clustering coefficient is the average "intensity" of
%   triangles around a node.
%
%   Input:      
%          W = weighted undirected connection matrix (assuming all weights
%          are positive)
%          method = 'W'      method Wang et al.
%                   'B'      method Barrat et al.
%                   'H'      method Holme et al.
%                   'O'      method Onnela et al. 
%                   'Op_am'  method Opsahl et al. (arithmetric mean)
%                   'Op_gm'  method Opsahl et al. (geometric mean)
%                   'Op_min' method Opsahl et al. (minimum)
%                   'Op_max' method Opsahl et al. (maximum)
%                   'Z'      method Zhang et al.
%                   'Mgm'    method Miyajima et al. (geometric mean)
%                   'Mhm'    method Miyajima et al. (harmonic mean)
%                   'Mmin'   method Miyajima et al. (minimum)
%
%          if method is not specified, 'Mhm' is taken as the default
%
%   Output:     
%          C = clustering coefficient vector (for each node)
%
%   References: 
%              method 'B' 
%              Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, 
%              and Alessandro Vespignani. The architecture of complex 
%              weighted networks. PNAS 101 (2004), 3747-3752
%              method 'H'
%              Holme, P., Park, S. M., Kim, B. J., Edling, C. R. (2007). 
%              Korean university life in a network perspective: Dynamics of 
%              a large affiliation network. Physica A: Statistical 
%              Mechanics and its Applications 373 (2007), 821–830.
%              method 'O' 
%              Onnela, J.-P., Saramaki, J., Kertesz, J., Kaski, K. 
%              Intensity and coherence of motifs in weighted complex 
%              networks. Phys. Rev E 71 (2005), 065103
%              method 'Op'
%              Opsahl, T., Panzarasa, P. Clustering in weighted networks. 
%              Social networks 31 (2009), 155–163.
%              method 'Z'
%              Zhang, Bin and Horvath, Steve, A general framework for 
%              weighted gene co-expression network analysis. Statistical 
%              applications in genetics and molecular biology 4 (2005),1544 
%              methods 'Mgm', 'Mhm' and 'Mmin'
%              Kent Miyajima and Takashi Sakuragawa. Continuous and robust 
%              clustering coefficients for weighted and directed networks 
%              arXiv preprint (2014), 1412.0059
%              method 'W'
%              PhD thesis Yu Wang (2015)
%
% authors: Yu Wang and Patrick Dupont
% date:   26/11/2015
% history: september 2016 avoided the use of nansum (so for this routine 
%                    the statistics toolbox in matlab is no longer needed).
%__________________________________________________________________________
% @(#)LCN_clustering_coef_wu.m        v0.14       last modified: 2016/09/27

N = size(W,1);
W(1:N+1:end) = 0;% set the diagonal to 0 (no self connections)
C = zeros(N,1); % initialize output
if nargin < 2
   method = 'Mhm';
end

if strcmp(method,'W') == 1 
   [index_x,index_y] = find(tril(ones(N-1),-1)); %pairwise index
   MAX_W = max(W(:));
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_min  = min(angles,[],2);

       % numerator
       index_tri = sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangles(:,1:2) = angles;
       triangles(:,3)   = W(index_tri);
       triangles_w_min  = min(triangles,[],2);
       
       angles_min(isnan(angles_min)) = 0;
       triangles_w_min(isnan(triangles_w_min)) = 0;
       C(i) = sum((triangles_w_min).^3)/sum(MAX_W.*(angles_min).^2); 
   end
elseif strcmp(method,'B') == 1
   K    = sum(W~=0,2); 
   S    = sum(W,2); 
   cyc2 = zeros(N,1);
   for i=1:N
       rep_row = repmat(W(i,:),N,1);
       rep_col = repmat(W(:,i),1,N);
       cyc2(i) = sum(sum((rep_row+rep_col).*(rep_row~=0).*(rep_col~=0).*(W~=0)))/2;
   end
   C = cyc2./(S.*(K-1));
elseif strcmp(method,'H') == 1
   denom_W1 = zeros(N,1);
   for i = 1:N
       sa_W = W(i,:);
       denom_W1(i)=sum(sa_W).^2; % denominator
   end
   cyc3 = diag((W^3))./max(W(:));
   denom_W1(cyc3==0) = inf; % if no 3-cycles exist, make C=0
   C = cyc3./ denom_W1; % clustering coefficient
elseif strcmp(method,'Mgm') == 1
   [index_x,index_y]=find(tril(ones(N-1),-1));%pairwise index
   for i=1:N
       W_t     = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_gm = sqrt(sqrt(angles(:,1).*angles(:,2))*max(W(:))); % Geometric mean
       % numerator
       index_tri=sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangles(:,1:2)= angles(:,1:2);
       triangles(:,3)  = W(index_tri);
       triangles_gm = sqrt(sqrt(angles(:,1).*angles(:,2)).*triangles(:,3)); % Geometric mean
       if sum(triangles_gm)==0
          C(i)=0;   % if no 3-cycles exist  make C=0
       else
          angles_gm(isnan(angles_gm)) = 0; 
          C(i)= sum(triangles_gm)./sum(angles_gm);
       end
   end
elseif strcmp(method,'Mhm') == 1
   [index_x,index_y]=find(tril(ones(N-1),-1));%pairwise index
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_hm = HarmonicMean(HarmonicMean(angles(:,1),angles(:,2)),max(W(:))); % Harmonic mean
       % numerator
       index_tri = sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangles(:,1:2) = angles(:,1:2);
       triangles(:,3)   = W(index_tri);
       triangles_hm = HarmonicMean(HarmonicMean(angles(:,1),angles(:,2)),triangles(:,3));%Geometric mean
       if sum(triangles_hm)==0
          C(i)=0;   %if no 3-cycles exist  make C=0
       else
          C(i) = sum(triangles_hm)./sum(angles_hm);
       end
   end
elseif strcmp(method,'Mmin') == 1
   [index_x,index_y]=find(tril(ones(N-1),-1));
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_min  = min(angles,[],2); % minmum of the triplet
       % numerator
       index_tri=sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangles(:,1:2) = angles(:,1:2);
       triangles(:,3)   = W(index_tri);
       triangles_min    = min(triangles,[],2); % minimum of the closed triplet
       if sum(triangles_min)==0
          C(i) = 0;   %if no 3-cycles exist  make C=0
       else
          C(i) = sum(triangles_min)./sum(angles_min);
       end
   end
elseif strcmp(method,'O') == 1
   K    = sum(W~=0,2);            	
   cyc3 = diag((W.^(1/3))^3)./max(W(:));           
   K(cyc3==0) = inf;             % if no 3-cycles exist, make C=0 (via K=inf)
   C    = cyc3./(K.*(K-1));      % clustering coefficient
elseif strcmp(method,'Op_min') == 1
   [index_x,index_y] = find(tril(ones(N-1),-1));
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
    
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_min  = min(angles,[],2);
  
       % numerator
       index_tri = sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangle_closed = (W(index_tri)>0)';
       
       C(i) = sum(angles_min.*triangle_closed)./sum(angles_min);
   end    
elseif strcmp(method,'Op_max') == 1
   [index_x,index_y] = find(tril(ones(N-1),-1));
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
    
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_max  = max(angles,[],2).*prod(angles>0,2);
  
       % numerator
       index_tri = sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangle_closed = (W(index_tri)>0)';
       
       C(i) = sum(angles_max.*triangle_closed)./sum(angles_max);
   end    
elseif strcmp(method,'Op_am') == 1
   [index_x,index_y] = find(tril(ones(N-1),-1));
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
    
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_mean = mean(angles,2).*prod(angles>0,2); % take only into account if both angles are nonzero
  
       % numerator
       index_tri = sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangle_closed = (W(index_tri)>0)';
       
       C(i) = sum(angles_mean.*triangle_closed)./sum(angles_mean);
   end    
elseif strcmp(method,'Op_gm') == 1
   [index_x,index_y] = find(tril(ones(N-1),-1));
   for i=1:N
       W_t = W(:,i);
       index_i = 1:N;
       index_i(i) = [];
    
       % denominator
       angles(:,1) = W_t(index_i(index_x));
       angles(:,2) = W_t(index_i(index_y));
       angles_gm   = sqrt(prod(angles,2));
  
       % numerator
       index_tri = sub2ind([N,N],index_i(index_x),index_i(index_y));
       triangle_closed = (W(index_tri)>0)';
       
       C(i) = sum(angles_gm.*triangle_closed)./sum(angles_gm);
   end    
elseif strcmp(method,'Z') == 1
   denom_W1 = zeros(N,1);
   for i = 1:N
       sa_W = W(i,:);
       denom_W1(i)=sum(sa_W).^2 - sum(sa_W.^2); % denominator
   end
   cyc3 = diag((W^3))./max(W(:));
   denom_W1(cyc3==0) = inf; % if no 3-cycles exist, make C=0
   C = cyc3./ denom_W1; % clustering coefficient
else
    C = [];
    disp(['ERROR: method ' method ' not implemented'])
end

function HM = HarmonicMean(x,y)
HM = 2./(1./x + 1./y);