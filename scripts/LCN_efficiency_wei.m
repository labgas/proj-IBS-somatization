function E = LCN_efficiency_wei(W,local,method)
% LCN_efficiency_wei     
% This function calculates different version for the global and local 
% efficiency, local efficiency.
%
%   The global efficiency is the average of inverse shortest path length,
%   and is inversely related to the characteristic path length.
%
%   The local efficiency is the global efficiency computed on the
%   neighborhood of the node, and is related to the clustering coefficient.
%
%   INPUT
%       W    = weighted undirected or directed connection matrix (all 
%              weights in W must be between 0 and 1).   
%       local = optional argument. 
%               local = 0 computes global efficiency (default)
%               local = 1 computes local efficiency
%       method = 'W'      method  Wang et al.
%                'Z'      method2 Wang et al.
%                'P'      method3 Wang et al.
%                'BCT'    method used in the brain connectivity toolbox
%                if no method is specified, 'P' is taken as default.
%
%   OUTPUT
%       E = global or local efficiency depending on the variable local
%           local = 0: E = Eglob, the global efficiency (scalar)
%           local = 1: E = Eloc, the local efficiency per node (vector)
%
%   Notes:
%   The  efficiency is computed using an auxiliary connection-length
%   matrix L, defined as L_ij = 1/W_ij for all nonzero L_ij; This has an
%   intuitive interpretation, as higher connection weights intuitively
%   correspond to shorter lengths.
%
% author: Yu Wang and Patrick Dupont
% date:   23/03/2015
% history: march 2015 - added different versions for the calculation of the
%                       local efficiency and corrected formula for method
%                       'W'
%          september 2016 avoided the use of nansum (so for this routine 
%                    the statistics toolbox in matlab is no longer needed).
%__________________________________________________________________________
% @(#)LCN_efficiency_wei.m          v0.22         last modified: 2016/09/27

if nargin < 3
   method = 'P';
end
% if (strcmp(method,'W') ~= 1 && strcmp(method,'Z') ~= 1 && strcmp(method,'BCT') ~= 1 && strcmp(method,'W_old') ~= 1)
%    disp(['WARNING: no valid method for LCN_efficiency_wei is specified. Taking the default method ''Z'''])
%    method = 'Z';
% end

n = length(W);                                    %number of nodes
L = W;
A = W~=0;
ind = L~=0;
L(ind) = 1./L(ind);                             %connection-length matrix
MAX_W  = max(W(:));
[index_x,index_y] = find(tril(ones(n-1),-1));
if exist('local','var') && local                %local efficiency
    E = zeros(n,1);
    if strcmp(method,'W') == 1
       for u = 1:n
%            disp(['calculation local efficiency: ' num2str(100*u/n) '% done']); % temporary added to see progress
           V = find(A(u,:)|A(:,u).'); % neighbors
           W_t = W(:,u);
           index_i = 1:n;
           index_i(u) = [];        
           angles(:,1) = W_t(index_i(index_x));
           angles(:,2) = W_t(index_i(index_y));
           angles_min  = min(angles,[],2);
           Len_new = L(V,V).*(1./repmat(W(V,u),1,length(V))).*(1./repmat(W(u,V),length(V),1)); % Len_new = 1./(w_ij.*w_ik.*w_jk)
           e = distance_inv_wei(Len_new); %inverse distance matrix
           se = e+e.';
           Dis = nan(n);      
           Dis(V,V) = se/2;
           Dis(u,:) = [];
           Dis(:,u) = [];
           Dis_v    = Dis((tril(ones(n-1),-1))==1);
           
           angles_min(isnan(angles_min)) = 0;
           Dis_v(isnan(Dis_v)) = 0;
           
           numer = sum((angles_min.^3).*Dis_v); % numerator        

           if max(angles_min)~=0
              denom = sum(MAX_W.*(angles_min.^2)); % denominator
              E(u) = numer/denom;
           end
       end
    elseif strcmp(method,'Z') == 1
       for u = 1:n
%            disp(['calculation local efficiency: ' num2str(100*u/n) '% done']); % temporary added to see progress
           V = find(A(u,:)|A(:,u).'); % neighbors
%            W_t = W(:,u); % older version
           W_t = W(:,u);
           index_i = 1:n;
           index_i(u) = [];        
           angles(:,1) = W_t(index_i(index_x));
           angles(:,2) = W_t(index_i(index_y));           
           angles_prod  = (angles(:,1).^(1/3)).*(angles(:,2).^(1/3));
           Len_new = L(V,V).^(1/3); % The distance is calculated by changing the connection matrix
           e = distance_inv_wei(Len_new); % inverse distance matrix
           se = e+e.';
           Dis = nan(n);      
           Dis(V,V) = se/2;
           Dis(u,:) = [];
           Dis(:,u) = [];
           Dis_v    = Dis((tril(ones(n-1),-1))==1);
           
           angles_prod(isnan(angles_prod)) = 0;
           Dis_v(isnan(Dis_v)) = 0;         
           numer = 2*sum(angles_prod.*Dis_v); % numerator        

           if max(angles_prod)~=0
              denom = sum(MAX_W).^(1/3)*2*sum(angles_prod); % denominator
              E(u) = numer/denom;
           end
       end   
    elseif strcmp(method,'P') == 1
       for u = 1:n
%            disp(['calculation local efficiency: ' num2str(100*u/n) '% done']); % temporary added to see progress
           V = find(A(u,:)|A(:,u).'); % neighbors
%            W_t = W(:,u); % older version
           W_t = W(:,u)./MAX_W;
%            node_strength = sum(W_t);
           index_i = 1:n;
           index_i(u) = [];        
           angles(:,1) = W_t(index_i(index_x));
           angles(:,2) = W_t(index_i(index_y));           
           angles_prod  = angles(:,1).*angles(:,2);
%            Len_new = L(V,V); % The distance is calculated by changing the connection matrix
           Len_new = L(V,V).*(1./repmat(W(V,u),1,length(V))).*(1./repmat(W(u,V),length(V),1)); % Len_new = 1./(w_ij.*w_ik.*w_jk)
           e = distance_inv_wei(Len_new); % inverse distance matrix
           se = e+e.';
           Dis = nan(n);      
           Dis(V,V) = se/2;
           Dis(u,:) = [];
           Dis(:,u) = [];
           Dis_v    = Dis((tril(ones(n-1),-1))==1);
           
           angles_prod(isnan(angles_prod)) = 0;
           Dis_v(isnan(Dis_v)) = 0;
           numer = 2*sum(angles_prod.*Dis_v); % numerator        

           if max(angles_prod)~=0
              denom = 2*sum(angles_prod); % denominator
              E(u) = numer/denom;
           end
       end
    elseif strcmp(method,'BCT') == 1
       E = efficiency_wei(W,local);
    end
else
    E = efficiency_wei(W,0); %global efficiency
end

% the function below is taken from the brain connectivity toolbox from the
% script efficiency_wei.m
function D=distance_inv_wei(W_)
    n_ = length(W_);
    D = inf(n_);                                      %distance matrix
    D(1:n_+1:end) = 0;

    for u=1:n_
        S = true(1,n_);                               %distance permanence (true is temporary)
        W1_ = W_;
        V = u;
        while 1
            S(V)     = 0;                                 %distance u->V is now permanent
            W1_(:,V) = 0;                             %no in-edges as already shortest
            for v = V
                T = find(W1_(v,:));                   %neighbours of shortest nodes
                D(u,T) = min([D(u,T);D(u,v)+W1_(v,T)]);%smallest of old/new path lengths
            end
            minD = min(D(u,S));
            if isempty(minD)||isinf(minD),          %isempty: all nodes reached;
               break,                              %isinf: some nodes cannot be reached
            end;
            V=find(D(u,:)==minD);
        end
    end
    D = 1./D;                                         %invert distance
    D(1:n_+1:end) = 0;