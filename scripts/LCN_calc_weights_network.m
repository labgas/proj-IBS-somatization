function w = LCN_calc_weights_network(Z,option)         
%
% calculates the weigth of a network described by the matrix of Z-scores
% coming from a Fisher r-to-Z transform of the (partial) correlations among
% timeseries of nodes.
%
% FORMAT w = LCN_calc_weights_network(Z,option)
%
% input: 
%   Z matrix of Z-scores
%   option
%       1 = weighted according to w = (2*(normcdf(Z,0,1)-0.5)).^4;
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	February, 2015
% history: 	september 2016 modified to work without the statistics toolbox
%                          as well
%__________________________________________________________________________
% @(#)LCN_calc_weights_network.m	0.2           last modified: 2016/09/24

if option == 1
   stats_toolbox = license('checkout','statistics_toolbox');
   if stats_toolbox
      w = (2*(normcdf(Z,0,1)-0.5)).^4;
   else
      w = (2*(LCN_normcdf(Z,0,1)-0.5)).^4;
   end
   w(1:size(w,1)+1:end) = 0;  %clear diagonal
else
   disp('not yet implemented')
end

end
