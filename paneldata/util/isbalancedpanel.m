function [ isBalanced, idx, n, T, Tid, Tmean, Thmean] = isbalancedpanel( id, time )
%ISBALANCEDPANEL Internal Function
%   Internal Function
%
%   Copyright 2013-2017 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, August, 2017
%
    % Get uniques
    id_uniq = unique(id);
    time_uniq = unique(time);
    idt = [id time]; 
    
    % Check no id and time is repeated
    if size(idt,1) ~= size(unique(idt,'rows'))
        error('Multiple observations with same id and time');
    end    
    
    % Sort panel
    [~, idx] = sortrows(idt,[1 2]);
    
    % Get number of units
    n = length(id_uniq);
    T = length(time_uniq);
    
    % Get number of observations
    N = size(idt,1);     
    
    % Check balanced
    if N == (n * T)
        isBalanced = 1;
        Tid = (T*ones(n,1));
        Tmean = T;
        Thmean = T;
    else
        isBalanced = 0;
        % Compute number of time periods per ID
        %Tid = nan(n,1);
        %{
        for i=1:n
            Tid(i,1) = sum(id==i);
        end
        %}
        Tid = accumarray(id,id,[],@(j) size(j,1),NaN) ;
        
        % Remove NaN's generated if there are gaps in ID's
        Tid(isnan(Tid)) = [];
        
        % Compute mean and harmonic mean
        Tmean = mean(Tid);
        Thmean = n ./ sum(1./Tid);
    end

end

