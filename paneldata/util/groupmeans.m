function [ means ] = groupmeans( id, X, varargin )
%GROUPMEANS Internal Function
%   Internal Function
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % Change: Added the option 'fillval' to 'NaN' to fill extra output with
    % NaN instead of 0. And then, removes the way. This prevents more group
    % means appearing if id has jumps. Ex: 1, 2, 4, 5.

    % Parse Additional options
    p = inputParser;
    if verLessThan('matlab', '8.2')
        addPar = @(v1,v2,v3,v4) addParamValue(v1,v2,v3,v4);
    else
        addPar = @(v1,v2,v3,v4) addParameter(v1,v2,v3,v4);
    end
    addPar(p,'replicate',0,@(x) isnumeric(x));
    p.parse(varargin{:})
    options = p.Results;
    
    if size(X,2) == 1
        means = accumarray(id,X,[],@mean,NaN); 
    else
        % X is a matrix
        means = cell2mat(arrayfun(@(j) accumarray(id,X(:,j),[],@mean,NaN),1:size(X,2),'UniformOutput', false));
    end
    
    % Remove NaNs rows if there are NaNs in all columns.
    means(all(isnan(means),2),:) = [];
    
    % Replicate elements
    if options.replicate
        % Get Tid for each id
        Tid = accumarray(id,id,[],@(j) size(j,1),NaN);
        
        % Remove NaNs that appear if there are jumps in id's
        Tid(isnan(Tid)) = [];
        
        % Replicate
        means = replicate_elements(means,Tid);
    end

end

