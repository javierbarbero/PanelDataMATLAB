function [ X, Xnames ] = extracttable( X )
%EXTRACTTABLE Internal Function
%   Internal Function
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % If MATLAB R2013b (8.2) or higher
    if ~verLessThan('matlab', '8.2')
        if istable(X)
            % Extract variable names only if the size is equal to the size
            % of X
            Xnames = X.Properties.VariableNames;
            if size(Xnames) ~= size(X,2)
                Xnames = [];
            end
            
            % Convert table to array
            X = table2array(X);
        else
            Xnames = [];
        end
    else
        Xnames = [];
    end

end

