function [ isti, diff ] = istinvariant( id,  X, tol )
%ISTINVARIANT Internal Function
%   Internal Function
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%
    % Default tolerance
    if nargin < 4
        tol = 1e-10;
    end

    % Compute group means and replicate for all observations
    Xbar = groupmeans(id,X,'replicate',1);
    
    % Substract means
    diff = X - Xbar;
    
    % Check if all are zero
    isti = all(abs(diff) < tol);

end

