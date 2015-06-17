function [ p ] = normaldist( x, mu, sigma )
%NORMALDIST Return the cumulative distribution function of the normal
%   Return the cumulative distribution function of the normal
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    % Normal(0,1)
    if nargin < 2
        mu = 0;
        sigma = 1;
    end
    
    % Normal(mu,1)
    if nargin < 3
        sigma = 1;
    end
    
    p = 1/2 * erfc(-1/sqrt(2) * ((x-mu)/sigma));

end

