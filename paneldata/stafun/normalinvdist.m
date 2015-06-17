function [ x ] = normalinvdist( p, mu, sigma )
%NORMALINVDIST Return the inverse cumulative distribution function of the normal
%   Return the inverse cumulative distribution function of the normal
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
    
    x = mu - sqrt(2)*sigma*erfcinv(2*p);

end

