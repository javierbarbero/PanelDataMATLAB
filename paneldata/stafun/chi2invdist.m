function [ x ] = chi2invdist( p, v )
%CHI2INVDIST Return the inverse cumulative distribution function of a Chi-Square
%   Return the inverse cumulative distribution function of a Chi-Square
%   with v degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    x = 2*gammaincinv(p,v/2);

end

