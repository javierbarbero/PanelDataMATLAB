function [ f ] = finvdist( p, v1, v2 )
%FINVDIST Return the inverse cumulative distribution function of the F
%   Return the inverse cumulative distribution function of the F with v1
%   and v2 degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    f =  (v2*betaincinv(p,1/2*v1,1/2*v2))./(v1*(1-betaincinv(p,1/2*v1,1/2*v2)));

end

