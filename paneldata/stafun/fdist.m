function [ p ] = fdist( f, v1, v2 )
%FDIST Return the cumulative distribution function of the F
%   Return the cumulative distribution function of the F with v1 and v2
%   degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    x = (v1*f)./(v2+v1*f);
    p = betainc(x, 1/2*v1, 1/2*v2);

end

