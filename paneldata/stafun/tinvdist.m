function [ t ] = tinvdist( p, v )
%TINVDIST Return the inverse cumulative distribution function of a student-t
%   Return the inverse cumulative distribution function of a student-t with
%   v degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    % x = betaincinv(2*min(p,1-p), 1/2*v, 1/2);
    % t = sqrt(v*(1-x)./x);
    x = 2*min(p,1-p);
    t = sqrt(v*(1-betaincinv(x, 1/2*v, 1/2))./betaincinv(x, 1/2*v, 1/2));

end

