function [f,df,d2f] = wey_approx(K,pi,pj,p_approx)
% WEY_APPROX Computes the function and the first and second order derivatives 
%   of gas flow through pipelines equations.    
%   [F,DF,D2F] = WEY_APPROX(K,PI,PJ,P_APPROX)
% 
%   This function also calculates the approximation of the equation in the 
%   point of non-existence of the first derivative, through a polynomial� 
%   approximation proposed in the formulation. Inputs are as follow:
%   
%   K - Weymouth constant vector for all pipelines.
% 
%   PI - Quadratic pressure vector at node from for all pipelines.
% 
%   PJ - Quadratic pressure vector at node for of all pipelines.
%   
%   P_APPROX (optional) - Scalar which indicates the bounds where the 
%       function gets non-differentiable.
%   
%   For more details on the derivations behind the derivative code used
%   in MPNG information, see the user's manual.
%
%   See also MPNG
%
%   Main author: Wilson Gonz�lez Vanegas M.Sc

%   MPNG Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio Garc�a-Mar�n - Universidad Nacional de Colombia - Sede Manizales
%   Wilson Gonz�lez-Vanegas - Universidad Tecnol�gica de Pereira
%   Carlos E. Murillo-S�nchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

%% default magins for the weymouth polynomial approximation
if nargin < 4
   p_approx = 0.005*ones(numel(K),1); 
end
%% inizialization
f = zeros(numel(K),1);
df = zeros(numel(K),1);
d2f = zeros(numel(K),1);
%% finding bounds for the approximation
p = pi-pj;
id_l = p < -p_approx;
id_u = p > p_approx;
id_lu = p >= -p_approx & p <= p_approx;
%% calcualation outside the bounds
f(id_l) = -sqrt(-K(id_l).*p(id_l));
df(id_l) = K(id_l)./(2*sqrt(-K(id_l).*p(id_l)));
d2f(id_l) = (K(id_l).^2)./(4*sqrt((-K(id_l).*p(id_l)).^3));

f(id_u) = sqrt(K(id_u).*p(id_u));
df(id_u) = K(id_u)./(2*sqrt(K(id_u).*p(id_u)));
d2f(id_u) = -(K(id_u).^2)./(4*sqrt((K(id_u).*p(id_u)).^3));

%% calculation of the approximation inside the bounds
if not(isempty(id_lu))
    k1 = sqrt(K(id_lu).*p_approx(id_lu)); 
    k2 = K(id_lu)./(2*sqrt(K(id_lu).*p_approx(id_lu)));
    k3 = -((K(id_lu).^2)./(4*(sqrt((K(id_lu).*p_approx(id_lu)).^3))));
    c = (k1 - p_approx(id_lu).*k2 + (k3/3).*(p_approx(id_lu).^2))./((8/3)*(p_approx(id_lu).^5));
    b = (k3 - (20*c).*(p_approx(id_lu).^3))./(6*p_approx(id_lu));
    a = k2 - 0.5*p_approx(id_lu).*(k3 - (20*c).*(p_approx(id_lu).^3)) - (5*c).*(p_approx(id_lu).^4);
    
    f(id_lu) = a.*p(id_lu) + b.*(p(id_lu).^3) + c.*(p(id_lu).^5);
    df(id_lu) = a + (3*b).*(p(id_lu).^2) + (5*c).*(p(id_lu).^4);
    d2f(id_lu) = (6*b).*p(id_lu) + (20*c).*(p(id_lu).^3);
end