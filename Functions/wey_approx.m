%% ========================================================================
%% Aproximación de la ecuación de Weymouth y sus derivadas de primer y
% segundo orden para un conjunto de gasoductos dados en función de la 
% variable p=pi-pj, empleando la aproximación polinomial propuesta
% ========================================================================
% Autor: Wilson González Vanegas M.Sc(c)
% Grupo de investigación Automática, Universidad Tencnológica de Pereira
% GIPEM, Universidad Nacional de Colombia - sede Manizales
% ========================================================================

function [f,df,d2f] = wey_approx(K,pi,pj,p_approx)
if nargin < 4
   p_approx = 0.005*ones(numel(K),1); 
end

f = zeros(numel(K),1);
df = zeros(numel(K),1);
d2f = zeros(numel(K),1);

p = pi-pj;
id_l = p < -p_approx;
id_u = p > p_approx;
id_lu = p >= -p_approx & p <= p_approx;

f(id_l) = -sqrt(-K(id_l).*p(id_l));
df(id_l) = K(id_l)./(2*sqrt(-K(id_l).*p(id_l)));
d2f(id_l) = (K(id_l).^2)./(4*sqrt((-K(id_l).*p(id_l)).^3));

f(id_u) = sqrt(K(id_u).*p(id_u));
df(id_u) = K(id_u)./(2*sqrt(K(id_u).*p(id_u)));
d2f(id_u) = -(K(id_u).^2)./(4*sqrt((K(id_u).*p(id_u)).^3));

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