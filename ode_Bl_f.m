function [ df ] = ode_Bl_f( eta,f )
df    = zeros(3,1); 
df(1) = f(2);
df(2) = f(3);
df(3) = -f(1)*f(3);
end