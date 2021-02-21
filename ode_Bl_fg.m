function [ df ] = ode_Bl_fg( Pr,eta,f )
df    = zeros(5,1); 
df(1) = f(2);
df(2) = f(3);
df(3) = -f(1)*f(3);
df(4) = f(5);
df(5) = -Pr*f(1)*f(5);
end