function [ quadrature ] = GLC_quadrature( xmin, xmax, wanted_level )
% Generates an array of point corresponding to the Gauss-Lobatto-Tchebishev
% quadrature comprised between user parameters xmin and xmax. The wanted
% level for the quadrature is set by the wanted_level argument
N = wanted_level;
 
Points = 2; % points extremes de lintervale
 
x = [xmin,xmax];
 
 
for n=1:N
    
    Points = Points + 2^(n-1);
    
    x = [];
    
    for i=2:2:Points-1
        
        aux = -cos(pi*(i-1)/(Points-1));
        
        x = [x xmin + (xmax-xmin)*(aux+1)/2];
        
    end    
    
    quadrature = x;
end

end

