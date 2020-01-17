% Purepose:
%   Generate neighbor list
function [nbr] = InitNlist(x, boxlx, A)
    na = length(x(1,:));
    % Separation between two nearest neighbors at equilibrium
    r0= sqrt(3.0)*A/4.0;
    
    % Initialization 
    nbr = -1; 
    
    for i0 = 1: na; 
        counter = 0;
        for i1 = 1: na;
            if(i1 ~= i0) % cannot be self
                dx = x(:,i0) - x(:,i1);
                dx = dx-round(dx./boxlx).*boxlx;
                dr = sqrt(dot(dx,dx));
                
                if(dr < r0*1.01)
                    counter = counter + 1;
                    nbr(i0, counter) = i1;
                else
                    continue;
                end
            else
                continue;
            end
        end
    end
end