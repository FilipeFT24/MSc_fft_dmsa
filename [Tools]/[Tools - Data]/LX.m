classdef LX
    methods(Static)
        % >> Cell/face L1,L2 and L_infinity norms.
        function [L] = n(E,V)
            if nargin == 1
                L(1,:) = LX.L1(E);
                L(2,:) = LX.L2(E);
            else
                L(1,:) = LX.L1(E,V);
                L(2,:) = LX.L2(E,V);
            end
            L(3,:) = LX.L3(E);
        end
        % >> Auxiliary functions.
        %  > L1.
        function [L1] = L1(E,V)
            if nargin == 1
                L1 = mean(E);
            else
                L1 = sum (E.*V)./sum(V);
            end
        end
        %  > L2.
        function [L2] = L2(E,V)
            if nargin == 1
                L2 = mean(sqrt(E.^2));
            else
                L2 = sum (sqrt((E.*V).^2))./sum(sqrt(V.^2));
            end
        end
        %  > L_infinity.
        function [L3] = L3(E)
            L3 = max(E);
        end
    end
end