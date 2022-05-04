classdef Tools_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> Set directories.
        % >> 1.1. ---------------------------------------------------------
        function [] = Set_Directories_1D()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('[Post-processing]'));
            addpath(genpath('../[Tools]'));
        end
        % >> 1.2. ---------------------------------------------------------
        
        %% > 2. -----------------------------------------------------------
        % >> Sort structures.
        % >> 2.1. ---------------------------------------------------------
        %  > Sort "msh" (2D) fields.
        function [msh] = Sort_msh_2D(msh)
            % >> msh.
            msh        = orderfields(msh       ,{'c','d','f','v','struct'});
            %  > c.
            msh.c      = orderfields(msh.c     ,{'c','f','h','logical','Nc','Volume'});
            msh.c.c    = orderfields(msh.c.c   ,{'nb','xy'});
            msh.c.c.nb = orderfields(msh.c.c.nb,{'f','v'});
            msh.c.c.xy = orderfields(msh.c.c.xy,{'c','v'});
            msh.c.f    = orderfields(msh.c.f   ,{'if','xy','Sf'});
            msh.c.f.xy = orderfields(msh.c.f.xy,{'c','v'});
            msh.c.h    = orderfields(msh.c.h   ,{'h','xy'});
            %  > d.
            msh.d      = orderfields(msh.d     ,{'h'});
            %  > f.
            msh.f      = orderfields(msh.f     ,{'ic','iv','logical','Nf','xy'});
            msh.f.xy   = orderfields(msh.f.xy  ,{'c','v'});
            %  > v.
            msh.v      = orderfields(msh.v     ,{'ic','if','logical'});
        end
        
        %% > 3. -----------------------------------------------------------
        % >> Other functions.
        % >> 3.1. ---------------------------------------------------------
        %  > Similar (faster) version of other built-in functions (to compute distance between 2 points).
        function [D] = dist(XY)
            D = sqrt((XY(1,1)-XY(2,1)).^2+(XY(1,2)-XY(2,2)).^2);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "mean".
        function [y] = mean(x,j)
            y = sum(x,j,'default','includenan')./size(x,j);
        end
    end
end