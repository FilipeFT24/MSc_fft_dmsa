classdef obj
    %% > Simulation variables.
    %  > Non-overridable variables.
    properties(Constant)
        %% > #1.
        % >> Mesh.
        Xv_i            =  0;
        Xv_f            =  1;
        %  > Mesh type.
        msh_Type        = 'Non-uniform';
        %  > if 'msh_Type'='Non-uniform'...
        Stretching_Type = 'Wall';
        Nf_0            =  0;
        Nf_Loc          =  0.5;
        K_SF            =  1.01;
        Wall_Select     = 'East';
        %% > #2.
        % >> Boundary conditions.
        West_Boundary   = 'Dirichlet';
        East_Boundary   = 'Dirichlet';  
        %% > #3.
        % >> Integration method.
        n               =  4;         
        f_Select        = 'Manufactured';
        %% > #4.
        % >> Post-processing.
        SaveData        = 'F';
        LoadData        = 'F';
        PP              = 'T';        
    end
    %% > Test variables(1).
    %  > Non-overridable variables.
    properties(Constant)
        %% > #1.
        % >> Flow conditions.
        V               =  0;
        Gamma           =  0.1; 
        %% > #2.
        % >> Boundary values/manufactured function.
        %  > if 'f_Select'='Analytic'...
        Phi_W           =  0;
        Phi_E           =  1;
        gradPhi_E       =  1;
        gradPhi_W       =  0;      
    end
    %% > Test variables(2).
    %  > Non-overridable variables.
    properties(Constant)
        %% > #1.
        % >> Flux reconstruction method.
        Sim_Type        = 'Implicit';
        %  > if DC...
        n_LO            =  4;
        n_HO            =  6;   
        %% > #2.
        % >> Run type.
        i_Test          = 'Standard';
        %  > if 'i_Test'='Error decay'...
        Numb_Loops      =  10;
        Add_To          =  5;
        %% > #3.
        % >> Refinement routine parameters.
        i_ref           =  'F';
        K_ref           =  1.5;
        Numb_Cycles     =  1;
    end
    %  > Overridable variables.
    properties
        %% > #1.
        % >> Mesh.
        NC              =  15;
    end
    %% > Methods...
    methods
        function [obj] = Set_NC(obj)
            %% > Error decay.
            if strcmpi(obj.i_Test,'Error decay')
                for i = 1:obj.Numb_Loops
                    if i == 1
                        obj.NC(i) = obj.NC;
                    else
                        obj.NC(i) = obj.NC(i-1)+obj.Add_To;
                    end
                end
            end
        end
    end
end