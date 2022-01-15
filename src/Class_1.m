classdef Class_1
    methods(Static)
        %% > Wrap-up Class 1.
        % >> --------------------------------------------------------------
        % >> 1. Wrap-up SubClass_1_1.
        % >> 2. Wrap-up SubClass_1_2.
        % >> 3. Wrap-up SubClass_1_3.
        % >> 4. Order fields...
        % >> --------------------------------------------------------------
        function [inp,msh] = WrapUp_1()
            %  > 1.
            inp = SubClass_1_1.WrapUp_1_1();
            %  > 2.
            msh = SubClass_1_2.WrapUp_1_2(inp);
            %  > 3.
            msh = SubClass_1_3.WrapUp_1_3(inp,msh);
            %  > 4.
            inp = Class_1.Order_Fields(inp);
            msh = Class_1.Order_Fields(msh);
        end
        % >> 4. -----------------------------------------------------------
        function [struct_new] = Order_Fields(struct)
            [~,Order]  = sort(lower(fieldnames(struct)));
            struct_new = orderfields(struct,Order);
        end
    end
end