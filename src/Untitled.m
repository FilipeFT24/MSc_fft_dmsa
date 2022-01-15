%         % >> 1.4. ---------------------------------------------------------
%         function [add_to] = StencilExt_1(msh)
%             % >> Stencil indices.
%             for i = 1:size(msh.s.st,2)
%                 for j = 1:size(msh.s.st,1)
%                     st_i{i}{j} = msh.s.st{j,i};
%                 end
%                 st_i{i} = cell2mat(st_i{i});
%             end
%             % >> Add elements to stencil...
%             %  > Initialize.
%             add_to = cell(1,msh.f.NF);
%             for i = 1:msh.f.NF
%                 k = 1;
%                 for j = 1:msh.c.NC
%                     if (msh.c.mean(1,j) >= msh.s.lim.x_min(i) && msh.c.mean(1,j) <= msh.s.lim.x_max(i)) && ...
%                             (msh.c.mean(2,j) >= msh.s.lim.y_min(i) && msh.c.mean(2,j) <= msh.s.lim.y_max(i)) && ~ismembc(j,sort(st_i{i}))
%                         add_to{i}(k) = j;
%                         k = k+1;
%                     end
%                 end
%             end
%         end