%**********************************************************************
% Purpose:
%  Add local two-body dynamic matrix to the global one
% Arguments:
%  list2 -- two-atom list
%  D2 -- local two-body dynamic matrix
%   SUBROUTINE Add2bdDmat(list2, D2)
function [dmatallD] = Add2bdDmat(dmatallD, list2, D2, nbr)
%     INTEGER,DIMENSION(2)	:: list2
%     REAL(MYREAL),DIMENSION(2,3,2,3) :: D2
% 
%     INTEGER			:: ii, jj, i, j, n

    for ii = 1: 1; % 2; % ?? How about just use ii =1 ,
      i = list2(ii);
%       if (not(special(i))) %CYCLE
%           continue;
%       end

      for jj = 1: 2;
        j = list2(jj);

%         for n = 1: 4; % 17;  %?? What's n here? n=1:17? or n=1:4? 
%             if (nbr(i,n)==j) % THEN

                dmatallD(:,:,i,j) = dmatallD(:,:,i,j) - D2(:,:,ii,jj);
%                 if ( j==1 )
%                     continue;
%                 else
                    


%                 break;
%             end % IF
%         end % DO % n
      end % DO % jj

    end % DO % ii

end % SUBROUTINE Add2bdDmat
