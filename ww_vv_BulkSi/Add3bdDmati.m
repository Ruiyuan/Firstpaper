%**********************************************************************
% Purpose:
%  Add local three-body dynamic matrix to the global one
% Arguments:
%  list3 -- three-atom list
%  D3 -- local three-body dynamic matrix
%   SUBROUTINE Add3bdDmat(list3, D3)
%     INTEGER,DIMENSION(3)	:: list3
%     REAL(MYREAL),DIMENSION(3,3,3,3) :: D3
% 
%     INTEGER			:: ii, jj, i, j, n
function [dmatallD] = Add3bdDmati(dmatallD, list3, D3, nbr, nt)
%     [aa,bb] = min(abs(list3'-nt));
%     if(bb==1)
%         i = list3(bb);
%     else if(bb=2)
%         end
%     end
    for ii = 1: 1;
      i = list3(ii);
%       if (not(special(i))) % CYCLE
%         continue;
%       end
      for jj = 1:3; %  2: 2;
        j = list3(jj);
%             if(j == nt)
%                 continue;
%             end

%         for n = 1: 4; % 17; %%!!!????????? not sure whether choose n=4 or n=17??
%             if (nbr(i,n)==j) ;%if (dmat(i)%nbr(n).EQ.j) %THEN
%             dmat(i)%D(n,:,:) = dmat(i)%D(n,:,:) + D3(ii,:,jj,:)
                dmatallD(:,:,nt,j) = dmatallD(:,:,nt,j) - D3(:,:,ii,jj);
%                 if(j<=2)
%                 dmatallD(:,:,j,i) = dmatallD(:,:,j,i) - D3(:,:,ii,jj);
%                 end
                %EXIT
%                 break;
%             end %  IF
%         end %  DO % n
      end %  DO % jj

    end %  DO % ii

end %  SUBROUTINE Add3bdDmat