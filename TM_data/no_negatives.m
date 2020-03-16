function [Aex] = no_negatives(TM,volb)
% Correct transport matrix so all elements are non-negative
% 
%  Syntax:
%    [Aex] = no_negatives(TM,volb)
%
%  TM is uncorrected transport matrix
%  volb is vector of grid box volumes

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

    I = speye(size(TM));                % create identity matrix [nb x nb]
    dindx = find(I);                    % get index of diagonals [nb x 1]
    
    % convert from concentration flux to mass flux (* sink volume) [nb x nb]
    [i,j,v]=find(TM);
    v=v.*volb(i);
    Aex=sparse(i,j,v,size(TM,1),size(TM,2));
    
    Aneg        = Aex.*(Aex<0);         % identify negatives [nb x nb]
    Aneg(dindx) = 0;                    % exclude diagonals  [nb x nb]
    
    Aex = Aex - Aneg + (-Aneg');        % subtract negative off-diagonals 
                                        % then add transpose of negative off-diagonals, changing sign [nb x nb]
    colneg      = sum(Aneg,1)';         % sum of negatives in columns (downstream negative mass) [nb x 1]
    rowneg      = sum(Aneg,2);          % sum of negatives in rows    (  upstream negative mass) [nb x 1]
    
    Aex(dindx)  = Aex(dindx) + colneg;  % add sum of column negatives to the diagonal [nb x nb]
    Aex(dindx)  = Aex(dindx) + rowneg;  % add sum of row negatives to the diagonal    [nb x nb]
    
    % convert back to concentration flux (/ snk volume) [nb x nb]    % convert from concentration flux to mass flux (* sink volume) [nb x nb]
    [i,j,v]=find(Aex);
    v=v./volb(i);
    Aex=sparse(i,j,v,size(TM,1),size(TM,2));
    
end


