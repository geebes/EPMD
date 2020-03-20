function [Aex] = conserve_mass(TM,volb)
    % The original transport matrix is given as a concentration flux. This
    % means that each element transfroms an upstream concentration into a
    % downstream concentration, implicitly accounting for changes in grid cell
    % size. We first convert the concentration flux matrix to a mass flux
    % matrix, with each element describing the actual water flux from one grid 
    % cell to another (irrespective of grid size).   
    
    disp('Correcting mass conservation')
    % convert from concentration flux to mass flux (* sink volume) [nb x nb]
    [i,j,v]=find(TM);
    v=v.*volb(i);
    tm=sparse(i,j,v,size(TM,1),size(TM,2));
    
    % This matrix *should* (ideally) have the following properties:
    %   columns should sum to the local water mass
    %       (because each column describes the sink (downstream) distribution for each location)
    %   rows should also sum to the local water mass
    %       (because each row describes the source (upstream) distribution for each location)
 
    % These properties are not exactly met, for various reasons, so we apply a mass correction.
    
    % First we get an index of all elements on the diagonal    
    I = speye(size(tm));                % create identity matrix [nb x nb]
    adindx      = find(I);              % get index of diagonals [nb x 1]
    
    % We correct them by subtracting the column sum of the off diagonals from the local volume
    tm(adindx)  = volb - (sum(tm,1)' - tm(adindx));
    % This ensures that the columns of the mass transport matrix sum to exactly
    % the local grid cell volume, such that mass is conserved during transport
    
    % Finally, we convert transport matrix back to concentration fluxes...
    [i,j,v]=find(tm);
    v=v./volb(i);
    Aex=sparse(i,j,v,size(tm,1),size(tm,2));
    
    
    
end


