function [B] = extract_surf(Af,Ib,Ii,volb)

    vertical='trapped';

    % convert from concentration flux to mass flux (* sink volume) [nb x nb]
    [i,j,v]=find(Af);
    v=v.*volb(i);
    B=sparse(i,j,v,size(Af,1),size(Af,2));

    Bdn = B(Ii,Ib); % downwelling fluxes (Ib to Ii) % Downwelling flux should loop back into surface (i.e. trapped particles)
    Bup = B(Ib,Ii); % upwelling fluxes (Ii to Ib)   % Upwelling flux is drawn from surface box (assumes upwelling matches surface conc.)
    B   = B(Ib,Ib); % fluxes within surface layer

    Fdn = sum(Bdn,1); % sum of downwelling fluxes out of each surface cell (prevented from leaving)
    Fup = sum(Bup,2); % sum of upwelling fluxes into each surface cell     (replaced with surface conc.)

    I        = speye(size(B));
    dindx    = find(eye(size(B))); % get index for diagonal
    
    switch vertical
        case 'trapped'
            B(dindx) = B(dindx) + Fdn'; % Tracer is trapped in surface
        case 'upwelled'
            B(dindx) = B(dindx) + Fup; % Upwelling conc. = surface conc.
    end

    % convert from concentration flux to mass flux (* sink volume) [nb x nb]
    volume=volb(Ib);
    [i,j,v]=find(B);
    v=v./volume(i);
    B=sparse(i,j,v,size(B,1),size(B,2));

end