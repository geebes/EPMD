function [ diag_fcns ] = diagnostics ( )
    
    diag_fcns.last_year      = @last_year;
    diag_fcns.load_meta      = @load_meta;
    diag_fcns.geo_refs       = @geo_refs;
    diag_fcns.load_variable  = @load_variable;
    diag_fcns.reshape_bio    = @reshape_bio;
    diag_fcns.reshape_forcing= @reshape_forcing;
    diag_fcns.unpack_biomass = @unpack_biomass;
    diag_fcns.make_tree      = @make_tree;
    diag_fcns.draw_tree      = @draw_tree;
    diag_fcns.tree_meta      = @tree_meta;
    diag_fcns.prune_tree     = @prune_tree;
    diag_fcns.cmap_2d        = @cmap_2d;
end



%%
function [ eco_pars, gen_pars, bgc_pars, ocn_pars, Ib ] = load_meta(matObj)
% load model metadata and georeferencing data

    % structural metadata
    eco_pars=matObj.eco_pars;
    gen_pars=matObj.gen_pars;
    bgc_pars=matObj.bgc_pars;
    ocn_pars=matObj.ocn_pars;

    Ib=ocn_pars.Ib; % index of surface layer grid

end

%%
function [ vi, vj, vk, ni, nj, nk, xe, ye, xm, ym, u, v ] = geo_refs(ocn_pars)

    Ib=ocn_pars.Ib;
    vi=ocn_pars.i(Ib);
    vj=ocn_pars.j(Ib);
    vk=ocn_pars.k(Ib);
    
    ni=numel(unique(vi));
    nj=numel(unique(vj));
    nk=numel(unique(vk));

    % get edge coordinates of surface grid
    Xe=ocn_pars.lon_edges([unique(vi);max(vi)+1]);
    Ye=ocn_pars.lat_edges([unique(vj);max(vj)+1]);    
    % mesh grid of edges
    [xe ye]=meshgrid(Xe,Ye);
    
    % midpoint locations at surface
    xm = ocn_pars.lon(vi);
    ym = ocn_pars.lat(vj);
    
    % calculate water flow vectors (defined relative to lat/lon coord system)
    % note this might not be perfectly accurate, as does not use great circle path
    B=ocn_pars.TMB{1};
    u=B'*xm;
    v=B'*ym;

end

%%
function [ data ] = load_variable(matObj,variable,nyr)
% extract output data from matObj file 
% alsop converst from cell to matrix form
% 'variable' is one of:
%     {'AGE'     }
%     {'DOP'     }
%     {'GENOME'  }
%     {'PHY'     }
%     {'PO4'     }

    data = cell2mat(matObj.(variable)(nyr,1));

end

%%
function [ matrix, index ] = reshape_bio(data,vi,vj,eco_pars,ocn_pars)

    ni=numel(ocn_pars.lon);
    nj=numel(ocn_pars.lat);
    nt=eco_pars.ntroph;
    nv=eco_pars.nsize;   
    
    [~,~,tr]=unique(-eco_pars.trophic);
    [~,~,sz]=unique(eco_pars.V);
    
    ii=repmat(vi,1,eco_pars.jpmax);
    jj=repmat(vj,1,eco_pars.jpmax);
    kk=repmat(tr',ocn_pars.ni,1);
    ll=repmat(sz',ocn_pars.ni,1);
    
    % reshape matrix into 4 dimensions (lon, lat, trait 1, trait 2)
    
    
    mindx  = sub2ind([nj ni nt nv] , jj(:) , ii(:) , kk(:) , ll(:) );
    matrix = zeros([ni nj nt nv]).*NaN;
    matrix(mindx) = data;
    
    % create standard index for data in reshaped array
    index=zeros(size(matrix(:,:,:,:)));
    index(1:end)=1:numel(index);
end

%%
function [ matrix, index ] = reshape_forcing(data,vi,vj,eco_pars,ocn_pars)

    ni=numel(ocn_pars.lon);
    nj=numel(ocn_pars.lat);
    
    % reshape matrix in2 dimensions (lon, lat)
    mindx  = sub2ind([nj ni] , vj(:) , vi(:) );
    matrix = zeros([ni nj]).*NaN;
    matrix(mindx) = data;
    
    % create standard index for data in reshaped array
    index=zeros(size(matrix));
    index(1:end)=1:numel(index);
end
%%
function [tree D t_index] = make_tree(matObj,iyear,threshold)

    biomass = load_variable(matObj,'PHY',iyear);
    t_index = find(biomass>threshold);    % Extract indices of non-negligible biomass
    
    genome  = load_variable(matObj,'GENOME',iyear);
    ngenes  = size(genome,2);
    
    nbit=53; % log2(FlIntMax) = 53
    str='';
    for i=1:ngenes
        bin_gene = dec2bin(genome(t_index,i),nbit); % convert double value to binary string
        str=[str bin_gene];
    end
    % convert to numeric values (1/0)
    binstr=str-'0';
    
    % get Hamming distance as weighted number of base differences
    % D = size(binstr,2) * pdist(binstr,'hamming');
    ww=nthroot(2,10).^(0:ngenes-1);
    w=reshape(repmat(ww,nbit,1),1,[]); % weight genes to have different mutation rates
    w=w./mean(w);
    binstr=binstr./w; % multiply bits by inverse mutation weights (less likely worth more)
    D = pdist(binstr,'CityBlock')./size(binstr,2);
    
    tree = linkage(D,'average');
end
%%


function [ taxID, cmap] = draw_tree(tree,t_index,eco_pars,ocn_pars,nsp)

    % get traits corresponding to tree nodes
	tree_data = tree_meta(t_index,eco_pars,ocn_pars);
    
    % find distance cutoff for specified number of clusters
    cutoff = median([tree(end-(nsp-1),3) tree(end-(nsp-2),3)]);

    % rank nodes by trophic strategy and size
    [~,rankorder] = sortrows([tree_data.trophic tree_data.size],'descend'); 
    % optimally sort terminal leaf nodes by phenotype
    [leafOrder]   = sort_branches(tree,rankorder);
    % draw tree
    [H, nodes, outperm] = dendrogram(tree,0,'ColorThreshold',cutoff,'Orientation','left','ReOrder',leafOrder);
    % reverse Y direction
    ax1=gca;
    set(ax1,'YTick',[],'YDir','reverse');
    xl=xlim;
    
    % get nsp clusters
    T = cluster(tree,'maxclust',nsp);
    taxID=T;
    % get unique xcoord-node-cluster triplets
    ync=unique([nodes T],'rows');
    ync=[(1:numel(t_index))' ync(outperm,:)];
    
    % colour by genetic cluster 
    cmap=lhsdesign(nsp,3); % generate LHC colormap
    for i=1:numel(H)
        if ~ismember(H(i).Color,[0 0 0],'rows') % if node is not black (i.e. above defined taxonomic cutoff)
            % find an integer node within range of branch i
            ii=round(median(H(i).YData));
            % find cluster/s corresponding to branch
            ci=T(outperm(ii));
            H(i).Color=cmap(unique(ci),:); % change colour to row of cmap
        end
        H(i).LineWidth=0.75;
    end
    xlim([0 xl(2)]);
    xl=xlim;

    % grey out all branches with no phenotypic diversity among descendents
    [H,spind] = prune_tree(tree,H,tree_data);
    
    % get 2d trait colorspace
    [ rgb, rgb3 ] = cmap_2d(eco_pars);
    
    % create index of traits
    phenmat=repmat(1:eco_pars.nphy,ocn_pars.ni,1);
    % sort used data by dendrogram permutation
    iphen=phenmat(t_index(outperm));
    
    % create image based on phenotypes and rgb colour map
    phenotype=ones(size(rgb,1),2,size(rgb,2));
    phenotype(:,1,:)=rgb;
    phenotype=phenotype(iphen,:,:);
    % create image based on taxa and lhc colour map
    taxon    =ones(size(phenotype));
    taxon(:,1,:)=cmap(T(outperm),:);
    
    % specify locations of color bar images
    wd=xl(2)./30;
    xlim([-2*wd xl(2)]);    
    % plot colorbars with boxes around them
    hold on
    image(taxon,     'XData', -wd.*[0.5 1.5], 'YData', [1 numel(outperm)])
    image(phenotype, 'XData', -wd.*[1.5 2.5], 'YData', [1 numel(outperm)])
    rectangle('Position',[-wd 0 wd numel(outperm)],'LineWidth',0.75)
    rectangle('Position',[-2.*wd 0 wd numel(outperm)],'LineWidth',0.75)
    
    set(ax1,'Color','none')
    set(ax1,'LineWidth',0.75)
    set(gcf,'Color','w')
    
    set(gca,'FontSize',14)
    xlabel('Genetic distance','FontSize',14)
    
    % draw colour square
    ax2 = axes('Position',[0.12 0.725 .2 .2]);
    image(linspace(0,1,eco_pars.ntroph),linspace(0,1,eco_pars.nsize),rgb3)
    set(ax2,'XTick',[0 0.5 1],'YTick',[0 0.5 1],...
            'XTickLabel',{'auto','mixo','hetero'},...
            'YTickLabel',{'small','medium','large'},...
            'FontSize',12,...
            'Layer','bottom');
    axis square
    axis xy
    hold on
    ax=axis;
    rectangle('Position',[ax([1 3]) diff(ax(1:2)) diff(ax(3:4))],'LineWidth',0.75);
    
    uistack(ax2,'bottom')
    
    
end

%%


function [ tree_data ] = tree_meta(t_index,eco_pars,ocn_pars)

    nocn=ocn_pars.ni;
    vi=ocn_pars.i(ocn_pars.Ib);
    vj=ocn_pars.j(ocn_pars.Ib);
    
    % define EiE matrix grid coordinates
    ESD     = (6.*eco_pars.V./pi)'.^(1/3);
    ESDgrid = repmat(ESD,nocn,1);
    TROgrid = 1-repmat(eco_pars.trophic',nocn,1);
    LATgrid=repmat(vj,1,eco_pars.ntroph*eco_pars.nsize);
    LONgrid=repmat(vi,1,eco_pars.ntroph*eco_pars.nsize);
    
    % get trait and geographic coordinates of useable data
    tree_data.size    = ESDgrid(t_index);
    tree_data.trophic = TROgrid(t_index);
    tree_data.lat     = LATgrid(t_index);
    tree_data.lon     = LONgrid(t_index);
end

%%
function [H,spstr] = prune_tree(tree,H,tree_data)

    % make strings describing phenotypes
    [~,~,sz_ind]=unique(tree_data.size);
    [~,~,tr_ind]=unique(tree_data.trophic);
    spstr=cellstr(strcat(num2str(sz_ind,'%03i'),num2str(tr_ind,'%03i')));
    
    m = size(spstr,1);
    for i=1:m-1 % loop through nodes
        % get index of children
        n1=tree(i,1);
        n2=tree(i,2);

        % if children have identical phenoytpe
        if strcmp(spstr{n1},spstr{n2})
            % set parent phenotype to that of two (identical) children
            spstr{m+i}=spstr{n1};
            % set alpha to low value
            clr=1-H(i).Color(1:3);
            clr=1-clr.*0.25;
            H(i).Color=clr;
            H(i).LineWidth=0.25;
        else
            % set phenotype to Nan (so always different) 
            spstr{m+i}=NaN;
            H(i).Color=H(i).Color(1:3);
        end
    end
end

%%
function [ rgb, rgb3 ] = cmap_2d(eco_pars)


    ns=eco_pars.nsize;
    nt=eco_pars.ntroph;
    
    % define colormap
    x=[0 0 0.5 1 1]'; % x cordinates of corners and centre
    y=[0 1 0.5 0 1]'; % y cordinates of corners and centre
    r=[0 1 1 0 1]';   % corresponding r values
    g=[1 1 1 0 0]';   % corresponding g values
    b=[0 0 1 1 1]';   % corresponding b values
    
    xl=linspace(0,1,nt); % unique values of trait x
    yl=linspace(0,1,ns); % unique values of trait y
    
    [xx yy]=meshgrid(xl,yl); % generate trait grid
    
    Fr=scatteredInterpolant(x,y,r); % r values across grid
    Fg=scatteredInterpolant(x,y,g); % g values across grid
    Fb=scatteredInterpolant(x,y,b); % b values across grid
    
    rgb=[Fr(xx(:),yy(:)) Fg(xx(:),yy(:)) Fb(xx(:),yy(:))]; % rgb array

    rgb3=cat(3,Fr(xx,yy),Fg(xx,yy),Fb(xx,yy));
end

%%


























