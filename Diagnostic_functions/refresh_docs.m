cd ~/GitHub


!rm -rf /Users/baw103/GitHub/EPMD/doc


PATH = getenv('PATH');
if ~contains(PATH,'/usr/local/Cellar/graphviz/2.42.2/bin');
    setenv('PATH', [PATH ':/usr/local/Cellar/graphviz/2.42.2/bin']);
end

m2html('mfiles','EPMD',...
       'htmldir','EPMD/doc',...
       'ignoredDir',{'Diagnostic_functions','Richter_2020','SphereSample'},...
       'recursive','on',...
       'global','on',...
       'globalHypertextLinks','on',...
       'template','frame',...
       'index','menu',...
       'graph','on');