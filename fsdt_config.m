function cfg = fsdt_config()
.
cfg.name = 'cross_ply_demo';
cfg.plate.a = 0.50;                 
cfg.plate.b = 0.40;                 
cfg.mesh.nx = 12;
cfg.mesh.ny = 10;
ply = struct('t',0.001,'angle',0,'E1',135e9,'E2',10e9, ...
    'nu12',0.30,'G12',5e9,'G13',5e9,'G23',3.8e9,'rho',1600);
cfg.plies = repmat(ply,1,4);
angles = [0 90 90 0];
for k=1:numel(angles), cfg.plies(k).angle=angles(k); end
cfg.shearCorrection = 5/6;
cfg.boundary.transverse = 'ssss';    
cfg.boundary.membrane = 'ferreira';  
cfg.pulse.type = 'friedlander_cubic'; 
                                    
cfg.pulse.pmax = 1000;            
cfg.pulse.td = 0.001;               
cfg.pulse.decay = 1.5;              
cfg.pulse.includeNegative = true;
cfg.pulse.pmin = 100;               
cfg.pulse.tm = 0.004;               
cfg.load.direction = -1;            
cfg.load.patch = [];               
cfg.time.dt = 2e-5;             
cfg.time.tEnd = 0.04;             
cfg.damping.ratios = [0 0];         
cfg.damping.modes = [1 3];          
cfg.modal.count = 8;
cfg.output.point = [cfg.plate.a/2 cfg.plate.b/2]; 
cfg.output.folder = fullfile(fileparts(mfilename('fullpath')),'results',cfg.name);
cfg.output.writeFiles = true;
cfg.output.makePlots = true;
cfg.output.figureVisible = 'on';
cfg.output.saveFullHistory = false;
end
