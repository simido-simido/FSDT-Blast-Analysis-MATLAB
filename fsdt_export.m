function fsdt_export(result)
%FSDT_EXPORT SI CSV histories, modes, sampled ply stresses, and MAT results.
folder=result.cfg.output.folder;
if ~exist(folder,'dir'), mkdir(folder); end
fsdt_write_csv(fullfile(folder,'response.csv'), ...
    {'time_s','gauge_pressure_Pa','probe_w_m','probe_w_over_h', ...
    'probe_velocity_m_per_s','probe_acceleration_m_per_s2','center_w_m'}, ...
    [result.time result.pressurePa result.probeDisplacementM ...
    result.normalizedProbeDisplacement result.probeVelocityMps ...
    result.probeAccelerationMps2 result.centerDisplacementM]);
fsdt_write_csv(fullfile(folder,'energy.csv'), ...
    {'time_s','kinetic_J','strain_J','external_work_J','dissipated_J', ...
    'energy_balance_error_J','relative_equilibrium_residual'}, ...
    [result.time result.energy result.equilibriumResidual]);
m=result.modal;
fsdt_write_csv(fullfile(folder,'frequencies.csv'), ...
    {'mode','omega_rad_per_s','frequency_Hz','transverse_mass_fraction','eigen_residual'}, ...
    [(1:numel(m.omega))' m.omega m.frequencyHz m.transverseMassFraction m.residual]);
np=numel(result.cfg.plies); nt=numel(result.time);
headers={'time_s'}; data=result.time;
for k=1:np
    for face=1:2
        labels={'bottom','top'};
        prefix=sprintf('ply%d_%s_',k,labels{face});
        comp={'sigma_x_Pa','sigma_y_Pa','tau_xy_Pa','sigma_1_Pa','sigma_2_Pa','tau_12_Pa'};
        for j=1:6, headers{end+1}=[prefix comp{j}]; end
        sg=reshape(result.stress.globalPa(:,face,k,:),3,nt)';
        sl=reshape(result.stress.localPa(:,face,k,:),3,nt)';
        data=[data sg sl]; %#ok<AGROW>
    end
end
fsdt_write_csv(fullfile(folder,'probe_ply_stress.csv'),headers,data);
n=result.mesh.nn;
fsdt_write_csv(fullfile(folder,'mesh_nodes.csv'),{'node','x_m','y_m'}, ...
    [(1:n)' result.mesh.nodes]);
fsdt_write_csv(fullfile(folder,'mesh_elements.csv'),{'element','n1','n2','n3','n4'}, ...
    [(1:result.mesh.ne)' result.mesh.elements]);
fsdt_write_csv(fullfile(folder,'peak_probe_time_field.csv'), ...
    {'node','x_m','y_m','w_m','theta_x_rad','theta_y_rad','u_m','v_m'}, ...
    [(1:n)' result.mesh.nodes reshape(result.peakDisplacement,n,5)]);
save(fullfile(folder,'result.mat'),'result','-v7');
fid=fopen(fullfile(folder,'summary.txt'),'w');
if fid<0, error('Cannot create summary.'); end
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'Ferreira-derived linear FSDT blast FEM\nCase: %s\n',result.cfg.name);
fprintf(fid,'Units: m, kg, s, Pa. Ply order: bottom to top.\n');
fprintf(fid,'Element: Q4; 2x2 membrane/bending; 1x1 shear; 2x2 mass.\n');
fprintf(fid,'Support: %s; membrane: %s\n',result.cfg.boundary.transverse,result.cfg.boundary.membrane);
fprintf(fid,'Nodes=%d; elements=%d; dt=%.12g s\n',n,result.mesh.ne,result.actualDt);
fprintf(fid,'Probe=(%.9g, %.9g) m; stress is one-sided in element %d\n', ...
    result.probe.point(1),result.probe.point(2),result.probe.element);
fprintf(fid,'Peak absolute probe displacement=%.12g m; t=%.12g s\n',result.peakProbeAbsM,result.peakTimeS);
fprintf(fid,'Maximum nodal |w|/h over full history=%.12g\n',result.maxAbsOverH);
fprintf(fid,'Maximum relative energy-balance error=%.12g\n',result.energyRelativeError);
fprintf(fid,'Maximum equilibrium residual=%.12g\n',max(result.equilibriumResidual));
fprintf(fid,'First computed frequency=%.12g Hz\n',m.frequencyHz(1));
fprintf(fid,'Static probe displacement at +pmax pressure=%.12g m\n',result.staticAtPmaxM);
fprintf(fid,'No nonlinear membrane stiffening, damage, delamination, or blast-air solution.\n');
fprintf(fid,'Demonstration pulse inputs are prescribed, not experimentally calibrated.\n');
end
