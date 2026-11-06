function fsdt_plot(r)
vis=r.cfg.output.figureVisible;
figure('Name','FSDT blast histories','Visible',vis,'Color','w','Position',[100 100 1200 780]);
subplot(2,2,1); plot(r.time*1000,r.pressurePa/1000,'LineWidth',1.3);
xlabel('Time [ms]'); ylabel('Gauge pressure [kPa]'); grid on;
subplot(2,2,2); plot(r.time*1000,r.normalizedProbeDisplacement,'LineWidth',1.3);
xlabel('Time [ms]'); ylabel('Probe w/h'); grid on;
subplot(2,2,3); hold on;
nt=numel(r.time); np=numel(r.cfg.plies);
for k=1:np
    plot(r.time*1000,reshape(r.stress.localPa(1,2,k,:),nt,1)/1e6,'LineWidth',1);
end
xlabel('Time [ms]'); ylabel('Ply top-face sigma_1 at probe [MPa]'); grid on;
labels=cell(1,np); for k=1:np, labels{k}=sprintf('Ply %d',k); end
legend(labels,'Location','best');
subplot(2,2,4); plot(r.time*1000,r.energy(:,1:4),'LineWidth',1);
xlabel('Time [ms]'); ylabel('Energy [J]'); grid on;
legend('Kinetic','Strain','External work','Dissipated','Location','best');
if r.cfg.output.writeFiles
    print(gcf,fullfile(r.cfg.output.folder,'histories.png'),'-dpng','-r150');
end
figure('Name','FSDT displacement field','Visible',vis,'Color','w','Position',[100 100 900 680]);
n=r.mesh.nn;
X=reshape(r.mesh.nodes(:,1),r.mesh.nx+1,r.mesh.ny+1);
Y=reshape(r.mesh.nodes(:,2),r.mesh.nx+1,r.mesh.ny+1);
W=reshape(r.peakDisplacement(1:n),r.mesh.nx+1,r.mesh.ny+1);
surf(X,Y,W,W*1000,'FaceColor','interp','EdgeColor',[.6 .6 .6]);
view(35,30); axis tight; colorbar; xlabel('x [m]'); ylabel('y [m]'); zlabel('w [m]');
title(sprintf('Field at peak probe time %.3f ms; colour: w [mm]',r.peakTimeS*1000));
if r.cfg.output.writeFiles
    print(gcf,fullfile(r.cfg.output.folder,'peak_field.png'),'-dpng','-r150');
end
end
