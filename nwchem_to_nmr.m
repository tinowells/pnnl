%
function nwchem_to_nmr(dfile_path)

% tdfread reads all columns as seperate values
% How absolutely useful lol
nwc_mat=tdfread(dfile_path,'\t');

% Flip atom matrix from nx1 char array to 1xn cell array
isos = reshape(cellstr(nwc_mat.atom),[1,length(nwc_mat.atom)]);

% Converting atom list to isotopes list; required for NMR sims
% This is gross, figure out how to one line this. 
for ii=1:length(isos)
    if isos{ii}=='C'
        isos{ii}='13C';
    else
        isos{ii}='1H';
    end
end

% Flip chemical shifts from nx1 to 1xn 
shifts = reshape(nwc_mat.shift,[1,length(nwc_mat.shift)]);
shifts = num2cell(shifts,[1 length(shifts)]);

sys.isotopes=isos;
sys.magnet=5.9;
inter.zeeman.scalar=shifts;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';
bas.sym_group={'S3','S3','S3'};%Methyl group symmetry
bas.sym_spins={[14 15 16],[17 18 19],[20 21 22]};

% Sequence parameters
parameters.angle=pi/4;
parameters.offset=1200;
parameters.sweep=2000;
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@cosy,parameters,'nmr');

% Apodization
fid=apodization(fid,'cosbell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)));

% Destreaking
spectrum=destreak(spectrum);

% Plotting
%tinos_fig = figure(); 
plot_2d(spin_system,abs(spectrum),parameters,20,[0.0025 0.05 0.0025 0.05],2,256,6,'positive');
saveas(tinos_fig,'~/Downloads/testsave.png')


end
