% Magnitude mode COSY-45 spectrum of rotenone.
% 
% http://dx.doi.org/10.1002/jhet.5570250160
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% gareth.charnock@oerc.ox.ac.uk

% This is an edited version of cosy45_rotenone
% By Tino Wells, Post-Bachelors' RA Deep Learning,
% Machine Learning, & Cheminformatics @ Pacific 
% Northwest National Laboritory. For the original, 
% see original file. 

function cosy45_rotenone_pytest(isotope_list,sys_magnet,chemshifts,zcs)

% Spin system
sys.isotopes=isotope_list;
%Command Window: iso_lst={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}

% Interactions
sys.magnet = sys_magnet;
%Command Window: mag=5.9

inter.zeeman.scalar = chemshifts;
%Command Wondow: shifts={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91 3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72 3.72 3.76 3.76 3.76}

%Command Wondow: zcs =[3 4 12.1; 4 5 3.1; 3 5 1.0; 3 8 1.0; 1 8 1.0; 6 7 8.6; 5 8 4.1; 7 9 0.7; 7 10 0.7; 9 10 15.8; 10 11 9.8; 9 11 8.1; 13 14 1.5; 12 14 0.9; 22 22 0];

% ASSUMES Zeeman coupling scalars are in 3x2xn matrix form
for iter=1:length(zcs)
    inter.coupling.scalar{zcs(iter,1),zcs(iter,2)}=zcs(iter,3); 
end

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