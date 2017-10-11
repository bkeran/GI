(* ::Package:: *)

clear;clc;close all;

a=AllComponents;

a.E=45.7;
disp(['Energy of the particles is ',num2str(a.E),' keV.']);
a.N=[800,324];
disp(['There is [',num2str(a.N),'] particles in the beam.']);
a.p=5e-6;
a.dc=0.5;
a.SourceToG1Distance=10000;
a.grid_spacing=[5e-6/200,5e-6/200];

a.source_position=[-a.N(2)/2*a.grid_spacing(2),-a.N(1)/2*a.grid_spacing(1),0];

Source=a.PlaneWaveSource();

% G1=a.Create1DGrating({'Si',2.336,500e-6}, {'Ni',8.908,8e-6});
G1 = a.Rectangular1DGrating();  % mora biti horizontalan!
a.phShift=pi/2;
% absorb = 1-value + (value == 1)*eps;
absorb = 1;
a.absorb_exp = log (absorb)/2;
G1=exp (a.absorb_exp.*G1).*exp(i*a.phShift.*G1);
disp(['Grating G1 has the phase shift ',num2str(a.phShift),' and absorption ',num2str(exp(a.absorb_exp)),'.']);
disp(['First fractional Talbot distance is ',num2str(a.talbot),' m.']);

psi=a.SourceAfterG1(Source, 0, G1);

D=zeros(800,324);

z = linspace(0.00001,2.2,800);

for Dz=1:length(z)
    a.G1ToG2Distance=z(Dz);
    psi_out=a.G1ToG2(psi);
    I=mean(abs (psi_out).^2,2);
    D(:,Dz)=I;
end
hold on;
imagesc(z,[],D);
line([a.talbot,a.talbot],[0,523],'LineWidth',1,'Color','black')
line([2*a.talbot,2*a.talbot],[0,523],'LineWidth',1,'Color','red')
line([3*a.talbot,3*a.talbot],[0,523],'LineWidth',1,'Color','black')
line([4*a.talbot,4*a.talbot],[0,523],'LineWidth',1,'Color','red')
xlabel('propagation distance / m')
ylabel('intensity') % / Wm\:207b\.b2')
xlim([0.00001 2.2]);
ylim([0 400]);






