(* ::Package:: *)

clear;clc;close all;

tic

a=AllComponents;

a.E=45.7;
a.N=[1800,4000];
a.p=5e-6;
a.dc=0.5;
a.detector_size=[250e-6,250e-6];
a.detector_pixel=[3,5];
a.sigma=[1e-6,1e-6];
a.SourceToG1Distance=1;

a.SetGridSpacing();
disp(['Grid spacing in x direction is ',num2str(a.grid_spacing(2)),' and in y direction is ',num2str(a.grid_spacing(1)),'.']);

a.source_position=[-a.N(2)/2*a.grid_spacing(2),-a.N(1)/2*a.grid_spacing(1),0];
% a.sample_position=[-a.N(2)/2*a.grid_spacing,-a.N(1)/2*a.grid_spacing,0];
% a.grating_position=[-a.N(2)/2*a.grid_spacing,-a.N(1)/2*a.grid_spacing,0];
a.phase_step=10;

Source=a.PlaneWaveSource();

a.G1ToG2Distance=a.talbot;

% G=a.Rectangular2DGrating();
% G1=a.Create1DGrating({'Si',2.336,500e-6}, {'Ni',8.908,8e-6});
% G2=a.Create1DGrating({'Si',2.336,500e-6}, {'Au',19.32,150e-6});

G1 = a.Rectangular1DGrating();
a.phShift=pi/2;
% absorb = 1-value + (value == 1)*eps;
absorb = 1;
a.absorb_exp = log (absorb)/2;
G1=exp (a.absorb_exp.*G1).*exp(i*a.phShift.*G1);

G2 = a.Rectangular1DGrating();
a.phShift=0;
absorb = eps;
a.absorb_exp = log (absorb)/2;
G2=exp (a.absorb_exp.*G2).*exp(i*a.phShift.*G2);

% a.absorb_exp=log(1/sqrt(2));
% a.phShift=0;
% a.absorb_exp=0;
% a.phShift=pi/2;
Sample=a.TestSample();

psi_sample=a.SourceAfterG1(Source, Sample, G1);
psi_flat=a.SourceAfterG1(Source, 0, G1);

% psi_sample=a.waveFieldProp(a.G1ToG2Distance,psi_sample);
% psi_flat=a.waveFieldProp(a.G1ToG2Distance,psi_flat);

psi_sample=a.G1ToG2(psi_sample);
psi_flat=a.G1ToG2(psi_flat);

I_sample=a.PhaseStepping(psi_sample,G2);
I_flat=a.PhaseStepping(psi_flat,G2);

I_detector _sample=a.Scale2Detector(I_sample);
I_detector _flat=a.Scale2Detector(I_flat);


x_stepping=[];
for jj=1:a.phase_step
    x_stepping=[x_stepping,(jj-1)*a.pixel_step*a.grid_spacing(2)];
end

for ii=1:2
    for jj=1:2
        pixel=[ii,jj];  % piksel detektora koji zelimo gledati 
        v1_flat=squeeze(I_detector _flat(pixel(1),pixel(2),:));     % racuna vrijednosti intenziteta na zeljenom pikselu detektora
        f1_flat=fit(x_stepping.',v1_flat,'fourier1');      % fit na zeljenom pikselu detektora 
        v1_sample=squeeze(I_detector _sample(pixel(1),pixel(2),:));     % racuna vrijednosti intenziteta na zeljenom pikselu detektora
        f1_sample=fit(x_stepping.',v1_sample,'fourier1');
        if ii==1
            subplot(2,2,ii+(jj-1))
        elseif ii==2
            subplot(2,2,ii+jj)%+1)
%         elseif ii==3
%             subplot(3,3,ii+jj+3)
        end
        hold on
        plot(x_stepping,v1_flat.','s');    % crta zeljeni piksel detektora za ff
        plot(f1_flat,'blue');    % crta fit na zeljenom pikselu za ff
        plot(x_stepping,v1_sample.','+');      % crta zeljeni piksel detektora za prop. s uzorkom
        plot(f1_sample,'green');     % crta fit na zeljenom pikselu uzorka za prop. s uzorkom
        xlim([0 5e-6])
        hold off
        legend('flat','flat','sample','sample')
        title(['Value on pixel [', num2str(ii), ',', num2str(jj),']'])
        xlabel('phase steps / m')
        ylabel('intensity') % / Wm\:207b\.b2')
        a0_flat=f1_flat.a0;
        a1_flat=f1_flat.a1;
        b1_flat=f1_flat.b1;
        a0_sample=f1_sample.a0;
        a1_sample=f1_sample.a1;
        b1_sample=f1_sample.b1;
        
        a0_s=a0_sample;
        a0_f=a0_flat;
        a1_s=sqrt(a1_sample^2+b1_sample^2);
        a1_f=sqrt(a1_flat^2+b1_flat^2);
        fi_s=atan(-b1_sample/a1_sample);
        fi_f=atan(-b1_flat/a1_flat);
        
        T=a0_s/a0_f;
        D=a1_s*a0_f/(a1_f*a0_s);
        PS=abs(fi_s-fi_f);
        V_flat=2*a1_f/a0_f;
        V_sample=2*a1_s/a0_s;
        V_relative=V_sample/V_flat;
        disp(['For pixel [', num2str(ii), ',', num2str(jj),'] T = ',num2str(T*100), '%, D = ',...
            num2str(D),', PS = ',num2str(PS*180/pi),'\[Degree], V_f = ',num2str(abs(V_flat)*100),'%, V_s = ',...
            num2str(abs(V_sample)*100),'% and V_rel is ',num2str(abs(V_relative)*100),'%.'])
    end
end

toc
