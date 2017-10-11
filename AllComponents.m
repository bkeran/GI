(* ::Package:: *)

classdef AllComponents < handle
% sources, samples, gratings, detectors
    properties
        c=299792458;                 % [m/s]         
        eV=1.6021764874e-19;       
        h_planck=6.6260689633e-34;   %[Js]
        r_e=2.81794032e-15;          %[m]
        u=1.660539040e-27;           %[kg]
    end
    properties
        absorb_exp         % [a] exponent in absorption term
        beta               % [1] absorption index
        dc                 % [1] duty cycle
        delta              % [1] refractive index
        E                  % [keV] energy
        G1ToG2Distance     % [m]
        grid_spacing       % [m] resolution, space between the particles
        k                  % [1/m] wave vector
        lambda             % [m] wavelength
        N                  % [1,1] number of particles
        p                  % [m] grating period 
        phShift            % [1] phase shift
        sigma              % [1] 2D Gauss filter standard deviation
        SourceToG1Distance % [m]
        source_position    
        sample_position
        grating_position
        detector_size      % [m,m] one detector pixel size
        detector_pixel     % [1,1] total number of pixels covered by simulation grid
        phase_step         % [1] total number of phase steps
        pixel_step         % [1] number of particles that are skipped in phase stepping
        talbot             % [m] fractional Talbot distance
    end
   
methods
    function Material(obj, name, material_density)
        % Method that uses external "mucal" library in order to calculate
        % complex refractive index constants beta and delta
        M=mucal(name, 0, obj.E, 0);  % mucal returns [Z,Ar,mu]
        obj.beta=M(3)*100*obj.lambda/(4*pi);
        obj.delta=M(1)*material_density*1000*obj.r_e*obj.lambda^2/(2*pi*M(2)*obj.u);
    end
    function SetGridSpacing(obj)
        % Calculating spacing between the particles for given detector size
        % and particle number.
        obj.grid_spacing=(obj.detector_pixel-1).*obj.detector_size./(obj.N-1);
    end
%--------------------------------------------------------------------------
% Sources
%--------------------------------------------------------------------------
    function psi=PlaneWaveSource(obj)
        % Function that creates plane wave source and calculates
        % wavelength, wave vector and fractional Talbot distance.
        obj.lambda=obj.c.*obj.h_planck ./ (obj.E.*obj.eV.*1000);
        obj.k=2.*pi./obj.lambda;
        obj.talbot=obj.p^2./(2*obj.lambda);
        psi=ones(obj.N(1),obj.N(2));
        psi=psi.*exp(i*obj.k*obj.source_position(3));%./ ( sqrt ((obj.N(1)-1)*(obj.N(2)-1)).*obj.grid_spacing ); 
    end
% -------------------------------------------------------------------------
% Amplitude Transmission Function
% -------------------------------------------------------------------------
    function setRefracAbsorb(obj,name,density,height)
        % Setting phase shift and exponent in absorption term for a
        % material of choice.
        obj.Material(name,density);
        obj.phShift=-obj.k*obj.delta.*height;
        obj.absorb_exp=-obj.k.*obj.beta.*height;
    end
    function tau=AmplitudeTransmission(obj)
        % Calculates amplitude transmission function for given phase shift
        % and apsorption.
        tau=exp(obj.absorb_exp+i*obj.phShift);
    end
%--------------------------------------------------------------------------
% Gratings
%--------------------------------------------------------------------------
    function G_xy=MakeGrating(obj)
        % Creates patterns of zeros and ones for 1D and 2D gratings.
        scale=[round(obj.p/obj.grid_spacing(1)),round(obj.p/obj.grid_spacing(2))];  % number of particles in one grating period
        periods=round(obj.N./scale);    % number of periods in the whole grating
        xx=round(scale*obj.dc);    % number of ones in one grating period
        ha_xx=round(xx./2);     % number of zeros in one grating period
        G_x=[ones(1,ha_xx(1)) zeros(1,scale(1)-xx(1)) ones(1,xx(1)-ha_xx(1))];
        G_y=[ones(1,ha_xx(2)) zeros(1,scale(2)-xx(2)) ones(1,xx(2)-ha_xx(2))];
        gra={G_x,G_y};        
        for j=1:2
            cornum(j)=obj.N(j)-length(gra{j})*periods(j);     % number of particles that we need to correct
            if cornum(j)==0     % if correction number is 0, we just repeat "gra" period times
                G_xy{j}=repmat(gra{j},1,periods(j));
            elseif cornum(j)<0
                erase=round(abs (cornum(j))/periods(j));     % number of extra particles that we need to erase from "gra"
                G_xy{j}=[];
                b=[gra{:,j}];
                b=b(:,1:end-erase);     % erasing excces
                gra{j}=b;
                G_xy{j}=repmat(gra{j},1,periods(j));    % repeating new "gra" periods puta and creating a grid
                cornum2=obj.N(j)-length(G_xy{j});     % additional check if there are any extra or missing particles
                if cornum2<0
                    b=[G_xy{:,j}];
                    b=b(:,1:end-abs(cornum2));      % erasing cornum2 -last- extra particles
                    G_xy{j}=b;
                elseif cornum2>0
                    b=[G_xy{:,j}];
                    for jj=1:cornum2
                        b=[b,1];     % adding cornum2 -last- missing particles
                    end
                    G_xy{j}=b;
                end
            elseif cornum(j)>0
                add=round(abs (cornum(j))/periods(j));       % number of missing particles that we need to add into "gra"
                G_xy{j}=[];
                b=[gra{:,j}];
                for ii=1:add
                    b=[b,1];        % adding missing particles
                end
                gra{j}=b;
                G_xy{j}=repmat(gra{j},1,periods(j));        % repeating new "gra" periods puta and creating a grid
                cornum2=obj.N(j)-length(G_xy{j});       % additional check if there are any extra or missing particles
                if cornum2<0
                    b=[G_xy{:,j}];
                    b=b(:,1:end-abs(cornum2));      % erasing cornum2 -last- extra particles
                    G_xy{j}=b;      % final grating
                elseif cornum2>0
                    b=[G_xy{:,j}];
                    for jj=1:cornum2
                        b=[b,1];        % adding cornum2 -last- missing particles
                    end
                    G_xy{j}=b;      % final grating
                end
            end
        end
    end
    function G=Rectangular1DGrating(obj)
        % Creates 1D rectangular grating with values of only ones and
        % zeros.
        G_xy=obj.MakeGrating();
        G=ones(length(G_xy{1}),length(G_xy{2}));
        for j=1:length(G_xy{1})
            for l=1:length(G_xy{2})
%                 G(j,l)=G (j,l).*G_xy{1}(j);      % horizontal grating
                G(j,l)=G (j,l).*G_xy{2}(l);      % vertical grating
            end
        end
     end
     function G=Rectangular2DGrating(obj)
        % Creates 2D rectangular grating with values of only ones and
        % zeros.
        G_xy=obj.MakeGrating();
        G=ones(length(G_xy{1}),length(G_xy{2}));
        for j=1:length(G_xy{1})
            for l=1:length(G_xy{2})
                G(j,l)=G (j,l).*G_xy{1}(j).*G_xy{2}(l);
            end
        end
     end
     function G=Create1DGrating(obj,supstrate,lines)
         % Creates 1D grating with arbitrary values, depending on the 
         % materials used, of it's supstrate and lines.
         G=obj.Rectangular1DGrating();
         obj.setRefracAbsorb(supstrate{1},supstrate{2},supstrate{3});
         tau_supstrate=obj.AmplitudeTransmission();
         obj.setRefracAbsorb(lines{1},lines{2},lines{3});
         tau_lines=obj.AmplitudeTransmission();
         for j=1:obj.N(1)
            for l=1:obj.N(2)
                if G(j,l)==0
                    G(j,l)=tau_supstrate;
                else
                    G(j,l)=tau_supstrate*tau_lines;
                end
            end
         end        
     end
     function G=Create2DGrating(obj,supstrate,lines)
         % Creates 2D grating with arbitrary values, depending on the 
         % materials used, of it's supstrate and lines.
         G=obj.Rectangular2DGrating();
         obj.setRefracAbsorb(supstrate{1},supstrate{2},supstrate{3});
         tau_supstrate=obj.AmplitudeTransmission();
         obj.setRefracAbsorb(lines{1},lines{2},lines{3});
         tau_lines=obj.AmplitudeTransmission();
         for j=1:obj.N(1)
            for l=1:obj.N(2)
                if G(j,l)==0
                    G(j,l)=tau_supstrate;
                else
                    G(j,l)=tau_supstrate*tau_lines;
                end
            end
         end        
     end
% -------------------------------------------------------------------------
% Samples
% -------------------------------------------------------------------------
     function S=TestSample(obj)
        % Creates test sample which has different values on first four
        % detector pixels.
        S=ones(obj.N(1),obj.N(2));
%         S(:,:)=obj.AmplitudeTransmission();
        x=linspace(0,obj.N(2)*obj.grid_spacing(2),obj.N(2));
%         x=repmat(x,[obj.N(1) 1])
        fi=pi/2*obj.p/(obj.lambda*obj.G1ToG2Distance)*x;
        for ii=1:obj.N(1)
            S(ii,:)=exp(i*fi);
        end
%         detector_grid=round(obj.N./obj.detector_pixel);        
%         S(1:detector_grid(1),1:detector_grid(2))=1/sqrt(2);
%         S(1:detector_grid(1),(detector_grid(2)+1):2*detector_grid(2))=1/sqrt(20);
%         S((detector_grid(1)+1):2*detector_grid(1),1:1*detector_grid(2))=exp(i*pi/2);
%         S((detector_grid(1)+1):2*detector_grid(1),(detector_grid(2)+1):2*detector_grid(2))=1/sqrt(2)*exp(i*pi/4);
     end
     function S=FlatSpheres(obj, material, d_spheres, p_spheres, n_spheres)
        % Creates spheres of given diameters and probabilities on random
        % positions and gives their amplitude transmission function for
        % projection approximation mode.
        proj_approx=RandomSpheres(obj.N, d_spheres, p_spheres, n_spheres, obj.grid_spacing);
%         proj_approx=RandomSpheres(obj.N, d_spheres, p_spheres, n_spheres, obj.grid_spacing);
        % RandomSpheres(dimensions of x and y,sphere diameter, diameter probability, number of spheres, distance between particles)
%         S=ones(obj.sample_size(1),obj.sample_size(2));
        S=ones(obj.N(1),obj.N(2));
        for ii=1:obj.N(1)
            for jj=1:obj.N(2)
                obj.setRefracAbsorb(material{1},material{2},proj_approx(ii,jj));
                S(ii,jj)=obj.AmplitudeTransmission();
            end
        end
    end
% -------------------------------------------------------------------------
% Propagation
% -------------------------------------------------------------------------
     function psi=SourceAfterG1(obj,source,sample,grating)
         % Function that gives values of the wavefront after passing the
         % (sample and the) grating G1. If there is no sample existing,
         % sample==0, otherwise it multiplies everything together.
         if sample==0
             psi=source.*grating;
         else
             psi=source.*sample.*grating;
         end
     end
     function G=ft2(obj,g,delta)
         % 2D numerical implementation of the Fourier transform.
         G=fftshift(fft2(fftshift(g)))*delta(1)*delta(2);
     end
     function g=ift2(obj,G,delta_f1,delta_f2)
         % 2D numerical implementation of the inverse Fourier transform.
         Ny=size(G,1);
         Nx=size(G,2);
         g=ifftshift(ifft2(ifftshift(G)))*(Nx*Ny*delta_f1*delta_f2);
     end
     function [x2, y2, Uout] = AngularSpectrumPropagation(obj, Uin, d1, d2, Dz)
         % Numerical implementation of the Fresnel propagation.
         [N_y,N_x]=size(Uin);
         [x1,y1]=meshgrid((-N_x/2 : 1: (N_x/2-1))*d1(2),(-N_y/2 : 1: (N_y/2-1))*d1(1));   % source-plane coordinates
         r1sq=x1.^2+y1.^2;
         df1_x=1/(N_x*d1(2));  % spatial frequencies in the source plane
         df1_y=1/(N_y*d1(1));
         [fX,fY]=meshgrid((-N_x/2 : 1 : (N_x/2-1))*df1_x,(-N_y/2 : 1 : (N_y/2-1))*df1_y);
         fsq=fX.^2+fY.^2;
         m=d2/d1;   % scaling parameter
         [x2,y2]=meshgrid((-N_x/2 : 1 : (N_x/2-1))*d2(2),(-N_y/2 : 1 : (N_y/2-1))*d2(1));     % observation-plane coordinates
         r2sq=x2.^2+y2.^2;
         Q1=exp(i*obj.k/2*(1-m)/Dz*r1sq);
         Q2=exp(-i*pi^2*2*Dz/m/obj.k*fsq);
         Q3=exp(i*obj.k/2*(m-1)/(m*Dz)*r2sq);      
         Uout=Q3.*obj.ift2(Q2.*obj.ft2(Q1.*Uin/m,d1),df1_x,df1_y);  % computing the propagated field
     end
     function Uout=waveFieldProp(obj,Dz,Uin)
        % calculates the itensity of a propagated wave field f along the
        % z-axis. "psi" is the absolute square root of Eq. (8).

        du=1./(obj.N.*obj.grid_spacing);  % [1/m] sampling distance in k-space
        [UU,VV]=meshgrid((-obj.N(2)/2 : 1 : (obj.N(2)/2-1))*du(2),(-obj.N(1)/2 : 1 : (obj.N(1)/2-1))*du(1));
       
        ff=obj.ft2(Uin,obj.grid_spacing);
%         ff  = fftshift(fft(ifftshift(Uin)));           % FFT of wave field
        H   = exp(-i.*pi.*obj.lambda.*Dz.*(UU.^2+VV.^2)); % Fresnel kernel
        C   = exp (i.*obj.k.*Dz)./(2*obj.SourceToG1Distance);           % const
        C   = C/abs(C);                              % normalized amplitude
%         Uout = C .* fftshift(ifft(ifftshift(ff.*H))); % convolution
        Uout=C.*obj.ift2(ff.*H,du(2),du(1));    % convolution
        
    end
     function psi_out=G1ToG2(obj,psi)
         [x,y,psi_out]=obj.AngularSpectrumPropagation(psi,obj.grid_spacing,obj.grid_spacing,obj.G1ToG2Distance);
     end
% -------------------------------------------------------------------------
% Phase Stepping
% -------------------------------------------------------------------------
     function intensity=PhaseStepping(obj,psi,grating)
         period_pixels=round(obj.p/obj.grid_spacing(2));        % number of particles in one period (x coordinate)
         obj.pixel_step=round(period_pixels/obj.phase_step);     % number of particles (x coordinate) that are skipped in phase stepping
         intensity=zeros(obj.N(1),obj.N(2));
         for jj=1:obj.phase_step
             grating_h=grating;
             psi_out=psi.*grating;
             intensity(:,:,jj)=abs (psi_out).^2;
             intensity(:,:,jj)=obj.GaussianSmoother(intensity(:,:,jj), obj.sigma);
             for kk=size(grating,2):(-1):1       % phase stepping (rearranging the grating matrix)
                if (kk-obj.pixel_step)>=1 
                    grating(:,kk)=grating_h(:,kk-obj.pixel_step);
                else
                    grating(:,kk)=grating_h(:,size(grating,2)-abs(kk-obj.pixel_step));
                end
             end
         end
     end
% -------------------------------------------------------------------------
% Detector
% -------------------------------------------------------------------------
     function result=Scale2Detector(obj,intensity)
        detector_grid=round(obj.N./obj.detector_pixel);     % number of particles that are covered by one detector pixel
        obj.detector_pixel=round(obj.detector_pixel);
        result=zeros(obj.detector_pixel(1),obj.detector_pixel(2));
        x=[0];
        y=[0];
        for kk=1:obj.detector_pixel(1)
            x=[x,kk*detector_grid(1)];
        end
        x(end)=obj.N(1);
        for ll=1:obj.detector_pixel(2)
            y=[y,ll*detector_grid(2)];
        end
        y(end)=obj.N(2);
        for jj=1:(obj.phase_step)
            intensity_tmp=zeros(obj.detector_pixel(1),obj.detector_pixel(2));
            for kk=1:(length(x)-1)
                for ll=1:(length(y)-1) 
                    delta_grid _ 1=[];
                    for mm=1:(x(kk+1)-x(kk))
                        delta_grid _ 1=[delta_grid _ 1,mm*obj.grid_spacing(1)];
                    end
                        delta_grid _ 2=[];
                    for mm=1:(y(ll+1)-y(ll))
                        delta_grid _ 2=[delta_grid _ 2,mm*obj.grid_spacing(2)];
                    end
                    intensity_crop=intensity((x(kk)+1):x(kk+1),(y(ll)+1):y(ll+1),jj);
                    tmp2=trapz(delta_grid _ 1,intensity_crop);
                    tmp3=trapz(delta_grid _ 2,tmp2);
                    intensity_tmp(kk,ll)=tmp3;
                end
            end
            result(:,:,jj)=intensity_tmp;
        end
     end
%--------------------------------------------------------------------------
% Smoother
%--------------------------------------------------------------------------
     function Iout=GaussianSmoother(obj, Iin, sigma)
         [XX,YY]=meshgrid((-obj.N(2)/2 : 1: obj.N(2)/2-1)*obj.grid_spacing(2),(-obj.N(1)/2 : 1: obj.N(1)/2-1)*obj.grid_spacing(1));
         
         sigma_proj=sigma.*obj.G1ToG2Distance./(obj.SourceToG1Distance);    % projected source size (FWHM)
         sigma_proj _sq=sigma_proj.^2/(8.*log(2));
         gaussX   = exp(-(1/2).*(XX.^2)./sigma_proj_sq(2));
         gaussY   = exp(-(1/2).*(YY.^2)./sigma_proj_sq(1));
         srcgauss = gaussX .* gaussY;
         srcgauss = srcgauss./sum(sum(srcgauss));       % normalized Gauss
         
%          gam=fftshift(fft2(ifftshift(srcgauss)));   % damping factor
%          If=fftshift(fft2(ifftshift(Iin)));          % intensity Fourier transform
%          Iout=abs(fftshift(ifft2(ifftshift(If.*gam))));
         
         gam=obj.ft2(srcgauss,obj.grid_spacing);   % damping factor
         If=obj.ft2(Iin,obj.grid_spacing);          % intensity Fourier transform
         Iout=abs(obj.ift2(If.*gam,1./(obj.N(2).*obj.grid_spacing(2)),1./(obj.N(1).*obj.grid_spacing(1))));
     end
end
end
