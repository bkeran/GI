function results=Material(varargin)
% complex refractive index
% Material(name,energy,material_density)

c=299792458;
eV=1.6021764874e-19;
h_planck = 6.6260689633e-34;
r_e=2.81794032e-15;  %[m]
u=1.660539040e-27; %[kg/m3]

if nargin ~= 3
    error('Incorrect input.');
    return;
end
name = varargin{1};
energy = varargin{2};
material_density = varargin{3};
M=mucal(name, 0, energy, 0);  % M=[Z,Ar,mu]
lambda=c*h_planck / (energy*eV*1000);
results.beta=M(3)*100*lambda/(4*pi);
results.delta=M(1)*material_density*1000*r_e*lambda^2/(2*pi*M(2)*u);
end