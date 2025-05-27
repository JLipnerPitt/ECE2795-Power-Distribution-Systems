clc
clear
j=1i;

% Given values
VN1 = 2500;
V1 = 2400;
V2 = 240;
pf = 0.95;
angle = acos(pf);
a = 2400/240;

% Bases
Vbase1 = 2400;
Vbase2 = 240;
Sbase1 = 25e3;
Sbase2 = 37.5e3;
Zbase1 = Vbase1^2/Sbase1;
Zbase2 = Vbase1^2/Sbase2;

% Impedances
Zline1 = (5000/5280)*(0.306+j*0.6272);
Zline2 = (2500/5280)*(0.306+j*0.6272);
ZT1pu = 0.018*exp(j*40*pi/180);
ZT2pu = 0.02*exp(j*50*pi/180);
ZT1 = ZT1pu*Zbase1;
ZT2 = ZT2pu*Zbase2;

% Max diversified real powers
P12 = 57.89e3;
P23 = 41.56e3;
PT1 = 22.71e3;
PT2 = 41.56e3;

% Complex powers
S12 = (P12/pf)*exp(j*angle);
S23 = (P23/pf)*exp(j*angle);
ST1 = (PT1/pf)*exp(j*angle);
ST2 = (PT2/pf)*exp(j*angle);

% Current into segment N1-N2
I12 = conj(S12/V1);
disp(['I12 = ', Phasor(I12)])

% Voltage at N2
VN2 = V1 - I12*Zline1;
disp(['VN2 = ', Phasor(VN2)])

% Current going into T1
IT1 = conj(ST1/VN2);
disp(['IT1 = ', Phasor(IT1)])

% Secondary voltage on T1
VT1sprime = VN2-(IT1*ZT1);
VT1s = VT1sprime/a;
disp(['VT1s = ', Phasor(VT1s)])
fprintf('\n')


%Current flowing into segment N2-N3
IN23 = conj(S23/VN2);
disp(['IN23 = ', Phasor(IN23)])

% Voltage at N3
VN3 = VN2 - IN23*Zline2;
disp(['VN3 = ', Phasor(VN3)])

% Current going into T2
IT2 = conj(ST2/VN2);
disp(['IT2 = ', Phasor(IT2)])

% Secondary voltage on T2
VT2sprime = VN3-(IT2*ZT2);
VT2s = VT2sprime/a;
disp(['VT2s = ', Phasor(VT2s)])
fprintf('\n')


function string = Phasor(z)
    %phasor: converts complex numbers to phasor form. 
    %displays phasor as a string and gives the angle in degrees.
    
    %Takes one paramter (z). z is the complex number you want to
    %convert to a phasor.
    
    magnitude = abs(z);
    %phase = round(angle(z), 3)*(180/pi);
    x = real(z);
    y = imag(z);
    phase = atand(y/x);
    if (phase < 10^-10 && phase > 0)
        phase = 0;
    end
    string = [num2str(magnitude, '%.3f'),'∠',num2str(phase, '%.6f'),'°'];
end