clc
clear
j=1i;

% Given values
Smax = 15.5e3 + j*7.5e3;
Zline1 = (0.4421+j*0.3213)*(300/1000);
Zline2 = (0.4421+j*0.3213)*(470/1000);
Zline3 = (0.4421+j*0.3213)*(750/1000);
Zline4 = (0.4421+j*0.3213)*(820/1000);
Vp = 2400;
Vs = 240;
a = Vp/Vs;
ZT1pu = 0.01 + j*0.03;
ZT2pu = 0.01 + j*0.03;
ZT3pu = 0.015 + j*0.035;
ZT4pu = 0.015 + j*0.035;

% Bases
S1base = 37.5e3;
S2base = 37.5e3;
S3base = 50e3;
S4base = 50e3;

Vpbase = 2400;
Vsbase = 240;

Z1base = Vpbase^2/S1base;
Z2base = Vpbase^2/S2base;
Z3base = Vpbase^2/S3base;
Z4base = Vpbase^2/S4base;

ZT1 = Z1base*ZT1pu;
ZT2 = Z2base*ZT2pu;
ZT3 = Z3base*ZT3pu;
ZT4 = Z4base*ZT4pu;

% (a)
ST1maxdiversified = 4*Smax;
ST2maxdiversified = 4*Smax;
ST3maxdiversified = 5*Smax;
ST4maxdiversified = 5*Smax;

disp('(a)')
disp(['ST1maxdiversified = ', Format(ST1maxdiversified, 1e3), ' kVA'])
disp(['ST2maxdiversified = ', Format(ST2maxdiversified, 1e3), ' kVA'])
disp(['ST3maxdiversified = ', Format(ST3maxdiversified, 1e3), ' kVA'])
disp(['ST4maxdiversified = ', Format(ST4maxdiversified, 1e3), ' kVA'])
fprintf('\n')


Smaxseg12 = ST1maxdiversified + ST2maxdiversified + ST3maxdiversified + ST4maxdiversified;
Smaxseg23 = ST2maxdiversified + ST3maxdiversified + ST4maxdiversified;
Smaxseg34 = ST3maxdiversified + ST4maxdiversified;
Smaxseg45 = ST4maxdiversified;

disp(['Smaxseg12 = ', Format(Smaxseg12, 1e3), ' kVA'])
disp(['Smaxseg23 = ', Format(Smaxseg23, 1e3), ' kVA'])
disp(['Smaxseg34 = ', Format(Smaxseg34, 1e3), ' kVA'])
disp(['Smaxseg45 = ', Format(Smaxseg45, 1e3), ' kVA'])
fprintf('\n')
fprintf('\n')

disp('(b)')

% T1
VN1 = 2600;
I12 = conj(Smaxseg12/VN1);
VN2 = VN1 - I12*Zline1;
IT1 = conj(ST1maxdiversified/VN2);
VT1sprime = VN2 - (IT1*ZT1);
VT1s = VT1sprime/a;

disp('T1')
disp(['I12 = ', Phasor(I12)])
disp(['VN2 = ', Phasor(VN2)])
disp(['IT1 = ', Phasor(IT1)])
disp(['VT1s = ', Phasor(VT1s)])
fprintf('\n')

% T2
I23 = conj(Smaxseg23/VN2);
VN3 = VN2 - I23*Zline2;
IT2 = conj(ST2maxdiversified/VN3);
VT2sprime = VN2 - (IT2*ZT2);
VT2s = VT2sprime/a;

disp('T2')
disp(['I23 = ', Phasor(I23)])
disp(['VN3 = ', Phasor(VN3)])
disp(['IT2 = ', Phasor(IT2)])
disp(['VT2s = ', Phasor(VT2s)])
fprintf('\n')

% T3
I34 = conj(Smaxseg34/VN3);
VN4 = VN3 - I34*Zline3;
IT3 = conj(ST3maxdiversified/VN3);
VT3sprime = VN3 - (IT3*ZT3);
VT3s = VT3sprime/a;

disp('T3')
disp(['I34 = ', Phasor(I34)])
disp(['VN4 = ', Phasor(VN4)])
disp(['IT3 = ', Phasor(IT3)])
disp(['VT3s = ', Phasor(VT3s)])
fprintf('\n')

% T4
I45 = conj(Smaxseg45/VN4);
VN5 = VN4 - I45*Zline4;
IT4 = conj(ST4maxdiversified/VN4);
VT4sprime = VN4 - (IT4*ZT4);
VT4s = VT4sprime/a;

disp('T4')
disp(['I34 = ', Phasor(I45)])
disp(['VN4 = ', Phasor(VN5)])
disp(['IT3 = ', Phasor(IT4)])
disp(['VT3s = ', Phasor(VT4s)])
fprintf('\n')
fprintf('\n')

disp('(c)')
Iload1 = Smaxseg12/18;
Iload2 = Smaxseg23/18;
Iload3 = Smaxseg34/18;
Iload4 = Smaxseg45/18;

disp(['Iload1 = ', Phasor(Iload1)])
disp(['Iload2 = ', Phasor(Iload2)])
disp(['Iload3 = ', Phasor(Iload3)])
disp(['Iload4 = ', Phasor(Iload4)])
fprintf('\n')
fprintf('\n')

disp('(d)')




function string = Format(z, power)
    %phasor: Formats complex numbers for usage in the 'disp' function.
    
    %Takes two paramters (z, power). z is the complex number you want to
    %format and power is used for unit formatting.
    
    x = real(z)/power;
    y = imag(z)/power;
    string = [num2str(x, '%.3f'),'+ j*',num2str(y, '%.3f')];
end


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


