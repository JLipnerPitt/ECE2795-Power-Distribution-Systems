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

% Calculating power each segment sees
Sload = Smaxseg12/18;
Sload12 = Sload*18;
Sload23 = Sload*14;
Sload34 = Sload*10;
Sload45 = Sload*5;
disp(['Sload = ', Phasor(Sload/1e3), ' kVA'])
disp(['Sload12 = ', Phasor(Sload12/1e3), ' kVA'])
disp(['Sload23 = ', Phasor(Sload23/1e3), ' kVA'])
disp(['Sload34 = ', Phasor(Sload34/1e3), ' kVA'])
disp(['Sload45 = ', Phasor(Sload45/1e3), ' kVA'])
fprintf('\n')

% T1
I12 = conj(Sload12/VN1);
VN2 = VN1 - (I12*Zline1);
IT1 = conj(Sload*4/VN2);
VT1sprime = VN2 - (IT1*ZT1);
VT1s = VT1sprime/a;
disp('T1')
disp(['I12 = ', Phasor(I12)])
disp(['VN2 = ', Phasor(VN2)])
disp(['IT1 = ', Phasor(IT1)])
disp(['VT1s = ', Phasor(VT1s)])
fprintf('\n')

% T2
I23 = conj(Sload23/VN2);
VN3 = VN2 - (I23*Zline2);
IT2 = conj(Sload*4/VN3);
VT2sprime = VN3 - (IT2*ZT2);
VT2s = VT2sprime/a;
disp('T2')
disp(['I23 = ', Phasor(I23)])
disp(['VN3 = ', Phasor(VN3)])
disp(['IT2 = ', Phasor(IT2)])
disp(['VT2s = ', Phasor(VT2s)])
fprintf('\n')

% T3
I34 = conj(Sload34/VN3);
VN4 = VN3 - (I34*Zline3);
IT3 = conj(Sload*5/VN3);
VT3sprime = VN3 - (IT3*ZT3);
VT3s = VT3sprime/a;
disp('T3')
disp(['I34 = ', Phasor(I34)])
disp(['VN4 = ', Phasor(VN4)])
disp(['IT3 = ', Phasor(IT3)])
disp(['VT3s = ', Phasor(VT3s)])
fprintf('\n')

% T4
I45 = conj(Sload45/VN4);
VN5 = VN4 - (I45*Zline4);
IT4 = conj(Sload*4/VN5);
VT4sprime = VN5 - (IT4*ZT4);
VT4s = VT4sprime/a;

disp('T4')
disp(['I34 = ', Phasor(I45)])
disp(['VN4 = ', Phasor(VN5)])
disp(['IT3 = ', Phasor(IT4)])
disp(['VT3s = ', Phasor(VT4s)])
fprintf('\n')
fprintf('\n')


disp('(d)')

% Calculating currents drawn from each transformer
Iload = I12/18;
Iload12 = Iload*18;
Iload23 = Iload*14;
Iload34 = Iload*10;
Iload45 = Iload*5;

% Calculating voltages at each node
VN2 = VN1 - (Iload12*Zline1);
VN3 = VN2 - (Iload23*Zline2);
VN4 = VN3 - (Iload34*Zline3);
VN5 = VN4 - (Iload45*Zline4);

% Calculating transformer secondary voltages
VT1sprime = VN2 - (Iload*4*ZT1);
VT1s = VT1sprime/a;

VT2sprime = VN3 - (Iload*4*ZT2);
VT2s = VT2sprime/a;

VT3sprime = VN4 - (Iload*5*ZT3);
VT3s = VT3sprime/a;

VT4sprime = VN5 - (Iload*5*ZT4);
VT4s = VT4sprime/a;

% Displaying results
disp(['Iload = ', Phasor(Iload), ' A'])
disp(['Iload12 = ', Phasor(Iload12), ' A'])
disp(['Iload23 = ', Phasor(Iload23), ' A'])
disp(['Iload34 = ', Phasor(Iload34), ' A'])
disp(['Iload45 = ', Phasor(Iload45), ' A'])
fprintf('\n')

disp(['VT1s = ', Phasor(VT1s)])
disp(['VT2s = ', Phasor(VT2s)])
disp(['VT3s = ', Phasor(VT3s)])
disp(['VT4s = ', Phasor(VT4s)])
fprintf('\n')
fprintf('\n')


disp('(e)')
% Calculating alocation factor
AF = Smaxseg12/175e3;
disp(['AF = ', Phasor(AF)])

% Calculating transformer powers
ST1 = AF*S1base;
ST2 = AF*S2base;
ST3 = AF*S3base;
ST4 = AF*S4base;
S12 = Smaxseg12;
S23 = ST2 + ST3 + ST4;
S34 = ST3 + ST4;
S45 = ST4;

% Current into segment N1-N2
I12 = conj(S12/VN1);
disp(['I12 = ', Phasor(I12)])

% Voltage at N2
VN2 = VN1 - I12*Zline1;
disp(['VN2 = ', Phasor(VN2)])

% Current going into T1
IT1 = conj(ST1/VN2);
disp(['IT1 = ', Phasor(IT1)])

% Secondary voltage on T1
VT1sprime = VN2-(IT1*ZT1);
VT1s = VT1sprime/a;
disp(['VT1s = ', Phasor(VT1s)])
fprintf('\n')


% Current flowing into segment N2-N3
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


% Current flowing into segment N3-N4
IN34 = conj(S34/VN3);
disp(['IN34 = ', Phasor(IN34)])

% Voltage at N4
VN4 = VN3 - IN34*Zline3;
disp(['VN4 = ', Phasor(VN4)])

% Current going into T3
IT3 = conj(ST3/VN4);
disp(['IT3 = ', Phasor(IT3)])

% Secondary voltage on T3
VT3sprime = VN4-(IT3*ZT3);
VT3s = VT3sprime/a;
disp(['VT3s = ', Phasor(VT3s)])
fprintf('\n')

% Current flowing into segment N4-N5
IN45 = conj(S45/VN4);
disp(['IN45 = ', Phasor(IN45)])

% Voltage at N5
VN5 = VN4 - IN45*Zline4;
disp(['VN5 = ', Phasor(VN5)])

% Current going into T4
IT4 = conj(ST4/VN5);
disp(['IT4 = ', Phasor(IT4)])

% Secondary voltage on T4
VT4sprime = VN5-(IT4*ZT4);
VT4s = VT4sprime/a;
disp(['VT4s = ', Phasor(VT4s)])
fprintf('\n')
fprintf('\n')


disp('(f)')


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


