clear all, close all, clc

%% Initial values

d = 20; %length or diameter
P = 70; %power  -- (Target value)
U0 = 10; %wind speed
B = 3; %number of blades 
rho = 1.2;

maxTSR = 6;  %tip speed ratio
TSR = zeros(1,20);
psy = zeros(1,20);
Ftip = zeros(1,20);
% Vp = zeros(1,20);
Vr = zeros(1,20);
theta = zeros(1,20);
dr = 0.5;

%% Airfoil parameters

AoA = 5.25;
Cl = 1.0518;
Cd = 0.00813;

x = 20; %number of blade elements

%% INDUCTION FACTORS

a = 1/3*ones(1,20);
b = zeros(1,20);

tic

while true
    
    
    for n = 1:20
    
        TSR(n) = maxTSR*(n/2)/(d/2); %Local Tip speed ratio -- *
        psy(n) = atan((1-a(n))/((1+b(n))*TSR(n))); %Inflow angle -- 
%         psy(n) = (2/3) * atan ( 1/TSR(n) );
        theta(n) = (psy(n) - AoA*pi/180)*180/pi; %twist angle
 
        Ftip(n) = real((2/pi)*acos(exp(-1*((B*(d-(n/2)))/(2*d*sin(psy(n))))))); %Tip loss factor -- *
        
        RPM = TSR(20)*60*U0/(pi*d); %RPM of turbine
        Vang = 2*pi*RPM/60; %angular velocity
%         Vp = [3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60];
        Vp(n) = Vang*(n/2); %peripheral velocity
        Vr(n) = sqrt(U0^2 + Vp(n)^2); %resultant velocity
        
        Copt4(n) = ((8*pi*(n/2))/(B*Cl))* (1-cos(psy(n))); %Optimized linearization *

        solidity(n) = (B*Copt4(n))/(2*pi*(n/2)); % *
        
        Cn(n) = Cl*cos(psy(n)) + Cd*sin(psy(n));
        Ct(n) = Cl*sin(psy(n)) - Cd*cos(psy(n));
        
        a_new1(n) = 1/(1+(4*(sin(psy(n)))^2)/(solidity(n)*Cn(n)));
        b_new1(n) = real(1/((4*Ftip(n)*(sin(psy(n))*cos(psy(n)))/(solidity(n)*Ct(n)))-1));
        
        dr = dr + 0.5;

       
    end
    
 
    
        if (abs(a_new1(1)-a(1)) < 0.0001) && (abs(a_new1(2)-a(2)) < 0.0001) && (abs(a_new1(3)-a(3)) < 0.0001) && (abs(a_new1(4)-a(4)) < 0.0001) && (abs(a_new1(5)-a(5)) < 0.0001) && (abs(a_new1(6)-a(6)) < 0.0001) && (abs(a_new1(7)-a(7)) < 0.0001) && (abs(a_new1(8)-a(8)) < 0.0001) && (abs(a_new1(9)-a(9)) < 0.0001) && (abs(a_new1(10)-a(10)) < 0.0001) && (abs(a_new1(11)-a(11)) < 0.0001) && (abs(a_new1(12)-a(12)) < 0.0001) && (abs(a_new1(13)-a(13)) < 0.0001) && (abs(a_new1(14)-a(14)) < 0.0001) && (abs(a_new1(15)-a(15)) < 0.0001) && (abs(a_new1(16)-a(16)) < 0.0001) && (abs(a_new1(17)-a(17)) < 0.0001) && (abs(a_new1(18)-a(18)) < 0.0001) && (abs(a_new1(19)-a(19)) < 0.0001) && (abs(a_new1(20)-a(20)) < 0.0001) && (abs(b_new1(1)-b(1)) < 0.0001) && (abs(b_new1(2)-b(2)) < 0.0001) && (abs(b_new1(3)-b(3)) < 0.0001) && (abs(b_new1(4)-b(4)) < 0.0001) && (abs(b_new1(5)-b(5)) < 0.0001) && (abs(b_new1(6)-b(6)) < 0.0001) && (abs(b_new1(7)-b(7)) < 0.0001) && (abs(b_new1(8)-b(8)) < 0.0001) && (abs(b_new1(9)-b(9)) < 0.0001) && (abs(b_new1(10)-b(10)) < 0.0001) && (abs(b_new1(11)-b(11)) < 0.0001) && (abs(b_new1(12)-b(12)) < 0.0001) && (abs(b_new1(13)-b(13)) < 0.0001) && (abs(b_new1(14)-b(14)) < 0.0001) && (abs(b_new1(15)-b(15)) < 0.0001) && (abs(b_new1(16)-b(16)) < 0.0001) && (abs(b_new1(17)-b(17)) < 0.0001) && (abs(b_new1(18)-b(18)) < 0.0001) && (abs(b_new1(19)-b(19)) < 0.0001) && (abs(b_new1(20)-b(20)) < 0.0001)  
         break;
        end
        
        y = abs(a_new1(1)-a(1));
        z = abs(b_new1(1)-b(1));
        
        for n = 1:20
            a(n) = a_new1(n);
            b(n) = b_new1(n);
        end

end


%% Torque and Thrust

for n = 1:20 
    
     dr = 263;
     dM(n) = (1/2)*rho*Vr(n)*B*(Cl*sin(psy(n))- Cd*cos(psy(n)))*Copt4(n)*(n/2)*dr; % From genetic algorithm paper
     M = sum(dM(n)); % Torque
 
     P = M*Vang; % Power
     Cp = P / ((1/2)*rho*pi*((d/2)^2)*(U0^3)); %% coefficient of power
 
 end

 Copt4(1) = 1.6;
 rad = 0.5:0.5:10;
 plot(rad,Copt4)
 axis([0 10.5 0 2])
 
 
 toc


% for n = 1:20 
% 
%     dM(n) = (1/2)*rho*Vr(n)*B*(Cl*sin(psy(n))- Cd*cos(psy(n)))*Copt4(n)*(n/2)*dr; %genetic algorithm
%     M = sum(dM(n)); % Torque
% 
%     P = M*Vang; % Power
%     Cp = P / ((1/2)*rho*pi*((d/2)^2)*(U0^3)); %% coefficient of power
% 
% end
% 
% rad = 0.5:0.5:10;
% plot(rad,Copt4)
% 
% 
% toc