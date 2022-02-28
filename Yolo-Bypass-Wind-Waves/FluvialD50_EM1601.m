%Estimate D50 from flow velocities in YB

V = 2.0; %ft/s from USACE DWSE report (580,000cfs/190,000ft^2)
a = 0.254815; %slope of levee (arctan(levee height/levee length from toe to top))
kl = sqrt(1-((sin(a)^2)/(sin(0.4)^2))); %eq V1-5-129
gamr = 160; %lbf/ft^3 (stone)
gamw = 62.4; %lbf/ft^3 (fresh water)
%params from EM 1110-2-1601 Eq. 3-3
g = 32.2; %ft/s^2
Sf =1.1;
Cs = 0.3;
Cv = 1;
Ct = 1;
h = 31.30; %ft, NAVD88 from USACE DWSE report

D30 = Sf*Cs*Cv*Ct*h*(sqrt(gamw/(gamr-gamw))*(V/sqrt(kl*g*h)))^2.5;
D50 = D30*(1.4^(1/3));
fprintf('Estimated D50 based off of flow velocities in YB of %0.2f cfs: %0.2f ft\n',V,D50) 
