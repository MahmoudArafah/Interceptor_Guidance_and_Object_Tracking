clc
close all
clear

dT = 0.1;
time=0:dT:100;

D0 = 37000;
VT = 200;
sgma=40;
sgma=(sgma./180).*pi;
epsilon0 = 35;
epsilon0 = (epsilon0/180)*pi;
thetaT = 0;
thetaT = (thetaT/180)*pi;
phiT=0;
phiT=(phiT./180).*pi;
g=10000;

XT0 = D0*cos(epsilon0).*cos(sgma);
YT0 = D0*sin(epsilon0);
ZT0 = D0*cos(epsilon0)*sin(sgma);
VT_hor = VT*cos(thetaT).*cos(phiT);
VT_ver = VT*sin(thetaT);
VT_z   = VT*cos(thetaT)*sin(phiT);
XT = XT0;
YT = YT0;
ZT = ZT0;


XM0 = 0;
YM0 = 0;
ZM0 = 0;
XM = XM0;
YM = YM0;
ZM = ZM0;

P = 4;
VM = P*VT;



phiM =sgma;
% thetaM = atan2((YT-YM),(XT-XM));
thetaM = epsilon0
VM_hor = VM*cos(thetaM)*cos(phiM);
VM_ver = VM*sin(thetaM);
VM_z  = VM*cos(thetaM)*sin(phiM);

K1=4;
K2=1;



epsilon = epsilon0;
D = D0;

figHandle = figure;
grid on
flag = 0;

for i=2:length(time)
    
    VT_hor = VT*cos(thetaT).*cos(phiT);
    VT_ver = VT*sin(thetaT);
    VT_z   = VT*cos(epsilon0)*sin(phiT);
    XT = XT + VT_hor*dT;
    YT = YT + VT_ver*dT;
    ZT = ZT + VT_z*dT;
    
    XT_vector(i) = XT; 
    YT_vector(i) = YT;
    ZT_vector(i) = ZT;
    
    D_dot = (VT*cos(sgma-phiT)*cos(epsilon - thetaT) - ...
            VM*cos(sgma-phiM)*cos(epsilon-thetaM ));%thetaM - epsilon +
    Dn = D + D_dot * dT; 
        D_vector(i) = Dn;

    if abs(Dn) > abs(D)
        flag = flag+1;
    end
    if flag == 2 
        break;
    end
 
    
    D = Dn;
 
    epsilon_dot = (-VT*cos(sgma-phiT)*sin(epsilon-thetaT) + ...
                    VM*cos(sgma-phiM)*sin(epsilon-thetaM ))/D; 
    epsilon = epsilon + epsilon_dot * dT; 
    thetaM =thetaM+ K1.*epsilon_dot.*dT;
    thetaM_vector(i)=thetaM;

    sgma_dot = ((-VT*sin(sgma - phiT)*cos(thetaT) + ...
                    VM*sin(sgma-phiM)*cos(thetaM)))/(D.*cos(epsilon));
    sgma = sgma + sgma_dot * dT; 
    phiM =phiM+ K2.*sgma_dot.*dT;
    phiM_vector(i)=phiM;

    
    VM_hor = VM*cos(thetaM)*cos(sgma);
    VM_ver = VM*sin(thetaM);
    VM_z =VM*cos(thetaM)*sin(sgma);
    XM = XM + VM_hor*dT;
    YM = YM + VM_ver*dT;
    ZM = ZM + VM_z*dT;

    JNmy(i) = VM .* (thetaM_vector(i)-thetaM_vector(i-1))./dT./g;
    JNmz(i) = VM .* cos(thetaM) .* (phiM_vector(i)-phiM_vector(i-1))./dT./g;
    



    hold on;
    box on;
    plot3(XT,ZT,YT,'r*',XM,ZM,YM,'k*');    
%          
%     
    xlabel('X-coordinate (meters)');
    zlabel('Y-coordinate (meters)');
    ylabel('Z-coordinate (meters)');
%    imgframe =0;
%     imgframe = imgframe+1;
%     images1(imgframe) = getframe(figHandle);
end


  le_temp=1:length(JNmz);
figure(2)
 plot(le_temp,JNmy,le_temp,JNmz)
grid on
xlabel('time [sec]')
ylabel('Normal acceleration  [g]')
legend('JNmy','JNmz')
% figure(3)
% plot(le_temp,thetaM_vector)
% figure(4)
% plot(le_temp,rad2deg(phiM_vector))
    axis([0 500 -1 1])

