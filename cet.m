% Sample Program
clc
clear
%naming convention 
%los is line of sight
% u for upper b for bottom 
% number is number of bounce
%xyz is xyz
%g is for 24GHZ
f = 2.4*10^9; % frequency (given)
fg= 2.4*10^10;
c = 3*10^8; % speed of light
beta = 2*pi*f/c; % wave number
betag=  2*pi*fg/c;
wuy=6; wby=0;
index=sqrt(3);
% to solve for the direct ray
tx = 0; ty = 4; tz = 2; % tx [x,y,z] location
rx = 10; ry = 3; rz = 2; % rx [x,y,z] location
rxv = rx:.01:22; % from p to q

%los calc
dlos = dist(tx,ty,tz,rxv,ry,rz); % distance
Edlos = ecalc(dlos,beta); % direct ray (line of sight)
Edflos = logcalc(Edlos); % direct ray alone
Edlosg = ecalc(dlos,betag);
Edflosg = logcalc(Edlosg);

% one bounce top
tyu=8;
order=1;
[Edu1,Edfu1,Edu1g,Edfu1g]=tcalc(index,tx,tyu,tz,rxv,ry,rz,wuy,beta,betag,order);

% one bounce bottom
tyb=-4;
[Edb1,Edfb1,Edb1g,Edfb1g]=tcalc(index,tx,tyb,tz,rxv,ry,rz,wby,beta,betag,order);

% two bounce u->b
tybb=-8;
order=2;
[Edub2,Edfub2,Edub2g,Edfub2g]=tcalc(index,tx,tybb,tz,rxv,ry,rz,wby,beta,betag,order);

% two bounce b->u
tyuu=16;
[Edbu2,Edfbu2,Edbu2g,Edfbu2g]=tcalc(index,tx,tyuu,tz,rxv,ry,rz,wuy,beta,betag,order);

%three bounce u->b->u
tyuuu=20;
order=3;
[Edubu3,Edfubu3,Edubu3g,Edfubu3g]=tcalc(index,tx,tyuuu,tz,rxv,ry,rz,wuy,beta,betag,order);

%three bounce b->u->b
tybbb=-16;
[Edbub3,Edfbub3,Edbub3g,Edfbub3g]=tcalc(index,tx,tybbb,tz,rxv,ry,rz,wby,beta,betag,order);

%four bounce u->b->u->b
order=4;
tybbbb=-20;
[Edubub4,Edfubub4,Edubub4g,Edfubub4g]=tcalc(index,tx,tybbbb,tz,rxv,ry,rz,wby,beta,betag,order);

%four bounce b->u->b->u
tyuuuu=22;
[Edbubu4,Edfbubu4,Edbubu4g,Edfbubu4g]=tcalc(index,tx,tyuuuu,tz,rxv,ry,rz,wby,beta,betag,order);

%adding lines tgt
Edsum=Edlos+Edu1+Edb1+Edub2+Edbu2+Edubu3+Edbub3+Edubub4+Edbubu4;
Edfsum=logcalc(Edsum);
Edsuml=Edlos+Edu1+Edb1+Edub2+Edbu2+Edubu3+Edbub3;
Edfsum1=logcalc(Edsuml);
Edsumg=Edlosg+Edu1g+Edb1g+Edub2g+Edbu2g+Edubu3g+Edbub3g+Edubub4g+Edbubu4g;
Edfsumg=logcalc(Edsumg);
Edsumg1=Edlosg+Edu1g+Edb1g+Edub2g+Edbu2g+Edubu3g+Edbub3g;
Edfsumg1=logcalc(Edsumg1);

% graph plotting part
plot(rxv,Edflos)
hold on
plot(rxv,Edfu1)
plot(rxv,Edfb1)
plot(rxv,Edfub2)
plot(rxv,Edfbu2)
plot(rxv,Edfubu3)
plot(rxv,Edfbub3)
plot(rxv,Edfubub4)
plot(rxv,Edfbubu4)
plot(rxv,Edfsum,'r')
plot(rxv,Edfsum1,'k')
xlabel('Line Segment "pq"/m'),ylabel('E-field/dB')
grid on
title('graph of E-field agianst Line segment')
legend('los','1 bounce upper','1 bounce bottom','2 bounce ub','2 bounce bu','3 bounce ubu','3 bounce bub', ...
    '4 bounce ubub','4 bounce bubu','sum4','sum3')

% fucntion defs
function [Ed,Edf,Edg,Edfg]=tcalc(index,tx,ty,tz,rxv,ry,rz,wall,beta,betag,order)
    dt=dist(tx,ty,tz,rxv,ry,rz);
    wallx=wcalc(ry,ty,rxv,tx,wall);
    dR=dist(wallx,wall,rz,rxv,ry,rz);
    refeff=effcalc(dR,index);
    Ed=refeff.^order.*ecalc(dt,beta);
    Edf=logcalc(Ed);
    Edg=refeff.^order.*ecalc(dt,betag);
    Edfg=logcalc(Edg);
end
function dist=dist(tx,ty,tz,rx,ry,rz)
    dist=sqrt(((rx-tx).^2)+((ry-ty)^2)+((rz-tz)^2));
end
function wallx=wcalc(ry,ty,rx,tx,wall)
    m=(ry-ty)./(rx-tx);
    wallx=(wall-ty)./m;
end
function l=logcalc(Ed)
    l=20*log10(abs(Ed));
end
function E=ecalc(d,beta)
    E=(1./d).*exp(-j*beta*d);
end
function eff=effcalc(dist,index)
    anglei=acosd(3./dist);
    anglet=asind(1/index.*sind(anglei));
    eff=(cosd(anglei)-index*cosd(anglet))./(cosd(anglei)+index*cosd(anglet));
end