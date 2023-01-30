//Scilab code for Computational Analysis of Orbital Trajectories of Marbles on a 3D Printed Surface.
//Code written by Pooja and Preety, from the
//Orbital Mechanics Through 3D Printing Group,
//Department of Physics, and D.S. Kothari Center for Research and Innovation in Science Education,
//Miranda House, University of Delhi.
//Undergraduate Research Project (June 2022). Code release date: January 30, 2023.
//We welcome all constructive suggestions/feedback. Please reach us at 
//orbitalmechanics3dprintinggroup_mh@googlegroups.com 

//----------------------------------CODE----------------------------------
clear
clc
xdel(winsid())

sheet0 = readxls("C:\Users\user\Downloads\Closed _Orbit_30_FPS.xls")
s0 = sheet0(1)
T_1 = s0.value()
T_2 = s0.text()
Track_Data = T_1(3:$, 2:3)
xT = Track_Data(:,1)
yT = Track_Data(:,2)
N = length(xT)

c1 = 1.5208506
c2 = -12.071687
c3 = 6.4208506

for i = 1:1:N
    rT(i) = sqrt(xT(i)^2 + yT(i)^2)
    zT(i) = c1 + c2/rT(i)
end

//Centroid:
xm = sum(xT)/N
ym = sum(yT)/N
zm = sum(zT)/N

//Centered Data:
rx = xT - xm
ry = yT - ym
rz = zT - zm

//-----------------Plane Fitting----------------------//
//~~~~~~~~Approach 1: PCA~~~~~~~~~~~//
xvar = sum(rx.*rx)/N
yvar = sum(ry.*ry)/N
zvar = sum(rz.*rz)/N
xycov = sum(rx.*ry)/N
yzcov = sum(ry.*rz)/N
xzcov = sum(rx.*rz)/N

varM = [xvar, xycov, xzcov; xycov, yvar, yzcov; xzcov, yzcov, zvar]
[vec_Cov, val_Cov] = spec(varM)
eigPCA = diag(val_Cov)
[val_P, index_P] = min(eigPCA)
n1 = vec_Cov(:, index_P)
disp(n1, "Normal Vector (PCA): ")

//~~~~~~~~Approach 2: SVD~~~~~~~~~~//
PD(1,:) = rx
PD(2,:) = ry
PD(3,:) = rz
[PU, PS, PV] = svd(PD)
po = min(size(PD))
//smallest singular value:
for j = 1:1:po
    for i = 1:1:length(PS(:,j))
        if PS(i,j) > 0 then
            psVal(j) = PS(i,j)
        end    
    end    
end
if max(size(psVal)) ~= po then
    disp("psVal does not contain all singular values!")
end
[psmin, pindex] = min(psVal)
n2 = PU(:, pindex)
disp(n2, "Normal Vector (SVD):")

//To compare results of PCA and SVD:
n_cross = cross(n1,n2)
mag_cross = norm(n_cross)
if mag_cross > 1D-8 then
    disp("Warning: PCA and SVD results are different !") 
end    

//For Kepler's First Law:-------------------//
for k = 1:1:N-2
    j =k
    rf1 = [xT(j) yT(j) zT(j)]
    rf2 = [xT(j+1) yT(j+1) zT(j+1)]
    rf3 = [xT(j+2) yT(j+2) zT(j+2)]
    
    rf13 = rf3 - rf1
    rf12 = rf2 - rf1
    nf = cross(rf13, rf12)
    unf = nf/norm(nf)
    N_Matrix(k,:) = unf(1,1:3)
    
    //Find Angle Between Planes:
    vc = cross(n2', unf)
    vcsign = sign(vc(3))
    dot = sum(n2'.*unf)
    angle(k) = acos(dot)
    if angle(k) > %pi/2 then
        angle(k) = %pi - angle(k)
    end    
    Angle(k) = vcsign*angle(k)*180/%pi
end
//-------------------------------------------------------//

//-----------Projection on Best Fitted Plane--------------//
//distance of points from best fitted plane:
for i = 1:1:N
    dist_1(i) = rx(i)*n1(1) + ry(i)*n1(2) + rz(i)*n1(3) 
    dist_2(i) = rx(i)*n2(1) + ry(i)*n2(2) + rz(i)*n2(3) 
end

//Projected Points: (w.r.t original coordinate system)
for i = 1:1:N
    imP_PCA(i,:) = [xT(i) yT(i) zT(i)] - dist_1(i)*n1'
    imP_SVD(i,:) = [xT(i) yT(i) zT(i)] - dist_2(i)*n2'
end

//----------------3D to 2D Conversion---------------------//
//Method 1: By Rotating the Coordinate Axes:
function [BasisR, data] = Rotation(xR, yR, zR)
    origin = [xm; ym; zm]
    xNew = xR - origin(1)
    yNew = yR - origin(2)
    zNew = zR - origin(3)
    
    i=1
    p1 = [origin(1) origin(2) origin(3)]
    p2 = [xR(i+4) yR(i+4) zR(i+4)]
    p3 = [xR(i+8) yR(i+8) zR(i+8)]
    p21 = p1 - p2
    mp21 = norm(p21)
    p13 = p3 - p1
    mp13 = norm(p13)
    
    //First Row of R Vector:
    ex1 = p21(1)/mp21
    ex2 = p21(2)/mp21
    ex3 = p21(3)/mp21
    ex = [ex1 ex2 ex3]
    
    //3rd Row of R vector:
    Ez1 = p13(3)*ex2 - p13(2)*ex3
    Ez2 = p13(1)*ex3 - p13(3)*ex1
    Ez3 = p13(2)*ex1 - p13(1)*ex2
    Ez = [Ez1 Ez2 Ez3]
    mEz = norm(Ez)
    
    ez1 = Ez1/mEz
    ez2 = Ez2/mEz
    ez3 = Ez3/mEz
    ez = [ez1 ez2 ez3]
    
    //2nd Row of R vector:
    Ey1 = ez3*ex2 - ez2*ex3
    Ey2 = ez1*ex3 - ez3*ex1
    Ey3 = ez2*ex1 - ez1*ex2
    Ey = [Ey1 Ey2 Ey3]
    mEy = norm(Ey)
    
    ey1 = Ey1/mEy
    ey2 = Ey2/mEy
    ey3 = Ey3/mEy
    ey = [ey1 ey2 ey3]
    
    err = 1D-8
    if norm(ex)-1>err | norm(ey)-1>err | norm(ez)-1>err then
        disp("Warning: norm is not equal to 1.")    
    end
    
    BasisR = [ex(1) ex(2) ex(3); ey(1) ey(2) ey(3); ez(1) ez(2) ez(3)]
    for i = 1:1:N
        data(:,i) = BasisR*[xNew(i); yNew(i); zNew(i)]
    end
    data = data'
endfunction
xR1 = imP_PCA(:,1)
yR1 = imP_PCA(:,2)
zR1 = imP_PCA(:,3)
[Basis_1, Data_1] = Rotation(xR1, yR1, zR1)
xR2 = imP_SVD(:,1)
yR2 = imP_SVD(:,2)
zR2 = imP_SVD(:,3)
[Basis_2, Data_2] = Rotation(xR2, yR2, zR2)

//Method 2: Using Basis Vectors from PCA and SVD:
eigP = eigPCA
s_val = psVal
for i = 1:1:3
    [val, in] = max(abs(eigP))
    eigP(in) = 0
    Basis_3(i,:) = vec_Cov(:,in)
    [val, in] = max(abs(s_val))
    s_val(in) = 0
    Basis_4(i,:) = PU(:,in)
end

for i = 1:1:N
    Data_3(:,i) = Basis_3*[xR1(i)-xm; yR1(i)-ym; zR1(i)-zm]
    Data_4(:,i) = Basis_3*[xR2(i)-xm; yR2(i)-ym; zR2(i)-zm]
    Data_5(:,i) = Basis_4*[xR1(i)-xm; yR1(i)-ym; zR1(i)-zm]
    Data_6(:,i) = Basis_4*[xR2(i)-xm; yR2(i)-ym; zR2(i)-zm]
end
Data_3 = Data_3'
Data_4 = Data_4'
Data_5 = Data_5'
Data_6 = Data_6'
//--------------------------------------------------------//

//-----------------Ellipse Fitting-------------------------//
xp = Data_2(:,1)
yp = Data_2(:,2)

xavg = sum(xp)/N
yavg = sum(yp)/N

xnew = xp - xavg
ynew = yp - yavg

//~~~~Ellipse Fit Method 1 - Direct Ellipse Fitting:~~~~~
D1 = [xnew.*xnew, xnew.*ynew, ynew.*ynew]
D2 = [xnew, ynew, ones(N,1)]
S1 = D1'*D1
S2 = D1'*D2
S3 = D2'*D2
T = -inv(S3)*S2'
M = S1 + S2*T
M = [M(3,:)./2; -M(2,:); M(1,:)./2]
[Evec, Eval] = spec(M)
Cond = 4*Evec(1,:).*Evec(3,:) - Evec(2,:).^2
p_1 = Evec(:, find(Cond>0))
q_1 = [p_1; T*p_1]

//Ellipse Parameters:
R(:,1) = q_1

//~~~~~~~~~Ellipse Fit Method 2 - FITZGIBBON'S APPROACH~~~~~~~~~~
D = [xnew.*xnew, xnew.*ynew, ynew.*ynew, xnew, ynew, ones(N,1)]
Sn = D'*D
C(6,6) = 0
C(1,3) = 2
C(2,2) = -1
C(3,1) = 2
[gevec, geval] = spec(inv(Sn)*C)
[PosR, PosC] = find(geval > 0 & ~isinf(geval))
q_2 = gevec(:, PosC)

R(:,2) = q_2

//~~~~~~~~~~~~~Ellipse Fit Method 3 - SVD:~~~~~~~~~~~~~~
QD(1,:) = xnew.*xnew
QD(2,:) = xnew.*ynew
QD(3,:) = ynew.*ynew
QD(4,:) = xnew
QD(5,:) = ynew
QD(6,:) = ones(1,N)
[U, S, V] = svd(QD)
o = min(size(QD))
//Smallest Singular Value:
for j = 1:1:o
    for i = 1:1:length(S(:,j))
        if S(i,j) > 0 then
            sVal(j) = S(i,j)
        end    
    end    
end

if max(size(sVal)) ~= o then
    disp("sVal does not conatin all singular values!")
end
[smin, index] = min(sVal)
q_3 = U(:, index)
R(:,3) = q_3

//~~~~~~~~~~~Ellipse Fit Method 4 - LSF:~~~~~~~~~~~~~~
QA = [xnew.*xnew, xnew.*ynew, ynew.*ynew, xnew, ynew]
w = ones(N,1)
q_4 = QA\w
q_4(6) = -1
R(:,4) = q_4

//To compare results of different ellipse fitting approaches:
for i = 1:1:4
    if R(6, i) > 0 | R(6,i) ~= -1 then
        Rnew(:,i) = -R(:,i)/R(6,i)
    else
        Rnew(:,i) = R(:,i)    
    end            
end
disp("Ellipse Parameters:")
disp([" Method 1 ", " Method 2 ", " Method 3", "  Method 4"])
disp(Rnew)

dR1 = Rnew(:,1)-Rnew(:,2)
dR2 = Rnew(:,1)-Rnew(:,3)
dR3 = Rnew(:,1)-Rnew(:,4)
for i = 1:1:length(dR1)
    if abs(dR1(i)) > 1D-03 | abs(dR2(i)) > 1D-03 | abs(dR3(i)) > 1D-03 then
        disp("Warning: Difference in ellipse parameters is greater than 1D-03: ")
        disp(dR1(i), "Diff 1: ")
        disp(dR2(i), "Diff 2: ")
        disp(dR3(i), "Diff 3: ")
    end    
end
//-----------------------------------------------------------//

//-----------Ellipse Parameters and Applications---------------//
i = 1
a0 = R(1,i)
b0 = R(2,i)
c0 = R(3,i)
d0 = R(4,i)
e0 = R(5,i)
f0 = R(6,i)

//Principal Axis Transformation:
tol = 1D-8
if abs(b0) > 1D-08 then
    theta = 0.5*acot((a0-c0)/b0)
    a = a0*cos(theta)^2 + b0*sin(theta)*cos(theta) + c0*sin(theta)^2
    b = b0*(cos(theta)^2 - sin(theta)^2) + 2*(c0-a0)*sin(theta)*cos(theta)
    c = a0*sin(theta)^2 - b0*sin(theta)*cos(theta) + c0*cos(theta)^2
    d = d0*cos(theta) + e0*sin(theta)
    e = -d0*sin(theta) + e0*cos(theta)
    f = f0
else
    theta = 0
    a = a0; b = b0; c = c0; d = d0; e = e0; f = f0    
end

//Ellipse Parameters Extraction:
k = (d^2)/(4*a) + (e^2)/(4*c) - f
a1 = sqrt(k/a)
b1 = sqrt(k/c)
xc = -d/(2*a)
yc = -e/(2*c)


if a1 > b1 then
    sa = a1
    sb = b1
    ec = sqrt(1-(sb/sa)^2)
    xf = xc + sa*ec
    yf = yc
    
    rp = sa*(1-ec)
    ra = sa*(1+ec)
    
    xPg = xf + rp
    yPg = yf
    xAg = xf - ra
    yAg = yf
    disp(" --- Major axis is along x-axis ! ---")
else
    sa = b1
    sb = a1
    ec = sqrt(1-(sb/sa)^2)
    xf = xc 
    yf = yc + sa*ec
    
    rp = sa*(1-ec)
    ra = sa*(1+ec)
    
    xPg = xf
    yPg = yf + rp
    xAg = xf
    yAg = yf - ra
    disp(" --- Major axis is along y-axis ! ---")        
end

//Time Period Estimation:
t_Tracker = T_1(3:$,1)
d_T = diff(t_Tracker)
dT = d_T(1)
TimeP_Tracker = t_Tracker(N)-t_Tracker(1)

xP1 = xR2 - mean(xR2)
yP1 = yR2 - mean(yR2)
zP1 = zR2 - mean(zR2)

for i = 1:1:N
    R_3D(i,:) = [xP1(i) yP1(i) zP1(i)]
    mR_3D(i) = norm(R_3D(i,:))
    R_3Dcap(i,:) = R_3D(i,:)/mR_3D(i)
end

for i = 1:1:N
    Rdot(i) = sum(R_3Dcap(i,:).*R_3Dcap(1,:))
    R_angle(i) = acos(Rdot(i))
end

for i = 1:1:N
    c_RP(i,:) = cross([xP1(i), yP1(i), zP1(i)], [xP1(1), yP1(1), zP1(1)]) 
end

for i = 2:1:N
    if c_RP(3,3) < 0 then
        if c_RP(i,3) > 0 | i == N then
            R_angle(i) = 2*%pi - R_angle(i)
        end
    else
        if c_RP(i,3) < 0 | i == N then
            R_angle(i) = 2*%pi - R_angle(i)
        end                  
    end 
end

for i = 1:1:N
    dev_90(i) = %pi/2 - R_angle(i)
    dev_180(i) = %pi - R_angle(i)
    dev_270(i) = 3*%pi/2 - R_angle(i)
end
cc = 180/%pi
[aP1, bP1] = min(abs(dev_90))
[aP2, bP2] = min(abs(dev_180))
[aP3, bP3] = min(abs(dev_270))
aP = [aP1; aP2; aP3]

Time_90 = 4*(bP1-1)*dT
Time_180 = 2*(bP2-1)*dT
Time_270 = (4/3)*(bP3-1)*dT
TimeP = [Time_90; Time_180; Time_270]

[ta, tb] = min(aP)
Time_Period = TimeP(tb)

//Mean Anomaly and Eccentric Anomaly:
err = 1D-8
N2 = 2001
t = linspace(0, Time_Period, N2)
dt = t(2) - t(1)

for i = 1:1:N2
    Me(i) = (2*%pi/Time_Period)*t(i)
    if Me(i) < %pi then
        E(i) = Me(i) + ec/2
    else
        E(i) = Me(i) - ec/2
    end
    ratioNR = 1
    while abs(ratioNR) > err  
        ratioNR = (E(i) - ec*sin(E(i)) - Me(i))/(1 - ec*cos(E(i)))
        E(i) = E(i) - ratioNR
    end        
end

//Ellipse Coordinates:
if a1 > b1 then
    x = sa*cos(E) + xc
    y = sb*sin(E) + yc
else
    x = sb*sin(E) + xc
    y = sa*cos(E) + yc    
end

//Gravitation Parameter u and Angular Momentum h:
for i = 1:1:N2-1
    vx(i) = (x(i+1) - x(i))/dt
    vy(i) = (y(i+1) - y(i))/dt
    vz(i) = 0
    v(i) = sqrt(vx(i)^2 + vy(i)^2 + vz(i)^2)
    
    xmid(i) = (x(i+1) + x(i))/2
    ymid(i) = (y(i+1) + y(i))/2
    xRel(i) = xmid(i) - xf
    yRel(i) = ymid(i) - yf
    zRel(i) = 0
    r(i) = sqrt(xRel(i)^2 + yRel(i)^2 + zRel(i)^2)
    
    u1(i) = 0.5*(v(i)^2)/(1/r(i) - 1/(2*sa))
    G = 6.674D-11
    M1(i) = (1/G)*u1(i)
    H(i,:) = cross([xRel(i) yRel(i) zRel(i)],[vx(i) vy(i) vz(i)])
    magH(i) = sqrt(H(i,1)^2 + H(i,2)^2 + H(i,3)^2)
end

u_check = diff(u1)
for i = 1:1:N2-2
    if abs(u_check(i)) > 1D-1 then
        disp("Warning: all u1 values are not same ! u_check = " +string(u_check(i)))
    end    
end

u = mean(u1)
Mass = mean(M1)
h1 = mean(magH)
h = sqrt(sa*u*(1 - ec^2))

//Other Quantities:
//Perigee and Apogee Velocity:
vp = h/rp
va = h/ra
//True Anomaly Averaged Radius:
ravg_True = sqrt(rp*ra)
//True Anomaly when r = r_True:
True_A = acos((1/ec)*(h*h/(u*ravg_True)-1)) 
//Speed of marble when r = r_True:
v_True = sqrt((2*u)*(1/ravg_True - 1/(2*sa)))
//Flight Path Angle:
Gamma = atan(ec*sin(True_A)/(1+ec*cos(True_A)))

//LRL Vector:
for i = 1:1:N2-1
    LRL(i,:) = cross([vx(i) vy(i) vz(i)],[H(i,1) H(i,2) H(i,3)]) - u*[xRel(i) yRel(i) zRel(i)]/r(i)
    mLRL(i) = sqrt(LRL(i,1)^2 + LRL(i,2)^2 + LRL(i,3)^2)
end

//Areal Velocity Calculation:
ThetaCal = 2*atan(sqrt((1+ec)/(1-ec))*tan(E./2))
for i = 1:1:N2
    Area_1(i) = 0.5*sa*sb*[E(i) - ec*sin(E(i))]
    Area_2(i) = 0.5*sa*sb*[2*atan(sqrt((1-ec)/(1+ec))*tan(ThetaCal(i)/2)) - (ec*sqrt(1-ec^2)*sin(ThetaCal(i))/(1+ec*cos(ThetaCal(i))))]
end

E_Area = %pi*sa*sb
for i = 2:1:N2
    if Area_2(i) < 0 | i == N2 then
        Area_2(i) = E_Area - abs(Area_2(i))
    end    
end

//Using Geometrical Argument:
function F1= AreaCal(a,b,th)
    F1 = (a*b/2)*[th - atan((b-a)*sin(2*th)/(b+a+(b-a)*cos(2*th)))]
endfunction

for i = 1:1:N2
    Cross_1(i,:) = cross([x(i)-xc y(i)-yc 0],[x(1)-xc y(1)-yc 0])
    Cross_2(i,:) = cross([x(i)-xf y(i)-yf 0],[x(1)-xf y(1)-yf 0])
end
for i = 1:1:N2
    Rc(i,:) = [x(i)-xc y(i)-yc]
    mRc(i) = norm(Rc(i,:))
    Rc(i,:) = Rc(i,:)/norm(Rc(i,:))
    Rf(i,:) = [x(i)-xf y(i)-yf]
    mRf(i) = norm(Rf(i,:))
    Rf(i,:) = Rf(i,:)/norm(Rf(i,:))
    
    dot_1(i) = sum(Rc(i,:).*Rc(1,:))
    dot_2(i) = sum(Rf(i,:).*Rf(1,:))
    Th_1(i) = acos(dot_1(i))
    Th_2(i) = acos(dot_2(i))
end
for i = 2:1:N2    
    if Cross_1(3,3) > 0 then
        if Cross_1(i,3) < 0 | i == N2 then
            Th_1(i) = 2*%pi - Th_1(i)
        end
    else
        if Cross_1(i,3) > 0 | i == N2 then
            Th_1(i) = 2*%pi - Th_1(i)
        end
    end
    if Cross_2(3,3) > 0 then
        if Cross_2(i,3) < 0 | i == N2 then
            Th_2(i) = 2*%pi - Th_2(i)
        end
    else
        if Cross_2(i,3) > 0 | i == N2 then
            Th_2(i) = 2*%pi - Th_2(i)
        end
    end      
end
for i = 1:1:N2
    A_TrigC(i) = 0.5*((mRc(i,:)*mRc(1,:)*sin(Th_1(i)-Th_1(1))))
    A_sector(i) = abs(AreaCal(sa,sb,Th_1(i))-AreaCal(sa,sb,Th_1(1)))
    A_chord(i) = A_sector(i)-A_TrigC(i)                        
    A_TrigF(i) = 0.5*((mRf(i,:)*mRf(1,:)*sin(Th_2(i)-Th_2(1))))
    A_focus(i) = A_TrigF(i)+ A_chord(i)
    Area_3(i) = A_focus(i)
end

Areal_v1 = diff(Area_1)/dt
Areal_v2 = diff(Area_2)/dt
Areal_v3 = diff(Area_3)/dt
Av1 = mean(Areal_v1)
Av2 = mean(Areal_v2)
Av3 = mean(Areal_v3)

Av1_check = diff(Areal_v1)
Av2_check = diff(Areal_v2)
Av3_check = diff(Areal_v3)
for i = 1:1:N2-2
    if abs(Av1_check(i)) > 1D-3 then
        disp("Warning: all {Areal_v1} values are not same ! Av1_check("+string(i)+") = "+string(Av1_check(i)))
    end
    if abs(Av2_check(i)) > 1D-3 then
        disp("Warning: all {Areal_v2} values are not same ! Av2_check("+string(i)+") = "+string(Av2_check(i)))
    end 
    if abs(Av3_check(i)) > 1D-3 then
        disp("Warning: all {Areal_v3} values are not same ! Av3_check("+string(i)+") = "+string(Av3_check(i)))
    end     
end

//Kepler's 3rd Law:
function [Ratio_1, Ratio_2] = k_3rd_Law(TimePeriod, sa, mu)
    Ratio_1 = (TimePeriod^2)/(sa^3)
    Ratio_2 = (mu*TimePeriod^2)/(sa^3)
endfunction
[Ratio_1, Ratio_2] = k_3rd_Law(Time_Period, sa, u)
//------------------------------------------------------------//

//------------Export Results-----------------//
cd("C:\Users\user\Downloads")
WriteFile = fullfile("C:\Users\user\Downloads\Orbit_Parameters_.csv")
csvWrite(["a", "b", "e", "T", "h", "\mu"; string(sa), string(sb), string(ec), string(Time_Period), string(h), string(u)], WriteFile)
//-------------------------------------------//

//Display Results:
disp("Semi-major Axis, a = "+string(sa))
disp("Eccentricity, e = "+string(ec))
disp("Time Period, T = "+string(Time_Period))
disp("Gravitational Parameter, u = "+string(u))
disp("Areal Velocity 1, Av1 = "+string(Av1))
disp("Areal Velocity 2, Av2 = "+string(Av2))
disp("Areal Velocity 3, Av3 = "+string(Av3))
disp("T^2/a^3 = "+string(Ratio_1))
disp("u*T^2/a^3 = "+string(Ratio_2))


//-------------------------FIGURES------------------------------//
xdel(winsid())
cd("C:\Users\user\Downloads")
//A static sketch containing the 3D points at different times, in space (3D figure).
second = 3
comet3d(xT, yT, zT)
fig = gca()
fig.box = "on"
fig.cube_scaling = "on"
fig.rotation_angles = [35,45]
fig.children(1).mark_size = 15
fig.children(1).mark_background = -2
fig.children(2).mark_mode = "on"
fig.children(2).mark_style = 0
fig.children(2).mark_size_unit = "point"
fig.children(2).mark_size = 4
fig.children(2).mark_foreground = 5
fig.children(2).mark_background = 5
fig.data_bounds=[-8,-9.5,-2;6,6,0.01];
fig.sub_ticks = [0,0,0]
sleep(second,"s")

clf
scatter3(xT, yT, zT,"markerEdgeColor", "red", "markerFaceColor", "red")
fig = gca()
fig.box = "on"
fig.cube_scaling = "on"
fig.rotation_angles = [35,45]
fig.children.mark_style = 9
fig.children.mark_size_unit = "point"
fig.children.mark_size = 4
fig.data_bounds=[-8,-9.5,-2;6,6,0.01];
fig.sub_ticks = [0,0,0]
xgrid(8,0,10)
sleep(second,"s")

// Distribution of angles of deviation from best fitted plane.
//(Indicating Kepler's law not strictly valid, but best fitted plane is a good representation).
clf
k = length(Angle)
plot(1:k, Angle', "black.")
plot(0:k+4, zeros(k+5,1)', "r")
xlabel(" k ")
ylabel("Angle (in degrees)")
sleep(second,"s")

clf
k = linspace(0,1,5)
for i = 1:1:N-2
    nP1 = N_Matrix(i,:)
    nP2 = n2'
    comet3d(nP1(1)*k, nP1(2)*k, nP1(3)*k, "colors", 2)    
end
comet3d(nP2(1)*k, nP2(2)*k, nP2(3)*k, "colors", 5)
fig = gca()
fig.cube_scaling = "on"
fig.rotation_angles = [64.75,378.25]
fig.sub_ticks = [0,0,0]

//Distribution of original data points about the best fitted plane identified through PCA/SVD.
clf
[xk, yk] = meshgrid(-10:1:10)
zk = -(n2(1).*(xk-xm) + n2(2).*(yk-ym))/n2(3) + zm
surf(xk,yk,zk,'facecol','direct','edgecol','lightblue')

scatter3(xT, yT, zT,"markerEdgeColor", "red", "markerFaceColor", "red")
fig = gca()
fig.box = "back_half"
fig.children(1).mark_style = 9
fig.children(1).mark_size_unit = "point"
fig.children(1).mark_size = 4
fig.children(2).mark_mode = "off"
fig.rotation_angles = [79.5,142.75]
xgrid(8,0,10)
sleep(second,"s")

//A: Projection of original data points on this best fitted plane.
//B: A plot of the results, having the actual 3D points, and their images on the plane indicated one through
// a solid image and other shaded differently.
clf
[xk, yk] = meshgrid(-10:1:10)
zk = -(n2(1).*(xk-xm) + n2(2).*(yk-ym))/n2(3) + zm
surf(xk,yk,zk,'facecol','direct','edgecol','lightblue')

xImage = imP_SVD(:,1)
yImage = imP_SVD(:,2)
zImage = imP_SVD(:,3)
param3d(xImage, yImage, zImage)
param3d(xT, yT, zT)
scatter3(xImage, yImage, zImage,"markerEdgeColor", "blue", "markerFaceColor", "blue")
scatter3(xT, yT, zT,"markerEdgeColor", "red", "markerFaceColor", "red")
fig = gca()
fig.box = "back_half"
fig.children(1).mark_style = 9
fig.children(1).mark_size_unit = "point"
fig.children(1).mark_size = 4
fig.children(2).mark_mode = "on"
fig.children(2).mark_style = 9
fig.children(2).mark_size_unit = "point"
fig.children(2).mark_size = 2
fig.rotation_angles = [83,106.75]
xgrid(8,0,10)
sleep(second,"s")

// Figure showing the original (x; y; z), and (xp; yp; zp = 0) coordinate systems 
clf
[xplot, yplot] = meshgrid(-7.75:1:7.75);
for i = 1:1:size(xplot,1)
    for j = 1:1:size(xplot,2)
        zplot(i,j) = c1 + c2/sqrt(xplot(i,j)^2 + yplot(i,j)^2);
    end    
end
surf(xplot, yplot, zplot,'facecol','white','edgecol','lightblue')
fig = gca()
fig.box = "back_half"
comet3d(xImage, yImage, zImage)
comet3d(xT, yT, zT, "colors",5)
sleep(second,"s")

// Figure showing ellipse fitted to projection data. 
//Showing that points are scattered about the best fitting ellipse.
clf
xRotate = x*cos(theta) - y*sin(theta)
yRotate = x*sin(theta) + y*cos(theta)   
plot(xRotate, yRotate,"red")
plot(xnew,ynew,"b.") 
legend("Best Fitted Ellipse","Projected data")
sleep(second,"s")

//A) Principal axis transformation, and
//B) Ellipse parameters marked out.
clf
plot(x,yc*ones(N2,1),'green')
plot(xc*ones(N2,1),y,'green')
plot(x,y,"red")
plot(xnew*cos(theta)+ynew*sin(theta), -xnew*sin(theta)+ynew*cos(theta),"b.")
plot(xc, yc, "o-m")
fig = gce()
fig.children(1).mark_style = 9
fig.children(1).mark_size = 8
fig.children(1).mark_foreground = -1
fig.children(1).mark_background = 6
plot(xf, yf, "o-y")
fig = gce()
fig.children(1).mark_style = 9
fig.children(1).mark_size = 8
fig.children(1).mark_foreground = -1
fig.children(1).mark_background = 7
plot(xPg, yPg, "o-g")
fig = gce()
fig.children(1).mark_style = 9
fig.children(1).mark_size = 8
fig.children(1).mark_foreground = -1
fig.children(1).mark_background = 3
plot(xAg, yAg, "o-r")
fig = gce()
fig.children(1).mark_style = 9
fig.children(1).mark_size = 8
fig.children(1).mark_foreground = -1
fig.children(1).mark_background = 5
legends(["Center","Focus","Perigee","Apogee"],[-1,-1,-1,-1], opt="ur")
fig = gce()
fig.children(2).mark_style = 9
fig.children(2).mark_foreground = -1
fig.children(2).mark_background = 5
fig.children(4).mark_style = 9
fig.children(4).mark_foreground = -1
fig.children(4).mark_background = 3
fig.children(6).mark_style = 9
fig.children(6).mark_foreground = -1
fig.children(6).mark_background = 7
fig.children(8).mark_style = 9
fig.children(8).mark_foreground = -1
fig.children(8).mark_background = 6

// LRL vector conservation
clf
plot(x,yc*ones(N2,1),'green')
plot(xc*ones(N2,1),y,'green')
plot(x,y,"red")
plot(xnew*cos(theta)+ynew*sin(theta), -xnew*sin(theta)+ynew*cos(theta),"b.")
plot(xc, yc, "o-m")
fig = gce()
fig.children(1).mark_style = 9
fig.children(1).mark_size = 8
fig.children(1).mark_foreground = -1
fig.children(1).mark_background = 6
plot(xf, yf, "o-y")
fig = gce()
fig.children(1).mark_style = 9
fig.children(1).mark_size = 8
fig.children(1).mark_foreground = -1
fig.children(1).mark_background = 7

[i(1),j(1)]= max(real(xmid))
[i(2),j(2)]= max(real(ymid))
[i(3),j(3)]= min(real(xmid))
[i(4),j(4)]= min(real(ymid))
interval = 125
j(5) = interval; j(6) = 5*interval
j(7) = 11*interval; j(8) = 15*interval

if a1 > b1 then
    for k = 1:1:length(j)
        plot([xmid(j(k)) xmid(j(k))+1*(LRL(j(k),1)/max(real(LRL(j(k),1))))],[ymid(j(k)) ymid(j(k))+LRL(j(k),2)],'->black')
    end    
else
    for k = 1:1:length(j)
        plot([xmid(j(k)) xmid(j(k))+LRL(j(k),1)],[ymid(j(k)) ymid(j(k))+1*(LRL(j(k),2)/max(real(LRL(j(k),2))))],'-^black')
    end        
end

// equal areal velocity sketch.
clf
plot(x,y,"red")
k=linspace(0,1,5)
for i = 1:200:2001
    plot([xf x(i)],[yf y(i)],'black')
end
plot(xPg, yPg, "o-g")
plot(xAg, yAg, "o-r")
scatter([xf xf],[yf yf],600,"yellow","fill")
//------------------------------------END--------------------------------------------
