% OBS: Recall your attention to change the extension from 
% flexible_rotor_modal_DTU_c.txt to flexible_rotor_modal_DTU_c.m
% before run Matlab !


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MACHINERY DYNAMICS LECTURES  (41514)                  %
% MEK - DEPARTMENT OF MECHANICAL ENGINEERING            %
% DTU - TECHNICAL UNIVERSITY OF DENMARK                 %
%                                                       %
%              Copenhagen, March 30th, 2017             %
%                                                       %
%                        Ilmar Ferreira Santos          %
%                                                       %
% ROTATING MACHINES -- NATURAL FREQUENCIES AND MODES    %
%                                                       %  
% EXPERIMENTAL RESULTS                                  %
% 13.0 (horizontal)                                     %
% 14.9 (vertical)                                       %
% 33.6 (horizontal)                                     %
% 43.0 (horizontal)                                     %
% 46.0 (vertical)                                       %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clear all;
 close all;
 for NEL = 14:2:80

 folder = 'Free_disc_disc';
 dataAquisition = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEFINITION OF THE STRUCTURE OF THE MODEL   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 NE=NEL;         % number of shaft elements (Minimum 14)
 GL = (NE+1)*4; % number of degree of freedom
 ND=2;          % number of discs
 NM=2;          % number of bearings
 
 penalty = 1.05;
 
 [l,Rext,Rint,nodePos] = DivideShaftElements(NE);
 
 %Disc position
 CD1=(65+40)/1000;            % node - disc 1 
 CD2=(CD1 + 90/1000);         % node - disc 2
 
 %Define span in which discs are positioned
 CD1pos = [65,105]./1000;
 CD2pos = [145,245]./1000;
 
 
 discNode = MountOnNode([CD1,CD2],nodePos)
 CD1 = discNode(1);
 CD2 = discNode(2);

 
 %Bearing Position 
 CMM1=(25+40+(1010-835)+(835-780)+(780-660)/2)/1000;    % node - bearing 1
 CMM2=(1150-52-23/2)/1000;                              % node - bearing 2
 bearingNode = MountOnNode([CMM1,CMM2],nodePos);
 CMM1 = bearingNode(1);
 CMM2 = bearingNode(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CONSTANTS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 E    = 2.0E11;  % {elasticity modulus [N/m^2}
 RAco = 7800;    % {steel density [kg/m^3]}
 RAl  = 7800;    % {aluminum density [kg/m^3]}
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         OPERATIONAL CONDITIONS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Omega=  0*2*pi;  % angular velocity [rad/s]
  Omegarpm = Omega*60/2/pi; 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   GEOMETRY OF THE ROTATING MACHINE           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%(A) DISCS    %
%%%%%%%%%%%%%%%
%%% DISC 1 %%%
 Rd1   = (295/2)/1000;                         % disc radius [m]
 espD1 = 80/1000;                              % disc thickness  [m]
 MasD1 = pi*Rd1^2*espD1*RAl;                   % disc mass [kg]
 Id1   = 1/4*MasD1*Rd1^2+1/12*MasD1*espD1^2;   % transversal mass moment of inertia of the disc [Kgm^2]
 Ip1   = 1/2*MasD1*Rd1*Rd1;                    % polar mass moment of inertia of the disc [Kgm^2]
 
%%% Disc 2 %%%
 Rd2   = (295/2)/1000;                         % disc radius [m]
 espD2= 100/1000;                              % disc thickness  [m]
 MasD2 = pi*Rd2^2*espD2*RAl;                 % disc mass [kg]
 Id2   = 1/4*MasD2*Rd2^2+1/12*MasD2*espD2^2;   % transversal mass moment of inertia of the disc [Kgm^2]
 Ip2   = 1/2*MasD2*Rd2*Rd2;                    % polar mass moment of inertia of the disc [Kgm^2]
 
 %%%%%%%%%%%%%%%
 %(B) BEARINGS %
 %%%%%%%%%%%%%%%
 
 MasM =   0*0.40698;                      % bearing mass [kg](housing + ball bearings)
 h=1/1000;                              % beam thickness [m]
 b=28.5/1000;                           % beam width [m]
 Area=b*h;                              % beam cross section area [m^2] 
 I=b*h^3/12;                            % beam moment of inertia of area [m^4]  
 lr=7.5/100;                            % beam length [m]
 Kty0=2*12*E*I/(lr^3)*0.000001;                  % equivalent beam flexural stiffness [N/m]
 Ktz0=2*E*Area/lr*0.0000000001;                      % equivalent bar stiffness [N/m]
% Bearing 1 - Damping 
 Dty1 = 0.0 ;                 
 Dtz1 = 0.0 ;                 
 Dry1 = 0.0 ;                
 Drz1 = 0.0 ;                  
% Bearing 2 - Damping
 Dty2 = 0.0 ;                 
 Dtz2 = 0.0 ;                 
 Dry2 = 0.0 ;                 
 Drz2 = 0.0 ;  
% Bearing 1 - Stiffness
 Kty1 = Kty0 ;        
 Ktz1 = Ktz0 ;            
 Kry1 = 0.0 ;                 
 Krz1 = 0.0 ;                 
% Bearing 2 - Stiffness
 Kty2 = Kty0 ;                
 Ktz2 = Ktz0 ;                
 Kry2 = 0.0 ;                 
 Krz2 = 0.0 ;  
 
 %%%%%%%%%%%%%%%
 %(C) SHAFT    %
 %%%%%%%%%%%%%%%
 


% external radius of shaft elements [m]
    for i=1:NE,
     rx(i)=Rext(i);
    end

% internal radius of shaft elements [m]
    for i=1:NE,
     ri(i)=Rint(i);
    end

% density of shaft elements [kg/m]
    for i=1:NE,
     ro(i) = RAco;
    end

% transversal area of the shaft elements [m^2]}
    for i=1:NE,
     St(i) = pi*(rx(i)^2-ri(i)^2);
    end

% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE,
     II(i)=pi*(rx(i)^4-ri(i)^4)/4;
    end
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MOUNTING THE GLOBAL MATRICES           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 disp('MOUNTING THE GLOBAL MATRICES - WAIT!')
 disp(' ')

% Defining the global matrices with zero elements
   M=zeros(GL);
   G=zeros(GL);
   K=zeros(GL);   
   
   
%%%%%%%%%%%%%%%%%%%%%%   
% GLOBAL MASS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%

disp('MOUNTING THE GLOBAL MASS MATRIX - WAIT!')
disp(' ')

% Mass matrices of shaft elements due to linear and angular movements

a=1; b=8;

for n=1:NE,   

  MteAux= [156       0         0          22*l(n)    54        0         0          -13*l(n)
           0         156       -22*l(n)   0          0         54        13*l(n)    0
           0         -22*l(n)  4*l(n)^2   0          0         -13*l(n)  -3*l(n)^2  0
           22*l(n)   0         0          4*l(n)^2   13*l(n)   0         0          -3*l(n)^2
           54        0         0          13*l(n)    156       0         0          -22*l(n)
           0         54        -13*l(n)   0          0         156       22*l(n)    0
           0         13*l(n)   -3*l(n)^2  0          0         22*l(n)   4*l(n)^2   0
           -13*l(n)  0         0          -3*l(n)^2  -22*l(n)  0         0          4*l(n)^2];


  Mte = ((ro(n)*St(n)*l(n))/420)*MteAux;
  

  MreAux= [36      0        0         3*l(n)    -36      0       0         3*l(n)
           0       36       -3*l(n)   0         0        -36     -3*l(n)   0
           0       -3*l(n)  4*l(n)^2  0         0        3*l(n)  -l(n)^2   0
           3*l(n)  0        0         4*l(n)^2  -3*l(n)  0       0         -l(n)^2
           -36     0        0         -3*l(n)   36       0       0         -3*l(n)
           0       -36      3*l(n)    0         0        36      3*l(n)    0
           0       -3*l(n)  -l(n)^2   0         0        3*l(n)  4*l(n)^2  0
           3*l(n)  0        0         -l(n)^2   -3*l(n)  0       0         4*l(n)^2];

  Mre = ((ro(n)*St(n)*(rx(n)^2-ri(n)^2))/(120*l(n)))*MreAux;
 
  MauxT=Mte+Mre;

   for f=a:b,
    for g=a:b,
     M(f,g)=M(f,g)+MauxT(f-(n-1)*4,g-(n-1)*4);
    end
   end

a=a+4; b=b+4;

end

% Adding the mass matrices of the disc elements 

   M((CD1-1)*4+1,(CD1-1)*4+1)=M((CD1-1)*4+1,(CD1-1)*4+1)+MasD1;
   M((CD1-1)*4+2,(CD1-1)*4+2)=M((CD1-1)*4+2,(CD1-1)*4+2)+MasD1;
   M((CD1-1)*4+3,(CD1-1)*4+3)=M((CD1-1)*4+3,(CD1-1)*4+3)+Id1;
   M((CD1-1)*4+4,(CD1-1)*4+4)=M((CD1-1)*4+4,(CD1-1)*4+4)+Id1;
   
   M((CD2-1)*4+1,(CD2-1)*4+1)=M((CD2-1)*4+1,(CD2-1)*4+1)+MasD2;
   M((CD2-1)*4+2,(CD2-1)*4+2)=M((CD2-1)*4+2,(CD2-1)*4+2)+MasD2;
   M((CD2-1)*4+3,(CD2-1)*4+3)=M((CD2-1)*4+3,(CD2-1)*4+3)+Id2;
   M((CD2-1)*4+4,(CD2-1)*4+4)=M((CD2-1)*4+4,(CD2-1)*4+4)+Id2;

% Adding the mass matrices of the bearing elements 

   M((CMM1-1)*4+1,(CMM1-1)*4+1)=M((CMM1-1)*4+1,(CMM1-1)*4+1)+MasM;
   M((CMM1-1)*4+2,(CMM1-1)*4+2)=M((CMM1-1)*4+2,(CMM1-1)*4+2)+MasM;
   
   M((CMM2-1)*4+1,(CMM2-1)*4+1)=M((CMM2-1)*4+1,(CMM2-1)*4+1)+MasM;
   M((CMM2-1)*4+2,(CMM2-1)*4+2)=M((CMM2-1)*4+2,(CMM2-1)*4+2)+MasM;
   
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL GYROSCOPIC MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('MOUNTING THE GLOBAL GYROSCOPIC MATRIX - WAIT!')
disp(' ')
           
% Gyroscopic matrix of shaft elements

a=1; b=8;

for n=1:NE,

   GeAux=[0        -36       3*l(n)   0          0        36       3*l(n)    0
          36       0        0         3*l(n)     -36      0        0         3*l(n)
          -3*l(n)  0        0         -4*l(n)^2  3*l(n)   0        0         l(n)^2
          0        -3*l(n)  4*l(n)^2  0          0        3*l(n)   -l(n)^2   0
          0        36       -3*l(n)   0          0        -36      -3*l(n)   0
          -36      0        0         -3*l(n)    36       0        0         -3*l(n)
          -3*l(n)  0        0         l(n)^2     3*l(n)   0        0         -4*l(n)^2
          0        -3*l(n)  -l(n)^2   0          0        3*l(n)   4*l(n)^2  0        ];

   Ge = 2*Omega*((ro(n)*St(n)*(rx(n)^2+ri(n)^2))/(120*l(n)))*GeAux;
   
   for f=a:b,
    for g=a:b,
     G(f,g)=G(f,g)+Ge(f-(n-1)*4,g-(n-1)*4);
    end
   end

a=a+4; b=b+4;

end

% Adding the gyroscopic matrices of the disc elements

   G((CD1-1)*4+3,(CD1-1)*4+4)=G((CD1-1)*4+3,(CD1-1)*4+4)-Omega*Ip1;
   G((CD1-1)*4+4,(CD1-1)*4+3)=G((CD1-1)*4+4,(CD1-1)*4+3)+Omega*Ip1;
   
   G((CD2-1)*4+3,(CD2-1)*4+4)=G((CD2-1)*4+3,(CD2-1)*4+4)-Omega*Ip2;
   G((CD2-1)*4+4,(CD2-1)*4+3)=G((CD2-1)*4+4,(CD2-1)*4+3)+Omega*Ip2;
   
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL STIFFNESS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('MOUNTING THE GLOBAL STIFFNESS MATRIX - WAIT!')
disp(' ')
           
% Stiffness matrix of shaft elements due to bending

a=1; b=8;

for n=1:NE,
    nodes = [n,n+1];
  KbeAux= [12      0        0         6*l(n)    -12      0       0         6*l(n)
           0       12       -6*l(n)   0         0        -12     -6*l(n)   0
           0       -6*l(n)  4*l(n)^2  0         0        6*l(n)  2*l(n)^2  0
           6*l(n)  0        0         4*l(n)^2  -6*l(n)  0       0         2*l(n)^2
           -12     0        0         -6*l(n)   12       0       0         -6*l(n)
           0       -12      6*l(n)    0         0        12      6*l(n)    0
           0       -6*l(n)  2*l(n)^2  0         0        6*l(n)  4*l(n)^2  0
           6*l(n)  0        0         2*l(n)^2  -6*l(n)  0       0         4*l(n)^2];
       if ND == 1
           if all([(nodePos(nodes(1))>=CD2pos(1)), (nodePos(nodes(1))<=CD2pos(2)),(nodePos(nodes(2))>=CD2pos(1)), (nodePos(nodes(2))<=CD2pos(2))])
               Kbe = ((penalty*E*II(n))/(l(n)^3))*KbeAux;
               disp('Element inside bearing')
               
           else
               Kbe = ((E*II(n))/(l(n)^3))*KbeAux;
           end
           
       elseif ND == 2
           if all([(nodePos(nodes(1))>=CD1pos(1)), (nodePos(nodes(1))<=CD1pos(2)),(nodePos(nodes(2))>=CD1pos(1)), (nodePos(nodes(2))<=CD1pos(2))])
               Kbe = ((penalty*E*II(n))/(l(n)^3))*KbeAux;
               disp('Element inside bearing 1')
               
           else
               Kbe = ((E*II(n))/(l(n)^3))*KbeAux;
           end
           
           if all([(nodePos(nodes(1))>=CD2pos(1)), (nodePos(nodes(1))<=CD2pos(2)),(nodePos(nodes(2))>=CD2pos(1)), (nodePos(nodes(2))<=CD2pos(2))])
               Kbe = ((penalty*E*II(n))/(l(n)^3))*KbeAux;
               disp('Element inside bearing 2')
               
           else
               Kbe = ((E*II(n))/(l(n)^3))*KbeAux;
           end
       end

  for f=a:b,
   for g=a:b,
    K(f,g)=K(f,g)+Kbe(f-(n-1)*4,g-(n-1)*4);
   end
  end

a=a+4; b=b+4;

end

% Adding the stiffness matrices of the bearing elements

   K((CMM1-1)*4+1,(CMM1-1)*4+1)=K((CMM1-1)*4+1,(CMM1-1)*4+1)+Ktz1;
   K((CMM1-1)*4+2,(CMM1-1)*4+2)=K((CMM1-1)*4+2,(CMM1-1)*4+2)+Kty1;
   
   K((CMM2-1)*4+1,(CMM2-1)*4+1)=K((CMM2-1)*4+1,(CMM2-1)*4+1)+Ktz2;
   K((CMM2-1)*4+2,(CMM2-1)*4+2)=K((CMM2-1)*4+2,(CMM2-1)*4+2)+Kty2;



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %    GLOBAL MATHEMATICAL MODEL                 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
  Mglob=[   M                zeros(size(M,1)) 
            zeros(size(M,1)) K             ];

  Kglob=[   G                K 
           -K                zeros(size(M,1))];
      
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %              MODAL ANALYSIS                  %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 disp('CALCULATING NATURAL FREQUENCIES AND MODE SHAPES - WAIT!')
 disp(' ')

 % Calculating Eigenvectors and Eigenvalues

 [U,lambda]=eig(-Kglob,Mglob);
 [lam,p]=sort(diag(lambda));
 U=U(:,p);
 disp('Frequency in rad/s:')
 imag(lam(10:4:28))
 disp('Frequency in Hz:')
 eigFreq = imag(lam(10:4:28))./(2*pi)

 if dataAquisition == 1
     name = strcat(folder,'/',strcat(folder,'_',sprintf('%d.mat',NE)));
     save(name,'eigFreq');   
 end
 end
 return

 % Number of divisions in time for plotting the mode shapes
  nn=99;

  N=size(U,1);
  maximo=num2str((N-2)/2);
  ModoVirt=N;

  ModoVirt=input([' Enter the number of the mode shape to be plotted, zero to esc, highest mode ',maximo,': ']);

  %For visualizing the mode shapes:

  ModoReal=2*ModoVirt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        LOOP TO PLOT THE MODES SHAPES            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 while ModoReal>0,

  % Natural frequencies
  wn=imag(lam(ModoReal));
  ttotal=8/abs(wn);
  dt=ttotal/nn;
  t=0:dt:ttotal;

  % Defining v and w real e imaginary for each node

  y=1:4:GL;
  z=2:4:GL;

  for i=1:(NE+1),
   vr(i)=real(U(y(i),ModoReal));
   vi(i)=imag(U(y(i),ModoReal));
   wr(i)=real(U(z(i),ModoReal));
   wi(i)=imag(U(z(i),ModoReal));
  end

  % Calculation the modal displacements v and w

  for i=1:(NE+1),
   v(i,:)=vr(i)*cos(wn*t)+vi(i)*sin(wn*t);
   w(i,:)=wr(i)*cos(wn*t)+wi(i)*sin(wn*t);
  end

  Zero=diag(zeros(length(t)))';
  Um=diag(eye(length(t)))';

  for i=1:(NE+1)
   pos(i,:)=Zero+(i-1)*Um;
  end
  
  clf
  hold on
  for cont=1:NE+1,
   plot3(pos(cont,:),w(cont,:),v(cont,:),'k','LineWidth',2.5);
  end

  nm=num2str(ModoVirt);
  fn=num2str(abs(wn/2/pi));
  dfi=num2str(Omegarpm);
  title(['Angular Velocity: ',dfi,' rpm       Mode: ',nm,'   Nat. Freq.:  ',fn,' Hz'],'FontSize',14)
  view(-25,20);
  grid on;
  ModoVirt=input([' Enter the number of the mode shape to be plotted, zero to esc, highest mode ',maximo,': ']);
  ModoReal=2*ModoVirt;
  figure(ModoVirt)
 
end