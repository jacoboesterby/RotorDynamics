function [l,Rext,Rint,nodePos] = DivideShaftElements(NE)
if NE<14
    error('To less elements defined. Minimum 14 elements required')
end
    
 % Length of initial the shaft sections [mm]
 l(1)  = 15;
 l(2)  = 10;
 l(3)  = 40;
 l(4) = 40;
 l(5) = (40+50);
 l(6) = ((1010-835)-(80+50));
 l(7)  = 835-780;
 l(8) = (780 - 660)/2;
 l(9) = (780 - 660)/2;
 l(10)  = 660-35;
 l(11)  = 35;
 l(12)  = 23/2;
 l(13) = 23/2;
 l(14)  = 52;
 %Convert to [m]
 l = l./1000;
 
% Inner and outer radius of initial sections
 Rextsec(1) = 30/1000;       % shaft external radius [m]
 Rintsec(1) = 7/1000;        % shaft internal radius [m]  
 
 Rextsec(2) = (62.5/2)/1000;
 Rintsec(2) = 0/1000;
 
 Rextsec(3) = (65/2)/1000;
 Rintsec(3) = 0/1000;
 
 Rextsec(4) = 35/1000;
 Rintsec(4) = 0/1000;
 
 Rextsec(5) = 35/1000;
 Rintsec(5) = 0/1000;
 
 Rextsec(6) = 35/1000;
 Rintsec(6) = 0/1000;
 
 Rextsec(7) = 45/1000;
 Rintsec(7) = 0/1000;
 
 Rextsec(8) = (99.6/2)/1000;
 Rintsec(8) = 0/1000;
 
 Rextsec(9) = (99.6/2)/1000;
 Rintsec(9) = 0/1000;
 
 Rextsec(10) = 45/1000;
 Rintsec(10) = 0/1000;
 
 Rextsec(11) = 35/1000;
 Rintsec(11) = 0/1000;
 
 Rextsec(12) = 25/1000;
 Rintsec(12) = 0/1000;
 
 Rextsec(13) = 25/1000;
 Rintsec(13) = 0/1000;
 
 Rextsec(14) = 20/1000;
 Rintsec(14) = 7/1000;


 
 %Length of each section [mm]
 sec(1) = 15;
 sec(2) = 10;
 sec(3) = 40;
 sec(4) = 40;
 sec(5) = (40+50);
 sec(6) = ((1010-835)-(80+50));
 sec(7) = 835-780;
 sec(8) = (780 - 660)/2;
 sec(9) = (780 - 660)/2;
 sec(10) = 660-35;
 sec(11) = 35;
 sec(12) = 23/2;
 sec(13) = 23/2;
 sec(14) = 52;
 sec = sec./1000;
 
 for k=1:length(sec)
     if k>1
         accLength(k) = accLength(k-1)+sec(k);
     else
         accLength(k) = sec(k);
     end
 end
 
 for k=1:(NE-14)
    ind = find(max(l)==l,1);
    if ind <length(l)
    l = [l(1:ind-1),l(ind)/2,l(ind)/2,l(ind+1:end)];
    else
        l = [l(1:ind-1),l(ind)/2,l(ind)/2];
    end
     
 end
 
 for k=1:length(l)
     if k==length(l)
         Rext(k) = Rextsec(end);
         Rint(k) = Rintsec(end);
     else
         ind = find(sum(l(1:k))<=accLength == 1,1)
         Rext(k) = Rextsec(ind);
         Rint(k) = Rintsec(ind);
     end
     
 end
     
 
  nodePos = zeros(1,length(l)+1);
 for k = 2:length(l)+1
     nodePos(k) = nodePos(k-1) + l(k-1)
 end
 