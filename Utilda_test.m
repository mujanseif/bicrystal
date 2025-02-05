function [Ux,Uy,Uz]=Utilda_test(rn,links,gnl,NU,xnodes,dx,dy,dz)

C=[dx/2,dy/2,dz/2];
nodenum=size(gnl,1);                                        %Number of FE nodes
segnum=size(links,1);                                       %Number of dislocation segments
Utilda=zeros(nodenum,3);                                    %Preallocate the displacement
vertices=[0,0,0;                                            %Vertices of cuboid
          dx,0,0;
          dx,dy,0;
          dx,dy,dz;
          dx,0,dz;
          0,0,dz;
          0,dy,dz;
          0,dy,0];
faces=[1,2,3,8;                                             %Faces of cuboid as defined by vertices
       1,2,5,6;
       1,6,7,8;
       2,3,4,5;
       3,4,7,8;
       4,5,6,7];
normals=[1,0,0;                                             %Normalised surface normal vectors for all faces
         0,1,0;
         0,0,1];

for i=1:segnum
    
    if rn(links(i,1),4) == 6 && rn(links(i,2),4) == 67      %Ignore all segments connecting surface nodes to virtual nodes
        continue
    elseif rn(links(i,2),4) == 6 && rn(links(i,1),4) == 67
        continue
    end
    
   A=rn(links(i,1),1:3);                                   %Coordinates of first node connected to current segment
    B=rn(links(i,2),1:3);                                   %Coordinates of second node connected to current segment
    b=links(i,3:5);                                         %Burgers vector of AB
    
    if rn(links(i,1),4) == 67 ||  rn(links(i,2),4) == 67                     %Finds virtual segments
        out=1;
        surfptsA=[];
        surfptsB=[];
        
        for v=1:size(normals,1)
            lineA=[A,normals(v,:)];
            lineB=[B,normals(v,:)];
            surfptsA=[surfptsA;intersectLineMesh3d(lineA,vertices,faces)];
            surfptsB=[surfptsB;intersectLineMesh3d(lineB,vertices,faces)];
        end
        
        [~,Aindex]=min((surfptsA(:,1)-A(1)).^2+(surfptsA(:,2)-A(2)).^2+(surfptsA(:,3)-A(3)).^2);
        [~,Bindex]=min((surfptsB(:,1)-B(1)).^2+(surfptsB(:,2)-B(2)).^2+(surfptsB(:,3)-B(3)).^2);
        Aprime=surfptsA(Aindex,:);
        Bprime=surfptsB(Bindex,:);
%         AB=B-A;
%         AAprime=Aprime-A;
%         planenormal=cross(AB,AAprime)/norm(cross(AB,AAprime));
%         AC=C-A;
%         planecheck=dot(planenormal,AC/norm(AC));
%         planecheck=sqrt(planecheck.^2);
        Ctemp=A+0.5.*(Bprime-A);
%         if planecheck>eps
%             pyramid=[A;B;C;Aprime;Bprime];                        %Designates original virtual nodes, closure point and projected nodes as vertices
%             tess=convhulln(pyramid,{'QJ','Pp'});                  %Finds the convex hull of the designated vertices
%             tol=1e1;
%             checker=inhull(xnodes(:,1:3),pyramid,tess,tol);       %Flags FE nodes inside the convex hull
%         else
%             checker=zeros(size(xnodes,1),1);
%         end
%         Aprime=A-10e9*b;
%         Aprime2=A+10e9*b;
%         Bprime=B-10e9*b;
%         Bprime2=B+10e9*b;
%         
%         if norm(Aprime-A)>norm(Aprime2-A)
%             Aprime=Aprime2;
%         end
%         
%         if norm(Bprime-B)>norm(Bprime2-B)
%             Bprime=Bprime2;
%         end
        
    else
        out=0;
    end
    
    for j=1:nodenum
        nodepoint=xnodes((gnl(j)),1:3);                           %Coordinates of the FE node
          p=nodepoint';                                             %Column vector for use in displacement_et
          
          if out==0
%           Utilda(j,:)=Utilda(j,:)+displacement_et(p,A',B',C',b',NU);
              Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',NU)+displacement_et_plas(p,A',B',C',b');
%           elseif checker(gnl(j))==0
%               Utilda(j,:)=Utilda(j,:)+displacement_et_plas(p,A',B',C',b')+displacement_et_plas(p,Aprime',A',C',b')+displacement_et_plas(p,B',Bprime',C',b');
          else
%               Utilda(j,:)=Utilda(j,:)+solang_correction(p,A',B',Bprime',Aprime',C',b');
              Utilda(j,:)=Utilda(j,:)+displacement_et_plas(p,Aprime',Bprime',C',b')+displacement_et_plas(p,Aprime',A',Ctemp',b')+displacement_et_plas(p,A',B',Ctemp',b')+displacement_et_plas(p,B',Bprime',Ctemp',b')+displacement_et_plas(p,Bprime',Aprime',Ctemp',b');
          end
    end
    
    Ux=Utilda(:,1);                                                     %Organises outputs
    Uy=Utilda(:,2);
    Uz=Utilda(:,3);
end

end

function [u] = displacement_et_el(p,A,B,b,nu)

%Calculates only elastic components of displacement for segment AB using
%code from displacement_et

con1 = (1-2*nu)/(8*pi*(1-nu));
con2 = 1/(8*pi*(1-nu));


% R vectors (from P to nodes)
RA = A - p;
RB = B - p;

lamA = safenorm(RA);
lamB = safenorm(RB);

vecAB = B - A;

  
tAB = safenorm(vecAB);

    

% f = fab
f = fab(b, tAB, lamA, lamB, RA, RB);

% g = gab
g = gab(b, lamA, lamB);

% update displacement inc. solid angle
u = - con1.*f + con2.*g;

u = u';

end

function [u] = displacement_et_plas(p,A,B,C,b)
% Calculates the plastic (solid angle) components of displacement for
% triangle ABC using code from displacement_et


% R vectors (from P to nodes)
RA = A - p;
RB = B - p;
RC = C - p;

lamA = safenorm(RA);
lamB = safenorm(RB);
lamC = safenorm(RC);

vecAB = B - A;
vecBC = C - B;
  
tAB = safenorm(vecAB);
tBC = safenorm(vecBC);

% calculate slip plane normal
% vec_a = segs(1,:)/norm(segs(1,:));
%     vec_b = -segs(end,:)/norm(segs(end,:));

    n = cross(tAB,tBC);
    
omega  = solang(lamA, lamB, lamC,n);
% update displacement inc. solid angle
u =  - b*omega/(4*pi);

u = u';

end

function [ out ] = fab( b, tAB, lamA, lamB, RA, RB )
%calculates f vector

numerator = norm(RB)*(1 + dot(lamB,tAB));
denominator = norm(RA)*(1 + dot(lamA,tAB));

   % this can happen if field point is on the dislocation segment
    if abs(denominator) < eps
        logarithm  = 0;
    elseif abs(numerator) < eps
        logarithm = 0;
    else
        logarithm = real(log(numerator/denominator)); % for large step calculations, log will produce an imaginary 
    end    

out = cross(b,tAB)*logarithm;

end

function [ out ] = gab( b, lamA, lamB )
%calculates g vector

numerator = dot(b,cross(lamA,lamB)) * (lamA + lamB);
denominator = 1 + dot(lamA,lamB);

if abs(denominator) < eps

    out = [0; 0; 0];

else
    out = numerator/denominator;
end

end

function [ omega ] = solang( lamA, lamB, lamC, n)
%calculates solid angle due to triangular loop ABC
% 1) this may not be robust, check different formulae for the 
% spherical excess (wiki spherical geometry)
% 2) it could be a lot faster to calculate a,b,c in the main loop 
% for each pair of lambdas, then use a,b,c as input arguments 

a = acos(dot(lamB,lamC));
b = acos(dot(lamC,lamA));
c = acos(dot(lamA,lamB));

s = (a + b + c)/2;

svec = [s, s-a, s-b, s-c]/2;

temp = tan(svec(1))*tan(svec(2))*tan(svec(3))*tan(svec(4));

if temp < eps
    temp = 0;
end
omega = 4*real(atan(sqrt(temp)));

sgn = sign(dot(lamA,n));

omega = -sgn*omega;


%-------- check sign explicitly for c code ------
% if dot(lamA,n) > eps
%     sgn2 = 1;
% elseif dot(lamA,n) < -eps    
%     sgn2 = -1;
% else
%     sgn2 = 0;
% end
% if sgn2 ~= sgn
%     disp('sign error')
%    % pause
% end
%--------------------------

end

function unitR = safenorm(R)

normR = sqrt(R(1)*R(1) + R(2)*R(2)+R(3)*R(3));

if normR > eps
    unitR = R/normR;
else
    unitR = zeros(3,1);
end

end