%%
% Script responsible for model the KPZ equation in 1D, using a
% correction function from https://doi.org/10.1103/PhysRevE.77.031134
%
% @author: Mariana P. M. A Baroni
% @last access: April 30, 2020
%
% Don't forget to cite it properly!

%****Sanitizing
clear all
clc
close all

%****Dimension
dt = 0.05;
dx = 1;
L = 512%64%128%256%512;
T = 10000;
x = 0:dx:L;
t = 0:dt:T;
lt = length(t);
lx = length(x);

%****Parameters
vj = 1;%0:1:10;
lj = 6.93;%0:.1:10;
c = 1;
D = 1;
d = 1;
r = sqrt( (2*D)/(dx^d) );

%****Random numbers
a=-1/2;b=1/2;
rand('seed',0);
ru=a+(b-a)*rand(lt,lx); %uniform distribution

for vi = 1:length(vj) %you could test a vector of parameters
    v=vj(vi);
    
    for li = 1:length(lj)
        lambda=lj(li);
        
        %****Initial condition
        h=zeros(1,lx);
        h0=zeros(1,lx);
        
        %****G
        g = (lambda*lambda*D)/(v*v*v);
        
        for n = 1:lt
            
            for i = 1:lx
                
                if i==1 %periodic boundary condition
                    
                    %discretization without correction function
                    %h(i) = h0(i) + (dt/(dx*dx))* (
                    %            v*( h0(lx) - 2*h0(i) + h0(i+1) ) +
                    %             (lambda/8)*( ( h0(i+1) - h0(lx) ).^2)
                    %             ) + r*sqrt(12*dt)*ru(n,i) ;
                    
                    %discretization with correction function
                    h(i) = h0(i) + dt*( (v/(dx*dx))*(h0(lx)-2*h0(i)+h0(i+1)) + lambda*(1+(1/2)*( (1-exp(-c*( ((h0(i+1)-h0(lx))/(2*dx)).^2 ) ))/c  ) )) + r*sqrt(12*dt)*ru(n,i) ;
                    
                elseif i==lx %periodic boundary condition
                    
                    %discretization without correction function
                    %h(i)=h0(i)+dt*( ((v/(dx*dx))*(h0(i-1)-2*h0(i)+h0(1))) + lambda*(1+(1/(4*dx))*(h0(1)-h0(i-1)).^2) ) +r*sqrt(12)*ru(n,i) ;
                    
                    %discretization with correction function
                    h(i)=h0(i) + dt*( (v/(dx*dx))*(h0(i-1)-2*h0(i)+h0(1)) + lambda*(1+(1/2)*( (1-exp(-c*( ((h0(1)-h0(i-1))/(2*dx)).^2 ) ))/c ) ) ) +r*sqrt(12*dt)*ru(n,i);
                    
                else
                    
                    %discretization without correction function
                    %h(i)=h0(i)+dt*( ((v/(dx*dx))*(h0(i-1)-2*h0(i)+h0(i+1))) + lambda*(1+(1/(4*dx))*(h0(i+1)-h0(i-1)).^2) )+r*sqrt(12)*ru(n,i) ;
                    
                    %discretization with correction function
                    h(i)=h0(i)+dt*( (v/(dx*dx))*(h0(i-1)-2*h0(i)+h0(i+1)) + lambda*(1+(1/2)*( (1-exp(-c*( ((h0(i+1)-h0(i-1))/(2*dx)).^2 ) ))/c  ) ) ) + r*sqrt(12*dt)*ru(n,i);
                
                end
                
                
            end
            
            h0=h;
            
            %*****Save the profile on a file
            %name=strcat('KPZprofile_v',num2str(v),'_la_',num2str(lambda),'.txt');
            %fid = fopen(name,'w');
            %fprintf(fid,'%12.8f\n',h);
            %fclose(fid);
            %clear pp name
            
            %***Interface width
            W(n) = sqrt( (1/L)*sum( (h - mean(h)).*(h - mean(h)) ) );
            
%             figure(1)
%             set(gca,'FontSize',12)
%             plot(x,h,'b','Linewidth',2)
%             title('KPZ 1D')
%             xlabel('x')
%             ylabel('h(x)')
%             hold on
%             axis([1 512 -1 max(h)])
        end
        
        figure(2)
        set(gca,'FontSize',12)
        plot(x,h,'b','Linewidth',2)
        title('profile KPZ 1D')
        xlabel('x')
        ylabel('h(x)')
        hold on
        axis([1 512 -1 max(h)])
        
        figure(3)
        set(gca,'FontSize',12)
        loglog(t,W,'b*','Linewidth',2)
        title('Interface Width')
        xlabel('t')
        ylabel('W(t)')
        
        %*****Finding the w_saturation
        n00 = find(t==1);
        nff = n;
        hm = mean(W(n00:nff));
        j = find(W>=(hm-0.01),1,'first');
        
        %*****Finging the growth exponent
        n0 = find(t == 1);%10;%ceil(size(t,2)*.01); 
        nf = j-10;%find(t == 100);%ceil(size(t,2)*.4); 
        [slope] = polyfit(log(t(n0:nf)), log(W(n0:nf)),1);
        
        figure(4)
        set(gca,'FontSize',12)
        loglog(t,W,'b*','Linewidth',2)
        title('Interface Width')
        xlabel('t')
        ylabel('W(t)')
        hold on
        loglog(t,hm*ones(1,size(W,2)),'-*r')
        text(t(n0)+2,W(j)+1,['w_{sat} = ',num2str((hm))]);
        %loglog(t(n00:nff),hm*ones(1,size(n00:nff,2)),'-r')
        hold on
        
        %y = polyval(slope,n0:nf);
        y = 10.^(slope(1)*log10(t) + slope(2)); %interpolation
       
        
        %*****Finding the tx
        txi = find(y>=(hm),1,'first');
               
         plot(t(n0:txi),y(n0:txi), '-*r');
        text(t(ceil(n0/2))+2,y(ceil(n0/2))+2,['\beta = ',num2str(slope(1))]);
        
        hold on
        loglog(t(txi)*ones(1,size(W,2)),W,'*-k')
        text(t(txi)+100,y(txi)+1,['t_x = ',num2str(t(txi))]);
       
    end
end

