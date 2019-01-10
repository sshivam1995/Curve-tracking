close all
dt=0.01; tf=160; t=0:dt:tf; tbar=0:dt:100;
N=size(t,2); Nbar=size(tbar,2); T=0.5; delta_T=0.005*T;
alpha=30;
r_g=zeros(3,2,N);

for counter=1:3
    if counter==1
    	vel_ini=(5/18)*15;
    elseif counter==2
        vel_ini=(5/18)*25;
    else
        vel_ini=(5/18)*35;
    end

udot=zeros(2,Nbar);

%% Car parameters 
lr=1.628; lf=1.218; m=1586.9; Iz=2315.3; Calpha=[-35000;-35000];      
accln_sat=50;
angle_sat=pi/2;

%% Z is xdot,ydot,phidot, X, Y, phi
Z=zeros(6,N);
%Z0=[25;0;0;0;0;1.13];
Z(:,1)=[vel_ini;0;0;-60;0;0];
%Z(:,1)=Z0;

%% define reference r(1) is x_pos, r(2) is y_pos;

V=zeros(2,N);M=N-1;

lim1=floor(120/(vel_ini*dt));
lim2=floor(pi*30/(vel_ini*dt));



for i=1:lim1
    V(:,i)=[vel_ini; 0];
    V(:,i+lim1+lim2)=[-vel_ini;0];
    V(:,i+2*(lim1+lim2))=[vel_ini; 0];
    V(:,i+3*(lim1+lim2))=[-vel_ini;0];
    V(:,i+4*(lim1+lim2))=[vel_ini; 0];
end

for i=1:lim2
    V(:,i+lim1)=[vel_ini*cos(pi*i/lim2);vel_ini*sin(pi*i/lim2)];
    V(:,i+1*(lim1+lim2)+lim1)=[-vel_ini*cos(pi*i/lim2);-vel_ini*sin(pi*i/lim2)];
    V(:,i+2*(lim1+lim2)+lim1)=[vel_ini*cos(pi*i/lim2);vel_ini*sin(pi*i/lim2)];
    V(:,i+3*(lim1+lim2)+lim1)=[-vel_ini*cos(pi*i/lim2);-vel_ini*sin(pi*i/lim2)];
    V(:,i+4*(lim1+lim2)+lim1)=[vel_ini*cos(pi*i/lim2);vel_ini*sin(pi*i/lim2)];
end



r0=[-60;0];
r=zeros(2,N);

r(:,1)=r0;

for i=2:N
    r(:,i)=r(:,i-1)+V(:,i)*dt;
end    


%% Define future reference, T time ahead
ref=r;
for i=1:N-T/dt
   ref(:,i)=r(:,i+T/dt); 
end    
for i=N-T/dt+1:N
    ref(:,i)=ref(:,N-T/dt);
end

%% u is delta_f (front wheel angle wrt body) and a (acceleration)
u=zeros(2,Nbar);
u(:,1)=[0;0];

% Track  vehicle motion
pos_tracker=zeros(2,N);
pos_tracker(:,1)=r0;

%% derivative of z at all time instances
zdot=zeros(6,Nbar);

g_u=zeros(2,Nbar);
g_u_prime=zeros(2,2,Nbar);
inv_matrix=zeros(2,2,Nbar);

%% Main time loop
for i=2:Nbar 
    if mod(i,2000)==0
        i
    end
    
    zdot(:,i)=dzdt(Z(:,i-1),u(:,i-1),lr,lf,m,Iz,Calpha);
  
%     for k=1:6  
%         if abs(zdot(k,i))>1
%             zdot(k,i)=zdot(k,i)/abs(zdot(k,i)); 
%         end
%     end    
     
    Z(:,i)=Z(:,i-1)+zdot(:,i)*dt;  
    pos_tracker(:,i)=[Z(4,i);Z(5,i)];
    
    %Future prediction
    [g_u(:,i),g_u_prime(:,:,i)]=g_rt(Z(:,i-1),u(:,i-1), T, delta_T, lr,lf,Calpha,m,Iz);
    
    if i==2
       g_u(:,1)=g_u(:,2); 
    end 
    
    r_g(counter,:,i)=ref(:,i)-g_u(:,i);
    
    inv_matrix(:,:,i)=inv(g_u_prime(:,:,i));
    
    u(:,i)=u(:,i-1)+alpha*(inv_matrix(:,:,i))*(ref(:,i)-g_u(:,i))*dt; 
    
    if i>N-T/(dt)
        %delta=u(:,N-T/(dt))/(T/(dt));
        %u(:,i)=u(:,i-1)-delta;
        u(:,i)=u(:,i-1);
    end    
    
    if abs(u(1,i))>angle_sat
        u(1,i)=u(1,i)*angle_sat/abs(u(1,i));
    end    
    
    if abs(u(2,i))>accln_sat
        u(2,i)=u(2,i)*accln_sat/abs(u(2,i));
    end
    
    udot(:,i-1)=(u(:,i)-u(:,i-1))/dt;
    
end

%% storing variables for graphing

lat_accl(counter,:)=u(2,:);
jerk(counter,:)=udot(2,:);


path_heading=mod(atan2(V(2,1:Nbar),V(1,1:Nbar)),2*pi);
actual_heading=mod(Z(6,1:Nbar),2*pi);

error_heading=mod(actual_heading-path_heading,2*pi);
for i=1:size(error_heading,2)
        if error_heading(i)>6
            error_heading(i)=error_heading(i)-2*pi;
        end 
end

head_error(counter,:)=180/pi*error_heading;


[close_point,lateral_error_norm(counter,:),arc]=distance2curve(r',pos_tracker');


error=pos_tracker-r;   
norm_error(counter,:)=vecnorm(error(:,1:Nbar));

position(counter,:,:)=pos_tracker(:,1:Nbar);

end


%% Plot of accln
figure (1)
plot(tbar,lat_accl(1,:),'LineWidth',1.5)
hold on
plot(tbar,lat_accl(2,:),'LineWidth',1.5)
plot(tbar,lat_accl(3,:),'LineWidth',1.5)

x1=xlabel('Time$~[s]$');
 y1=ylabel('$a_l~[m/s^2]$');
 set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
  leg1=legend('$15~km/h$','$25~km/h$','$35~km/h$');
 set(leg1,'Interpreter','latex')
 
 set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('accln','-dsvg','-r0')

% plot(tbar,lat_accl(1,:))
% hold on
% plot(tbar,lat_accl(2,:))
% plot(tbar,lat_accl(3,:))
% 
% xlabel('time (s)')
%  ylabel('Longitudinal acceleration (m/s^2)')
%  legend('15 km/h','25 km/h','35 km/h')
 
%% Plot of jerk
% figure (2)
% plot(tbar,jerk(1,:))
% hold on
% plot(tbar,jerk(2,:))
% plot(tbar,jerk(3,:))
% %title('Longitudinal jerk vs time');
% xlabel('time (s)')
%  ylabel('Jerk (m/s^3)')
%  legend('15 km/h','25 km/h','35 km/h')
 
%% PLot of tracking error norm
%figure (3)


% plot (tbar,norm_error(1,:));
% hold on
% plot (tbar,norm_error(2,:));
% plot (tbar,norm_error(3,:));
% %title('tracking error norm vs time');
% xlabel('time (s)')
%  ylabel('Total Error (m)')
%  legend('15 km/h','25 km/h','35 km/h')
 

%% plot of heading error 
figure (4)

plot(tbar,head_error(1,:),'LineWidth',1.5)
hold on
plot(tbar,head_error(2,:),'LineWidth',1.5)
plot(tbar,head_error(3,:),'LineWidth',1.5)
%title('Heading error vs time')
x1=xlabel('Time$~[s]$');
y1=ylabel('$\Delta\psi~[^\circ]$');
 set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')

 leg1=legend('$15~km/h$','$25~km/h$','$35~km/h$');
 set(leg1,'Interpreter','latex')
     set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('head_error','-dsvg','-r0')
     

%% PLot of lateral error
figure (5)

plot(tbar,lateral_error_norm(1,1:Nbar),'LineWidth',1.5)
hold on
plot(tbar,lateral_error_norm(2,1:Nbar),'LineWidth',1.5)
plot(tbar,lateral_error_norm(3,1:Nbar),'LineWidth',1.5)
%title ('Normal Error')
 x1=xlabel('Time$~[s]$');
 y1=ylabel('Lateral Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
 leg1=legend('$15~km/h$','$25~km/h$','$35~km/h$');
 set(leg1,'Interpreter','latex')
     
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('lat_error','-dsvg','-r0')

 %% plot of track
 figure (6)
 
 pos_15_x(1,:)=position(1,1,1:Nbar); pos_15_y(1,:)=position(1,2,1:Nbar);
 pos_25_x(1,:)=position(2,1,1:Nbar); pos_25_y(1,:)=position(2,2,1:Nbar);
 pos_35_x(1,:)=position(3,1,1:Nbar); pos_35_y(1,:)=position(3,2,1:Nbar);
 

 plot (r(1,1:Nbar),r(2,1:Nbar),'LineWidth',1.5);
 hold on

 plot(pos_15_x,pos_15_y,'LineWidth',1.5);
 plot(pos_25_x,pos_25_y,'LineWidth',1.5);
 plot(pos_35_x,pos_35_y,'LineWidth',1.5);

 ylim([-10 70])
 
x1=xlabel('$z_1~[m]$');
 y1=ylabel('$z_2~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
 leg1=legend('Reference','$15~km/h$','$25~km/h$','$35~km/h$');
 set(leg1,'Interpreter','latex')
 
 set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('path','-dsvg','-r0')


%title ('Tracking path')
%  xlim([-100 100]) 
%  ylim([-10 70])
 
% s_pos =[-87 9 -82 14];
% % location of the zoom-in plot  
% t_pos = [-40 10 25 40];    
% ah = gca; 
% % generate a zoom-in plot.  
% zoomPlot(ah, s_pos, t_pos); 
%  
% 
% s = plot(pos_tracker(1,1),pos_tracker(2,1),'o','MarkerFaceColor','black');      % bot 1
% q = plot(r(1,1),r(2,1),'o','MarkerFaceColor','green');                          % ref
% %r = 
% 
% for k = 2:Nbar
%     s.XData = pos_tracker(1,k);
%     s.YData = pos_tracker(2,k);
%       
%     q.XData = r(1,k);
%     q.YData = r(2,k);
%     
%     drawnow
% end
    
%% plot of r-g
figure (8)
for i=1:N
    e1(1,i)=norm(r_g(1,:,i));
    e2(1,i)=norm(r_g(2,:,i));
    e3(1,i)=norm(r_g(3,:,i));
end

plot (tbar,e1(1:Nbar),'LineWidth',1.5);
hold on
plot (tbar,e2(1:Nbar),'LineWidth',1.5);
plot (tbar,e3(1:Nbar),'LineWidth',1.5);

x1=xlabel('Time$~[s]$');
 y1=ylabel('Control Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
 leg1=legend('$15~km/h$','$25~km/h$','$35~km/h$');
 set(leg1,'Interpreter','latex')
 
 set(gcf, 'color', 'none');
set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('Control_error','-dsvg','-r0')