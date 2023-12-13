function orbit_propagator(semi_major,ecc,inclination_int,RAAN,AOP,n_orbits,nu_0)
clc 
 close all
%Constants
   mu=3.986004418e5; % km^3/s^2
   RE=6371;
   inclination_int=deg2rad(inclination_int);
   RAAN=deg2rad(RAAN);
   AOP=deg2rad(AOP);
   nu_0=deg2rad(nu_0);
   period=2*pi*sqrt(semi_major^3/mu);
   we=(2*pi+2*pi/365.26)/(24*3600);
   %% Function handles
    
   nu=@(E) 2*atan(tan(E/2)*sqrt((1+ecc)/(1-ecc)));
   Eccentric_anomaly=@(nu) (2*atan(tan(nu/2)*sqrt((1-ecc)/(1+ecc))));
   mean_anomaly=@(t) 2*pi*t/period;
   Kepler_eqn=@(E) E-ecc*sin(E);
   Rx=[];
   Ry=[];
   Rz=[];
   longitude=[];
   latitude=[];
   theta=[];
   rot_ex=[];
   rot_ey=[];
   rot_ez=[];
   x_plane=[];
   y_plane=[];
   iter=1;
   %% Initial and final time
        E_0=Eccentric_anomaly(nu_0);
        if(nu_0>pi)
        mean_anomaly_0=Kepler_eqn(E_0)+2*pi;
        else
             mean_anomaly_0=Kepler_eqn(E_0);
        end
        Init_time=mean_anomaly_0*period/(2*pi);
        final_time=n_orbits*period;
        time_domain=linspace(Init_time,final_time,10000);

%% Time loop
        for t=1:length(time_domain) 
            time=time_domain(t);
             Kepler_eqn_iter=@(E) E-ecc*sin(E) - mean_anomaly(time);
                E=fzero(Kepler_eqn_iter,E_0);
                E_0=E;
                r=(semi_major.*(1-ecc.^2))./(1+ecc*cos(nu(E)));
                coord(1,1)=r*cos(nu(E));
                coord(2,1)=r*sin(nu(E));
                coord(3,1)=0;
                x_plane=[x_plane;coord(1,1)];
                y_plane=[y_plane;coord(2,1)];
                coord=z_rot(coord,-AOP);
                coord=x_rot(coord,-inclination_int);
                coord=z_rot(coord,-RAAN);
                %save r in a vector
                rot_ex=[rot_ex;coord(1)];
                rot_ey=[rot_ey;coord(2)];
                rot_ez=[rot_ez;coord(3)];
                %create earth rotation
                theta(iter)=we*(time-Init_time);
                R_rel=z_rot(coord,theta(iter));
                %save the relative R in a vector
                Rx=[Rx;R_rel(1)];
                Ry=[Ry;R_rel(2)];
                Rz=[Rz;R_rel(3)];
                %calculating the longitude and the latitude
                l=R_rel(1)/norm(R_rel);
                m=R_rel(2)/norm(R_rel);
                n=R_rel(3)/norm(R_rel);
            
            
                latitude_temp=asind(n);
                %scaling the longitude to -180:180
                if m>0
                    longitude_temp=acosd(l/cosd(latitude_temp));
                else
                    longitude_temp=-acosd(l/cosd(latitude_temp));
                end
                %saving it into a vector
                longitude=[longitude;longitude_temp];
                latitude=[latitude;latitude_temp];
                iter=iter+1;
        end

        %% Plotting
        figure('Color','K')
        world_map=imread("land_ocean_ice_2048.jpg","jpg");
        ax=axes('XLim',[-2*semi_major 2*semi_major],'YLim',[-2*semi_major 2*semi_major],'ZLim',[-2*semi_major 2*semi_major],'Color','K');
        [x,y,z]=ellipsoid(0,0,0,RE,RE,RE,50);
        hold on;
        earth=surface(x,y,-z,'FaceColor','texturemap','CData',world_map,'EdgeColor','none','FaceAlpha',0.9);
        hg = hgtransform('Parent',ax);
        set(earth,'Parent',hg);
        axis equal
        axis off
        view(-42,14);
        u=0;
        %The equtorial plane
        x_vec=linspace(-2*semi_major,2*semi_major,40);
        y_vec=linspace(-2*semi_major,2*semi_major,40);
        z_mat=zeros(40,40);
        [x_mat,y_mat]=meshgrid(x_vec,y_vec);
        surf(x_mat,y_mat,z_mat);
        for vv=1:10:length(theta)
            earth_rot=makehgtform('zrotate',theta(vv));
            pl(vv)=plot3(rot_ex(vv),rot_ey(vv),rot_ez(vv),'v','MarkerFaceColor','auto','MarkerSize',10);
            set(hg,'Matrix',earth_rot);
            drawnow
            u=u+1;
            if u>1
                 delete(pl(vv));
                pl(vv)=plot3(rot_ex(vv),rot_ey(vv),rot_ez(vv),'--mx','LineWidth',0.8,'MarkerSize',2);
            end
        end
        figure(2)
        world_map=imread("land_ocean_ice_2048.jpg","jpg");
        image(-180:180,90:-1:-90,world_map);
        set(gca,'YDir','normal');
        axis normal
        ax=gca;
        axis([-180 180 -90 90]);
        ax.YTick=-90:15:90;
        ax.XTick=-180:20:180;
        hold on 
        axis equal
        axis tight
        grid on 
        for k=1:10:length(latitude)
        plot(longitude(k),latitude(k),'.r');
        pause(10e-6)
        end
        figure (3)
        plot(x_plane,y_plane);
end