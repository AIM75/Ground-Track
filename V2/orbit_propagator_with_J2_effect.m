%function orbit_propagator(semi_major,ecc,inclination_int,RAAN,AOP,n_orbits,nu_0)
clc
clearvars
close all
%% COEs
    semi_major= input("add the orbit semi major axies: ");
    inclination_int= input("add the orbit's inclination: ");
    ecc= input("add the orbit's eccentricity: ");
      
     if ecc>0 && ecc<1
            RAAN= input("add the Right Assension of the acending node: ");
            AOP= input("add the argument of perigee: ");
            elseif ecc==0
            AOP=0;
            RAAN=0;
            else 
            error("this programme is made to draw the ground track for elliptic & circular orbits only.")
     end
    n_orbits= input("add the number of orbits: ");
    nu_0= input("add the true anomaly: ");

%% Constants
       mu=3.986004418e5; % km^3/s^2
       RE=6371;% km
       period=2*pi*sqrt(semi_major^3/mu);
       we=2*pi/24/3600;%(2*pi/86164);%(2*pi+2*pi/365.26)/(24*3600);
   %Convert Degree 2 radian
       inclination_int=deg2rad(inclination_int);
       RAAN=deg2rad(RAAN);
       AOP=deg2rad(AOP);
       nu_0=deg2rad(nu_0);
%J2 effect 'earth Oblateness'
       J2=1.08262668e-3;
       fac=-3/2*sqrt(mu)*J2*RE^2/(1-ecc^2)^2/semi_major^(7/2);
       nodal_regression=fac*cos(inclination_int);                   % ' dRAAN/dt '
       perigee_rotation=fac*(5/2*sin(inclination_int)^2-2);         % ' dAOP/dt '
   %% Function handles
       nu=@(E) 2*atan(tan(E/2)*sqrt((1+ecc)/(1-ecc)));
       Eccentric_anomaly=@(nu) (2*atan(tan(nu/2)*sqrt((1-ecc)/(1+ecc))));
       mean_anomaly=@(t) 2*pi*t/period;
       Kepler_eqn=@(E) E-ecc*sin(E);
%% Initialize Plot variables
        Rx0=[];
        Ry0=[];
        Rz0=[];
        Rx=[];
        Ry=[];
        Rz=[];
        longitude0=[]; 
        latitude0=[];
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
                E_0=E; %Updating the inintial condition
                %Orbital plane calculation
                    r=(semi_major.*(1-ecc.^2))./(1+ecc*cos(nu(E)));
                    coord(1,1)=r*cos(nu(E));
                    coord(2,1)=r*sin(nu(E));
                    coord(3,1)=0;
                %Saving the Orbital plane coordinates in a vector for plotting
                    x_plane=[x_plane;coord(1,1)];
                    y_plane=[y_plane;coord(2,1)];
                %Convert from the orbital plane frame of reference to the The Geocentric-equatorial Coordinate System.
                    coord0=z_rot(coord,-AOP);
                    coord0=x_rot(coord0,-inclination_int);
                    coord0=z_rot(coord0,-RAAN);
                %Apply the Earth rotation
                    theta(iter)=we*(time-Init_time);
                    R_rel0=z_rot(coord0,theta(iter));
                %Save the relative R in a vector for plotting the ground track
                    Rx0=[Rx0;R_rel0(1)];
                    Ry0=[Ry0;R_rel0(2)];
                    Rz0=[Rz0;R_rel0(3)];
                %J2 effect Calculations
                    RAAN_new=nodal_regression*(time-Init_time)+RAAN;
                    AOP_new=perigee_rotation*(time-Init_time)+AOP;
                %Convert from the orbital plane frame of reference to the The Geocentric-equatorial Coordinate System + J2 effect
                    coord1=z_rot(coord,-AOP_new);
                    coord1=x_rot(coord1,-inclination_int);
                    coord1=z_rot(coord1,-RAAN_new);
                %Save the orbit coordinates in a vector for plotting the 3D orbit
                    rot_ex=[rot_ex;coord1(1)];
                    rot_ey=[rot_ey;coord1(2)];
                    rot_ez=[rot_ez;coord1(3)];
                %Apply earth rotation
                    R_rel=z_rot(coord1,theta(iter));
                %save the relative R in a vector
                    Rx=[Rx;R_rel(1)];
                    Ry=[Ry;R_rel(2)];
                    Rz=[Rz;R_rel(3)];
                %calculating the longitude and the latitude
                    l0=R_rel0(1)/norm(R_rel0);
                    m0=R_rel0(2)/norm(R_rel0);
                    n0=R_rel0(3)/norm(R_rel0);

                    latitude_temp0=asind(n0);
                    %scaling the longitude to -180:180
                        if m0>0
                            longitude_temp0=acosd(l0/cosd(latitude_temp0));
                        else
                            longitude_temp0=-acosd(l0/cosd(latitude_temp0));
                        end
                    %saving it into a vector
                        longitude0=[longitude0;longitude_temp0];
                        latitude0=[latitude0;latitude_temp0];

                  %With J2 effect
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
        plot(longitude,latitude,'.r');
        plot(longitude0,latitude0,'.k')
        legend("J2 effect","No J2 effect");
        figure (3)
        plot(x_plane,y_plane);
     