clc
clearvars
close all
% COEs
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
        %% Initialize Plot variables
        Ecc_anomaly=nan(1,length(time_domain));
        True_anomaly=nan(1,length(time_domain));
        r=nan(1,length(time_domain));
        longitude=nan(1,length(time_domain));
        latitude=nan(1,length(time_domain));
        theta=nan(1,length(time_domain));
        R=[];
%% Time loop
        for t=1:length(time_domain) 
                time=time_domain(t);
                Kepler_eqn_iter=@(E) E-ecc*sin(E) - mean_anomaly(time);
                Ecc_anomaly(t)=fzero(Kepler_eqn_iter,E_0);
                E_0=Ecc_anomaly(t); %Updating the inintial condition
        end
        True_anomaly=nu(Ecc_anomaly);
        r=(semi_major.*(1-ecc.^2))./(1+ecc*cos(True_anomaly));
        coord(1,:)=r.*cos(True_anomaly);
        coord(2,:)=r.*sin(True_anomaly);
        coord(3,:)=0;
        RAAN_new=nodal_regression*(time_domain-Init_time)+RAAN;
        AOP_new=perigee_rotation*(time_domain-Init_time)+AOP;
        theta=we*(time_domain-Init_time);
        %Coordinates with J2 effect
        X=sum([-sin(RAAN_new).*cos(inclination_int).*sin(AOP_new)+cos(RAAN_new).*cos(AOP_new);-sin(RAAN_new).*cos(inclination_int).*cos(AOP_new)-cos(RAAN_new).*sin(AOP_new);sin(RAAN_new).*sin(inclination_int)].*coord,1);
        Y=sum([cos(RAAN_new).*cos(inclination_int).*sin(AOP_new)+sin(RAAN_new).*cos(AOP_new);cos(RAAN_new).*cos(inclination_int).*cos(AOP_new)-sin(RAAN_new).*sin(AOP_new);-cos(RAAN_new).*sin(AOP_new)].*coord,1);
        Z=sum([sin(inclination_int).*sin(AOP_new);sin(inclination_int).*cos(AOP_new);repelem(cos(inclination_int),length(time_domain))].*coord,1);
        X_earth_rot=sum([cos(theta);sin(theta);zeros(1,length(time_domain))].*[X;Y;Z],1);
        Y_earth_rot=sum([-sin(theta);cos(theta);zeros(1,length(time_domain))].*[X;Y;Z],1);
        Z_earth_rot=sum([zeros(1,length(time_domain));zeros(1,length(time_domain)); ones(1,length(time_domain))].*[X;Y;Z],1);
        z=vecnorm([X_earth_rot;Y_earth_rot;Z_earth_rot]);
        l=X_earth_rot./z;
        m=Y_earth_rot./z;
        n=Z_earth_rot./z;
        latitude=asind(n);
        q=m>0;
        longitude=acosd(l./cosd(latitude));
        longitude(~q)=-acosd(l(~q)./cosd(latitude(~q)));
        z=ones(40);

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
            pl(vv)=plot3(X(vv),Y(vv),Z(vv),'v','MarkerFaceColor','auto','MarkerSize',10);
            set(hg,'Matrix',earth_rot);
            drawnow
            u=u+1;
            if u>1
                 delete(pl(vv));
                pl(vv)=plot3(X(vv),Y(vv),Z(vv),'--mx','LineWidth',0.8,'MarkerSize',2);
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
        plot(longitude,latitude,'.','ColorMode','auto');
        text(longitude(1),latitude(1),'\leftarrow start point','Color','red','FontSize',20);

        xlabel('Longitude')
        ylabel('Latitude')
        title('Ground Track')