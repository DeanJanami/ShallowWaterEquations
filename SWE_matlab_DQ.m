%% Dean Pakravan: 757389
% Qijie Li: 927249 
% The Shallow Water Equations
% September 21st 2018
%% Main function
function SWE_matlab_DQ()

    clear all; close all; clc;
    
    global Delta_x N_x Delta_y N_y g;

    % Simulation parameters
    x_min           =  0.00;
    x_max           = 100.00;
    y_min           =  0.00;
    y_max           = 100.00;
    t_min           =  0.00;
    t_max           = 20.00;
    N_x             = 101;
    N_y             = 101;
    N_t             = 101;
    Delta_x         =  (x_max-x_min)/(N_x-1);
    Delta_y         =  (y_max-y_min)/(N_y-1);
    Delta_t         =  (t_max-t_min)/(N_t-1);      
    [x, y]           = meshgrid(x_min:Delta_x:x_max, y_min:Delta_y:y_max);
    t               = t_min:Delta_t:t_max;
    g               = 9.81; % Our gravitational constant

    
    % Allocate arrays
    phi              = zeros(N_x, N_y, 3);
    % Vx(i,j) = phi(i,j,1)
    % Vy(i,j) = phi(i,j,2)
    %  h(i,j) = phi(i,j,3)
    
    
%% Initial Conditions
    % Vx = 0 and Vy = 0 for all x,y
    phi(:,:,1)       = 0;
    phi(:,:,2)       = 0;

    % Inital drops
    phi(:,:,3)       = 1 + 0.5*exp((-1/25)*((x-30).^2 + ((y-30).^2)))...
                     + 0.5*exp((-1/25)*((x-80).^2 + ((y-20).^2)))...
                     + 0.5*exp((-1/25)*((x-50).^2 + ((y-70).^2)));
    
    %% Plotting preallocation
    h=figure('WindowStyle', 'normal');
    Solution        = surf(x,y,phi(:,:,3), 'Linestyle', 'none');
    axis([x_min x_max y_min y_max 0 4],'on');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view([45 25]);
    colormap('winter');
    drawnow;
    
    %% Time marching for-loop

    
    filename = 'testAnimated.gif';
    
    for l=1:N_t-1
        
        k1          = f(phi);
        k2          = f(phi + Delta_t/2*k1);
        k3          = f(phi + Delta_t/2*k2);
        k4          = f(phi + Delta_t*k3);
        
        phi     =   phi + Delta_t  *(k1/6 + k2/3 + k3/3 + k4/6);

        % Plot the solution
        set(Solution,'Zdata', phi(:,:,3));
        title(['t = ' num2str(t(l+1))]);
        drawnow;
        
           %Capture the plot as an image 
          frame = getframe(h); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 

          %Write to the GIF File 
          if l == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
        
        
    end
        
return

%% RK4 function
function k = f(phi)
    global Delta_x N_x Delta_y N_y g;

    k    = zeros(N_x,N_y,3);
    for i=1:N_x
        for j = 1:N_y
            % Periodic boundary conditions on i
            n_1 = i; n_2 = i; n_3 = i; n_4 = i; n_5 = i; n_6 = i;
            if i==3
                n_1 = N_x+3;
            elseif i==2
                n_1 = N_x+2;
                n_2 = N_x+2;
            elseif i==1
                n_1 = N_x+1;
                n_2 = N_x+1;
                n_3 = N_x+1;
            elseif i==N_x-2
                n_6 = -2;
            elseif i==N_x-1
                n_5 = -1;
                n_6 = -1;
            elseif i==N_x
                n_4 = 0;
                n_5 = 0;
                n_6 = 0;
            end
            % Periodic boundary conditions on j
            m_1 = j; m_2 = j; m_3 = j; m_4 = j; m_5 = j; m_6 = j;
            if j==3
                m_1 = N_y+3;
            elseif j==2
                m_1 = N_y+2;
                m_2 = N_y+2;
            elseif j==1
                m_1 = N_y+1;
                m_2 = N_y+1;
                m_3 = N_y+1;
            elseif j==N_y-2
                m_6 = -2;
            elseif j==N_y-1
                m_5 = -1;
                m_6 = -1;
            elseif j==N_y
                m_4 = 0;
                m_5 = 0;
                m_6 = 0;
            end
            % k for Vx
            k(i,j,1)    = (-g/Delta_x)*((-1/60)*phi(n_1-3,j,3) + (3/20)*phi(n_2-2,j,3) + (-3/4)*phi(n_3-1,j,3) + (3/4)*phi(n_4+1,j,3) + (-3/20)*phi(n_5+2,j,3) + (1/60)*phi(n_6+3,j,3))...
                            + (-phi(i,j,1)/Delta_x)*((-1/60)*phi(n_1-3,j,1) + (3/20)*phi(n_2-2,j,1) + (-3/4)*phi(n_3-1,j,1) + (3/4)*phi(n_4+1,j,1) + (-3/20)*phi(n_5+2,j,1) + (1/60)*phi(n_6+3,j,1)) ...
                            + (-phi(i,j,2)/Delta_y)*((-1/60)*phi(i,m_1-3,1) + (3/20)*phi(i,m_2-2,1) + (-3/4)*phi(i,m_3-1,1) + (3/4)*phi(i,m_4+1,1) + (-3/20)*phi(i,m_5+2,1) + (1/60)*phi(i,m_6+3,1));
            
            % k for Vy
            k(i,j,2)    = (-g/Delta_y)*((-1/60)*phi(i,m_1-3,3) + (3/20)*phi(i,m_2-2,3) + (-3/4)*phi(i,m_3-1,3) + (3/4)*phi(i,m_4+1,3) + (-3/20)*phi(i,m_5+2,3) + (1/60)*phi(i,m_6+3,3))...
                            + (-phi(i,j,1)/Delta_x)*((-1/60)*phi(n_1-3,j,2) + (3/20)*phi(n_2-2,j,2) + (-3/4)*phi(n_3-1,j,2) + (3/4)*phi(n_4+1,j,2) + (-3/20)*phi(n_5+2,j,2) + (1/60)*phi(n_6+3,j,2))...
                            + (-phi(i,j,2)/Delta_y)*((-1/60)*phi(i,m_1-3,2) + (3/20)*phi(i,m_2-2,2) + (-3/4)*phi(i,m_3-1,2) + (3/4)*phi(i,m_4+1,2) + (-3/20)*phi(i,m_5+2,2) + (1/60)*phi(i,m_6+3,2));
          

            % k for h
            k(i,j,3)    = (-1/Delta_x)*((-1/60)*phi(n_1-3,j,1)*phi(n_1-3,j,3) + (3/20)*phi(n_2-2,j,1)*phi(n_2-2,j,3) + (-3/4)*phi(n_3-1,j,1)*phi(n_3-1,j,3) + (3/4)*phi(n_4+1,j,1)*phi(n_4+1,j,3) + (-3/20)*phi(n_5+2,j,1)*phi(n_5+2,j,3) + (1/60)*phi(n_6+3,j,1)*phi(n_6+3,j,3)) ...
                         + (-1/Delta_y)*((-1/60)*phi(i,m_1-3,2)*phi(i,m_1-3,3) + (3/20)*phi(i,m_2-2,2)*phi(i,m_2-2,3) + (-3/4)*phi(i,m_3-1,2)*phi(i,m_3-1,3) + (3/4)*phi(i,m_4+1,2)*phi(i,m_4+1,3) + (-3/20)*phi(i,m_5+2,2)*phi(i,m_5+2,3) + (1/60)*phi(i,m_6+3,2)*phi(i,m_6+3,3));
        end
    end
    
return
    
    