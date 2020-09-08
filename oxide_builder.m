function oxide_builder
    clc; 
    clear all;
    
    base_time = -4:0.1:3.0;
    time = 10.^base_time; %0.0:0.1:10.0;
    Tk = 298.0;
    kf = 504;
    ap = 0.32; %0.3; %0.43; %
    Ef = 0.6; %0.8; %0.4; %
    i0 = 4.0e-10; %1.0e-10;
    Vapp = 0.6;
    imono = monolayer_model(time);
    i_growth = interface_model(kf,ap,Ef,i0,Vapp,time,Tk);
  
    i_pass = passive_model(time);
    
    i_tot = i_growth + i_pass + imono;
%     i_check = 10.^-time;    
    
    figure(7)
    hold on
    box on
    plot(time, i_growth,'bo')
    plot(time, i_pass, 'r+')
    plot(time, imono, 'g^')
    plot(time, i_tot,'-k')
%     plot(time, i_check,'r+')
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    ylim([1.0e-15 1.0])
    hold off
    
    figure(8)
    hold on
    box on
    plot(time, i_growth,'bo')
    plot(time, i_pass, 'r+')
    plot(time, imono, 'g^')
    plot(time, i_tot,'-k')
%     plot(time, i_check,'r+')
    ax = gca;
%     ax.XScale = 'log';
%     ax.YScale = 'log';
    ylim([1.0e-6 0.5])
    xlim([0.0 10.0])
    hold off    
    
end

function ig = interface_model(kf, ap, Ef, i0, U, dt, T)
    F = 96485;
    R = 8.314;
    gp = (ap*F)/(R*T);
    num = i0 * exp(gp * U);
    
    denom = 1.0 + (gp * Ef * kf * num).*dt;
    
    ig = num./denom;
end

function im = monolayer_model(dt)
    i0 = 0.3e-3; %3.0e-1;
    k = 2.0e3; %5.0e3;
    im = i0.*exp(-k.*dt);
end

function ip = passive_model(dt)
    i0_pass = 1.0e-2/(80*60);
    ip = zeros(1,length(dt));
    ip(1,:) = i0_pass;
end