function [FM,OscPar,FMtype] = OscillatorType(OsID)

OscType={'Rossler' ,'Lorenz','Chen','HR','Duffing','vanderPol'};
FMtype = char(OscType(OsID));

syms a b c real
syms x(t) y(t) z(t)

switch FMtype
    case 'Rossler'
        FM=[ -y(t)-z(t); x(t)+a*y(t); b+z(t)*(x(t)-c)];
        ap=0.2;
        bp=0.2;
        cp=5.7;
    case 'Lorenz'
        FM=expand([ a*(y(t)-x(t)); x(t)*(b-z(t))-y(t); x(t)*y(t)-c*z(t)]);
        ap=10;
        bp=28;
        cp=2;
    case 'Chen'
        % Reference: xdot=a(y-x), ydot=(c-a-z)x+cy, zdot=xy-Bz, with
        % a=35, c=28, B=8/3 [paper's own symbols skip 'b']. Renamed here
        % in order of first appearance (a->a, c->b since it's the 2nd
        % distinct symbol, B->c since it's the 3rd) -- so eq2 and eq3
        % below use b/c, NOT the paper's own c/B labels directly.
        FM=[ a*(y(t)-x(t)); x(t)*(b-a-z(t))+b*y(t); x(t)*y(t)-c*z(t)];
        ap=35;
        bp=28;
        cp=8/3;
    case 'HR'
        FM=[ y(t)+3*(x(t))^2-(x(t))^3-z(t)+a;1-5*(x(t))^2-y(t);-b*z(t)+b*c*(x(t)+1.6)];
        ap=3.2;
        bp=0.006;
        cp=4;
    case 'Duffing'
        % Reference (forced Duffing): xdot=y, ydot=-h*y-x^3+q*sin(eta*t),
        % with eta=1, h=0.1, q=5.6. Non-autonomous (explicit t via the
        % forcing term) -- converted to autonomous 3-state form via the
        % standard "z=t" trick: add zdot=1 (so z(t)=t given z(0)=0), and
        % replace t with z inside sin(). Renamed in order of first
        % appearance: h->a, q->b, eta->c.
        FM=[ y(t); -a*y(t)-x(t)^3+b*sin(c*z(t)); 1];
        ap=0.1;
        bp=5.6;
        cp=1;
    case 'vanderPol'
        % Reference (forced van der Pol): xdot=y,
        % ydot=-x+d*(1-x^2)*y+F*sin(eta*t), with d=3, F=15, eta=4.065.
        % Same z=t autonomization trick as Duffing above. Renamed in
        % order of first appearance: d->a, F->b, eta->c.
        FM=[ y(t); -x(t)+a*(1-x(t)^2)*y(t)+b*sin(c*z(t)); 1];
        ap=3;
        bp=15;
        cp=4.065;
    otherwise
        FM=[ -y(t)-z(t); x(t)+a*y(t); b+z(t)*(x(t)-c)];
end

OscPar=[ap,bp,cp];
end