% ID_Ortogonalna_Korelacija_Izracun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Debelina linij
deb_lin = 2;

%% Show the properties of the noise

if s > 0
    white_noise=s*randn(n,N);

    % Lastnosti suma. Simulirali bomo N potekov suma in povprecili
    % spektralno mocnostno gostoto.
    t=linspace(0,0.1*(n-1),n);
    for i=1:N
        v(:,i)=lsim(1,[1/omega_g 1],white_noise(:,i),t);
        avtokorelacija(:,i) = xcorrp(v(:,i),v(:,i))/t(2);
    end
    % frekvencna os
    df = 1/t(end);
    f_os = 0:df:(n-1)*df;
    % frekvencni spekter
    fv=abs(fft(v));
    % Spektralna mocnostna gostota
    PHI =t(2)* abs(fv).^2/t(end);
    PHI = mean(PHI.');
    % Teoreticno
    PHI_T = s^2./(1+(2*pi*f_os/omega_g).^2);
    % prikaz rezultatov  
    figure
    % set(gcf,'Units','normalized', 'Outerposition', [0 0.04 1 0.96]); % poveèamo graf èez cel ekran
    subplot(211)
    hLu = semilogx(f_os(2:end/2),PHI(2:end/2),'b',f_os(2:end/2),PHI_T(2:end/2),'r');
    set(hLu,'linewidth',deb_lin);
    title('Power spectral density of the noise','Fontsize',16)
    xlabel('$f$ [Hz]','Fontsize',16,'Interpreter','Latex');
    ylabel('$\Phi_{nn}(f)$','Fontsize',16,'Interpreter','Latex');
    l = legend('measured','true');
    set(l,'Fontsize',16,'Interpreter','Latex');
    besedilo2 = ['$\Phi_{nn}(\omega)=\frac{\Phi_0}{1+',...
        '\left(\frac{\omega}{\omega_g}\right)^2}$'];
    hh2 = text(0.02, 0.02, besedilo2,'Units','normalized');
    set(hh2, 'Interpreter','Latex','Color',[0 0 1],'Fontsize',20, ...
        'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    %
    % Se avtokorelacijska funkcija suma
    Avtok = mean(avtokorelacija.');
    tau = t-(t(end)+t(2))/2;
    % teoreticno
    Avtok_T = s^2*omega_g/2*exp(-abs(tau*omega_g));
    
    subplot(212)
    hLu = plot(tau,fftshift(Avtok),'b',tau,Avtok_T,'r');
    set(hLu,'linewidth',deb_lin);
    title('Autocorrelation of the noise','Fontsize',16)
    xlabel('$\tau$ [s]','Fontsize',16,'Interpreter','Latex');
    ylabel('$\phi_{nn}(\tau)$','Fontsize',16,'Interpreter','Latex');
    l = legend('measured','true');
    set(l,'Fontsize',16,'Interpreter','Latex');
    besedilo3 = ['$\phi_{nn}(\tau)=', ...
        '\frac{\Phi_0 e^{-\mid\tau\mid\omega_g}}{2}\omega_g$'];
    hh2 = text(0.02, 0.98, besedilo3,'Units','normalized');
    set(hh2, 'Interpreter','Latex','Color',[0 0 1],'Fontsize',20, ...
        'VerticalAlignment','top','HorizontalAlignment','left');
    %
    % pocakamo da pritisnemo tipko
    gca_pos = get(gca,'Position');
    hh = text(1+(1-gca_pos(1)-gca_pos(3)), -gca_pos(2), 'Press any key', ...
        'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',14, 'Color',[0 0.5 0]);
    pause
    delete(hh);

end


%% zanka izracunov pri posameznih frekvencah
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prealokacija
ampl = nan(1,N); faza = nan(1,N);
% Stevilo period opazovanja izberemo tako, da se pri prikazu kaj vidi
Stevilo_period_opazovanja = round(linspace(2,30,N));

for i=1:N
    % definiramo vhodni signal; casovna os se prilagaja frekvenci
    Perioda = 1/f(i);
    % zacetni cas je vsaj 300 sek, da ninejo prehodni pojavi
    t_zacetni = 300;
    Zacetnih_period = ceil(t_zacetni/Perioda);
    % zacetni cas mora biti na zacetku periode
    t_zacetni = Zacetnih_period * Perioda;
    t_koncni = t_zacetni + Stevilo_period_opazovanja(i) * Perioda;
    t=linspace(0,t_koncni,n)';
    %
    % vhodni signal je cosinus
    x=cos(2*pi*f(i)*t);
    %
    % Za izracun bomo rabili se sinus
    x1=sin(2*pi*f(i)*t);
    %
    % Indeksi zacetka in konca prikazanega signala
    ind1 = round(t_zacetni/t(2));
    ind3 = round(t_koncni/t(2));
    % izracunamo izhodni signal in ga zmotimo s sumom
    y=lsim(b,a,x,t);
    v=0;
    if s > 0
        if omega_g < 100
            v = v+lsim(1,[1/omega_g 1],white_noise(:,i),t);
        else
            v = v+white_noise(:,i);
        end
    end
    if A_sin > 0
        v = v+A_sin* sin(omega_g*t);
    end
    if Drift >0
        v = v + Drift * (t-t(ind1));
    end
    y=y+v;
    %
    % prikaz vhodnega in izhodnega signala na zaslonu le  za 5 izbranih
    % frekvenc
    if sum([1  13  26  38 N] == i) ~= 0
        %
        % Narisemo vhodni in izhodni signal pri trenutni frekvenci
        figure
        % set(gcf,'Units','normalized', 'Outerposition', [0 0.04 1 0.96]); % poveèamo graf èez cel ekran
        subplot(211)
        hLu = plot(t(ind1:ind3),x(ind1:ind3),'g',t(ind1:ind3),y(ind1:ind3),'b');
        %plot(t(ind1:ind3),x(ind1:ind3),'g',t(ind1:ind3),y(ind1:ind3),'b');
        set(hLu,'linewidth',deb_lin);
        if flag_detrend
             ydetr = detrend(y);
             hold on;
             hLu1 = plot(t(ind1:ind3),ydetr(ind1:ind3),'b--');
             set(hLu1,'linewidth',deb_lin);
        end
        set(gca,'xlim',[t_zacetni t_koncni]);%,'ylim',[-3 3])
        title(['frequency = ',num2str(f(i)),' Hz'],'Fontsize',16)
        xlabel('$t$ [s]','Fontsize',16,'Interpreter','Latex');
        ylabel('$u(t)$, $y(t)$','Fontsize',16,'Interpreter','Latex');
        l = legend(hLu,'$u(t)$','$y(t)$');
        set(l,'Fontsize',16,'Interpreter','Latex', 'AutoUpdate','off');
        %
        % Ceprav za izracun ni potrebno izracunati celotnih korelacijskih
        % funkcij, jih tu za ilustracijo eliminacije suma prikazimo:
        xx=xcorrp(x(ind1:ind3),x(ind1:ind3));
        xy=xcorrp(x(ind1:ind3),y(ind1:ind3));
        %
        % Ilustrativni prikaz avtokorelacijske funkcije vhoda
        % in krizne korelacije med vhodom in izhodom
        %
        %      Vpliv suma postane zanemarljiv.
        n4=ind3-ind1;
        subplot(212);
        hLu =  plot(t(1:n4),xx(1:n4),'g',t(1:n4),xy(1:n4),'b');
        %plot(t(ind1:ind3),x(ind1:ind3),'g',t(ind1:ind3),y(ind1:ind3),'b');
        set(hLu,'linewidth',deb_lin);
        set(gca,'xlim',[0 t_koncni-t_zacetni]);%,'ylim',[-1.5 1.5])
        title(['frequency = ',num2str(f(i)),' Hz'],'Fontsize',16)
        xlabel('$\tau$ [s]','Fontsize',16,'Interpreter','Latex');
        ylabel('$\phi_{uu}(\tau)$, $\phi_{uy}(\tau)$','Fontsize',16,'Interpreter','Latex');
        l = legend(hLu,'$\phi_{uu}(\tau)$','$\phi_{uy}(\tau)$');
        set(l,'Fontsize',16,'Interpreter','Latex', 'AutoUpdate','off');
        %
        % vzorcenje realne komponente korelacijskih funkcij (samo za prikaz na
        % zaslonu, racuna se drugace
        Stevilo_period = floor(t_koncni/Perioda)-1;
        t_real = (Stevilo_period - Zacetnih_period)*Perioda;
        % izracun indeksa az vorcenje realne komponente - pri polni periodi
        ind_real = round(t_real/t(2))+1;
        % izracun indeksa az vorcenje imaginarne komponente - pri 3/4
        % periode
        t_imag = (Stevilo_period - Zacetnih_period - 1/4)*Perioda;
        ind_imag = round(t_imag/t(2))+1;
        %
        % Prikazemo tocko na koncu in na 3/4 periode, ki sta realni in
        % imaginarni del frekvenènega odziva
        hold on
        h_Re = plot(t(ind_real),xy(ind_real),'ro');
        set(h_Re,'Markersize',24,'linewidth',3);
        h_Re = plot(t(ind_real),xy(ind_real),'r+');
        set(h_Re,'Markersize',24,'linewidth',3);
        %
        h_Im = plot(t(ind_imag),xy(ind_imag),'mo');
        set(h_Im,'Markersize',24,'linewidth',3);
        h_Im = plot(t(ind_imag),xy(ind_imag),'m+');
        set(h_Im,'Markersize',24,'linewidth',3);
        %
        besedilo4a = ['$\phi_{uu}(\tau)  = \lim_{T_i\rightarrow \infty}',...
            '\frac{1}{T_i} \int_{-\frac{T_i}{2}}^{\frac{T_i}{2}} u(t)u(t+\tau)dt$'];
        besedilo4b = ['$\phi_{uy}(\tau)  = \lim_{T_i\rightarrow \infty}',...
            '\frac{1}{T_i} \int_{-\frac{T_i}{2}}^{\frac{T_i}{2}} u(t)y(t+\tau)dt$'];
        hh2 = text(0.5, 0.02, [besedilo4a,'  ;  ',besedilo4b],'Units','normalized');
        set(hh2, 'Interpreter','Latex','Color',[0 0 1],'Fontsize',20, ...
            'VerticalAlignment','bottom','HorizontalAlignment','center');
        % pocakamo da si ogledamo eno frekvenco
        gca_pos = get(gca,'Position');
        hh = text(1+(1-gca_pos(1)-gca_pos(3)), -2*gca_pos(2), 'Press any key', ...
            'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',14, 'Color',[0 0.5 0]);
        pause
        delete(hh);
    end
    % ali naredimo detrend meritve izhodnega signala
    if flag_detrend
        y=detrend(y);
    end
    %
    % izracun re. in imag. komponente frekvencnega odziva
    %
    % Formule: 3.49 in 3.51 pomenijo 2 krat srednja vrednost produktov
    re=y(ind1:ind3)'*x(ind1:ind3)*2/(ind3-ind1);
    im=-y(ind1:ind3)'*x1(ind1:ind3)*2/(ind3-ind1);
    ampl(i)=sqrt(re^2+im^2);
    faza(i)=atan2(im,re);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      'AMPLITUDNI odziv pri N frekvencah dobljen z metodo ortogonalne'
%      'korelacije je prikazan na zgornji sliki. Poleg so prikazani se    '
%      'pravi odziv in pas moznih amplitudnih odzivov.                    '
%      '                                                                  '
%      'Ocenjena in prava vrednost FAZNEGA odziva sta prikazani na        '
%      'sliki spodaj.                                                     ';
%
% glajenje faze
faza=unwrap(faza);
%
% izracun pravega frekvencnega odziva
[ampl1,faza1]=bode(b,a,2*pi*f);
%
% in pasu dveh stand. deviacij
% n T_p, to je cas integracije
nTp = t_koncni -t_zacetni;
napaka = 0;
if s + A_sin + Drift == 0
    besedilo5 = '';
else
    besedilo5 = ['$\sigma_{G}^2=$  '];
end
if s > 0
    if omega_g >=100
        % izracun pasu napake za beli sum po furmuli 3.71:
        napaka = napaka + 2*2*s./(1*sqrt(nTp));
        besedilo5 = [besedilo5,'$\frac{4\Phi_0}{U_0^2nT_p}$'];
    else
        % izracun pasu napake za barvni sum po furmuli 3.74:
        % Vzamemo teoreticni potek spektralne moènostne gostote barvnega
        % suma
        PHI_T = s^2./(1+(2*pi*f/omega_g).^2);
        napaka = napaka+2*PHI_T./(1*sqrt(nTp));
        besedilo5 = [besedilo5,'$\frac{4\Phi_{nn}(\omega_0)}{U_0^2nT_p}$'];
    end
end
if A_sin > 0
    % izracun pasu napake za sinusno motnjo po furmulah 3.90, 3.91:
    % najprej izracun P -3.91
    normomega = omega_g./(2*pi*f);
    cs=cos(pi*normomega*Stevilo_period_opazovanja(i));
    stevec=2*normomega.*sqrt(1-(1-normomega.*normomega).*(cs.*cs));
    imen=abs(1-normomega.*normomega);
    P=(stevec./imen).*(abs(sin(pi*normomega*Stevilo_period_opazovanja(i) ...
        ))./(pi*normomega*Stevilo_period_opazovanja(i)));
    % nato pas 2 krat stand. deviacija  3.90
    napaka = sqrt(napaka^2+(2*A_sin/1 * P).^2);
    besedilo5 = [besedilo5,' + $\frac{N_0}{U_0}P$'];
end
if Drift >0
    if ~ flag_detrend
        % izracun pasu napake za drift po furmuli 3.99,
        napaka = sqrt(napaka.^2+(2*Drift./(1*2*pi*f).^2));
        besedilo5 = [besedilo5,' + $\frac{4D_0^2}{U_0^2\omega_0^2}$'];
    else
        napaka = nan(size(napaka));
        besedilo5 = '';
    end
end
napaka1=ampl1+napaka';
napaka2=ampl1-napaka';
%
% Omejimo napako na spodnji rob slike
napaka2 = max(napaka2,0.0001);
%
figure
% set(gcf,'Units','normalized', 'Outerposition', [0 0.04 1 0.96]); % poveèamo graf èez cel ekran
subplot(211);
hLu = loglog(f,ampl1,'r',f,ampl,'.k',f,napaka1,'--m',f,napaka2,'--m');
set(hLu,'linewidth',deb_lin,'Markersize',18);
set(gca,'xlim',[f(1) f(end)],'ylim',[0.001 100])
title('Amplitude response','Fontsize',16)
xlabel('$f$ [Hz]','Fontsize',16,'Interpreter','Latex');
ylabel('$|G(\jmath 2 \pi f)|$','Fontsize',16,'Interpreter','Latex');
l = legend('true','identified');
set(l,'Fontsize',16,'Interpreter','Latex');

hh2 = text(0.02, 0.02, besedilo5,'Units','normalized');
set(hh2, 'Interpreter','Latex','Color',[0 0 1],'Fontsize',20, ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');

besedilo6 = ['$\Re [G(\jmath\omega_0)]=', ...
    '\frac{2}{U_0nT_p}\int_0^{nTp}y(t)\cos(\omega_0t)dt$'];
besedilo7 = ['$\Im [G(\jmath\omega_0)]=\frac{-2}', ...
    '{U_0nT_p}\int^{nT_p}_0y(t)\sin(\omega_0t)dt$'];

hh2 = text(0.02, 0.98, [{besedilo6},{besedilo7}],'Units','normalized');
set(hh2, 'Interpreter','Latex','Color',[0 0 1],'Fontsize',20, ...
    'VerticalAlignment','top','HorizontalAlignment','left');
%
% Faza
subplot(212)
hLu = semilogx(f,faza1,'r',f,faza*180/pi,'.k');
set(hLu,'linewidth',deb_lin,'Markersize',18);
set(gca,'xlim',[f(1) f(end)],'ylim',[-300 0])
title('Phase response','Fontsize',16)
xlabel('$f$ [Hz]','Fontsize',16,'Interpreter','Latex');
ylabel('$\angle \bigl ( G(\jmath 2 \pi f)\bigr )$','Fontsize',16,'Interpreter','Latex');
l = legend('true','identified');
set(l,'Fontsize',16,'Interpreter','Latex');
