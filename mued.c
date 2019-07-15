model 'mued'/1001.
option chepPDWidth=236.
alias 
	nmax=6.

keys 
	maxkk=1.

do_if maxkk==1.

        let
                cos1=cos(1),
                cos2=0,
                cos3=0,
                cos4=0,
                cos5=0,
                cos6=0,
                sin1=sin(1),
                sin2=0,
                sin3=0,
                sin4=0,
                sin5=0,
                sin6=0.

do_else_if maxkk==2.
        let
                cos1=cos(1),
                cos2=cos(2),
                cos3=0,
                cos4=0,
                cos5=0,
                cos6=0,
                sin1=sin(1),
                sin2=sin(2),
                sin3=0,
                sin4=0,
                sin5=0,
                sin6=0.
do_else_if maxkk==4.
        let
                cos1=cos(1),
                cos2=cos(2),
                cos3=cos(3),
                cos4=cos(4),
                cos5=0,
                cos6=0,
                sin1=sin(1),
                sin2=sin(2),
                sin3=sin(3),
                sin4=sin(4),
                sin5=0,
                sin6=0.
do_else_if maxkk==6.
        let
                cos1=cos(1),
                cos2=cos(2),
                cos3=cos(3),
                cos4=cos(4),
                cos5=cos(5),
                cos6=cos(6),
                sin1=sin(1),
                sin2=sin(2),
                sin3=sin(3),
                sin4=sin(4),
                sin5=sin(5),
                sin6=sin(6).
end_if.

option ReduceGamma5=0.

/* parameters */
/**************/

/* free parameters */

parameter
	alphae=1/128,
	pi=acos(-1),
	ge=sqrt(4*pi*alphae),
	gc=1.21978,
	sw=0.471813,
	CKMs12=0.22506,
	CKMs23=0.0410788,
	CKMs13=0.00357472,
	MZ0=91.1876,
	Mh0=125.

SetTexName([
	alphae='\\alpha',
	pi='\\pi',
	ge='e',
	gc='g_s',
	sw='s_W',
	CKMs12='s_{12}',
	CKMs13='s_{13}',
	CKMs23='s_{23}',
	MZ0='m_Z',
	Mh0='m_h'
]).
	
/* derived parameters */

parameter
	cw=sqrt(1-sw**2),
	MW0=MZ0*cw,
	muH=Mh0/Sqrt2,
	lamH=((ge/sw)*Mh0/8/(MZ0*cw))**2.

SetTexName([
	cw='c_W',
	MW0='m_W',
	muH='\\mu',
	lamH='\\lambda'
]).
		
/* masses */

parameter
	Mu0=0.0022,
	Md0=0.0047,
	Mc0=1.28,
	Ms0=0.096,
	Mt0=175,
	Mb0=4.18.

SetTexName([
	Mu0='m_u',
	Md0='m_{d^i}',
	Mc0='m_c',
	Ms0='m_s',
	Mt0='m_t',
	Mb0='m_b'
]).

_x=[h,Z,W,u,d,c,s,t,b] in parameter
	M_x02=M_x0**2.

SetTexName([
	Mh02='m_h^2',
	MZ02='m_Z^2',
	MW02='m_W^2',
	Mu02='m_u^2',
	Md02='m_{d^i}^2',
	Mc02='m_c^2',
	Ms02='m_s^2',
	Mt02='m_t^2',
	Mb02='m_b^2'
]).

parameter
	n=1,
	KKmn=1,
	invR=500.

SetTexName([
	n='n',
	KKmn='m_{n}',
	invR='(n/KKmn)^{-1}'
]).

%SetTexName([
%	invR='(n/KKmn)^{-1}',
%	(n/KKmn)='(n/KKmn)'
%]).

_n=1-nmax in parameter
	Mhn_n=sqrt(Mh0**2+(_n/(n/KKmn))**2),
	Ma0n_n=sqrt(MZ0**2+(_n/(n/KKmn))**2),
	Macn_n=sqrt(MW0**2+(_n/(n/KKmn))**2),
	Mu1n_n=sqrt(Mu0**2+(_n/(n/KKmn))**2),
	Mu2n_n=sqrt(Mu0**2+(_n/(n/KKmn))**2),
	Md1n_n=sqrt(Md0**2+(_n/(n/KKmn))**2),
	Md2n_n=sqrt(Md0**2+(_n/(n/KKmn))**2),
	Mc1n_n=sqrt(Mc0**2+(_n/(n/KKmn))**2),
	Mc2n_n=sqrt(Mc0**2+(_n/(n/KKmn))**2),
	Ms1n_n=sqrt(Ms0**2+(_n/(n/KKmn))**2),
	Ms2n_n=sqrt(Ms0**2+(_n/(n/KKmn))**2),
	Mt1n_n=sqrt(Mt0**2+(_n/(n/KKmn))**2),
	Mt2n_n=sqrt(Mt0**2+(_n/(n/KKmn))**2),
	Mb1n_n=sqrt(Mb0**2+(_n/(n/KKmn))**2),
	Mb2n_n=sqrt(Mb0**2+(_n/(n/KKmn))**2),
	MGn_n=(_n/(n/KKmn)),
	MAn_n=(_n/(n/KKmn)),
	MZn_n=sqrt(MZ0**2+(_n/(n/KKmn))**2),
	MWn_n=sqrt(MW0**2+(_n/(n/KKmn))**2).

SetTexName([
	Mhn1='m_{h_n}',
	Ma0n1='m_{a_n}',
	Macn1='m_{a^\\pm_n}',
	Mu1n1='m_{u_n^1}',
	Mu2n1='m_{u_n^2}',
	Md1n1='m_{d_n^{1i}}',
	Md2n1='m_{d_n^{2i}}',
	Mc1n1='m_{c_n^1}',
	Mc2n1='m_{c_n^2}',
	Ms1n1='m_{s_n^1}',
	Ms2n1='m_{s_n^2}',
	Mt1n1='m_{t_n^1}',
	Mt2n1='m_{t_n^2}',
	Mb1n1='m_{b_n^1}',
	Mb2n1='m_{b_n^2}',
	MGn1='m_{g_n}',
	MAn1='m_{A_n}',
	MZn1='m_{Z_n}',
	MWn1='m_{W_n}'
]).

_n=1-nmax, _x=[h,a0,ac,u1,u2,d1,d2,c1,c2,s1,s2,t1,t2,b1,b2,G,A,Z,W] in parameter
	M_xn_n2=M_xn_n**2.

SetTexName([
	Mhn2='m_{h_n}^2',
	Ma0n2='m_{a_n}^2',
	Macn2='m_{a^\\pm_n}^2',
	Mu1n2='m_{u^1_n}^2',
	Mu2n2='m_{u^2_n}^2',
	Md1n2='m_{d^{1i}_n}^2',
	Md2n2='m_{d^{2i}_n}^2',
	Mc1n2='m_{c^1_n}^2',
	Mc2n2='m_{c^2_n}^2',
	Ms1n2='m_{s^1_n}^2',
	Ms2n2='m_{s^2_n}^2',
	Mt1n2='m_{t^1_n}^2',
	Mt2n2='m_{t^2_n}^2',
	Mb1n2='m_{b^1_n}^2',
	Mb2n2='m_{b^2_n}^2',	
	MGn2='m_{g_n}^2',
	MAn2='m_{A_n}^2',
	MZn2='m_{Z_n}^2',
	MWn2='m_{W_n}^2'
]).

_n=1-nmax in parameter
	MBBn_n=sqrt((MZ0*sw)**2+(_n/(n/KKmn))**2).

SetTexName([
	MBBn1='m_{B_n}'
]).

/* projection operators */

let
	PL=(1-gamma5)/2,
	PR=(1+gamma5)/2.

/* ckm matrix */

parameter
	CKMc12=sqrt(1-CKMs12**2),
	CKMc23=sqrt(1-CKMs23**2),
	CKMc13=sqrt(1-CKMs13**2).

SetTexName([
	CKMc12='c_{12}',
	CKMc13='c_{13}',
	CKMc23='c_{23}'
]).

parameter
	Vud=CKMc12*CKMc13,
	Vus=CKMs12*CKMc13,
	Vub=CKMs13,
	Vcd=(-CKMs12*CKMc23-CKMc12*CKMs23*CKMs13),
	Vcs=(CKMc12*CKMc23-CKMs12*CKMs23*CKMs13),
	Vcb=CKMs23*CKMc13,
	Vtd=(CKMs12*CKMs23-CKMc12*CKMc23*CKMs13),
	Vts=(-CKMc12*CKMs23-CKMs12*CKMc23*CKMs13),
	Vtb=CKMc23*CKMc13.

SetTexName([
	Vud='V_{ud}',
	Vus='V_{us}',
	Vub='V_{ub}',
	Vcd='V_{ci}',
	Vcs='V_{cs}',
	Vcb='V_{cb}',
	Vtd='V_{ti}',
	Vts='V_{ts}',
	Vtb='V_{tb}'
]).

OrthMatrix({{Vud,Vus,Vub},{Vcd,Vcs,Vcb},{Vtd,Vts,Vtb}}).

/* fields */
/**********/

/* sm */

scalar
	h0/h0:		(Higgs, mass Mh0, texname 'h_0', atexname 'h_0').

spinor
	u0/U0:		(up, color c3, mass Mu0, texname 'u_0', atexname '\\bar u_0'),
	d0/D0:		(down, color c3, mass Md0, texname 'd_0^i', atexname '\\bar d_0^i'),
	c0/C0:		(charm, color c3, mass Mc0, texname 'c_0', atexname '\\bar c _0'),
	s0/S0:		(strange, color c3, mass Ms0, texname 's_0', atexname '\\bar s _0'),
	t0/T0:		(top, color c3, mass Mt0, texname 't_0', atexname '\\bar t_0'),
	b0/B0:		(bottom, color c3, mass Mb0, texname 'b_0', atexname '\\bar b_0').

vector
	G0/G0:		(gluon, color c8, gauge, texname 'g_0', atexname 'g_0'),
	A0/A0:		(photon, gauge, texname '\\gamma_0', atexname '\\gamma_0'),
	Z0/Z0:		(Z, mass MZ0, gauge, texname 'Z_0', atexname 'Z_0'),
	Wp0/Wm0:	(W, mass MW0, gauge, texname 'W_0^+', atexname 'W_0^-').

/* kk */

_n=1-nmax in scalar
	hn_n/hn_n:	(Higgs_n, mass Mhn_n, texname 'h_{n}', atexname 'h_{n}'),
	a0n_n/a0n_n:	(neutral_n, mass Ma0n_n, texname 'a_{n}', atexname 'a_{n}'),
	apn_n/amn_n:	(charged_n, mass Macn_n, texname 'a_{n}^+', atexname 'a_{n}^-').

_n=1-nmax in spinor
	u1n_n/U1n_n:	(up1_n, color c3, mass Mu1n_n, texname 'u_{n}^1', atexname '\\bar u_{n}^1'),
	u2n_n/U2n_n:	(up2_n, color c3, mass Mu2n_n, texname 'u_{n}^2', atexname '\\bar u_{n}^2'),
	d1n_n/D1n_n:	(down1_n, color c3, mass Md1n_n, texname 'd_{n}^{1i}', atexname '\\bar d_{n}^{1i}'),
	d2n_n/D2n_n:	(down2_n, color c3, mass Md2n_n, texname 'd_{n}^{2i}', atexname '\\bar d_{n}^{2i}'),
	c1n_n/C1n_n:	(charm1_n, color c3, mass Mc1n_n, texname 'c_{n}^1', atexname '\\bar c_{n}^1'),
	c2n_n/C2n_n:	(charm2_n, color c3, mass Mc2n_n, texname 'c_{n}^2', atexname '\\bar c_{n}^2'),
	s1n_n/S1n_n:	(strange1_n, color c3, mass Ms1n_n, texname 's_{n}^1', atexname '\\bar s_{n}^1'),
	s2n_n/S2n_n:	(strange2_n, color c3, mass Ms2n_n, texname 's_{n}^2', atexname '\\bar s_{n}^2'),
	t1n_n/T1n_n:	(top1_n, color c3, mass Mt1n_n, texname 't_{n}^1', atexname '\\bar t_{n}^1'),
	t2n_n/T2n_n:	(top2_n, color c3, mass Mt2n_n, texname 't_{n}^2', atexname '\\bar t_{n}^2'),
	b1n_n/B1n_n:	(bottom1_n, color c3, mass Mb1n_n, texname 'b_{n}^1', atexname '\\bar b_{n}^1'),
	b2n_n/B2n_n:	(bottom2_n, color c3, mass Mb2n_n, texname 'b_{n}^2', atexname '\\bar b_{n}^2').

_n=1-nmax in vector
	Gn_n/Gn_n:	(gluon_n, color c8, mass MGn_n, gauge, texname 'g_{n}', atexname 'g_{n}'),
	An_n/An_n:	(photon_n, mass MAn_n, gauge, texname '\\gamma_{n}', atexname '\\gamma_{n}'),
	Zn_n/Zn_n:	(Z_n, mass MZn_n, gauge, texname 'Z_{n}', atexname 'Z_{n}'),
	Wpn_n/Wmn_n:	(W_n, mass MWn_n, gauge, texname 'W_{n}^+', atexname 'W_{n}^-').

/* physical fields */
/******************/

/* scalars */

let
	phip0='Wp0.f',
	phim0='Wm0.f',
	phi30='Z0.f'.

SetTexName([
	'Wp0.f'='\\phi_0^+',
	'Wm0.f'='\\phi_0^-',
	'Z0.f'='\\phi_0^3'
]).

_n=1-nmax in let
	G5n_n='Gn_n.f',
	GAn_n='An_n.f',
	GZn_n='Zn_n.f',
	Gpn_n='Wpn_n.f',
	Gmn_n='Wmn_n.f'.

SetTexName([
	'Gn1.f'='G_{5n}',
	'An1.f'='G_{An}',
	'Zn1.f'='G_{Zn}',
	'Wpn1.f'='G_n^+',
	'Wmn1.f'='G_n^-'
]).

_n=1-nmax in let
	Wp5n_n=((_n/(n/KKmn))*Gpn_n+(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*apn_n)/MWn_n,
	Wm5n_n=((_n/(n/KKmn))*Gmn_n+(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*amn_n)/MWn_n,
	phipn_n=((((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*Gpn_n-(_n/(n/KKmn))*apn_n)/MWn_n,
	phimn_n=((((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*Gmn_n-(_n/(n/KKmn))*amn_n)/MWn_n,
	W15n_n=(Wp5n_n+Wm5n_n)/Sqrt2,
	W25n_n=(-Wp5n_n+Wm5n_n)/Sqrt2/i.

_n=1-nmax in let
	G1n_n=-(MZ0*sw)*(_n/(n/KKmn))/(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)/MWn_n*GAn_n-(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*MZn_n/(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)/MWn_n*GZn_n,
	G2n_n=(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*MZn_n/(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)/MWn_n*GAn_n-(MZ0*sw)*(_n/(n/KKmn))/(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)/MWn_n*GZn_n,
	W35n_n=-(_n/(n/KKmn))/MWn_n*G1n_n+(MZ0*sw)*(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw/MWn_n/MZn_n*G2n_n+(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw/MZn_n*a0n_n,
	BB5n_n=MWn_n/MZn_n*G2n_n-(MZ0*sw)/MZn_n*a0n_n,
	phi3n_n=(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw/MWn_n*G1n_n+(MZ0*sw)*(_n/(n/KKmn))/MWn_n/MZn_n*G2n_n+(_n/(n/KKmn))/MZn_n*a0n_n.

_x=[h,phip,phim,phi3] in let
	_x=(_x0)*cos(0)+Sqrt2*((_xn1)*cos1+(_xn2)*cos2+(_xn3)*cos3+(_xn4)*cos4+(_xn5)*cos5+(_xn6)*cos6).

_x=[W1,W2,W3,Wp,Wm,BB,G] in let
	_x5=Sqrt2*((_x5n1)*sin1+(_x5n2)*sin2+(_x5n3)*sin3+(_x5n4)*sin4+(_x5n5)*sin5+(_x5n6)*sin6).
	
let
	WW5={W15,W25,W35}.

/* spinors */

_n=1-nmax, _x=[u,d,c,s,t,b] in parameter
        f_x_n=atan(M_x0/(_n/(n/KKmn)))/2,
        cf_x_n=cos(f_x_n),
        sf_x_n=sin(f_x_n).

SetTexName([
	cfu1='\\cos (\\varphi_{u1})',
	cfd1='\\cos (\\varphi_{d1})',
	cfc1='\\cos (\\varphi_{c1})',
	cfs1='\\cos (\\varphi_{s1})',
	cft1='\\cos (\\varphi_{t1})',
	cfb1='\\cos (\\varphi_{b1})',
	sfu1='\\sin (\\varphi_{u1})',
	sfd1='\\sin (\\varphi_{d1})',
	sfc1='\\sin (\\varphi_{c1})',
	sfs1='\\sin (\\varphi_{s1})',
	sft1='\\sin (\\varphi_{t1})',
	sfb1='\\sin (\\varphi_{b1})'
]).

_n=1-nmax, _x=[u,d,c,s,t,b] in angle
        sin=sf_x_n,
        cos=cf_x_n.

_n=1-nmax, _x=[u,d,c,s,t,b] in let
	_xRn_n=-gamma5*cf_x_n*(_x1n_n)+sf_x_n*(_x2n_n),
	_xLn_n=gamma5*sf_x_n*(_x1n_n)+cf_x_n*(_x2n_n).

_x=[u,d,c,s,t,b] in let
        _xLL=PL*((_xLn1)*cos1+(_xLn2)*cos2+(_xLn3)*cos3+(_xLn4)*cos4+(_xLn5)*cos5+(_xLn6)*cos6),
	_xLR=PR*((_xLn1)*sin1+(_xLn2)*sin2+(_xLn3)*sin3+(_xLn4)*sin4+(_xLn5)*sin5+(_xLn6)*sin6),
	_xRR=PR*((_xRn1)*cos1+(_xRn2)*cos2+(_xRn3)*cos3+(_xRn4)*cos4+(_xRn5)*cos5+(_xRn6)*cos6),
        _xRL=PL*((_xRn1)*sin1+(_xRn2)*sin2+(_xRn3)*sin3+(_xRn4)*sin4+(_xRn5)*sin5+(_xRn6)*sin6).

_x=[u,d,c,s,t,b] in let
        _xL=PL*(_x0)*cos(0)+Sqrt2*((_xLL)+(_xLR)),
        _xR=PR*(_x0)*cos(0)+Sqrt2*((_xRR)+(_xRL)).

/* vectors */

let
	W10=(Wp0+Wm0)/Sqrt2,
	W20=(-Wp0+Wm0)/Sqrt2/i,
	W30=sw*A0+cw*Z0,
	BB0=cw*A0-sw*Z0.

_n=1-nmax in let
	W1n_n=(Wpn_n+Wmn_n)/Sqrt2,
	W2n_n=(-Wpn_n+Wmn_n)/Sqrt2/i,
	W3n_n=sw*An_n+cw*Zn_n,
	BBn_n=cw*An_n-sw*Zn_n.

_x=[G,W1,W2,W3,BB,Wp,Wm] in let
	_x=(_x0)*cos(0)+Sqrt2*((_xn1)*cos1+(_xn2)*cos2+(_xn3)*cos3+(_xn4)*cos4+(_xn5)*cos5+(_xn6)*cos6).

let
	WW={W1,W2,W3}.

/* lagrangian */
/**************/

/* higgs sector */

let
	phi={i*phip,(h+(2*(MZ0*cw)/(ge/sw))+i*phi3)/Sqrt2},
	Phi={-i*phim,(h+(2*(MZ0*cw)/(ge/sw))-i*phi3)/Sqrt2},
	phi2=Phi*phi.

let
	Dphi={i*deriv*phip+(i*(ge/sw)/2*W3+i*(ge/cw)/2*BB)*(i*phip)+i*(ge/sw)/Sqrt2*Wp*(h+(2*(MZ0*cw)/(ge/sw))+i*phi3)/Sqrt2,i*(ge/sw)/Sqrt2*Wm*(i*phip)+deriv*h/Sqrt2+i*deriv*phi3/Sqrt2+(-i*(ge/sw)/2*W3+i*(ge/cw)/2*BB)*(h+(2*(MZ0*cw)/(ge/sw))+i*phi3)/Sqrt2},
	D5phi={i*deriv5/(n/KKmn)*phip+(i*(ge/sw)/2*W35+i*(ge/cw)/2*BB5)*(i*phip)+i*(ge/sw)/Sqrt2*Wp5*(h+(2*(MZ0*cw)/(ge/sw))+i*phi3)/Sqrt2,i*(ge/sw)/Sqrt2*Wm5*(i*phip)+deriv5/(n/KKmn)*h/Sqrt2+i*deriv5/(n/KKmn)*phi3/Sqrt2+(-i*(ge/sw)/2*W35+i*(ge/cw)/2*BB5)*(h+(2*(MZ0*cw)/(ge/sw))+i*phi3)/Sqrt2},
	DPhi={-i*deriv*phim+(-i*(ge/sw)/2*W3-i*(ge/cw)/2*BB)*(-i*phim)-i*(ge/sw)/Sqrt2*Wm*(h+(2*(MZ0*cw)/(ge/sw))-i*phi3)/Sqrt2,-i*(ge/sw)/Sqrt2*Wp*(-i*phim)+deriv*h/Sqrt2-i*deriv*phi3/Sqrt2+(i*(ge/sw)/2*W3-i*(ge/cw)/2*BB)*(h+(2*(MZ0*cw)/(ge/sw))-i*phi3)/Sqrt2},
	D5Phi={-i*deriv5/(n/KKmn)*phim+(-i*(ge/sw)/2*W35-i*(ge/cw)/2*BB5)*(-i*phim)-i*(ge/sw)/Sqrt2*Wm5*(h+(2*(MZ0*cw)/(ge/sw))-i*phi3)/Sqrt2,-i*(ge/sw)/Sqrt2*Wp5*(-i*phim)+deriv5/(n/KKmn)*h/Sqrt2-i*deriv5/(n/KKmn)*phi3/Sqrt2+(i*(ge/sw)/2*W35-i*(ge/cw)/2*BB5)*(h+(2*(MZ0*cw)/(ge/sw))-i*phi3)/Sqrt2}.

lterm
	DPhi*Dphi-D5Phi*D5phi+muH**2*phi2-lamH*phi2**2.

/* gauge sector */

let
	FB^mu^nu=deriv^mu*BB^nu-deriv^nu*BB^mu,
	FB5^mu=deriv^mu*BB5-deriv5/(n/KKmn)*BB^mu,
	FW^mu^nu^a=deriv^mu*WW^nu^a-deriv^nu*WW^mu^a-(ge/sw)*eps^a^b^c*WW^mu^b*WW^nu^c,
	FW5^mu^a=deriv^mu*WW5^a-deriv5/(n/KKmn)*WW^mu^a-(ge/sw)*eps^a^b^c*WW^mu^b*WW5^c,
	FG^mu^nu^a=deriv^mu*G^nu^a-deriv^nu*G^mu^a-gc*f_SU3^a^b^c*G^mu^b*G^nu^c,
	FG5^mu^a=deriv^mu*G5^a-deriv5/(n/KKmn)*G^mu^a-gc*f_SU3^a^b^c*G^mu^b*G5^c.

lterm
	-1/4*FG**2+1/2*FG5**2-1/4*FW**2+1/2*FW5**2-1/4*FB**2+1/2*FB5**2.

lterm
	-1/2*(deriv*G-deriv5/(n/KKmn)*G5)**2-(deriv*Wp-deriv5/(n/KKmn)*Wp5-(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*phip)*(deriv*Wm-deriv5/(n/KKmn)*Wm5-(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*phim)-1/2*(deriv*W3-deriv5/(n/KKmn)*W35+(((ge/sw)*(2*(MZ0*cw)/(ge/sw))/2)/cw)*cw*phi3)**2-1/2*(deriv*BB-deriv5/(n/KKmn)*BB5-(MZ0*sw)*phi3)**2.

/* ghost sector */

let
	ghG='G0.c'*cos(0),
	GhG='G0.C'*cos(0).

lterm
	i*gc*f_SU3*deriv*GhG*G*ghG.

/* fermion sector */

_x=[L,R] in let 
	da_x=Vud*d_x+Vus*s_x+Vub*b_x,
	sa_x=Vcd*d_x+Vcs*s_x+Vcb*b_x,
	ba_x=Vtd*d_x+Vts*s_x+Vtb*b_x.

lterm % kinetic
	anti(psiL)*i*gamma*deriv*psiL+anti(psiR)*i*gamma*deriv*psiR
		where
			psiL=uL, psiR=uR;
			psiL=dL, psiR=dR;
			psiL=cL, psiR=cR;
			psiL=sL, psiR=sR;
			psiL=tL, psiR=tR;
			psiL=bL, psiR=bR.

lterm % kinetic-5
	-anti(psiL)*gamma5*deriv5/(n/KKmn)*psiL-anti(psiR)*gamma5*deriv5/(n/KKmn)*psiR
		where
			psiL=uL, psiR=uR;
			psiL=dL, psiR=dR;
			psiL=cL, psiR=cR;
			psiL=sL, psiR=sR;
			psiL=tL, psiR=tR;
			psiL=bL, psiR=bR.

lterm % u(1)
	anti(psiL)*i*gamma*i*(ge/cw)/2*YL*BB*psiL+anti(psiR)*i*gamma*i*(ge/cw)/2*YR*BB*psiR
		where
			psiL=uL, YL=1/3, psiR=uR, YR=4/3;
			psiL=cL, YL=1/3, psiR=cR, YR=4/3;
			psiL=tL, YL=1/3, psiR=tR, YR=4/3;
			psiL=dL, YL=1/3, psiR=dR, YR=-2/3;
			psiL=sL, YL=1/3, psiR=sR, YR=-2/3;
			psiL=bL, YL=1/3, psiR=bR, YR=-2/3.

lterm % u(1)-5
	-anti(psiL)*gamma5*i*(ge/cw)/2*YL*BB5*psiL-anti(psiR)*gamma5*i*(ge/cw)/2*YR*BB5*psiR
		where
			psiL=uL, YL=1/3, psiR=uR, YR=4/3;
			psiL=cL, YL=1/3, psiR=cR, YR=4/3;
			psiL=tL, YL=1/3, psiR=tR, YR=4/3;
			psiL=dL, YL=1/3, psiR=dR, YR=-2/3;
			psiL=sL, YL=1/3, psiR=sR, YR=-2/3;
			psiL=bL, YL=1/3, psiR=bR, YR=-2/3.

lterm % neutral component of su(2)
	anti(psi1L)*i*gamma*i*(ge/sw)/2*W3*psi1L-anti(psi2L)*i*gamma*i*(ge/sw)/2*W3*psi2L
		where
			psi1L=uL, psi2L=dL;
			psi1L=cL, psi2L=sL;
			psi1L=tL, psi2L=bL.

lterm % neutral component of su(2)-5
	-anti(psi1L)*gamma5*i*(ge/sw)/2*W35*psi1L+anti(psi2L)*gamma5*i*(ge/sw)/2*W35*psi2L
		where
			psi1L=uL, psi2L=dL;
			psi1L=cL, psi2L=sL;
			psi1L=tL, psi2L=bL.

lterm % charged component of su(2)
	anti(psi1L)*i*gamma*i*(ge/sw)/Sqrt2*Wp*psi2L+anti(psi2L)*i*gamma*i*(ge/sw)/Sqrt2*Wm*psi1L
		where
			psi1L=uL, psi2L=daL;
			psi1L=cL, psi2L=saL;
			psi1L=tL, psi2L=baL.

lterm % charged component of su(2)-5
	-anti(psi1L)*gamma5*i*(ge/sw)/Sqrt2*Wp5*psi2L-anti(psi2L)*gamma5*i*(ge/sw)/Sqrt2*Wm5*psi1L
		where
			psi1L=uL, psi2L=daL;
			psi1L=cL, psi2L=saL;
			psi1L=tL, psi2L=baL.

lterm % su(3)
	anti(psiL)*i*gamma*i*gc*lambda*G*psiL+anti(psiR)*i*gamma*i*gc*lambda*G*psiR
		where
			psiL=uL, psiR=uR;
			psiL=dL, psiR=dR;
			psiL=cL, psiR=cR;
			psiL=sL, psiR=sR;
			psiL=tL, psiR=tR;
			psiL=bL, psiR=bR.

lterm % su(3)-5
	-anti(psiL)*gamma5*i*gc*lambda*G5*psiL-anti(psiR)*gamma5*i*gc*lambda*G5*psiR
		where
			psiL=uL, psiR=uR;
			psiL=dL, psiR=dR;
			psiL=cL, psiR=cR;
			psiL=sL, psiR=sR;
			psiL=tL, psiR=tR;
			psiL=bL, psiR=bR.

/* yukawa sector */

let
	q1a={uL,daL},
	q2a={cL,saL},
	q3a={tL,baL}.

let
	H={i*phip,(h+(2*(MZ0*cw)/(ge/sw))+i*phi3)/Sqrt2},
	Hcc={-i*phim,(h+(2*(MZ0*cw)/(ge/sw))-i*phi3)/Sqrt2}.

lterm   
	-M2/MW0/Sqrt2*(ge/sw)*(anti(pl)*pr*H + anti(pr)*pl*Hcc )
                where
                        M2=Vud*Md0, pl=q1a, pr=dR;
                        M2=Vus*Ms0, pl=q1a, pr=sR;
                        M2=Vub*Mb0, pl=q1a, pr=bR;
                        M2=Vcd*Md0, pl=q2a, pr=dR;
                        M2=Vcs*Ms0, pl=q2a, pr=sR;
                        M2=Vcb*Mb0, pl=q2a, pr=bR;
                        M2=Vtd*Md0, pl=q3a, pr=dR;
                        M2=Vts*Ms0, pl=q3a, pr=sR;
                        M2=Vtb*Mb0, pl=q3a, pr=bR.

lterm   
	-M1/MW0/Sqrt2*(ge/sw)*(anti(pl)*i*tau2*pr*Hcc + anti(pr)*i*pl*tau2*H )
                where
                        M1=Mu0, pl=q1a, pr=uR;
                        M1=Mc0, pl=q2a, pr=cR;
                        M1=Mt0, pl=q3a, pr=tR.

/* conclusion */
/**************/

SetAngle(1-cw**2=sw**2).

SetAngle(1-CKMc12**2=CKMs12**2).
SetAngle(1-CKMc13**2=CKMs13**2).
SetAngle(1-CKMc23**2=CKMs23**2).

SetAngle(1-cfu1**2=sfu1**2).
SetAngle(1-cfu2**2=sfu2**2).
SetAngle(1-cfu3**2=sfu3**2).
SetAngle(1-cfu4**2=sfu4**2).
SetAngle(1-cfu5**2=sfu5**2).
SetAngle(1-cfu6**2=sfu6**2).

SetAngle(1-cfd1**2=sfd1**2).
SetAngle(1-cfd2**2=sfd2**2).
SetAngle(1-cfd3**2=sfd3**2).
SetAngle(1-cfd4**2=sfd4**2).
SetAngle(1-cfd5**2=sfd5**2).
SetAngle(1-cfd6**2=sfd6**2).

SetAngle(1-cfc1**2=sfc1**2).
SetAngle(1-cfc2**2=sfc2**2).
SetAngle(1-cfc3**2=sfc3**2).
SetAngle(1-cfc4**2=sfc4**2).
SetAngle(1-cfc5**2=sfc5**2).
SetAngle(1-cfc6**2=sfc6**2).

SetAngle(1-cfs1**2=sfs1**2).
SetAngle(1-cfs2**2=sfs2**2).
SetAngle(1-cfs3**2=sfs3**2).
SetAngle(1-cfs4**2=sfs4**2).
SetAngle(1-cfs5**2=sfs5**2).
SetAngle(1-cfs6**2=sfs6**2).

SetAngle(1-cft1**2=sft1**2).
SetAngle(1-cft2**2=sft2**2).
SetAngle(1-cft3**2=sft3**2).
SetAngle(1-cft4**2=sft4**2).
SetAngle(1-cft5**2=sft5**2).
SetAngle(1-cft6**2=sft6**2).

SetAngle(1-cfb1**2=sfb1**2).
SetAngle(1-cfb2**2=sfb2**2).
SetAngle(1-cfb3**2=sfb3**2).
SetAngle(1-cfb4**2=sfb4**2).
SetAngle(1-cfb5**2=sfb5**2).
SetAngle(1-cfb6**2=sfb6**2).

CheckHerm.
