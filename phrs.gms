$ontext
    File name: phrs.gms
    Author: Mohammed Alfaki, December, 2011.
    GAMS model for 't' heuristic with P-formulation.
$offtext

#===============================================================================
# Declare options
#===============================================================================
options optcr=0.001, limrow=0, limcol=0;

#===============================================================================
# Declare sets and parameters
#===============================================================================
set sl(i), lt(i), l(i);
lt(i) = i(i)-s(i); sl(i) = i(i)-t(i);
l(i)  = i(i)-s(i)-t(i);

$include terminal.gms

# Is there a path between (s,l)
parameter p(i,j);
#===============================================================================
# Compute p(s,i) using the Breadth-first-search (BFS)
#===============================================================================
alias (i,ii,iii), (j,jj), (l,h);
scalar u,path,n,mn;
set openlst(i);
set closedlst(i);
parameter val(i);
loop((ii,jj)$(ord(ii)<=ord(jj)),
  n = 1;
  path = 0;
  val(i) = +inf;
  openlst(i) = no;
  closedlst(i) = no;
  openlst(i)$(ord(i) eq ord(ii)) = yes;
  val(i)$(ord(i) eq ord(ii)) = n;
  while((card(openlst)>0 and path eq 0),
    mn = smin(openlst(i), val(i));
    loop(i$openlst(i),
      if((mn eq val(i)),
        u = ord(i);
        val(i) = +inf;
      );
    );
    openlst(i)$(ord(i) eq u) = no;
    if((ord(jj) eq u),
      path = 1;
    else
      closedlst(i)$(ord(i) eq u) = yes;
      loop((iii,j)$(ord(iii) eq u and a(iii,j)>0 and
                                         not closedlst(j) and not openlst(j)),
        n = n+1;
        val(j) = n;
        openlst(i)$(ord(i) eq ord(j)) = yes;
      );
    );
  );
  p(ii,jj) = path;
);

#===============================================================================
# Declare variables/bounds
#===============================================================================
variable cost, w(i,k);
positive variables f(i,j);

w.lo(l,k) = smin(s$(p(s,l)>0),q(s,k));
w.up(l,k) = smax(s$(p(s,l)>0),q(s,k));
f.up(i,j) = min(bu(i),bu(j))*a(i,j);

#===============================================================================
# Declare constraints
#===============================================================================
equations obj(t), sflowcaplb(s,t), sflowcapub(s,t), tflowcaplb(t),
          ptflowcapub(i,t), flowmassblnc(i,t), qualblnc(i,k,t), qualub(t,k),
          qualubTp(t,k);

#===============================================================================
# Define constraints
#===============================================================================
#-----------------------------Objective function--------------------------------
obj(tau).. cost =e= sum((i,j)$(a(i,j)*p(j,tau)>0), c(i,j)*f(i,j));
#-------------------------Raw material availabilities---------------------------
sflowcaplb(s,tau)$(bl(s)*p(s,tau)>0)..
                           sum(j$(a(s,j)*p(j,tau)>0), Fp(s,j)+f(s,j)) =g= bl(s);
sflowcapub(s,tau)$(bu(s)<+inf and p(s,tau)>0)..
                           sum(j$(a(s,j)*p(j,tau)>0), Fp(s,j)+f(s,j)) =l= bu(s);
#----------------Pool capacities and product demand restrictions----------------
tflowcaplb(tau)$(bl(tau)>0)..
                            sum(j$(a(j,tau)>0), Fp(j,tau)+f(j,tau)) =g= bl(tau);
ptflowcapub(lt,tau)$(bu(lt)<+inf and p(lt,tau)>0)..
                                sum(j$(a(j,lt)>0), Fp(j,lt)+f(j,lt)) =l= bu(lt);
#------------------------Flow mass balance around pools-------------------------
flowmassblnc(l,tau)$(p(l,tau)>0)..   sum(s$(a(s,l)>0), f(s,l)) - f(l,tau) =e= 0;
#-----------------------Commodity flow balance around pools---------------------
qualblnc(l,k,tau)$(a(l,tau)>0)..
                            w(l,k)*(f(l,tau) + sum(Tp$(a(l,Tp)>0), Fp(l,Tp))) -
                               sum(s$(a(s,l)>0), q(s,k)*(Fp(s,l)+f(s,l))) =e= 0;
#------------------------Product quality specifications-------------------------
qualub(tau,k)$(abs(q(tau,k))<+inf)..
                              sum(s$(a(s,tau)>0), (q(s,k)-q(tau,k))*f(s,tau)) +
                          sum(l$(a(l,tau)>0), (w(l,k)-q(tau,k))*f(l,tau)) =l= 0;
qualubTp(Tp,k)$(abs(q(Tp,k))<+inf)..
                                sum(s$(a(s,Tp)>0), (q(s,k)-q(Tp,k))*Fp(s,Tp)) +
                            sum(l$(a(l,Tp)>0), (w(l,k)-q(Tp,k))*Fp(l,Tp)) =l= 0;

#===============================================================================
# Solve the model
#===============================================================================
option nlp = baron;
model hrsnlp /obj, sflowcaplb, sflowcapub, tflowcaplb, ptflowcapub,
              flowmassblnc, qualblnc, qualub, qualubTp/;
hrsnlp.solprint = 2;
hrsnlp.tolinfeas = 1.e-3;
solve hrsnlp minimizing cost using nlp;
set v /1*4/;
parameter mdlstat, f_p(i,j), psiz(v);
mdlstat = hrsnlp.modelstat;
f_p(i,j) = 0.0;
loop(tau, f_p(i,j)$(a(i,j)*p(j,tau)>0) = f.l(i,j););
psiz('1') = hrsnlp.numvar-1;
psiz('2') = 0;
loop((s,l,tau)$(p(s,l)*p(l,tau)>0), psiz('2') = psiz('2')+1;);
psiz('3') = psiz('2');
loop((l,tau)$(a(l,tau)>0), psiz('3') = psiz('3')+1;);
psiz('3') = psiz('3')+card(k);
psiz('4') = hrsnlp.numequ-psiz('3')-1;
execute_unload 'results.gdx',cost.l=z_p,mdlstat,f_p,psiz;
#========================================END====================================