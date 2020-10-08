$ontext
    File name: pmodel.gms
    Author: Mohammed Alfaki, December, 2010.
    GAMS model for 't' heuristic algorithm to the pooling problem.
$offtext

#===============================================================================
# Declare options
#===============================================================================
options optcr=1.e-9, optca=1.e-6, limrow=0, limcol=0, reslim=3600;
 
#===============================================================================
# Declare sets and parameters
#===============================================================================
set sl(i), lt(i), l(i);
lt(i) = i(i)-s(i); sl(i) = i(i)-t(i);
l(i)  = i(i)-s(i)-t(i);
 
$include terminal.gms

set uT(t); uT(t) = Tp(t)+tau(t);
alias (uT,Tu)
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
variable cost;
positive variables y(i,t), f(i,j);

y.up(l,t) = p(l,t);
f.up(i,j) = min(bu(i),bu(j))*a(i,j);

#===============================================================================
# Declare constraints
#===============================================================================
equations obj(t), sflowcaplb(s,t), sflowcapub(s,t), tflowcaplb(t), 
          ptflowcapub(i,t), flowpathblnc(i,t), flowpathblncTp(i,t), qualub(t,k), 
          qualubTp(t,k), propblnc(i,t), rlt1(j,i,t);

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
#-----------------------Commodity flow balance around pools---------------------
flowpathblnc(l,tau)$(p(l,tau)>0).. 
                                 sum(j$(a(j,l)>0), (Fp(j,l)+f(j,l))*y(l,tau)) - 
                        sum(h$(a(l,h)*p(h,tau)>0), (Fp(l,h)+f(l,h))*y(h,tau)) -
                                        (Fp(l,tau)+f(l,tau))$(a(l,tau)>0) =e= 0;
flowpathblncTp(l,Tp)$(p(l,Tp)>0)..         sum(j$(a(j,l)>0), Fp(j,l)*y(l,Tp)) - 
        sum(h$(a(l,h)*p(h,Tp)>0), Fp(l,h)*y(h,Tp)) - Fp(l,Tp)$(a(l,Tp)>0) =e= 0;
#------------------------Product quality specifications-------------------------
qualub(tau,k)$(abs(q(tau,k))>0 and abs(q(tau,k))<+inf)..
                              sum(s$(a(s,tau)>0), (q(s,k)-q(tau,k))*f(s,tau)) + 
        sum((s,l)$(a(s,l)*p(l,tau)>0), (q(s,k)-q(tau,k))*f(s,l)*y(l,tau)) =l= 0;
qualubTp(Tp,k)$(abs(q(Tp,k))>0 and abs(q(Tp,k))<+inf)..
                               sum(s$(a(s,Tp)>0), (q(s,k) - q(Tp,k))*f(s,Tp)) +
           sum((s,l)$(a(s,l)*p(l,Tp)>0), (q(s,k)-q(Tp,k))*f(s,l)*y(l,Tp)) =l= 0;
#-------------------------Commodity proportion balances-------------------------
propblnc(l,uT)$(a(l,uT)>0)..             sum(Tu$(p(l,Tu)>0), y(l,Tu)) - 1 =e= 0;
#---------------------RLT for commodity proportion balances---------------------
rlt1(j,l,tau)$(a(j,l)*p(l,tau)>0).. 
                             sum(Tu$(p(l,Tu)>0), f(j,l)*y(l,Tu)) - f(j,l) =e= 0;
    
#===============================================================================
# Solve the model
#===============================================================================
option nlp = baron;
model hrsnlp /obj, sflowcaplb, sflowcapub, tflowcaplb, ptflowcapub, 
              flowpathblnc, flowpathblncTp, qualub, qualubTp, propblnc, rlt1/;
hrsnlp.solprint = 2;
hrsnlp.tolinfeas = 1.e-6;
solve hrsnlp minimizing cost using nlp;
scalars z_p, mdlstat;
z_p = cost.l;
mdlstat = hrsnlp.modelstat;
execute_unload 'results.gdx',z_p,mdlstat, f.l; 
#========================================END====================================