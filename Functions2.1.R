library(pracma)

# Functions
bSplineA2 <- function(x,df = NULL, knots = NULL, degree = 4L, intercept = FALSE,
               Boundary.knots = range(x, na.rm = TRUE))
{
  x1 = min(Boundary.knots)
  x2 = max(Boundary.knots)
  xtemp = seq(x1,x2,length.out=1e2)
  bstemp = bSpline(xtemp, df = df, knots = knots, degree = degree, intercept = intercept,
                  Boundary.knots = Boundary.knots)
  m = ncol(bstemp)
  xt = seq(x1,x2,length.out=q^degree*1e3)
  dx = xt[2]-xt[1]
  
  bsMat = bSpline(xt, df = df, knots = knots, degree = degree, intercept = intercept,
                  Boundary.knots = Boundary.knots)
  bsMat2 = deriv(bsMat,derivs = 2)

  A2 = diag(rep(1,m))
  for (j in 1:m)
  {
    for (k in j:m)
    {
      bsMat2jk = bsMat2[,j]*bsMat2[,k]
      A2[j,k] = dx*(sum(bsMat2jk)-(bsMat2jk[1]+bsMat2jk[1e5])/2)
      A2[k,j] = A2[j,k]
    }
  }
  return(A2)
}

bSplineVA2 <- function(x,df = NULL, knots = NULL, degree = 4L, intercept = FALSE,
               Boundary.knots = range(x, na.rm = TRUE))
{
  x1 = min(Boundary.knots)
  x2 = max(Boundary.knots)
  xt = seq(x1,x2,length.out=1e6)
  dx = xt[2]-x1[1]
  
  bsMat = bSpline(c(x1,x2), df = df, knots = knots, degree = degree, intercept = intercept,
                  Boundary.knots = Boundary.knots)
  m = ncol(bsMat)
  A2 = diag(rep(1,m))
  for (j in 1:m)
  {
    for (k in j:m)
    {
      fun_A2jk <- function(t)
      {
      bsMat = bSpline(t, df = df, knots = knots, degree = degree, intercept = intercept,
                  Boundary.knots = Boundary.knots)
      bsMat2 = deriv(bsMat,derivs = 2)
      bsMat2jk = bsMat2[,j]*bsMat2[,k]
      return(bsMat2jk)
      }
      #A2jk = integrate(fun_A2jk,lower=x1,upper=x2,rel.tol = 1e-6,
      #                 stop.on.error=F,subdivisions = 1000)
      #A2[j,k] = A2jk$value
      A2jk = integral(fun_A2jk,xmin=x1,xmax=x2,method='Kron')
      A2[j,k] = A2jk
      A2[k,j] = A2[j,k]
    }
  }
  return(A2)
}

bSplineVA0 <- function(x,df = NULL, knots = NULL, degree = 4L, intercept = FALSE,
                       Boundary.knots = range(x, na.rm = TRUE))
{
  x1 = min(Boundary.knots)
  x2 = max(Boundary.knots)
  xt = seq(x1,x2,length.out=1e6)
  dx = xt[2]-x1[1]
  
  bsMat = bSpline(c(x1,x2), df = df, knots = knots, degree = degree, intercept = intercept,
                  Boundary.knots = Boundary.knots)
  m = ncol(bsMat)
  A0 = diag(rep(1,m))
  for (j in 1:m)
  {
    for (k in j:m)
    {
      fun_A0jk <- function(t)
      {
        bsMat = bSpline(t, df = df, knots = knots, degree = degree, intercept = intercept,
                        Boundary.knots = Boundary.knots)
        bsMatjk = bsMat[,j]*bsMat[,k]
        return(bsMatjk)
      }
      A0jk = integrate(fun_A0jk,lower=x1,upper=x2)
      A0[j,k] = A0jk$value
      A0[k,j] = A0[j,k]
    }
  }
  return(A0)
}

Vec.chol2mat <- function(vec1)
{
  q = length(vec1)
  p = (sqrt(1+8*q)-1)/2
  cholmat = matrix(rep(0,p^2),p,p)
  ij = 1
  for (j in 1:p)
  {
    cholmat[j:p,j] = vec1[(ij):(ij+p-j)]
  }
  return(cholmat)
}

mplot <- function(x,y,data,id,col,xlim=NULL,ylim=NULL,...){
  data.id = data[,deparse(substitute(id))]
  data.x = data[,deparse(substitute(x))]
  data.y = data[,deparse(substitute(y))]
  
  col.try = tryCatch(col, error=function(cond){return(NULL)})
  
  if (is.null(col.try)){
    col = data[,deparse(substitute(col))]
  }
  
  data.idc = cbind(data.id,col)
  
  id.uni = unique(data.id)
  n.id = length(id.uni)
  
  #print(is.null(xlim))
  
  if (is.null(xlim) && is.null(ylim)){
    plot(NULL, xlim = range(data.x,na.rm = T), ylim = range(data.y,na.rm = T),...)
  }else if (is.null(xlim)){
    plot(NULL, xlim = range(data.x,na.rm = T),ylim=ylim,...)
  }else if (is.null(ylim)){
    plot(NULL, xlim = xlim, ylim = range(data.y,na.rm = T),...)
  }else{
    plot(NULL, xlim = xlim, ylim=ylim,...)
  }
  
  
  for (i in 1:n.id){
    Ii = (data.id == id.uni[i])
    xi = data.x[Ii]
    yi = data.y[Ii]
    col.i = data.idc[,2][Ii][1]
    
    lines(xi,yi,col=col.i)
  }
}


LCR_FixC_BLUP <- function(data,lcr.fit,yhat=NULL,AGE0=60,align.fun=align.fun.BS, pq = c(1,0), c.cov=2){
  DataFit = data$Data_Long
  DataFit_Uni = data$Data_Uni
  n_Uni = nrow(DataFit_Uni)
  
  # Refine the BLUE estimates
  lcr.fit.fixed = lcr.fit$coefficients$fixed
  
  sig_e = lcr.fit$sigma
  lcr.ci = intervals(lcr.fit)
  sig_temp = lcr.ci$reStruct$projid[,2]
  sig_a = sig_temp[1]
  sig_b = sig_temp[2]
  if (length(sig_temp)==3){
    cor_ab = lcr.ci$reStruct$projid[3,2]
  }else{
    cor_ab = 0
  }
  
  Sig_ab = matrix(c(sig_a^2, sig_a*sig_b*cor_ab, 
                    sig_a*sig_b*cor_ab, sig_b^2),2,2)*c.cov
  #Sig_ab = matrix(c(sig_a^2*2, sig_a*sig_b*cor_ab, 
  #                  sig_a*sig_b*cor_ab, sig_b^2*2),2,2)
  
  if (!is.null(lcr.ci$corStruct)){
    cor.fit = Initialize(corARMA(lcr.ci$corStruct[,2], form = ~ t_IntAGE0 | projid, p = pq[1], q = pq[2]),
                         data=DataFit)
    cor.all = corMatrix(cor.fit)
  }else{
    cor.all = NULL
  }
  
  
  lcr.formula = all.vars(lcr.fit$call$fixed[[2]])
  n.f = length(lcr.formula)
  
  XFit = matrix(rep(1,n_Uni),ncol=1)
  if (n.f > 1){
    c.names = c('c.(Intercept)')
    for (j in 1:(n.f-1)){
      c.names[j+1] = paste0('c.',lcr.formula[j+1])
    }
    
    X.names = lcr.formula[2:n.f]
    #print(X.names)
    XFit = cbind(XFit,DataFit_Uni[,X.names])
  }else{
    c.names = c('c')
  }
  
  chat = as.numeric(lcr.fit.fixed[c.names])
  
  t_Clk = seq(30,115,by=0.01)
  
  ab.est = cbind(rep(0,n_Uni),rep(0,n_Uni))
  for (i in 1:n_Uni){
    Ii = (DataFit$projid == DataFit_Uni$projid[i])
    yi = DataFit$y[Ii]
    ti = DataFit$t[Ii] - AGE0
    Xi = as.numeric(XFit[i,])
    
    ti = ti[!is.na(yi)]
    yi = yi[!is.na(yi)]
    
    if (isempty(yi)){
      next
    }
    
    ci = as.numeric(chat %*% Xi)
    CogClk.i = align.fun(t_Clk-AGE0,0,0,ci)
    
    yihat = yi#yhat[Ii]
    nui.1 = max(55,t_Clk[which.min(abs(CogClk.i - yihat[which.min(ti)]))])
    if (length(ti) >= 2){
      nui.2 = max(55,t_Clk[which.min(abs(CogClk.i - yihat[which.max(ti)]))])
      if (nui.1 < nui.2){
        bhat.i = -log((nui.2-nui.1)/(max(ti)-min(ti)))
      }else{
        bhat.i = 2
      }
      ahat.i = nui.1 - AGE0 - exp(-bhat.i)*min(ti)
    }else{
      bhat.i = 0
      ahat.i = nui.1 - AGE0 - min(ti)
    }
    
    if (is.null(cor.all)){
      cor.e.i = eye(length(ti))#cor.all[[i]]
    }else{
      cor.e.i = cor.all[[i]]
    }
    
    abhat.i = c(ahat.i, bhat.i)#c(0,0)#
    #print(abhat.i)
    #print(cbind(ti,yi))
    ab.i1 = optim(abhat.i, Omega_ab, 
                  c = ci, 
                  y = yi,t = ti,
                  sig_e = sig_e,Sig_ab = Sig_ab, 
                  align.fun = align.fun, cor.e = cor.e.i);
    ab.est[i,] = ab.i1$par
    #ab.i1 = optim(abhat.i, Omega_ab, 
    #              method = 'L-BFGS-B',
    #              lower = c(-100,-5),upper = c(100,5),
    #              c = ci, 
    #              y = yi,t = ti,
    #              sig_e = sig_e,Sig_ab = Sig_ab, 
    #              align.fun = align.fun, cor.e = cor.e.i);
    #ab.est[i,] = ab.i1$par
    #ab.i0 = optim(c(0,0), Omega_ab, c = ci, 
    #              y = yi,t = ti,
    #              sig_e = sig_e,Sig_ab = Sig_ab, 
    #              align.fun = align.fun, cor.e = cor.e.i);
    #print(c(ab.i1$value , ab.i0$value))
    
    #if (ab.i1$value <= ab.i0$value){
    #  ab.est[i,] = ab.i1$par
    #}else{
    #  ab.est[i,] = ab.i0$par
    #}
    
  }
  return(ab.est)
}




Omega_ab <- function(ab, c, y, t, sig_e, Sig_ab, align.fun, cor.e)
{
  a = ab[1];
  b = ab[2];
  
  m = length(y);
  
  yhat = align.fun(t,a,b,c);
  
  Omega0 = 1/(2*sig_e^2)*(t(y - yhat) %*% solve(cor.e) %*% (y - yhat)) + 
    1/2*t(ab)%*%solve(Sig_ab)%*%ab;
  return(c(Omega0))
}



CogAge_Est <- function(data,ab.est,AGE0=60){
  DataFit = data$Data_Long
  DataFit_Uni = data$Data_Uni
  n_Uni = nrow(DataFit_Uni)
  
  dt.AGE0 = 0
  nu_stack = rep(NA,nrow(DataFit))
  nu_bl = rep(NA,n_Uni)
  nu_lv = rep(NA,n_Uni)
  age_ad = rep(NA,n_Uni)
  age_mci = rep(NA,n_Uni)
  age_lv = rep(NA,n_Uni)
  cogage_ad = rep(NA,n_Uni)
  cogage_mci = rep(NA,n_Uni)
  #cogage_die = rep(NA,n_Uni)
  cogage_lv = rep(NA,n_Uni)
  y_bl = rep(NA,n_Uni)
  y_lv = rep(NA,n_Uni)
  y_ad = rep(NA,n_Uni)
  y_mci = rep(NA,n_Uni)
  for (i in 1:n_Uni)
  {
    id.i = DataFit_Uni$projid[i]
    Ii = (DataFit$projid==id.i)
    ti = DataFit$t[Ii] - AGE0 #DataFit$t_FitAGE0[Ii]
    yi = DataFit$y[Ii]
    
    
    ai = ab.est[i,1]
    bi = ab.est[i,2]
    nui = (ti)*exp(-bi)+ai+AGE0-dt.AGE0 ###################
    nu_stack[Ii] = nui
    nu_bl[i] = nui[which.min(ti)]#nui[1]
    nu_lv[i] = nui[which.max(ti)]#nui[length(nui)]
    
    if (nu_bl[i] < 30)
    {
      ki1 = (30-nu_lv[i])/(nu_lv[i]-nu_bl[i]);
      ki0 = 30 - ki1*nu_bl[i];
      bi = bi - log(ki1)
      ai = ki1*(ai+AGE0-dt.AGE0)+ki0-AGE0+dt.AGE0
      nui = (ti)*exp(-bi)+ai+AGE0-dt.AGE0 ###################
      nu_stack[Ii] = nui
      nu_bl[i] = nui[which.min(ti)]#nui[1]#
      nu_lv[i] = nui[which.max(ti)]#nui[length(nui)]#
    }
    
    if (nu_bl[i] > 110)
    {
      ki1 = (nu_lv[i]-110)/(nu_lv[i]-nu_bl[i]);
      ki0 = 110 - ki1*nu_bl[i];
      bi = bi - log(ki1)
      ai = ki1*(ai+AGE0-dt.AGE0)+ki0-AGE0+dt.AGE0
      nui = (ti)*exp(-bi)+ai+AGE0-dt.AGE0 ###################
      nu_stack[Ii] = nui
      nu_bl[i] = nui[which.min(ti)]#nui[1]#
      nu_lv[i] = nui[which.max(ti)]#nui[length(nui)]#
    }
    
    yi.na.rm = yi[!is.na(yi)]
    y_bl[i] = yi.na.rm[1]
    y_lv[i] = yi.na.rm[length(yi.na.rm)]
    ti.na.rm = ti[!is.na(yi)] 
    ti_lv = ti.na.rm[length(yi.na.rm)]
    
    cogage_ad[i] = (DataFit_Uni$age_bl[i]+DataFit_Uni$time2ad[i]-AGE0)*exp(-bi) + ai + AGE0-dt.AGE0
    cogage_mci[i] = (DataFit_Uni$age_bl[i]+DataFit_Uni$time2mci_i[i]-AGE0)*exp(-bi) + ai + AGE0-dt.AGE0
    cogage_lv[i] = ti_lv*exp(-bi) + ai + AGE0-dt.AGE0
    age_lv[i] = ti_lv + AGE0
    age_ad[i] = DataFit_Uni$age_bl[i]+DataFit_Uni$time2ad[i]
    age_mci[i] = DataFit_Uni$age_bl[i]+DataFit_Uni$time2mci_i[i]
    
    if (!is.na(age_ad[i])){
      y_ad[i] = yi.na.rm[which.min(abs(ti.na.rm - age_ad[i] + AGE0))]
    }
    if (!is.na(age_mci[i])){
      y_mci[i] = yi.na.rm[which.min(abs(ti.na.rm - age_mci[i] + AGE0))]
    }
    
    # age.die.i = DataFit_Uni$age_bl[i]+DataFit_Uni$time2die[i]-AGE0
    # cogage_die[i] = age.die.i*exp(-bi) + ai + AGE0-dt.AGE0
    # if (!is.na(cogage_die[i]) && (cogage_die[i] > 110))
    # {
    #   eb = (110 + log(cogage_die[i]-110+1) - nu_bl[i]) / (age.die.i + AGE0 - ti[which.min(ti)])
    #   if (is.na(log(eb)))
    #   {
    #     print(c(eb, cogage_die[i], age.die.i, nu_bl[i], age.die.i + AGE0 - ti[which.min(ti)]))
    #   }
    #   bi = -log(eb)
    #   ai = nu_bl[i] - eb*(ti[which.min(ti)]-AGE0) - AGE0 + dt.AGE0
    #   cogage_die[i] = age.die.i*exp(-bi) + ai + AGE0-dt.AGE0
    # }
  }
  return(list(cogage_bl = nu_bl,
              cogage_lv = cogage_lv,
              cogAge = nu_stack,
              age_ad = age_ad,
              age_mci = age_mci,
              cogage_ad = cogage_ad,
              cogage_mci = cogage_mci,
              #cogage_death = cogage_die,
              age_lv = age_lv,
              y_lv = y_lv,
              y_bl = y_bl,
              y_ad = y_ad,
              y_mci = y_mci))
}
