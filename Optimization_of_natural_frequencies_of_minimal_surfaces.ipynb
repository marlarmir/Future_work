##Future work
#Optimization of natural frequencies of minimal surfaces
N=len(list); #number of basis functions
u = [u_{min}+du* u for u in range (nu+1)];
v = [v_{min}+dv* v for v in range (nv+1)];
def integrateLHS(j, k, jj, kk): 
	cst = (-jj* 2-kk* 2) * (2/\pi)**2 * du* dv;
	fct_1= [\sin(j* x) * \sin(jj* x) for x in u];
	fct_2= [\sin(k* x) * \sin(kk* x) for x in v];
	f_1xf_2 = [y_1* y_2 for y_1 in fct_1 for y_2 in fct_2];
	return cst* \sum(f_1xf_2);
def compute_{ds}(a,b):
	z_t= [x+1j* y for x in u for y in v];
	nu = [((((a[0]* x) + a[1]) * x) + a[2]* x +a[3] for x in z_t];
	de = [((((b[0]* x) + b[1]) * x) + b[2])* x +b[3] for x in z_t];
	g = [a/b if abs(b)>1.^{-13} else 0.; for (a,b) in zip(nu, de)];
	ds = [(1+abs(g*)**2)**2* (abs(f)**2) for g* in g];
	return(ds);
def integrateRHS(j, k, jj, kk, ds):
	cst = (2/\pi) * (2/\pi) * du* dv;
	fct_1= [\sin(j* x) * \sin(jj* x) for x in u];
	fct_2= [\sin(k* x) * \sin(kk* x) for x in v];
	f_1xf_2 = [y_1* y_2 for y_1 in fct_1 for y_2 in fct_2];
	assert len(ds) == len(f_1xf_2);
	return cst* sum([a* b for (a,b) in zip(f1xf2, ds)]);
def assemble_{matrix}(a, b, ifLHS=True):
	lstMtrx=[];
	IF not ifLHS:
		ds = compute_{ds}(a,b);
	for i in range(N):
		for ii in range(N):
			j=list[i][0];
			k=list[i][1];
			jj=list[ii][0];
			kk=list[ii][1];
			IF ifLHS:
			lstMtrx.append(integrateLHS(j, k ,jj, kk));
		else:
			lstMtrx.append(integrateRHS(j, k ,jj, kk, ds));
	#stack entries in matrix format
	Mat = [[lstMtrx[i+j* N] for j in range(N)] for i in range(N)];
	return Mat;
def compute_{maxEig}(a,b):
	M_1 = assemble_{matrix}(a, b, ifLHS=True); #assemble LHS matrix
	M_2 = assemble_{matrix}(a, b, ifLHS=False);#assemble RHS matrix
	w, v = lg.eig(M_1, M_2);
	lam_{idx} = w.argmax();
	lam = w[lam_{idx}];
	ev = v[:,lam_{idx}];
	return lam, ev, np.array(M_2);
def compute_{d\lambda}(\lambda , ev, a, b, B):
	prefactor = \lambda/(ev.\cdot (B.\cdot (ev)));
	u* = [0 for i in u];
	for idx, indeces in enumerate(list):
		j = indeces[0];
		k = indeces[1];
		fct= [\sin(j* x) * \sin(k* y) for (x,y) in zip(u,v)];
		u* = [a + ev[idx]* b for (a,b) in zip(u_, fct)];
	u_2 = [a* a for a in u*];
	cst = \lambda * (2/\pi) * (2/\pi) * du* dv;
	zt= [x+1j* y for x in u for y in v];
	de = [((((b[0]* x) + b[1]) * x) + b[2])* x +b[3] for x in zt];
	dgda1 = [x**3/y for (x,y) in zip(zt, de)];
	dgda2 = [x**2/y for (x,y) in zip(zt, de)];
	dgda3 = [x/y for (x,y) in zip(zt, de)];
	dgda4 = [1/y for (x,y) in zip(zt, de)];
	de2 = [y* y for y in de];
	nu_ = [-((((a[0]* x) + a[1]) * x) + a[2])* x +a[3] for x in zt];
	nu_ = [a/b for (a,b) in zip(nu_,de2)];
	dgdb1 = [a* x**3 for (a,x) in zip(nu_,zt)];
	dgdb2 = [a* x**2 for (a,x) in zip(nu_,zt)];
	dgdb3 = [a* x for (a,x) in zip(nu_,zt)];
	dgdb4 = [a for (a,x) in zip(nu_,zt)];
	dlada1 = cst* \sum([a* b for (a,b) in zip(u2, dgda1)]);
	dlada2 = cst* \sum([a* b for (a,b) in zip(u2, dgda2)]);
	dlada3 = cst* \sum([a* b for (a,b) in zip(u2, dgda3)]);
	dlada4 = cst* \sum([a* b for (a,b) in zip(u2, dgda4)]);
	dladb1 = cst* \sum([a* b for (a,b) in zip(u2, dgdb1)]);
	dladb2 = cst* \sum([a* b for (a,b) in zip(u2, dgdb2)]);
	dladb3 = cst* \sum([a* b for (a,b) in zip(u2, dgdb3)]);
	dladb4 = cst* \sum([a* b for (a,b) in zip(u2, dgdb4)]);
	d\lambda_a = [dlada1, dlada2, dlada3, dlada4];
	d\lambda_b = [dladb1, dladb2, dladb3, dladb4];
	return d\lambda_a, d\lambda_b;
for step in [i*1*10**{-2} for i in range(1,7)]:
	a = [0.1, 0, 1, 0];
	b = [0, 0, 0, 1];
	lam1, ev, B = compute_{maxEig}(a,b);
	d\lambda_a, d\lambda_b= compute_{d\lambda}(lam1, ev, a, b, B);
	a = [x - step* dx for (x,dx) in zip(a, d\lambda_a)];
	b = [x - step* dx for (x,dx) in zip(b, d\lambda_b)];
	lam2, ev, B = compute_{maxEig}(a,b);
	print([-lam1, -lam2])
