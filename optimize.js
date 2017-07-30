//  the brent function here is a port of
//  Parabolic Interpolation and Brent’s Method in One Dimension, from the following url
//  http://www.aip.de/groups/soe/local/numres/bookcpdf/c10-2.pdf
//  i believe this is the same method used by R's optimize function



function optimize(f,lower,upper,tol,maximum,max_iterations) {
//   f = function to optimize
//   lower = lower limit of range to search
//   upper = upper limit of range to search
//   tol = tolerance (no point in making smaller than sqrt(epsilon)
// , or appr 10e-8 for double data type)
//   maximum = true for maximize, false for minimize
//   max_iterations = maximum number of iterations to run brent algorithm before giving up

//    lower = lower || 0;
//    upper = upper || 1;
    tol = tol || Math.pow(10,-8);
    maximum = maximum || false;
    max_iterations = max_iterations || 1000;

    var	lower_eps = 0
        , upper_eps = 0;
    	;
//    function func(x) {
//        return f(x);
//    }

    if(maximum) {
        function func(x) { return (-1*f(x)); }
    } else {
        function func(x) { return (f(x)); }
    }

    if(f(lower)) { lower_eps = 0; }
    else { lower_eps = Math.pow(10,-8); }

    if(f(upper)) { upper_eps = 0; }
    else { upper_eps = Math.pow(10,-8); }

	lower = lower + lower_eps;
    upper = upper - upper_eps;

	[x,y] =  brent(lower,upper,func,tol,max_iterations);

	return [x,y];
}


function brent(a,b,f,tol,max_iterations){
    var phi = (Math.sqrt(5) - 1)/2
        ,Cphi = 1-phi
        ,eps=Math.pow(10,-15)
        ;
    
    var e = 0
//        ,d
//       ,etemp
//        ,p
//        ,q
//        ,r
//        ,u
        ,v = b
        ,w = b
        ,x = b
//        ,fu = f(u)
        ,fv = f(v)
        ,fw = f(w)
        ,fx = f(x)
        ,xm = 0.5*(a+b)
        ,tol1 = tol*Math.abs(x)+eps
        ,tol2 = 2*tol1
        ,phi = (Math.sqrt(5) - 1)/2
	    ,Cphi = 1-phi
	    ,eps=Math.pow(10,-15)
        ;

    for (i = 0; i < max_iterations; i++) {
        xm = 0.5*(a+b);
        tol1 = tol*Math.abs(x)+eps;
        tol2 = 2*tol1;
        if (Math.abs(x-xm) <= (tol2 - 0.5*(b-a))) {     // check for convergence of solution
        //    console.log(i,"number of iterations to converge");
            xmin=x;
            return [xmin,fx];
        }

        if (Math.abs(e) > tol1) {
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if (q > 0.0) { p = -p; }
            q = Math.abs(q);
            etemp = e;
            e = d;
            if (Math.abs(p) >= Math.abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
                if(x >= xm) { e = a-x; }
                else { e = b-x; }
                d = Cphi*e;

            } else {
                d = p/q;
                u = x+d;
                if (u-a < tol2 || b-u < tol2) {
                    d = tol1 * Math.sign(xm-x);
                }
            }
        } else {
            if(x >= xm) { e = a-x; }
            else { e = b-x; }
            d = Cphi*e;
        }

        if(Math.abs(d) >= tol1) {  u = x+d; }
        else { u = x + tol1 * Math.sign(d); }

        fu=f(u);

        if(fu <= fx) {
            if(u >= x) { a = x; }
            else { b = x; }
            [v,w,x] = shift(v,w,x,u);
            [fv,fw,fx] = shift(fv,fw,fx,fu);
        } else {
            if(u < x) { a=u; }
            else { b = u; }
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }
    console.log("Too many iterations in brent");
    xmin=x;  //convergence not reached in max_iterations, return final value
    return [xmin,fx];
}


function shift(a,b,c,d) {
    a=b;
    b=c;
    c=d;
    return [a,b,c];
}

