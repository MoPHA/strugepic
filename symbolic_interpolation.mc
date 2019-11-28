/*
Eight order interpolation function, used for testing the numerical
c++ implementation.

W(X)= W1(x)*W1(y)*W1(z)

This function should be:

1. normalized 
2. symmetric W1(-x) = W1(x)
3. compact  W1(x) =0 |x| > 2
4. Smooth

Here it is a piecewise polynominal consisting of four parts.

*/

W1(x):=block(
    if x >2 then
        0
    elseif x<= 2 and x > 1 then
        float( 15/1024*x^8-15/128*x^7+49/128*x^6-21/32*x^5+35/64*x^4-x+1 )
    elseif x <= 1 and x > 0 then
        float(-15/1024*x^8-15/128*x^7+7/16*x^6-21/32*x^5+175/256*x^4-105/128*x^2+337/512)
    elseif x <= 0 and x > -1 then
        float(-15/1024*x^8+15/128*x^7+7/16*x^6+21/32*x^5+175/256*x^4-105/128*x^2+337/512)
    elseif x <= -1 and x > -2 then
        float( 15/1024*x^8+15/128*x^7+49/128*x^6+21/32*x^5+35/64*x^4+x+1  )
    else
        0
)$
/*
Derivative of W1
*/
W1d(q):=block(
    if q >2 then
        0
    elseif q<= 2 and q > 1 then(
        es:diff(15/1024*hh^8-15/128*hh^7+49/128*hh^6-21/32*hh^5+35/64*hh^4-hh+1,hh),
        float(ev(es,hh=q))
        )
    elseif q <= 1 and q > 0 then(
        es:diff(-15/1024*hh^8-15/128*hh^7+7/16*hh^6-21/32*hh^5+175/256*hh^4-105/128*hh^2+337/512,hh),
        float(ev(es,hh=q))
        )
    elseif q <= 0 and q > -1 then(
        es:diff(-15/1024*hh^8+15/128*hh^7+7/16*hh^6+21/32*hh^5+175/256*hh^4-105/128*hh^2+337/512,hh),
        float(ev(es,hh=q))
        )
    elseif q <= -1 and q > -2 then(
        es:diff(15/1024*hh^8+15/128*hh^7+49/128*hh^6+21/32*hh^5+35/64*hh^4+hh+1,hh),
        float(ev(es,hh=q))
        )
    else
        0
        
)$

W12(q):=block(
    if q > 2 then
        0
    elseif q < -1 then
        0
    else
        W1d(q)+W1d(q+1)+W1d(q+2)
)$

W(x,y,z):=W1(x)*W1(y)*W1(z)$

W_basis(id,x,y,z):=W(x-id[1],y-id[2],z-id[3])$
IW12(a,b):=quad_qag( W12(x),x,a,b,3,'epsrel=1d-6);

/*Index combinations for nonzero basis functions*/
/*
vals:[-1,0,1,2];
*/
array (com, fixnum, 64)$
idx:1$
for i:1 thru 4 step 1 do
(
    for j:1 thru 4 step 1 do
    (
        for k:1 thru 4 step 1 do(
            com[idx]:[i-2,j-2,k-2],
            idx:idx+1
            )
        )
            
)$
