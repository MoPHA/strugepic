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
    if abs(x) < 1 then
        1-abs(x)
    else
        0
)$

I_W1(a,b):=quad_qag( W1,x,a,b,3,'epsrel=1d-12);

W1d(x):=block(
       if x < 1 and x > 0 then
       -1
       else if x<0 and x> -1  then
       1
       else
       0
)$


W12(q):=block(
    if q > 0 and q <1 then
    1
    else
    0
)$

W(x,y,z):=W1(x)*W1(y)*W1(z)$

W_basis(id,x,y,z):=W(x-id[1],y-id[2],z-id[3])$
I_W12(a,b):=quad_qag( W12(x),x,a,b,3,'epsrel=1d-6);

/*Index combinations for nonzero basis functions*/
/*
vals:[0,1];
*/
array (com, fixnum, 8)$
idx:1$
for i:1 thru 2 step 1 do
(
    for j:1 thru 2 step 1 do
    (
        for k:1 thru 2 step 1 do(
            com[idx]:[i-1,j-1,k-1],
            idx:idx+1
            )
        )
            
)$
