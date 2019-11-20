/*
Eight order interpolation function, used for testing the numerical
c++ implementation.

W(X)= W1(x)*W2(y)*W3(z)

This function should be:

1. normalized 
2. symmetric W1(-x) = W1(x)
3. compact  W1(x) =0 x > 2
4. Smooth

We will test for 1-3 in this script

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

W(x,y,z):=W1(x)*W1(y)*W1(z)$

W_basis(id,x,y,z):=W(x-id[1],y-id[2],z-id[3])$

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





filename:"test_res.out"$
/* Lets generate some random points for testing*/
s1:make_random_state(true)$set_random_state(s1)$
num_samples:1000$
array (x, fixnum, num_samples)$
array (y, fixnum, num_samples)$
array (z, fixnum, num_samples)$
array (res, fixnum, num_samples)$

for i:0 thru num_samples step 1 do
    block(
        x[i]:random(5.0)-2.5,
        y[i]:random(5.0)-2.5,
        z[i]:random(5.0)-2.5
    )$
write_data(x,"data_x.out");
write_data(y,"data_y.out");
write_data(z,"data_z.out");
for i:0 thru num_samples step 1 do
    (
        res[i]:W(x[i],y[i],z[i])
    );
write_data(res,"data_val.out");


/*
Testing that the maxima implementation for the polynominal
is correct

Test are passing, so we don't run them again.
*/
/*
for i:1 thru num_samples step 1 do
    block(
        res:W(x[i],y[i],z[i]),
        neg_res:W(-x[i],-y[i],-z[i]),
        if not is(equal(res,neg_res)) then
            block(
            print("Function is not symmetric!!"),
            print("Coordinates: ",x[i],y[i],z[i]," Result: ",res," Result with negative args: ",neg_res," Iteration: ",i),
            quit()
            )
    )$
print("Function is symmetric")$

for i:1 thru num_samples step 1 do
    block(
        res:W(x[i],y[i],z[i]),
        neg_res:W(-x[i],-y[i],-z[i]),
        if ((x[i] > 2 or x[i] < -2 ) or (y[i] > 2 or y[i] < -2) or (z[i] > 2 or z[i] < -2) ) and not is(equal(res,0)) then
            block(
            print("Function is not compact!!"),
            print("Coordinates: ",x[i],y[i],z[i]," Result: ",res," Iteration: ",i),
            quit()
            )
 
    )$
print("Function is compact")$
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





for iter:1 thru num_samples step 1 do( 
    
    total:0,
    px:random(1.0),
    py:random(1.0),
    pz:random(1.0),
    
    for i:1 thru 64 step 1 do(
        total:total+W_basis(com[i],px,py,pz) 
    ),
    if abs(total -1)>0.000000001 then (
        print("Function is not normalized!!"),
        print("Coordinates: ",px,py,pz," Sum: ",total),
        quit()
    )
)$
print("Function is normalized")$
*/
