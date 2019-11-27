load("./symbolic_interpolation.mc");

s1:make_random_state(true)$set_random_state(s1)$
num_samples:1000$
array (x, fixnum, num_samples)$
array (y, fixnum, num_samples)$
array (z, fixnum, num_samples)$


for i:0 thru num_samples step 1 do
    block(
        x[i]:random(5.0)-2.5,
        y[i]:random(5.0)-2.5,
        z[i]:random(5.0)-2.5
    )$



for i:1 thru num_samples step 1 do
    block(
        res:W(x[i],y[i],z[i]),
        neg_res:W(-x[i],-y[i],-z[i]),
        if not is(equal(res,neg_res)) then
            block(
            print("Function is not symmetric!!"),
            print("Coordinates: ",x[i],y[i],z[i]," Result: ",res," Result with negative args: ",neg_res," Iteration: ",i),
            exit
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
            exit
            )
 
    )$
print("Function is compact")$
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
        exit
    )
)$
print("Function is normalized")$
exit(1);
