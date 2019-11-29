load("symbolic_interpolation.mc");

s1:make_random_state(true)$set_random_state(s1)$
num_samples:1000$

array (x, fixnum, num_samples)$
array (y, fixnum, num_samples)$
array (z, fixnum, num_samples)$
array (W_res, fixnum, num_samples)$
array (W1d_res,fixnum,num_samples)$
array (W12_res,fixnum,num_samples)$
array (W1_res,fixnum,num_samples)$
array (IW12_res,fixnum,num_samples)$
array (IW1_res,fixnum,num_samples)$

for i:0 thru num_samples step 1 do
    block(
        x[i]:random(5.0)-2.5,
        y[i]:random(5.0)-2.5,
        z[i]:random(5.0)-2.5
    )$
for i:0 thru num_samples step 1 do
    (
        W_res[i]:W(x[i],y[i],z[i]),
        W1_res[i]:W1(x[i]),
        W1d_res[i]:W1d(x[i]),
        W12_res[i]:W12(x[i]),
        IW12_res[i]:I_W12(x[i],y[i]),
        IW1_res[i]:I_W1(x[i],y[i])
    );

/*
Output is 1D comma separated, cause I could not figure out /bother how to do it properly in maxima
*/
write_data(x,"data_x.out");
write_data(y,"data_y.out");
write_data(z,"data_z.out");
write_data(W_res,"data_W.out");
write_data(W1d_res,"data_W1d.out");
write_data(W12_res,"data_W12.out");
write_data(W1_res,"data_W1.out");
write_data(IW12_res,"data_IW12.out");
write_data(IW1_res,"data_IW1.out");
