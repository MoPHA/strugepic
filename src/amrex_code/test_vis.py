import yt
print(yt.__version__)
ds=yt.load("plt_P00/test_particle")

print(ds.field_list)
