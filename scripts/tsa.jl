using tele
using NCDatasets
data_dir = "/home/brynn/Documents/JP/860/tsa_data/"
ersst_p = "sst.mnmean.nc"
nc = NCDataset(data_dir * ersst_p)
