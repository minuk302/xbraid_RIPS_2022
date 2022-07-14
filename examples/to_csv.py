import pandas as pd

path_to_u = '../examples/trischur-ex-14-adj.out.u.000'
read_file = pd.read_csv (path_to_u, header=None)
read_file.columns = ['u0', 'u1']
for k in range(read_file.shape[0]):
  s = read_file.iloc[k, 0]
  i = s.index(":")
  new_s = s[i+2:]
  read_file.iloc[k, 0] = new_s
read_file.to_csv (path_to_u + ".csv", index=None)
