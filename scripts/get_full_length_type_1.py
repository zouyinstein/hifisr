import pandas as pd
import sys

df_type_1 = pd.read_table(sys.argv[1], header=None)
df_type_1_fl = df_type_1[(df_type_1[5]-df_type_1[4])/df_type_1[1] > 0.99]
print("\n".join(list(df_type_1_fl[0])))
