import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mathtext
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants


plt.rcParams.update({'mathtext.default': 'default',
                     'mathtext.fontset': 'stix',
                     'font.family': 'Times New Roman',
                     'font.size': 12,
                     'figure.figsize': (3, 3)})


def outlier_iqr(df):

    for i in range(len(df.columns)):

        # 列を抽出する
        col = df.iloc[:, i]

        # 四分位数
        q1 = col.describe()['25%']
        q3 = col.describe()['75%']
        iqr = q3 - q1  # 四分位範囲

        # 外れ値の基準点
        outlier_min = q1 - (iqr) * 1.5
        outlier_max = q3 + (iqr) * 1.5

        # 範囲から外れている値を除く
        col[col < outlier_min] = None
        col[col > outlier_max] = None

    return df


df = pd.read_csv("../build/output.csv", header=0, names=['no_pivot', 'pivot'])
# print(df.head(5))
n = len(df['pivot']) + 1
df = np.log10(df)
plt.figure(figsize=(11.69, 8.27))
plt.plot(range(1, n), df['no_pivot'], label='no pivot')
plt.plot(range(1, n), df['pivot'], label='pivot')
plt.title('log-transformed no pivot or pivot $norm_2$')
plt.xlabel('count')
plt.ylabel('$log_{10}$  ', rotation=0)
plt.grid()
plt.legend()
plt.show()
