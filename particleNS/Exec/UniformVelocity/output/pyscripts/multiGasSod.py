import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

fig = plt.figure(figsize=(10,6))

ax = plt.axes((0.1,0.1,0.9,0.9))


df = pd.read_csv('output/txt/test9/time0.txt', sep=" ", header=None)
df.columns = ["x","rho","v","p","eps","yO2","yN2","Tg","gamma"]

print(df)

# fig, ax = plt.subplots()

ax.scatter(df['x'],df['p'], label='Line1')
# ax.plot(df['x'],df['y2'], label='Line2')
# ax.plot(df['x'],df['y3'], label='Line3')
# ax.plot(df['x'],df['y4'], label='Line4')
tick_spacing = 1
ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
plt.xlabel('x')
plt.ylabel('Density')
plt.title('Test plot')
plt.legend()
plt.show()
