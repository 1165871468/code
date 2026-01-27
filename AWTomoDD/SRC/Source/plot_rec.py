import numpy as np
import time
import matplotlib.pyplot as plt

def Max_min_normal(arr):
    min_v = np.min(arr)

    max_v = np.max(arr)

    normal_arr = (arr - min_v) / (max_v - min_v)

    return normal_arr

def Max_min_normal2(min_v, max_v, arr):
    normal_arr = (arr - min_v) / (max_v - min_v)

    return normal_arr

if __name__ == '__main__':
    Ns = 5; num_rec = 5
    
    for j in range(5, Ns+1):
        path1 = f'./D-o.txt'
        path2 = f'./D-o.txt-t'
        
        path1 = f'./D5-o'
        path2 = f'./D5-o-t'
        # path2 = path1

        for i in range(5, num_rec + 1):

                filnm1 = f'{path1}'
                filnm2 = f'{path2}'

                data1 = np.loadtxt(filnm1)
                data2 = np.loadtxt(filnm2)
        #        data3 = np.loadtxt(filnm3)

                x1 = data1[:, 0]
                y1 = data1[:, 1]
                z1 = data1[:, 2]

                x2 = data2[:, 0]
                y2 = data2[:, 1]
                z2 = data2[:, 2]
                
                numData = len(data1)

                plt.figure(i)
                plt.plot(x1, z1, color = 'green', label='Obs', linewidth = 1, linestyle='-')
                plt.plot(x2, z2, color = 'red', label='Syn', linewidth = 1, linestyle='--')
        #        plt.plot(x3, y3, color = 'black', linewidth = 1, linestyle='--')
                xx0 = 0
                xx1 = x1[numData-1]
        
                plt.xticks(np.arange(xx0, xx1, 1))
                plt.legend()

                p_name = f's{j:0d}_r{i:0d}.png'
                plt.savefig(p_name, dpi=150)
                plt.close(i)
