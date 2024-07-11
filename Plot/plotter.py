import matplotlib.pyplot as plt

def plot_results(yd):

    plt.plot(yd, '-b', label='Desired Output')
    # plt.plot(y[0][:-1], '--r', label='y1')
    plt.grid(True)
    plt.show()