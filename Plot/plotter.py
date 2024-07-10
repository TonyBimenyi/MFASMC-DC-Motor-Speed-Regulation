import matplotlib.pyplot as plt

def plot_results(yd):

    plt.plot(yd, '-b', label='Desired Output')
    plt.grid(True)
    plt.show()