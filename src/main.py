from random import seed
from time import process_time

from test import distance_k_n
from test import distance_t_k


def main():
    # For reproducibility
    seed(0)

    img_path = '../img/'
    distance_k_n.run(save_path=img_path, show_plot=False) # Takes about 80 seconds
    distance_t_k.run(save_path=img_path, show_plot=False) # Takes about 100 seconds


if __name__ == '__main__':
    main()