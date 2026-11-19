from read_file import read_inversion_config3d
from read_file import read_forward_config


def main():

    # config = read_inversion_config3d("../../Input/parameters_inversion3d.config")
    forward_config = read_forward_config("../../Input/parameters_forward.config")

    # print(config)
    print(forward_config)

if __name__ == "__main__":
    main()