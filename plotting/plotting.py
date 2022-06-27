import matplotlib.pyplot as plt
import typing


def plot_graph(
    x_label: str,
    y_label: str,
    points: {str: typing.Tuple[typing.List[float], typing.List[float], bool]},
    title: typing.Optional[str] = ''
):
    """
    Plots a 2D the graph in a pre-defined manner
    :param x_label: The name of the X axis;
    :param y_label: The name of the Y axis;
    :param points:
    if bool is true - line is plotted;
    else - scattered points are plotted;
    :return: shows the plot.
    """

    legend = []

    for key, point_set in points.items():

        if point_set[2]:
            plt.plot(point_set[0], point_set[1])
        else:
            plt.scatter(point_set[0], point_set[1])

        legend.append(key)

    plt.legend(legend)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.suptitle(title)
    plt.show()
