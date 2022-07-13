import typing

import matplotlib.pyplot as plt
import numpy


def plot_graph(
    x_label: str,
    y_label: str,
    points: {str: typing.Tuple[typing.List[float], typing.List[float], bool]},
    title: typing.Optional[str] = "",
):
    """
    Plots a 2D the graph in a pre-defined manner
    :param x_label: The name of the X axis;
    :param y_label: The name of the Y axis;
    :param points:
    if bool is true - line is plotted;
    else - scattered points are plotted;
    :param title: Title of the plot
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


def plot_surface(
    condition: bool,
    function: typing.Callable,
    x: typing.List[float],
    t: typing.List[float],
    p: typing.List[float],
    t_min: float,
    t_max: float,
    x_v: numpy.array,
):
    """
    Plots a 3D surface in a pre-defined manner
    :param condition: indicates if experimental points should be scattered
    :param function: Function to be plotted
    :param x: list of experimental values
    :param t: list of experimental values
    :param p: list of experimental values
    :param t_min: minimal value of t parameter to be displayed
    :param t_max: maximal value of t parameter to be displayed
    :param x_v: list of x values to be used for plotting
    :return: plots a 3D graph with or without scattered points
    """
    # TODO: docstring
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    if condition:
        ax.scatter(x, t, p, marker="o")

    t_v = numpy.linspace((t_min - 10), (t_max + 10), num=50)
    x_fit, t_fit = numpy.meshgrid(x_v, t_v)
    p_fit = numpy.array([function(x_fit[i], t_fit[i]) for i in range(len(x_fit))])
    ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

    ax.set_xlabel("First component fraction")
    ax.set_ylabel("Temperature K")
    ax.set_zlabel("Permeance")
    fig.suptitle("Fit Illustration", fontsize=10)
    plt.show()
