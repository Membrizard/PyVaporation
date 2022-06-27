import matplotlib.pyplot as plt
import typing
import numpy

# from optimizer import Measurements, PervaporationFunction


def plot_graph(
    x_label: str,
    y_label: str,
    points: {str: typing.Tuple[typing.List[float], typing.List[float], bool]},
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
    plt.show()


def plot_surface(
    x_label: typing.Optional[str] = None,
    t_label: typing.Optional[str] = None,
    p_label: typing.Optional[str] = None,
    suptitle: typing.Optional[str] = None,
):
    """
    Plots a 3D surface fit of the PervaporationFunction with or without the corresponding measurements
    :param pervaporation_function: Function to display
    :param measurements: experimental points
    :param x_label: composition axis label
    :param t_label: temperature axis label
    :param p_label: permeance axis label
    :param suptitle: title of the figure
    :return: shows the plot
    """

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    if measurements is not None:
        x = [m.x for m in measurements]
        t = [m.t for m in measurements]
        p = [m.p for m in measurements]
        ax.scatter(x, t, p, marker="o")
    else:
        x = [0, 1]
        t = [273.15, 373.15]

    x_v = numpy.linspace(min(x), max(x), num=50)
    t_v = numpy.linspace((min(t) - 10), (max(t) + 10), num=50)
    x_fit, t_fit = numpy.meshgrid(x_v, t_v)
    p_fit = numpy.array(
        [pervaporation_function(x_fit[i], t_fit[i]) for i in range(len(x_fit))]
    )
    ax.plot_surface(x_fit, t_fit, p_fit, alpha=0.2)

    ax.set_xlabel(x_label)
    ax.set_ylabel(t_label)
    ax.set_zlabel(p_label)
    fig.suptitle(suptitle, fontsize=10)
    plt.show()
